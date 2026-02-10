#!/usr/bin/env python3
import argparse
import math
import sys
from pathlib import Path

import pandas as pd
import numpy as np


REQ_COLS = [
    "xq2_bin", "ybin", "z_bin", "pt2_bin", "phi_bin",
    "phi_center_deg", "content", "error"
]

KEYS_FULL = ["ybin", "z_bin", "pt2_bin", "phi_bin"]
KEYS_NOPHI = ["ybin", "z_bin", "pt2_bin"]


def read_hist_csv(path: str, label: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"{label} file not found: {path}")
    df = pd.read_csv(p)

    missing = [c for c in REQ_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"{label} CSV missing columns {missing}. Found columns: {list(df.columns)}")

    # Ensure numeric dtypes
    for c in ["xq2_bin", "ybin", "z_bin", "pt2_bin", "phi_bin"]:
        df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
    for c in ["phi_center_deg", "content", "error"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Drop rows with invalid keys
    df = df.dropna(subset=["ybin", "z_bin", "pt2_bin", "phi_bin"]).copy()

    # Collapse potential duplicates by summing content and summing errors in quadrature
    # (safe even if there are no duplicates)
    g = df.groupby(KEYS_FULL, as_index=False)
    df2 = g.agg(
        xq2_bin=("xq2_bin", "first"),
        phi_center_deg=("phi_center_deg", "first"),
        content=("content", "sum"),
        error=("error", lambda s: float(np.sqrt(np.sum(np.square(s.fillna(0.0))))))
    )

    # Cast ints back to normal int (pandas Int64 -> python int where possible)
    for c in ["xq2_bin", "ybin", "z_bin", "pt2_bin", "phi_bin"]:
        df2[c] = df2[c].astype(int)

    return df2


def build_full_grid(df_list, nphi: int) -> pd.DataFrame:
    # Union of all (y,z,pt2) combinations present in any input
    base = pd.concat([d[KEYS_NOPHI] for d in df_list], ignore_index=True).drop_duplicates()
    base["__k"] = 1
    phi = pd.DataFrame({"phi_bin": np.arange(1, nphi + 1, dtype=int), "__k": 1})
    grid = base.merge(phi, on="__k").drop(columns="__k")
    return grid


def attach_xq2_mapping(grid: pd.DataFrame, df_list) -> pd.DataFrame:
    # Map ybin -> xq2_bin (take first available from any df)
    maps = []
    for d in df_list:
        maps.append(d[["ybin", "xq2_bin"]].drop_duplicates())
    m = pd.concat(maps, ignore_index=True).dropna().drop_duplicates(subset=["ybin"])
    out = grid.merge(m, on="ybin", how="left")
    return out


def merge_hist(grid: pd.DataFrame, df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    # Merge content/error onto full grid; fill missing with 0
    tmp = df[KEYS_FULL + ["content", "error"]].rename(
        columns={"content": f"{prefix}_content", "error": f"{prefix}_error"}
    )
    out = grid.merge(tmp, on=KEYS_FULL, how="left")
    out[f"{prefix}_content"] = out[f"{prefix}_content"].fillna(0.0)
    out[f"{prefix}_error"] = out[f"{prefix}_error"].fillna(0.0)
    return out


def ratio_and_error(D, eD, T, eT, M, eM, include_mc_stat: bool, zero_when_div0: bool):
    # Vectorized safe ratio
    M_is_zero = (np.abs(M) <= 0.0)

    val = np.zeros_like(D, dtype=float)
    err = np.zeros_like(D, dtype=float)

    good = ~M_is_zero
    if np.any(good):
        SF = T[good] / M[good]
        val_good = D[good] * SF

        if not include_mc_stat:
            err_good = np.abs(SF) * eD[good]
        else:
            # Full propagation for f = D*T/M
            termD = (T[good] / M[good]) * eD[good]
            termT = (D[good] / M[good]) * eT[good]
            termM = (T[good] * D[good] / (M[good] * M[good])) * eM[good]
            err_good = np.sqrt(termD * termD + termT * termT + termM * termM)

        val[good] = val_good
        err[good] = err_good

    if not zero_when_div0:
        val[M_is_zero] = np.nan
        err[M_is_zero] = np.nan

    return val, err


def method1_sum_phi_then_ratio(df: pd.DataFrame, nphi: int,
                              include_mc_stat: bool, zero_when_div0: bool) -> pd.DataFrame:
    """
    NEW Method 1:
      1) For each phi bin compute R_phi = D*T/M (with same error propagation as ratio_and_error)
         using only rows where D>0, T>0, M>0.
      2) Sum R_phi over phi within each (ybin,z_bin,pt2_bin).
      3) Error on the sum: sqrt(sum(eR_phi^2)).
    """

    # keep only phi rows where all three are strictly > 0
    ok = (df["D_content"] > 0.0) & (df["T_content"] > 0.0) & (df["M_content"] > 0.0)
    df_ok = df.loc[ok].copy()

    if df_ok.empty:
        return pd.DataFrame(columns=["xq2_bin", "ybin", "z_bin", "pt2_bin", "nphi", "content", "error"])

    # per-phi ratio
    R_phi, eR_phi = ratio_and_error(
        df_ok["D_content"].to_numpy(), df_ok["D_error"].to_numpy(),
        df_ok["T_content"].to_numpy(), df_ok["T_error"].to_numpy(),
        df_ok["M_content"].to_numpy(), df_ok["M_error"].to_numpy(),
        include_mc_stat=include_mc_stat,
        zero_when_div0=zero_when_div0
    )
    df_ok["R_phi"] = R_phi
    df_ok["eR_phi"] = eR_phi

    # if user requested NaN when div0, then drop NaNs before summing
    # (when zero_when_div0=True, div0 already gave 0 so it’s harmless)
    df_ok = df_ok[np.isfinite(df_ok["R_phi"]) & np.isfinite(df_ok["eR_phi"])].copy()
    if df_ok.empty:
        return pd.DataFrame(columns=["xq2_bin", "ybin", "z_bin", "pt2_bin", "nphi", "content", "error"])

    # sum over phi of the per-phi ratios
    g = df_ok.groupby(KEYS_NOPHI, as_index=False)
    summed = g.agg(
        xq2_bin=("xq2_bin", "first"),
        n_used=("phi_bin", "count"),
        content=("R_phi", "sum"),
        error=("eR_phi", lambda s: float(np.sqrt(np.sum(np.square(s.to_numpy(dtype=float)))))),
    )

    out = summed[["xq2_bin", "ybin", "z_bin", "pt2_bin"]].copy()
    out["nphi"] = summed["n_used"]   # how many phi bins actually contributed
    out["content"] = summed["content"]
    out["error"] = summed["error"]
    return out


def method2_ratio_then_avg_phi(df: pd.DataFrame, nphi: int, avg: str,
                              include_mc_stat: bool, zero_when_div0: bool) -> pd.DataFrame:
    # Per-phi ratio
    val, err = ratio_and_error(
        df["D_content"].to_numpy(), df["D_error"].to_numpy(),
        df["T_content"].to_numpy(), df["T_error"].to_numpy(),
        df["M_content"].to_numpy(), df["M_error"].to_numpy(),
        include_mc_stat=include_mc_stat,
        zero_when_div0=zero_when_div0
    )
    df2 = df.copy()
    df2["R"] = val
    df2["eR"] = err

    # Average over phi:
    # - unweighted: mean(R), err_mean = sqrt(sum(eR^2))/Nphi
    # - weighted: weights=1/eR^2, Rw = sum(wR)/sum(w), eRw = sqrt(1/sum(w))
    g = df2.groupby(KEYS_NOPHI, as_index=False)

    if avg == "weighted":
        # Protect against eR=0 -> infinite weight. Use 0 weight if eR<=0.
        df2["w"] = 0.0
        mask = np.isfinite(df2["eR"]) & (df2["eR"] > 0.0) & np.isfinite(df2["R"])
        df2.loc[mask, "w"] = 1.0 / np.square(df2.loc[mask, "eR"])

        agg = g.agg(
            xq2_bin=("xq2_bin", "first"),
            sumw=("w", "sum"),
            sumwR=("w", lambda s: float(np.sum(s * df2.loc[s.index, "R"]))),
        )
        out = agg[["xq2_bin", "ybin", "z_bin", "pt2_bin"]].copy()
        out["nphi"] = nphi
        out["content"] = np.where(agg["sumw"] > 0.0, agg["sumwR"] / agg["sumw"], 0.0 if zero_when_div0 else np.nan)
        out["error"] = np.where(agg["sumw"] > 0.0, np.sqrt(1.0 / agg["sumw"]), 0.0 if zero_when_div0 else np.nan)
        return out

    # default: unweighted
    agg = g.agg(
        xq2_bin=("xq2_bin", "first"),
        meanR=("R", "mean"),
        sum_eR2=("eR", lambda s: float(np.nansum(np.square(s))))
    )
    out = agg[["xq2_bin", "ybin", "z_bin", "pt2_bin"]].copy()
    out["nphi"] = nphi
    out["content"] = agg["meanR"].to_numpy()
    out["error"] = np.sqrt(agg["sum_eR2"].to_numpy()) / float(nphi)
    return out

def method3_sum_phi_all_then_ratio(df: pd.DataFrame, nphi: int,
                                  include_mc_stat: bool, zero_when_div0: bool) -> pd.DataFrame:
    """
    Sum over phi for D/T/M (including zeros and missing-filled-as-0),
    then compute ratio: (Σφ D) * (Σφ T) / (Σφ M).
    No additional division by nphi or n_used.
    """
    g = df.groupby(KEYS_NOPHI, as_index=False)

    summed = g.agg(
        xq2_bin=("xq2_bin", "first"),
        D=("D_content", "sum"),
        eD=("D_error", lambda s: float(np.sqrt(np.sum(np.square(s))))),
        T=("T_content", "sum"),
        eT=("T_error", lambda s: float(np.sqrt(np.sum(np.square(s))))),
        M=("M_content", "sum"),
        eM=("M_error", lambda s: float(np.sqrt(np.sum(np.square(s)))))
    )

    val, err = ratio_and_error(
        summed["D"].to_numpy(), summed["eD"].to_numpy(),
        summed["T"].to_numpy(), summed["eT"].to_numpy(),
        summed["M"].to_numpy(), summed["eM"].to_numpy(),
        include_mc_stat=include_mc_stat,
        zero_when_div0=zero_when_div0
    )

    out = summed[["xq2_bin", "ybin", "z_bin", "pt2_bin"]].copy()
    out["nphi"] = nphi
    out["content"] = val
    out["error"] = err
    return out



def main():
    ap = argparse.ArgumentParser(
        description="Build data*gen/rec from split CSVs in two phi-handling ways."
    )
    ap.add_argument("--data", default="tables_h_meas_data.csv", help="CSV for h_meas_data")
    ap.add_argument("--gen", default="tables_h_true_mc.csv", help="CSV for h_true_mc")
    ap.add_argument("--rec", default="tables_h_meas_mc.csv", help="CSV for h_meas_mc")
    ap.add_argument("--out_prefix", default="dataGenOverRec", help="Output prefix for CSVs")
    ap.add_argument("--nphi", type=int, default=0, help="Number of phi bins; 0 => infer from max(phi_bin)")
    ap.add_argument("--avg", choices=["unweighted", "weighted"], default="unweighted",
                    help="Averaging type for method2")
    ap.add_argument("--include_mc_stat", action="store_true",
                    help="If set: propagate GEN+REC errors too (full propagation). Default: data-only like your h_out.")
    ap.add_argument("--nan_when_div0", action="store_true",
                    help="If set: when REC=0 output NaN (instead of 0). Default: 0 like your macro.")
    args = ap.parse_args()

    dfD = read_hist_csv(args.data, "DATA")
    dfT = read_hist_csv(args.gen,  "GEN/TRUE")
    dfM = read_hist_csv(args.rec,  "REC/MEAS")

    nphi = args.nphi
    if nphi <= 0:
        nphi = int(max(dfD["phi_bin"].max(), dfT["phi_bin"].max(), dfM["phi_bin"].max()))
        if nphi <= 0:
            raise ValueError("Could not infer nphi. Please pass --nphi explicitly.")

    grid = build_full_grid([dfD, dfT, dfM], nphi=nphi)
    grid = attach_xq2_mapping(grid, [dfD, dfT, dfM])

    # Compute phi centers deterministically (optional but nice to have in merged df)
    grid["phi_center_deg"] = (grid["phi_bin"].astype(float) - 0.5) * (360.0 / float(nphi))

    merged = merge_hist(grid, dfD, "D")
    merged = merge_hist(merged, dfT, "T")
    merged = merge_hist(merged, dfM, "M")

    include_mc_stat = bool(args.include_mc_stat)
    zero_when_div0 = not bool(args.nan_when_div0)

    out1 = method1_sum_phi_then_ratio(merged, nphi, include_mc_stat, zero_when_div0)
    out2 = method2_ratio_then_avg_phi(merged, nphi, args.avg, include_mc_stat, zero_when_div0)
    out3 = method3_sum_phi_all_then_ratio(merged, nphi, include_mc_stat, zero_when_div0)
    
    p1 = f"{args.out_prefix}_method1_sumPhi_then_ratio.csv"
    p2 = f"{args.out_prefix}_method2_ratio_then_avgPhi_{args.avg}.csv"
    p3 = f"{args.out_prefix}_method3_sumPhi_all_then_ratio.csv"
    
    out1.to_csv(p1, index=False)
    out2.to_csv(p2, index=False)
    out3.to_csv(p3, index=False)
    
    print(f"Wrote:\n  {p1}\n  {p2}\n  {p3}")

if __name__ == "__main__":
    main()
