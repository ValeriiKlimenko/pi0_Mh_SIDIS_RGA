#!/usr/bin/env python3
"""

"""

import argparse
import os
import shutil
import subprocess
import sys

# Fixed settings
BASE_DATA_DIR = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_data/"
MACRO_PATH    = "source/define_bin_migr_hists_data.cxx"
COMPILE_MACRO = False  # use ACLiC (+). Set False to interpret.

# Map MX cut modes to directory names
MX_DIR_MAP = {
    0: "data_noMxcut/",
    1: "data_Mxcut_1/",
    2: "rec_data/",
}

def _q(s: str) -> str:
    """Escape for embedding in ROOT command lines."""
    return s.replace("\\", "\\\\").replace('"', '\\"')

def run_with_root_cli_link_first(dir_, outdir):
    if shutil.which("root") is None:
        raise SystemExit("ERROR: 'root' executable not found in PATH.")
    plus = "+" if COMPILE_MACRO else ""
    cmds = [
        f'.L {MACRO_PATH}{plus}',
        f'define_bin_migr_hists_data("{_q(dir_)}","{_q(outdir)}")',
        '.q',
    ]
    cmd_stream = "\n".join(cmds) + "\n"

    print("ROOT command stream:\n" + "\n".join("  " + c for c in cmds))

    proc = subprocess.run(
        ["root", "-l", "-b"],
        input=cmd_stream,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise SystemExit(f"ERROR: ROOT returned non-zero exit status {proc.returncode}")

def main():
    p = argparse.ArgumentParser(
        description="Run define_bin_migr_hists_data (ROOT CLI; link first)."
    )

    p.add_argument(
        "--mx-cut-mode",
        type=int,
        choices=[0, 1, 2],
        default=1,
        help=(
            "Missing mass cut mode / default directory:\n"
            "  0 = no cut       -> data_noMxcut/\n"
            "  1 = cut at 1.0   -> data_Mxcut_1/\n"
            "  2 = cut at 1.5   -> rec_data/\n"
            "Ignored if --dir is explicitly provided."
        ),
    )

    p.add_argument(
        "--dir",
        default=None,
        help=(
            "Optional: explicit subdirectory under base data dir "
            "(e.g. data_noMxcut/, data_Mxcut_1/, rec_data/). "
            "If omitted, directory is chosen from --mx-cut-mode."
        ),
    )

    args = p.parse_args()

    # Decide directory: explicit --dir wins, otherwise use mx-cut-mode mapping
    if args.dir is not None:
        dir_name = (args.dir or "").strip().strip("/") + "/"
        mx_info = f"(overridden by --dir, mx_cut_mode={args.mx_cut_mode})"
    else:
        dir_name = MX_DIR_MAP[args.mx_cut_mode]
        mx_info = f"(from mx_cut_mode={args.mx_cut_mode})"

    full_dir = os.path.join(BASE_DATA_DIR, dir_name)
    outdir = f"unfolding_{dir_name}"  # e.g. unfolding_data_noMxcut/

    if not os.path.isdir(os.path.dirname(full_dir)) and not os.path.isdir(full_dir):
        print(f"Warning: constructed data dir does not exist: {full_dir}", file=sys.stderr)

    print(f"Base data dir : {BASE_DATA_DIR}")
    print(f"Directory     : {dir_name} {mx_info}")
    print(f"Full input dir: {full_dir}")
    print(f"Output folder : {outdir}")
    print(f"Macro         : {MACRO_PATH}  ({'compile +' if COMPILE_MACRO else 'interpret'})")

    run_with_root_cli_link_first(full_dir, outdir)
    print("Done.")

if __name__ == "__main__":
    main()
