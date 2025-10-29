// zpt2phit_8x8x9_binning.cxx
// - Keeps the historical name "zpt2phit_8x8x9" (φ now 8 bins).
// - Restores original pT^2 binning: 10 analysis bins up to 1.0 + overflow 1.0–1.5 (=11 total).
// - Derives N_xq2bins at runtime from makeTH2PolyMars_Valerii().

#include <memory>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <iostream>

#include "TH2Poly.h"
#include "TH1D.h"
#include "TMath.h"

// -----------------------------------------------------------------------------
// Global binning constants / globals
// -----------------------------------------------------------------------------

static constexpr int N_Zbins  = 8;          // z bins
static constexpr int N_pTbins = 10;         // analysis bins up to 1.0
static constexpr int N_pTbins_with_overflow = N_pTbins + 1; // +1 overflow (1.0–1.5) => 11 total
const int N_phiTrbins = 8;                  // φ bins (Trento) – 8 here
int N_xq2bins = 0;                          // set at runtime from makeTH2PolyMars_Valerii()

// -----------------------------------------------------------------------------
// Forward declarations
// -----------------------------------------------------------------------------
struct bin_xBQ2_Valerii;
struct MarsZpt2;
struct Zpt2Phit_8x8x9;
struct IsEventInBins;

struct BinningContext;

TH2Poly* makeTH2PolyMars_Valerii();
TH2Poly* makeTH2PolyZpt2_Valerii(int xq2, const std::vector<double>& z);
std::vector<double> zbins_Valerii(int xqbin, int nXQ2);

BinningContext make_binning_context();
BinningContext& get_binning_context();

// -----------------------------------------------------------------------------
// Functors
// -----------------------------------------------------------------------------

struct bin_xBQ2_Valerii {
  std::shared_ptr<TH2Poly> poly_xQ2;
  int operator()(double x, double Q2) const {
    return static_cast<int>(poly_xQ2->FindBin(x, Q2));
  }
};

struct MarsZpt2 {
  std::shared_ptr<TH2Poly> poly_xQ2;
  std::shared_ptr<std::vector<std::shared_ptr<TH2Poly>>> poly_zpt2;
  int operator()(double x, double Q2, double z, double pt2) const {
    const int b = static_cast<int>(poly_xQ2->FindBin(x, Q2));
    if (b < 1 || b >= static_cast<int>(poly_zpt2->size())) return -1;
    const auto& h = (*poly_zpt2)[b];
    if (!h) return -1;
    return static_cast<int>(h->FindBin(z, pt2));   // 0 => underflow/out-of-range
  }
};

inline double wrap_phi_deg(double a){
  a = std::fmod(a, 360.0);
  if (a < 0) a += 360.0;
  // Map exact 0 to a tiny positive so 0° goes to first φ-bin instead of underflow
  if (a == 0.0) a = std::nextafter(0.0, 1.0);
  return a;
}

struct Zpt2Phit_8x8x9 {
  MarsZpt2 mars_zpt2;
  std::shared_ptr<TH1D> reft9; // name kept; nbins = N_phiTrbins (8)
  int operator()(double x, double Q2, double z, double pt2, double phit) const {
    const int zpt2 = mars_zpt2(x, Q2, z, pt2);
    if (zpt2 < 1) return -1;

    const double phi_wrapped = wrap_phi_deg(phit);
    const int    nphi = reft9->GetNbinsX();          // 8
    const int    phi  = reft9->FindBin(phi_wrapped); // 1..nphi
    // Composite: φ is fastest
    return (zpt2 - 1)*nphi + phi;                    // 1-based contiguous index
  }
};

struct IsEventInBins {
  std::shared_ptr<TH2Poly> poly_xQ2;
  std::shared_ptr<std::vector<std::shared_ptr<TH2Poly>>> poly_zpt2;
  std::shared_ptr<TH1D> reft9;
  bool operator()(double x, double Q2, double z, double pt2, double phit) const {
    const int b = static_cast<int>(poly_xQ2->FindBin(x, Q2));
    if (b < 1 || b >= static_cast<int>(poly_zpt2->size())) return false;

    const auto& h = (*poly_zpt2)[b];
    if (!h) return false;

    const double phi_wrapped = wrap_phi_deg(phit);
    const int    nphi = reft9->GetNbinsX();
    const int    phi  = reft9->FindBin(phi_wrapped);
    if (phi < 1 || phi > nphi) return false;

    const int zpt2bin = h->FindBin(z, pt2);
    return zpt2bin >= 1;
  }
};

// -----------------------------------------------------------------------------
// Bundle
// -----------------------------------------------------------------------------

struct BinningContext {
  std::shared_ptr<TH2Poly> poly_xQ2;                               // x–Q² mosaic
  std::shared_ptr<std::vector<std::shared_ptr<TH2Poly>>> poly_zpt2; // per-(x,Q²) z⊗pT² mosaics
  std::shared_ptr<TH1D> reft9;                                     // φ bins (8)

  bin_xBQ2_Valerii bin_xBQ2;
  MarsZpt2         mars_zpt2;
  Zpt2Phit_8x8x9   zpt2phit_8x8x9;
  IsEventInBins    isEventInBins;
};

// -----------------------------------------------------------------------------
// Builders
// -----------------------------------------------------------------------------

BinningContext make_binning_context() {
  BinningContext ctx;

  // (x,Q²) poly
  ctx.poly_xQ2 = std::shared_ptr<TH2Poly>(makeTH2PolyMars_Valerii());
  ctx.poly_xQ2->SetDirectory(nullptr);
  (void)ctx.poly_xQ2->FindBin(-1e9, -1e9); // warm-up

  // Derive nXQ2 from the poly and keep it globally available
  const int nXQ2_local = ctx.poly_xQ2->GetNumberOfBins();
  N_xq2bins = nXQ2_local;

  // Per-(x,Q²) z⊗pT² mosaics, allow index 0..nXQ2
  ctx.poly_zpt2 = std::make_shared<std::vector<std::shared_ptr<TH2Poly>>>(nXQ2_local + 1);
  for (int i = 0; i <= nXQ2_local; ++i) {
    auto h = std::shared_ptr<TH2Poly>(makeTH2PolyZpt2_Valerii(i, zbins_Valerii(i, nXQ2_local)));
    h->SetDirectory(nullptr);
    (void)h->FindBin(-1e9, -1e9);
    (*ctx.poly_zpt2)[i] = std::move(h);
  }

  // φ bins (8). Name kept as "reft9" for compatibility.
  ctx.reft9 = std::shared_ptr<TH1D>(new TH1D("reft","",N_phiTrbins,0,360));
  ctx.reft9->SetDirectory(nullptr);

  // Functors
  ctx.bin_xBQ2       = bin_xBQ2_Valerii{ctx.poly_xQ2};
  ctx.mars_zpt2      = MarsZpt2{ctx.poly_xQ2, ctx.poly_zpt2};
  ctx.zpt2phit_8x8x9 = Zpt2Phit_8x8x9{ctx.mars_zpt2, ctx.reft9};
  ctx.isEventInBins  = IsEventInBins{ctx.poly_xQ2, ctx.poly_zpt2, ctx.reft9};
  return ctx;
}

BinningContext& get_binning_context() {
  static BinningContext ctx = make_binning_context();
  return ctx;
}

// -----------------------------------------------------------------------------
// Bins geometry
// -----------------------------------------------------------------------------

// (x,Q²) mosaic (unchanged geometry)
TH2Poly* makeTH2PolyMars_Valerii(){
  auto poly = new TH2Poly();

  // 1st row
  double x1[] = {0.126,  0.18,    0.18,     0.126};
  double y1[] = {2,      2.658,   2,        2};

  double x2[] = {0.18,   0.18,    0.21,     0.21,   0.18};
  double y2[] = {2,      2.658,   2.658,    2,      2};

  double x3[] = {0.21,   0.21,    0.24,     0.24,   0.21};
  double y3[] = {2,      2.658,   2.658,    2,      2};

  double x4[] = {0.24,   0.28,    0.28,     0.24};
  double y4[] = {2.658,  2.658,   2.00,     2};

  // Split 0.28–0.45 into three bins
  double x5_a[] = {0.28,   0.28,    0.34,   0.34};
  double y5_a[] = {2.00,   2.658,   2.658,  2};

  double x5_b[] = {0.34,   0.34,    0.39,   0.39};
  double y5_b[] = {2.00,   2.658,   2.658,  2};

  double x5_c[] = {0.39,   0.39,    0.45,   0.39};
  double y5_c[] = {2.00,   2.658,   2.658,  2};

  // 2nd row
  double x6[] = {0.18,   0.24,    0.24,     0.18};
  double y6[] = {2.658,  3.625,   2.658,    2.658};

  double x7[] = {0.24,   0.24,    0.28,     0.28,   0.24};
  double y7[] = {2.658,  3.625,   3.625,    2.658,  2.658};

  double x8[] = {0.28,   0.28,    0.34,     0.34,   0.28};
  double y8[] = {2.658,  3.625,   3.625,    2.658,  2.658};

  double x9_a[] = {0.34,   0.34,    0.39,  0.39};
  double y9_a[] = {2.658,  3.625,   3.625,    2.658};

  double x9_b[] = {0.39,   0.39,    0.54025,  0.45};
  double y9_b[] = {2.658,  3.625,   3.625,    2.658};

  // 3rd row
  double x10[] = {0.24,  0.34,    0.34,     0.24};
  double y10[] = {3.625, 5.12,    3.625,    3.625};

  double x11[] = {0.34,  0.34,    0.39,     0.39,   0.34};
  double y11[] = {3.625, 5.12,    5.12,     3.625,  3.625};

  double x12[] = {0.39,  0.39,    0.6234,   0.57,   0.54025,  0.39};
  double y12[] = {3.625, 5.12,    5.12,     4.05,   3.625,    3.625};

  // 4th row
  double x13[] = {0.34,  0.677,   0.7896,   0.75,   0.708,    0.64,  0.6234, 0.34};
  double y13[] = {5.12,  10.185,  11.351,   9.52,   7.42,     5.4,   5.12,   5.12};

  // 0-row for bin migration (low Q²)
  double x0_1[] = {0.,   0.,    0.1,  0.1};
  double y0_1[] = {1.5,  2.0,   2.0,  1.5};

  double x0_2[] = {0.1,   0.1,    0.2,  0.2};
  double y0_2[] = {1.5,  2.0,   2.0,   1.5};

  double x0_3[] = {0.2,   0.2,    0.3,  0.3};
  double y0_3[] = {1.5,  2.0,   2.0,   1.5};

  double x0_4[] = {0.3,   0.3,    0.4,  0.4};
  double y0_4[] = {1.5,  2.0,   2.0,   1.5};

  poly->AddBin(4, x1, y1);
  poly->AddBin(5, x2, y2);
  poly->AddBin(5, x3, y3);
  poly->AddBin(4, x4, y4);
  poly->AddBin(4, x5_a, y5_a);
  poly->AddBin(4, x5_b, y5_b);
  poly->AddBin(4, x5_c, y5_c);
  poly->AddBin(4, x6, y6);
  poly->AddBin(5, x7, y7);
  poly->AddBin(5, x8, y8);
  poly->AddBin(4, x9_a, y9_a);
  poly->AddBin(4, x9_b, y9_b);
  poly->AddBin(4, x10, y10);
  poly->AddBin(5, x11, y11);
  poly->AddBin(6, x12, y12);
  poly->AddBin(8, x13, y13);

  // Low Q² bin migration (some may be intentionally empty)
  poly->AddBin(4, x0_1, y0_1);
  poly->AddBin(4, x0_2, y0_2);
  poly->AddBin(4, x0_3, y0_3);
  poly->AddBin(4, x0_4, y0_4);

  poly->SetLineStyle(1);
  poly->SetLineWidth(10);

  return poly;
}

// pT^2 edges per z-bin: 10 analysis bins up to 1.0 + overflow to 1.5 (=> 11 total)
TH2Poly* makeTH2PolyZpt2_Valerii(int /*xq2*/, const std::vector<double>& z){
  // p[i] holds the pT^2 edges for z-bin i (0..N_Zbins-1).
  std::vector<std::vector<double>> p(N_Zbins + 1);

  // Original 10-bin analysis edges up to 1.0:
  const std::vector<double> pt2_base = {
    0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.65, 0.80, 1.00
  }; // 11 edges => 10 bins

  // Fill for all z bins
  for (int z_bin_i = 0; z_bin_i < N_Zbins; ++z_bin_i) p[z_bin_i] = pt2_base;
  p[N_Zbins] = pt2_base;

  // Append overflow edge to 1.5 for each z bin (now 12 edges => 11 bins total)
  for (auto& pi : p) pi.push_back(1.50);

  auto poly = new TH2Poly(); poly->Sumw2();
  poly->SetLineStyle(1);
  poly->SetLineWidth(2);

  // Build z × pT^2 rectangles
  const int NzEdges = static_cast<int>(z.size());
  for (int iz = 0; iz < NzEdges - 1; ++iz) {
    const auto& pz = p[iz]; // pT^2 edges for this z-bin
    for (int j = 0; j < static_cast<int>(pz.size()) - 1; ++j) {
      poly->AddBin(z[iz], pz[j], z[iz+1], pz[j+1]);
    }
  }

  return poly;
}

// Uniform z edges (8 bins => 9 edges)
std::vector<double> zbins_Valerii(int xqbin, int nXQ2){
  if (xqbin < 0 || xqbin > nXQ2)
    throw std::out_of_range("zbins_Valerii: xqbin out of range");
  return {0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0};
}
