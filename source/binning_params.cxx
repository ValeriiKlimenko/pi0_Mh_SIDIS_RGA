/////////////////////////////////////////////////////////////////
/////////////////////////////NUMBER OF BINS
/////////////////////////////////////////////////////////////////

const int N_xq2bins = 16 + 3; //4 is low Q2 migration 1.5-2 GeV2, four bins in xB
//PT binning 8 -> 13 ->10
//analysis binnning, overflow is not included
const int N_pTbins = 10;
// includes 1 - 1.5 pt2 overflow bin
// it is actual total bin # when zpt2phit_8x8x9 is called because it includes overflow bin.
const int N_pTbins_with_overflow = N_pTbins + 1; 
const int N_phiTrbins = 9;
// 10 -> 8, and one overflow
// z binnning already includes overflow for whatever reason (0.8 - 1.0 )
const int N_Zbins = 8;

/// for pi0 mass fit:
const int N_pi0mm_bins = 60;
const double min_pi0mass = 0, max_pi0mass = 0.4; 

// number of sigmas in pi0mass fits for Integral Err. Calculation
// It is predefined parameter in extractPi0 - pi02D functions
const double n_sigma_pi0mass = 5.5;


//////////////////////

// functor types
struct bin_xBQ2_Valerii {
  std::shared_ptr<TH2Poly> poly_xQ2;
  int operator()(double x, double Q2) const;
};
struct MarsZpt2 {
  std::shared_ptr<TH2Poly> poly_xQ2;
  std::shared_ptr<std::vector<std::shared_ptr<TH2Poly>>> poly_zpt2;
  int operator()(double x, double Q2, double z, double pt2) const;
};
struct Zpt2Phit_8x8x9 {
  MarsZpt2 mars_zpt2;
  std::shared_ptr<TH1D> reft9;
  int operator()(double x, double Q2, double z, double pt2, double phit) const;
};
struct IsEventInBins {
  std::shared_ptr<TH2Poly> poly_xQ2;
  std::shared_ptr<std::vector<std::shared_ptr<TH2Poly>>> poly_zpt2;
  std::shared_ptr<TH1D> reft9;
  bool operator()(double x, double Q2, double z, double pt2, double phit) const;
};

// bundle
struct BinningContext {
  std::shared_ptr<TH2Poly> poly_xQ2;
  std::shared_ptr<std::vector<std::shared_ptr<TH2Poly>>> poly_zpt2;
  std::shared_ptr<TH1D> reft9;

  bin_xBQ2_Valerii bin_xBQ2;
  MarsZpt2         mars_zpt2;
  Zpt2Phit_8x8x9   zpt2phit_8x8x9;
  IsEventInBins    isEventInBins;
};

// creates all objects and functor *instances*
BinningContext make_binning_context();

// convenient singleton-style accessor
BinningContext& get_binning_context();


/////////////////////


// extern existing makers and constants
extern TH2Poly* makeTH2PolyMars_Valerii();
extern TH2Poly* makeTH2PolyZpt2_Valerii(int xq2, std::vector<double> z);
extern std::vector<double> zbins_Valerii(int xqbin);
extern const int N_xq2bins, N_phiTrbins;

int bin_xBQ2_Valerii::operator()(double x, double Q2) const {
  return static_cast<int>(poly_xQ2->FindBin(x, Q2));
}
int MarsZpt2::operator()(double x, double Q2, double z, double pt2) const {
  int b = static_cast<int>(poly_xQ2->FindBin(x, Q2));
  if (b < 1) b = 0;
  if (b >= static_cast<int>(poly_zpt2->size())) return -1;
  auto& h = (*poly_zpt2)[b];
  return h ? static_cast<int>(h->FindBin(z, pt2)) : -1;
}
int Zpt2Phit_8x8x9::operator()(double x, double Q2, double z, double pt2, double phit) const {
  int zpt2 = mars_zpt2(x, Q2, z, pt2);
  if (zpt2 < 0) return -1;
  int phi = reft9->FindBin(phit), nphi = reft9->GetNbinsX();
  if (phi < 1 || phi > nphi) return -1;
  return nphi * zpt2 + phi;
}
bool IsEventInBins::operator()(double x, double Q2, double z, double pt2, double phit) const {
  //if (!poly_xQ2 || !poly_zpt2 || !reft9) return true;
  int b = static_cast<int>(poly_xQ2->FindBin(x, Q2));
  if (b < 1) return false;
  
  if (b >= static_cast<int>(poly_zpt2->size())){ 
    // cout << b << endl; 
    return false;
  }
  
  auto& h = (*poly_zpt2)[b];
  //if (!h) return true;
  int phi = reft9->FindBin(phit), nphi = reft9->GetNbinsX();
  if (phi < 1 || phi > nphi) return false;
  int zpt2bin = h->FindBin(z, pt2);
  return zpt2bin >= 1;
}

BinningContext make_binning_context() {
  BinningContext ctx;

  ctx.poly_xQ2 = std::shared_ptr<TH2Poly>(makeTH2PolyMars_Valerii());
  ctx.poly_xQ2->SetDirectory(nullptr);
  (void)ctx.poly_xQ2->FindBin(-1e9, -1e9);

  ctx.poly_zpt2 = std::make_shared<std::vector<std::shared_ptr<TH2Poly>>>(N_xq2bins + 1);
  for (int i = 0; i <= N_xq2bins; ++i) {
    auto h = std::shared_ptr<TH2Poly>(makeTH2PolyZpt2_Valerii(i, zbins_Valerii(i)));
    h->SetDirectory(nullptr);
    (void)h->FindBin(-1e9, -1e9);
    (*ctx.poly_zpt2)[i] = std::move(h);
  }

  ctx.reft9 = std::shared_ptr<TH1D>(new TH1D("reft9","",N_phiTrbins,0,360));
  ctx.reft9->SetDirectory(nullptr);

  ctx.bin_xBQ2       = bin_xBQ2_Valerii{ctx.poly_xQ2};
  ctx.mars_zpt2      = MarsZpt2{ctx.poly_xQ2, ctx.poly_zpt2};
  ctx.zpt2phit_8x8x9 = Zpt2Phit_8x8x9{ctx.mars_zpt2, ctx.reft9};
  ctx.isEventInBins  = IsEventInBins{ctx.poly_xQ2, ctx.poly_zpt2, ctx.reft9};
  return ctx;
}

BinningContext& get_binning_context() {
  static BinningContext ctx = make_binning_context(); // created once, thread-safe
  return ctx;
}



/////////////////////////////////////////////////////////////////
/////////////////////////////BINS EDGES:
/////////////////////////////////////////////////////////////////

//Fuction generates a 2D histogram of the current xQ2 binning
//three additionall bins were added

//13 bins -> 16 bins now
 TH2Poly* makeTH2PolyMars_Valerii(){
    auto poly = new TH2Poly();
   
    //1st row
    double x1[] = {0.126,  0.18,    0.18,     0.126};//0.18
    double y1[] = {2,      2.658,   2,        2};

    double x2[] = {0.18,   0.18,    0.21,     0.21,   0.18};//.18
    double y2[] = {2,      2.658,   2.658,    2,      2};

    //It looks like one point appears two times (1st and last) for some reason it still rect.
    double x3[] = {0.21,   0.21,    0.24,     0.24,   0.21};//.18
    double y3[] = {2,      2.658,   2.658,    2,      2};

    double x4[] = {0.24,   0.28,    0.28,     0.24};
    double y4[] = {2.658,  2.658,   2.00,     2};

    //double x5[] = {0.28,   0.28,    0.45,     0.45,   0.387826, 0.28};
    //double y5[] = {2.00,   2.658,   2.658,    2.658,  2,        2};

    //Split by Valerii into two bins:
    double x5_a[] = {0.28,   0.28,    0.34,   0.34};
    double y5_a[] = {2.00,   2.658,   2.658,  2};

    double x5_b[] = {0.34,   0.34,    0.39,   0.39};
    double y5_b[] = {2.00,   2.658,   2.658,  2};

    double x5_c[] = {0.39,   0.39,    0.45,   0.39};
    double y5_c[] = {2.00,   2.658,   2.658,  2};
  
    //2nd row
    double x6[] = {0.18,   0.24,    0.24,     0.18};//.18
    double y6[] = {2.658,  3.625,   2.658,    2.658};

    double x7[] = {0.24,   0.24,    0.28,     0.28,   0.24};
    double y7[] = {2.658,  3.625,   3.625,    2.658,  2.658};
   
    double x8[] = {0.28,   0.28,    0.34,     0.34,   0.28};
    double y8[] = {2.658,  3.625,   3.625,    2.658,  2.658};

    //double x9[] = {0.34,   0.34,    0.54025,  0.45,   0.34};
    //double y9[] = {2.658,  3.625,   3.625,    2.658,  2.658};

    //Split by Valerii into two bins:
    double x9_a[] = {0.34,   0.34,    0.39,  0.39};
    double y9_a[] = {2.658,  3.625,   3.625,    2.658};

    double x9_b[] = {0.39,   0.39,    0.54025,  0.45};
    double y9_b[] = {2.658,  3.625,   3.625,    2.658};
    //3rd row
    double x10[] = {0.24,  0.34,    0.34,     0.24};
    double y10[] = {3.625, 5.12,    3.625,    3.625};

    double x11[] = {0.34,  0.34,    0.39,     0.39,   0.34};
    double y11[] = {3.625, 5.12,    5.12,     3.625,  3.625};

    double x12[] = {0.39,  0.39,    0.6234,   0.57,   0.54025,  0.39};
    double y12[] = {3.625, 5.12,    5.12,     4.05,   3.625,    3.625};
    //4th row
    double x13[] = {0.34,  0.677,   0.7896,   0.75,   0.708,    0.64,  0.6234, 0.34};
    double y13[] = {5.12,  10.185,  11.351,   9.52,   7.42,     5.4,   5.12,   5.12};

    // 0-row for bin migration
    double x0_1[] = {0.,   0.,    0.1,  0.1};
    double y0_1[] = {1.5,  2.0,   2.0,   1.5};

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
    //poly->AddBin(6, x5, y5);
    poly->AddBin(4, x5_a, y5_a);
    poly->AddBin(4, x5_b, y5_b);
    poly->AddBin(4, x5_c, y5_c);
    poly->AddBin(4, x6, y6);
    poly->AddBin(5, x7, y7);
    poly->AddBin(5, x8, y8);
    //poly->AddBin(5, x9, y9);
    poly->AddBin(4, x9_a, y9_a);
    poly->AddBin(4, x9_b, y9_b);
    poly->AddBin(4, x10, y10);
    poly->AddBin(5, x11, y11);
    poly->AddBin(6, x12, y12);
    poly->AddBin(8, x13, y13);

    // low Q2 bin migration:

    // this bin is always empty
    // remove at the next full chain rerun:
    poly->AddBin(4, x0_1, y0_1);

   
    poly->AddBin(4, x0_2, y0_2);
    poly->AddBin(4, x0_3, y0_3);
    poly->AddBin(4, x0_4, y0_4);
   

    poly->SetLineStyle(1);
    poly->SetLineWidth(10);

    return poly;
 }

// pt2bins Valerii's Version
vector<vector<double>> pbins(int &xq2){
   vector<vector<double>> p(N_Zbins + 1);
   // bins were deleted to acount for resolution as pt2 dependence
   //p[N_Zbins] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};//z overflow
   p[N_Zbins] = {0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.65,0.8,1};//z overflow

   if (xq2>= 0 && xq2 <= N_xq2bins){
     for (int z_bin_i = 0; z_bin_i < N_Zbins; z_bin_i ++)
       p[z_bin_i] = {0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.65,0.8,1};
   }else {
     throw out_of_range("xq2 bin is out of range in makeTH2PolyZpt2_Valerii function, check binnnig");
   }
   
   //pt overflow bin; pt > ptmx in z bin
   for (int i = 0; i < p.size(); i++){
     p[i].push_back(1.5);
   }

   vector<vector<double>> pf(p.size() + 1);
   for(int i = 0; i < p.size(); i++){
     pf[i+1] = p[i]; 
   }
   pf[0] = p[8];
   return pf;
}

// Valerii's version of pT binnig. 
//The bins borders did not match in pT variable. Now it does match. 

TH2Poly* makeTH2PolyZpt2_Valerii(int xq2, vector<double> z){
   // it means 10 Z bins + 1 additional
   // 13 pt bins + one additional for pt > 1
   vector<vector<double>> p(N_Zbins + 1);
   p[N_Zbins] = {0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.65,0.8,1};

   if (xq2>= 0 && xq2 <= N_xq2bins){
     for (int z_bin_i = 0; z_bin_i < N_Zbins; z_bin_i ++)
       p[z_bin_i] = {0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.65,0.8,1};
   }else {
     throw out_of_range("xq2 bin is out of range in makeTH2PolyZpt2_Valerii function, check binnnig");
   }
   //pt overflow bin; pt > ptmx in z bin
   for (int i = 0; i < p.size(); i++){
     p[i].push_back(1.5);
   }

   auto poly = new TH2Poly(); poly->Sumw2();
   poly->SetLineStyle(1);
   poly->SetLineWidth(2);

   for (int i = 0; i < z.size()-1; i++){
     if (i == 0){
       for (int j = 0; j < p[i].size()-1; j++){
         poly->AddBin(z[i], p[0][j], z[i+1], p[0][j+1]);
       }
     }
     else{
       for (int j = 0; j < p[i-1].size()-1; j++){
         //p first index is z bin, the second is 
         poly->AddBin(z[i], p[i-1][j], z[i+1], p[i-1][j+1]);
       }
     }
   }

    return poly;
 }

//Valerii, z bins are uniform:
vector<double> zbins_Valerii(int xqbin){
  vector<double> res;

  //if (xqbin==1) {res = {0,0.2,0.275,0.350,0.425,0.500,0.575,0.650,0.725,0.8,1};}

  if (xqbin>=0 && xqbin < N_xq2bins + 1){
    res = {0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1};
    return res;
  }
  else throw out_of_range("zbins_Valerii is out of range of xq2 bins");

  return {};
}