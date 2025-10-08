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
  if (!poly_xQ2 || !poly_zpt2 || !reft9) return true;
  int b = static_cast<int>(poly_xQ2->FindBin(x, Q2));
  if (b < 1) return false;
  if (b >= static_cast<int>(poly_zpt2->size())) return true;
  auto& h = (*poly_zpt2)[b];
  if (!h) return true;
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

/////////////////////////////////////////////////////////////////
/////////////////////////////INITIALIZATION://///////////////////
/////////////////////////////////////////////////////////////////




//Returns bin number in the xB vs. Q^2 binning developed by Marshall but added couple of bins by Valerii
  // For Mar's one see Table 4.2 in clas12 Sidis analysis note
  // Valerii's not present in the note yet
  //For whatever reason I was geting pointer to int instead of int while cern root doc says it returns int_t

  // now  (10/2025) it is a functor function in the main programm

  //int bin_xBQ2_Valerii(const double &x, const double &Q2)


//Returns the xQ2zpt2 bin number within the th2poly vec
// Valerii: zpPolyVec was updated to be uniform, poly_valerii replaced

// changed to functor 10/2025
/*

auto mars_zpt2(double &x, double &Q2, double &z, double &pt2){
   int mars_bin = poly_valerii->FindBin(x, Q2);
   if (mars_bin < 1){mars_bin = 0;}
  // retrun value should be 8*14 to match "int zpt_bins" from run_unf...cxx
  //if (mars_bin == NULL) cout <<"NULL in mars_zpt"<<x <<' ' << Q2 << endl;
  if (zpPolyVec.size() <= mars_bin) cout << mars_bin;
   return zpPolyVec[mars_bin]->FindBin(z, pt2); 
};

*/

//Returns zpt2phit bin number N_phiTrbins bins
  // Valerii: here is the binning that is in the Marshall analysis on taking over moment 

  // Updated by Valerii

  // now is functor 10/2025
/*
 auto zpt2phit_8x8x9(double &x, double &Q2, double &z, double &pt2, double &phit){
   int mars_bin = poly_valerii->FindBin(x, Q2);
   if (mars_bin < 1){mars_bin = 0;}
    
   //9 phi bins, 9 pt bins(8 real and 1 overflow), 10 z bins (8 real and 1 underflow and 1 overflow)
   // Z bins N was decreased by 2, 8->6, so it is 6 + 2 under and overflow
   
   // Valerii '9'* mars because there are 9 phit bins
   // mars_zpt2 is already updated by Valerii above
   return N_phiTrbins*mars_zpt2(x,Q2,z,pt2) + reft9->FindBin(phit);
};

*/

// check if the bin is covered by existing grid.
// for all the 

// converted into a functor
/*
bool GetisEventINbins(const double &x, const double &Q2,
                      const double &z, const double &pt2,
                      const double &phit)
{


  /*
  // --- configuration guards (requested behavior: return true on failure) ---
  if (!poly_valerii || !reft9) return true;

  auto bins2D = [](TH2 *h) -> int {
    if (!h) return 0;
    // TH2Poly has GetNumberOfBins(); fall back to X*Y for regular TH2
    if (auto *hp = dynamic_cast<TH2Poly*>(h)) return hp->GetNumberOfBins();
    return h->GetNbinsX() * h->GetNbinsY();
  };

  const int nBinsMars = bins2D(static_cast<TH2*>(poly_valerii));
  const int nBinsPhi  = reft9->GetNbinsX();

  if (nBinsMars <= 2 || nBinsPhi <= 2) return true;
  */
/*
  // --- normal logic ---
  int mars_bin = poly_valerii->FindBin(x, Q2);
  if (mars_bin < 1) return false;

  // zpPolyVec[mars_bin] must exist and have >2 bins; else return true
  if (mars_bin >= static_cast<int>(zpPolyVec.size()) || !zpPolyVec[mars_bin])
    return true;

  //const int nBinsZpt2 = bins2D(static_cast<TH2*>(zpPolyVec[mars_bin]));
  //if (nBinsZpt2 <= 2) return true;

  // phi bin check
  int phi_bin_tmp = reft9->FindBin(phit);
  if (phi_bin_tmp < 1 || phi_bin_tmp > reft9->GetNbinsX()) return false;

  // z-pt2 bin check
  int zpt2_bin_tmp = zpPolyVec[mars_bin]->FindBin(z, pt2);
  if (zpt2_bin_tmp < 1) return false;

  return true;
}
*/

 //Function generates a vector of TLines to draw the outline for the current binning
 vector<TLine> makeValeriiOutline(){
    vector<TLine> tvec(52);

vector<string> cvec = {
  "",
  "Q2 < 2.658 && xB < 0.18",
  "Q2 < 2.658 && xB > 0.18 && xB < 0.21",
  "Q2 < 2.658 && xB > 0.21 && xB < 0.24",
  "Q2 < 2.658 && xB > 0.24 && xB < 0.28",
  "Q2 < 2.658 && xB > 0.28",

  "Q2 > 2.658 && Q2 < 3.625 && xB < 0.24",
  "Q2 > 2.658 && Q2 < 3.625 && xB > 0.24 && xB < 0.28",
  "Q2 > 2.658 && Q2 < 3.625 && xB > 0.28 && xB < 0.34",
  "Q2 > 2.658 && Q2 < 3.625 && xB > 0.34",

  "Q2 > 3.625 && Q2 < 5.12 && xB < 0.34",
  "Q2 > 3.625 && Q2 < 5.12 && xB > 0.34 && xB < 0.40",
  "Q2 > 3.625 && Q2 < 5.12 && xB > 0.40",

  "Q2 > 5.12"
};

    //bin 1
    tvec[0] = TLine(0.126, 2, 0.18, 2.658); //.18
    tvec[1] = TLine(0.18, 2.658, 0.18, 2); //.18
    tvec[2] = TLine(0.18, 2., 0.126, 2);     //.18

    //bin 2
    tvec[3] = TLine(0.18, 2, 0.18, 2.658); //.18
    tvec[4] = TLine(0.18, 2.658, 0.21, 2.658); //.18
    tvec[5] = TLine(0.21, 2.658, 0.21, 2.); 
    tvec[6] = TLine(0.21, 2., 0.18, 2.); //.18

    //bin 3
    tvec[7] = TLine(0.21, 2, 0.21, 2.658); //.18
    tvec[8] = TLine(0.21, 2.658, 0.24, 2.658); //.18
    tvec[9] = TLine(0.24, 2.658, 0.24, 2.); 
    tvec[10] = TLine(0.24, 2., 0.21, 2.); //.18

    //bin 4
    tvec[11] = TLine(0.24, 2., 0.24, 2.658);//.18 
    tvec[12] = TLine(0.24, 2.658, 0.28, 2.658); 
    tvec[13] = TLine(0.28, 2.658, 0.28, 2.); //.18
    tvec[14] = TLine(0.28, 2., 0.24, 2.); //.18

    //bin 5
    tvec[15] = TLine(0.28, 2.658, 0.45, 2.658); 
    tvec[16] = TLine(0.45, 2.658, 0.387826, 2.); 
    tvec[17] = TLine(0.387826, 2., 0.28, 2); 




    //double x5[] = {0.28,   0.28,    0.45,     0.45,   0.387826, 0.28};
    //double y5[] = {2.00,   2.658,   2.658,    2.658,  2,        2};

    //Split by Valerii into two bins:
    double x5_a[] = {0.28,   0.28,    0.34,   0.34};
    double y5_a[] = {2.00,   2.658,   2.658,  2};

    double x5_b[] = {0.34,   0.34,    0.39,   0.387826};
    double y5_b[] = {2.00,   2.658,   2.658,  2};

    double x5_c[] = {0.39,   0.39,    0.45,   0.387826};
    double y5_c[] = {2.00,   2.658,   2.658,  2};

    //bin 6
    tvec[18] = TLine(0.18, 2.658, 0.24, 3.625);//.18 
    tvec[19] = TLine(0.24, 3.625, 0.24, 2.658); 
    tvec[20] = TLine(0.24, 2.658, 0.18, 2.658); //.18

    //bin 7 
    tvec[21] = TLine(0.24, 2.658, 0.24, 3.625); 
    tvec[22] = TLine(0.24, 3.625, 0.28, 3.625); 
    tvec[23] = TLine(0.28, 3.625, 0.28, 2.658); 
    tvec[24] = TLine(0.28, 2.658, 0.24, 2.658);  

    //bin 8
    tvec[25] = TLine(0.28, 2.658, 0.28, 3.625); 
    tvec[26] = TLine(0.28, 3.625, 0.34, 3.625); 
    //tvec[27] = TLine(0.34, 3.625, 0.34, 2.658); 
    tvec[27] = TLine(0.34, 3.625, 0.34, 2.0);// 5 bin divided 
    tvec[28] = TLine(0.34, 2.658, 0.28, 2.658);  

    //bin 9 
    tvec[29] = TLine(0.34, 2.658, 0.34, 3.625); 
    tvec[30] = TLine(0.34, 3.625, 0.54025, 3.625); 
    tvec[31] = TLine(0.54025, 3.625, 0.45, 2.658);    
    tvec[32] = TLine(0.45, 2.658, 0.34, 2.658);  

    //bin 10
    tvec[33] = TLine(0.24, 3.625, 0.34, 5.12); 
    tvec[34] = TLine(0.34, 5.12, 0.34, 3.625); 
    tvec[35] = TLine(0.34, 3.625, 0.24, 3.625); 

    //bin 11 
    tvec[36] = TLine(0.34, 3.625, 0.34, 5.12); 
    tvec[37] = TLine(0.34, 5.12, 0.40, 5.12); 
    //tvec[38] = TLine(0.40, 5.12, 0.40, 3.625);
    // 5 and 9 divided
    tvec[38] = TLine(0.39, 5.12, 0.39, 2.0);
    tvec[39] = TLine(0.40, 3.625, 0.34, 3.625);

    //bin 12 
    tvec[40] = TLine(0.40, 3.625, 0.40, 5.12); 
    tvec[41] = TLine(0.40, 5.12, 0.6234, 5.12); 
    tvec[42] = TLine(0.6234, 5.12, 0.57, 4.05);    
    tvec[43] = TLine(0.57, 4.05, 0.54025, 3.625);
    tvec[44] = TLine(0.54025, 3.625, 0.40, 3.625);

    //bin 13
    tvec[45] = TLine(0.34, 5.12, 0.677, 10.185); 
    tvec[46] = TLine(0.677, 10.185, 0.7896, 11.351);
    tvec[47] = TLine(0.7896, 11.351, 0.75, 9.52); 
    tvec[48] = TLine(0.75, 9.52, 0.708, 7.42); 
    tvec[49] = TLine(0.708, 7.42, 0.64, 5.4); 
    tvec[50] = TLine(0.64, 5.4, 0.6234, 5.12);
    tvec[51] = TLine(0.6234, 5.12, 0.34, 5.12);

    for (int i = 0; i < tvec.size(); i++){
      tvec[i].SetLineWidth(2);
    }

    return tvec;
 }