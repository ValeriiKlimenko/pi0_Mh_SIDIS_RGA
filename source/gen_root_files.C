#include "DC_fiducial_cuts_mars.cxx"
#include "QADB.h"

// clas12root headers
#include "hipo4/reader.h"
#include "clas12reader.h"

//#include <iguana/algorithms/.h>


using namespace clas12;


//Root file and plot path 
string filePath = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/";
string plotPath = "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/plots/";

//Mass and Energy parameters
const double m_p  = 0.93827208816;          //proton mass
const double b_E  = 10.6041;                //beam energy fall 2018
const double m_e  = 0.5109989461 * 0.001;   //electron mass
const double m_pi = 0.134977;               //pi0 mass

//Check that all the particles in FD
bool goodEvent(clas12::region_particle *e,
               clas12::region_particle *g1,
               clas12::region_particle *g2){
     bool res = 1;
     if (e->getStatus() > -2000 || e->getStatus() < -4000){res = 0;}
     if (g1->getStatus() < 2000 || g1->getStatus() > 4000){res = 0;}
     if (g2->getStatus() < 2000 || g2->getStatus() > 4000){res = 0;}
     if (e->getRegion() != FD || g1->getRegion() != FD || g2->getRegion() != FD){res = 0;}
     return res;
}

//Electron selection FD + Vz cut
bool goodEventElectron(clas12::region_particle *e){
     if (e->getStatus() > -2000 || e->getStatus() < -4000){return 0;}
     if (e->getRegion() != FD){return 0;}
     if (e->par()->getVz() <=-13 || e->par()->getVz() >= 12){return 0;}//Changed to RGA -13 < vz < 12 for inbending
     return 1;
}

// Gamma in FD
bool goodEventGamma(clas12::region_particle *g){
     if (g->getStatus() < 2000 || g->getStatus() > 4000){return 0;}
     if (g->getRegion() != FD){return 0;}
     return 1;
}

//Get Sector, format unknown
int  findSector(int &a, int &b, int &c){
     int res = 0;
     TString str = Form("%d%d%d", a, b, c);
     TString ref = ":111:222:333:444:555:666:";
     str.ReplaceAll("0", "");
     if (str.Length() > 0){res = (int)str[0] - 48;}
     return res;
}

//Electron and photon selection for rec events
bool goodEventRec(clas12::region_particle *e,
                  clas12::region_particle *g1,
                  clas12::region_particle *g2){
     
     if (e->par()->getStatus() > -2000 || e->par()->getStatus() < -4000){return 0;}
     if (g1->par()->getStatus() < 2000 || g1->par()->getStatus() > 4000){return 0;}
     if (g2->par()->getStatus() < 2000 || g2->par()->getStatus() > 4000){return 0;}
     if (e->getRegion() != FD || g1->getRegion() != FD || g2->getRegion() != FD){return 0;}
     if (e->par()->getVz() <=-13 || e->par()->getVz() >= 12){return 0;}//Changed to RGA -13 < vz < 12 for inbending

     TLorentzVector evec, provec, bvec;
     bvec.SetXYZM(0, 0, sqrt(b_E*b_E - m_e*m_e), m_e);
     evec.SetXYZM(e->par()->getPx(), e->par()->getPy(), e->par()->getPz(), m_e);
     if (evec.E() <= 2 || evec.E() >= 8){return 0;}
     provec.SetXYZM(0, 0, 0, m_p);

     if (-(bvec - evec).M2() <=1){return 0;}//Q^2 > 1 GeV
     if ((bvec - evec + provec).M() <= 2){return 0;}// W > 2GeV

     return 1; 
}

//Possible electron like canidates
bool eleCan(int pid){
    if (pid == 11){return 1;}//electron
    if (pid == 13){return 1;}//muon
    if (pid == 15){return 1;}//tau
    if (pid == -211){return 1;}//neg pion
    if (pid == -221){return 1;}//neg eta
    if (pid == -321){return 1;}//K-
    if (pid == -2212){return 1;}//neg proton
    return 0;
} 

//Possible gamma like canidates
bool gamCan(int pid){
    if (pid == 22){return 1;}//gamma
    if (pid == 111){return 1;}//pi0
    if (pid == 130){return 1;}//K0l
    if (pid == 310){return 1;}//K0s
    if (pid == 311){return 1;}//K0
    if (pid == 2112){return 1;}//neutron
    if (pid == -2112){return 1;}//antineutron
    return 0;
}


