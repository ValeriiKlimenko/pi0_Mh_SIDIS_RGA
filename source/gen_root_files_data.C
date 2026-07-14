#include "gen_root_files.C"
#include "Kinem_Funcs.C"
#include <iguana/algorithms/clas12/PhotonGBTFilter/Algorithm.h>


void gen_root_files_data(string inputPath, string hipoFile, string outputPath, bool is_greg_ai = false){
  //Graph styles
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1);
  gStyle->SetLineWidth(2);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.14);

  HipoChain chain;
  chain.Add((inputPath + hipoFile).c_str());
  TString outputFile = Form("%s%s.root", outputPath.c_str(), hipoFile.c_str());
  outputFile.ReplaceAll("nSidis_","");
  outputFile.ReplaceAll(".hipo", "");
 
  //create particles before looping to be more efficient
  auto db = TDatabasePDG::Instance();
  TLorentzVector p4_gamma1;
  TLorentzVector p4_gamma2;
  TLorentzVector p4_ele(0, 0, 0, m_e);

  // create iguana algorithms
  iguana::clas12::PhotonGBTFilter  algo_photon;   // filter the z-vertex (a filter algorithm)
  algo_photon.SetOption("pass", 2);  
  algo_photon.SetOption("o_threshold", 0.78);
  algo_photon.Start();

  auto iguana_action = [&algo_photon](clas12::clas12reader* cr)
  {
    algo_photon.Run(
      cr->getRECParticle(),
      cr->getRECCalorimeter(),
      cr->getRUNconfig()
    );
    return true;
  };

  cout<<"greg AI is:"<<is_greg_ai<<endl;

  //Loop over files
  for(int ifile=0;ifile<chain.GetNFiles();++ifile){

    clas12reader c12{chain.GetFileName(ifile).Data()};
    if (is_greg_ai) c12.SetReadAction(iguana_action);

    //Golden runs
    clas12databases db;
    c12.connectDataBases(&db);

    c12.applyQA("pass2");
    c12.db()->qadb_addQARequirement("MarginalOutlier");
    c12.db()->qadb_addQARequirement("TotalOutlier");
    c12.db()->qadb_addQARequirement("TerminalOutlier");
    c12.db()->qadb_addQARequirement("SectorLoss");

    //Particle Selection 
    c12.addExactPid(11,1);    //exactly 1 ele
    c12.addAtLeastPid(22,2);  //at least 2 gammas

    //Root file, tree, and variable initialization 
    auto tf = new TFile(outputFile.Data(), "RECREATE");
    auto tt = new TTree("h22", "h22");
    vector <string> bList = {"ex","ey","ez","esec","ev","ew",
                             "g1x","g1y","g1z","g1sec","g1v","g1w","g1pcal","g1b",
                             "g2x","g2y","g2z","g2sec","g2v","g2w","g2pcal","g2b",
                             "ePID", "eDC6x","eDC6y", "eDC6z","eDC18x", "eDC18y", "eDC18z", 
                             "eDC36x", "eDC36y", "eDC36z",
                              //Valerii:
                              //electron Cal
                              "e_pcal_Lu","e_pcal_Lv","e_pcal_Lw",
                              "e_ecin_Lu","e_ecin_Lv","e_ecin_Lw",
                              "e_ecout_Lu","e_ecout_Lv","e_ecout_Lw",
                              //photons cal:
                              "g1_pcal_Lu","g1_pcal_Lv","g1_pcal_Lw",
                              "g1_ecin_Lu","g1_ecin_Lv","g1_ecin_Lw",
                              "g1_ecout_Lu","g1_ecout_Lv","g1_ecout_Lw", 
                              "g2_pcal_Lu","g2_pcal_Lv","g2_pcal_Lw",
                              "g2_ecin_Lu","g2_ecin_Lv","g2_ecin_Lw",
                              "g2_ecout_Lu","g2_ecout_Lv","g2_ecout_Lw", 
                              //deposited energy:
                              "e_pcalE","e_ecinE","e_ecoutE",
                              "g1_pcalE","g1_ecinE","g1_ecoutE",
                              "g2_pcalE","g2_ecinE","g2_ecoutE",
                              //electrons Vz:
                               "e_vz", "e_sec_DC", 
                               "e_edge_R1", "e_edge_R2", "e_edge_R3"
                                //iguana:
                                //,"PhotonGBTF"
    
                               };
    vector <double> add(bList.size());
    for (int i = 0; i < bList.size(); i++){
      tt->Branch(bList[i].c_str(), &add[i], (bList[i] + "/D").c_str());
    }


    
    
    //Loop over all events in the file
    int countEvt=0;
    while(c12.next()==true){ 
      if(c12.getDetParticles().empty()){continue;}
      
      auto parts=c12.getDetParticles();
      auto ele=c12.getByID(11);
      auto gammas = is_greg_ai ? c12.getByID(22, true) : c12.getByID(22);

      if (!goodEventElectron(ele[0])){continue;}

      //Digamma loop
      TString str;
      for (int i = 0; i < gammas.size(); i++){
        auto gamma1 = gammas[i];
        if (!goodEventGamma(gamma1)){continue;} 
        for(int j = i+1; j < gammas.size(); j++){
          auto gamma2 = gammas[j];
          if (!goodEventGamma(gamma2)){continue;}
      
          p4_gamma1.SetXYZM(gamma1->par()->getPx(),gamma1->par()->getPy(),gamma1->par()->getPz(),0);
          p4_gamma2.SetXYZM(gamma2->par()->getPx(),gamma2->par()->getPy(),gamma2->par()->getPz(),0);
          p4_ele.SetXYZM(ele[0]->par()->getPx(), ele[0]->par()->getPy(), ele[0]->par()->getPz(), m_e);
          auto pi0 = p4_gamma1 + p4_gamma2;
      
          //Fill ttree if gammas are in FD
          if (goodEvent(ele[0], gamma1, gamma2))
          {
             int e_s, ga_s1, ga_s2;
             int g1_s[3], g2_s[3];
             e_s  = ele[0]->cal(PCAL)->getSector();

             //Iterating over calos
             vector<ushort> cal = {PCAL, ECIN, ECOUT};
             double ev = 0, ew = 0, 
                   g1v = 0, g1w = 0,
                   g2v = 0, g2w = 0;
             for (int k = 0; k < 3; k++){
               if (ev == 0){ev = ele[0]->cal(cal[k])->getLv();}
               if (ew == 0){ew = ele[0]->cal(cal[k])->getLw();}
               if (g1v == 0){g1v = gammas[i]->cal(cal[k])->getLv();}
               if (g1w == 0){g1w = gammas[i]->cal(cal[k])->getLw();}
               if (g2v == 0){g2v = gammas[j]->cal(cal[k])->getLv();}
               if (g2w == 0){g2w = gammas[j]->cal(cal[k])->getLw();}
               
               g1_s[k] = gamma1->cal(cal[k])->getSector();
               g2_s[k] = gamma2->cal(cal[k])->getSector();
             }
             
             ga_s1 = findSector(g1_s[0], g1_s[1], g1_s[2]);
             ga_s2 = findSector(g2_s[0], g2_s[1], g2_s[2]);
              
             //Betas
             double g1b = gammas[i]->par()->getBeta();
             double g2b = gammas[j]->par()->getBeta();

             add = {p4_ele.Px(), p4_ele.Py(), p4_ele.Pz(), (double)e_s, ev, ew,
                    p4_gamma1.Px(), p4_gamma1.Py(), p4_gamma1.Pz(), (double)ga_s1, g1v, g1w, (double)g1_s[0], g1b,
                    p4_gamma2.Px(), p4_gamma2.Py(), p4_gamma2.Pz(), (double)ga_s2, g2v, g2w, (double)g2_s[0], g2b,
                    (double)ele[0]->getPid(), ele[0]->traj(DC,6)->getX(), ele[0]->traj(DC,6)->getY(), ele[0]->traj(DC,6)->getZ(),
                                      ele[0]->traj(DC,18)->getX(),ele[0]->traj(DC,18)->getY(),ele[0]->traj(DC,18)->getZ(),
                                      ele[0]->traj(DC,36)->getX(),ele[0]->traj(DC,36)->getY(),ele[0]->traj(DC,36)->getZ(),
                     //Added by Valerii:
                    // for calor. Fid Cuts
                 	ele[0]->cal(PCAL)->getLu(),ele[0]->cal(PCAL)->getLv(),ele[0]->cal(PCAL)->getLw(),
                 	ele[0]->cal(ECIN)->getLu(),ele[0]->cal(ECIN)->getLv(),ele[0]->cal(ECIN)->getLw(),
                 	ele[0]->cal(ECOUT)->getLu(),ele[0]->cal(ECOUT)->getLv(),ele[0]->cal(ECOUT)->getLw(),
                 
                    gammas[i]->cal(PCAL)->getLu(),gammas[i]->cal(PCAL)->getLv(),gammas[i]->cal(PCAL)->getLw(),
                 	gammas[i]->cal(ECIN)->getLu(),gammas[i]->cal(ECIN)->getLv(),gammas[i]->cal(ECIN)->getLw(),
                 	gammas[i]->cal(ECOUT)->getLu(),gammas[i]->cal(ECOUT)->getLv(),gammas[i]->cal(ECOUT)->getLw(),

                    gammas[j]->cal(PCAL)->getLu(),gammas[j]->cal(PCAL)->getLv(),gammas[j]->cal(PCAL)->getLw(),
                 	gammas[j]->cal(ECIN)->getLu(),gammas[j]->cal(ECIN)->getLv(),gammas[j]->cal(ECIN)->getLw(),
                 	gammas[j]->cal(ECOUT)->getLu(),gammas[j]->cal(ECOUT)->getLv(),gammas[j]->cal(ECOUT)->getLw(),

                    // deposited Ener cal for SF cuts:
                    ele[0]->cal(PCAL)->getEnergy(),ele[0]->cal(ECIN)->getEnergy(),ele[0]->cal(ECOUT)->getEnergy(),
                    gammas[i]->cal(PCAL)->getEnergy(),gammas[i]->cal(ECIN)->getEnergy(),gammas[i]->cal(ECOUT)->getEnergy(),
                    gammas[j]->cal(PCAL)->getEnergy(),gammas[j]->cal(ECIN)->getEnergy(),gammas[j]->cal(ECOUT)->getEnergy(),

                    //Vz:s
                    ele[0]->par()->getVz(),
                    (double) ele[0]->getSector(),
                    ele[0]->traj(DC,6)->getEdge(),ele[0]->traj(DC,18)->getEdge(),ele[0]->traj(DC,36)->getEdge()
                 };

            tt->Fill();
          }//Good event 
       }//Gamma j
     }//Gamma i
     countEvt++;
    }//Event
    tt->Write();
    //tf->Print();
    //tf->Close();
  }//File loop

  algo_photon.Stop();
}//Function 
