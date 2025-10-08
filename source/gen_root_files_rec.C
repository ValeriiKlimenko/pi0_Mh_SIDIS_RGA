#include "gen_root_files.C"
#include <cmath>


//Tuple that contains particle information for use in calculating angular
//distance squared between rec and matched gen events
typedef tuple<TLorentzVector, int, int, double, double, double, int, int, int> lTuple;
//TLorentzVector, pid,parentpid, delta px, delta py, delta pz, index, parent index, match

//Function that compares two lTuples
//Note that phi has twice the possible range as 
bool comparelTuple(lTuple lt1, lTuple lt2){
  double res1 = pow(get<3>(lt1), 2) + pow(get<4>(lt1), 2) + 0*pow(get<5>(lt1), 2);
  double res2 = pow(get<3>(lt2), 2) + pow(get<4>(lt2), 2) + 0*pow(get<5>(lt2), 2);
  return res1 < res2;
}

//Function returns the angular distance squared between the rec and gen events
double getChisq(lTuple lt1){
  return pow(get<3>(lt1), 2) + pow(get<4>(lt1), 2) + 0*pow(get<5>(lt1), 2);
}


bool q2Cut(double e_px, double e_py, double e_pz, double Ebeam){
    
    double lowCut = 1.5;
    double highCut = 11;
    double e_E = std::sqrt(e_px*e_px + e_py*e_py + e_pz*e_pz);
    
    TLorentzVector ele(e_px,e_py,e_pz,e_E);
	TLorentzVector beam(0,0,Ebeam,Ebeam);
	TLorentzVector target(0,0,0,0.93827);
	TLorentzVector fGamma = beam - ele;
	const double Q2 =  -fGamma.M2();
	return (Q2 > lowCut) && (Q2 < highCut);
}

//Function to fix issues with two or more rec events being matched with the 
//same gen event. The rec event with the lowest angular distance, 
//i.e. r^2 = (phi_r - phi_g)^2 + (theta_r -theta-g)^2, gets the match and 
//if there are more rec events than gen, the events with the higher
//angular distance squared are matched and labeled with a 0.
void fixDup(vector<vector<lTuple>> &inMat){
     int dupCount = 2, whileInt = 0; int hjk = 0;
     while (dupCount!=0 && whileInt < inMat.size()){
       for (int i = 1; i < inMat.size(); i++){
         for (int j = 0; j < i; j++){
           if (get<0>(inMat[i][0]) == get<0>(inMat[j][0]))// Compares TLorentz Vectors
           {
             if (comparelTuple(inMat[i][0], inMat[j][0]))//Compares the r^2 value
             {
               inMat[j].push_back(inMat[j][0]);
               inMat[j].erase(inMat[j].begin());
             }
             else{
               inMat[i].push_back(inMat[i][0]);
               inMat[i].erase(inMat[i].begin());
             }
           }//if
        }//j
      }//i

      whileInt++;
     }
      for (int i = 1; i < inMat.size(); i++){
        for (int j = 0; j < i; j++){
          if (get<0>(inMat[i][0]) == get<0>(inMat[j][0])){
            if (comparelTuple(inMat[i][0],inMat[j][0])){
              get<8>(inMat[j][0]) = 0;
            }
            else{
              get<8>(inMat[i][0]) = 0;
            }
          }//if
        }//j
      }//i
}


void gen_root_files_rec(string inputPath, string hipoFile, string outputPath){

  //Gathering the file names from nSidis_files_mc_rec.txt
  HipoChain chain;
  chain.Add((inputPath + '/' + hipoFile).c_str());
  TString outputFile = Form("%s%s.root", outputPath.c_str(), hipoFile.c_str());
  outputFile.ReplaceAll(".hipo", "");

  //Returns pid
  auto getByFDID=[](int id, std::vector<region_part_ptr> parts)
       {return clas12::container_filter(parts, [id](region_part_ptr dr)
       {return dr->par()->getPid()==id;});
  };

  //create particles before looping to be more efficient
  auto db = TDatabasePDG::Instance();
  TLorentzVector p4_g1;
  TLorentzVector p4_g2;
  TLorentzVector p4_ele(0, 0, 0, m_e);
  double pi = TMath::Pi();

  //Loop over files
  for(int ifile=0;ifile<chain.GetNFiles();++ifile){
    clas12reader c12{chain.GetFileName(ifile).Data()};

    //Particle selection
    c12.addExactPid(11,1);   //exactly 1 electron
    c12.addAtLeastPid(22,2); //at least 2 gamma

    //Root file, tree, and variable initialization 
    cout<<outputFile<<endl;
    auto tf = new TFile(outputFile.Data(), "RECREATE");
    auto tt = new TTree("h22", "h22");
    vector <string> bList = {"ex","ey","ez","esec","ev","ew",
                             "g1x","g1y","g1z","g1sec","g1v","g1w","g1pcal","g1b",
                             "g2x","g2y","g2z","g2sec","g2v","g2w","g2pcal","g2b",
                             "emx", "emy", "emz", "emPID", "eVal", "emIndex", "emPIndex",
                             "g1mx", "g1my", "g1mz", "g1mPID", "g1mPPID", "g1mVal", "g1mIndex", "g1mPIndex", "g1match",
                             "g2mx", "g2my", "g2mz", "g2mPID", "g2mPPID", "g2mVal", "g2mIndex", "g2mPIndex", "g2match",
                             "ePID", "eDC6x", "eDC6y", "eDC6z", "eDC18x", "eDC18y", "eDC18z", "eDC36x", "eDC36y", "eDC36z",
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
                               "e_edge_R1", "e_edge_R2", "e_edge_R3",
                                "isMatched"
      
                             };
    vector <double> add(bList.size());
    for (int i = 0; i < bList.size(); i++){
      tt->Branch(bList[i].c_str(), &add[i], (bList[i] + "/D").c_str());
    }
  
    //Loop over all events in the file
    int countEvt = 0;
    while(c12.next()==true){
      if (c12.getDetParticles().empty()){continue;}

      //Accessing the mc bank 
      auto mcBank = c12.mcparts();
      int nRows = mcBank->getRows();
      vector<lTuple> evec, gvec;

      //Looping over particles in the event within the mc bank
      bool isThereElecGEN = 0, isTherePi0GEN = 0;

      double ex_gen = 0, ey_gen = 0, ez_gen= 0,
          pi0x_gen = 0, pi0y_gen = 0, pi0z_gen= 0;
      
      for (int i = 0; i < nRows; i++){
        mcBank->setEntry(i);
        
        //Parameters
        int pid         = mcBank->getPid();
        int index       = mcBank->getIndex();
        int parentIndex = mcBank->getParent();
        int dauIndex    = mcBank->getDaughter();

        

        
        // can be wrong pi0 so should not be used
        /*
        if (isTherePi0GEN == 0 && pid == 111){
          pi0x_gen = mcBank->getPx();
          pi0y_gen = mcBank->getPy();
          pi0z_gen = mcBank->getPz();
          isTherePi0GEN = 1;
        }
        */
        
        if (dauIndex > 0){continue;}

        mcBank->setEntry(parentIndex - 1);
        int parentPid   = mcBank->getPid();//Yields correct parentPid
        
        mcBank->setEntry(i);
        double px       = mcBank->getPx();
        double py       = mcBank->getPy();
        double pz       = mcBank->getPz();
        double energy   = sqrt(pow(mcBank->getP(),2) + pow(mcBank->getMass(),2));
        double theta    = mcBank->getTheta();
        double phi      = mcBank->getPhi();

        if (eleCan(pid)){
          evec.push_back(make_tuple((TLorentzVector){px, py, pz, energy}, 
                                     pid, 
                                     parentPid,
                                     theta,
                                     phi,
                                     pz,
                                     index,
                                     parentIndex,
                                     1)
                                     );
        }
        if (gamCan(pid)){
          gvec.push_back(make_tuple((TLorentzVector){px, py, pz, energy}, 
                                     pid, 
                                     parentPid,
                                     theta,
                                     phi,
                                     pz,
                                     index,
                                     parentIndex,
                                     1)
                                     );
        } 
      }//i

      auto ele=c12.getByID(11);
      p4_ele.SetXYZM(ele[0]->par()->getPx(), ele[0]->par()->getPy(),
                     ele[0]->par()->getPz(), m_e);
            
      //MC electron theta and phi difference loop
      for (int i = 0; i < evec.size(); i++){
        get<3>(evec[i]) = (p4_ele.Theta() - get<3>(evec[i]));
        double dphi = p4_ele.Phi() - get<4>(evec[i]);
        if (dphi > pi) {dphi -= 2.0*pi;}
        if (dphi < -pi){dphi += 2.0*pi;}
        get<4>(evec[i]) = dphi;//0.2
      } 
      sort(evec.begin(), evec.end(), comparelTuple);
            
      auto gs=c12.getByID(22);

      //Gathers and sorts the gamma vector
      vector<vector<lTuple>> gMat(gs.size());
      for (int i = 0; i < gMat.size(); i++){
        gMat[i] = gvec;

        double theta = gs[i]->getTheta();
        double phi   = gs[i]->getPhi();
        for (int j = 0; j < gMat[i].size(); j++){
          get<3>(gMat[i][j]) = (theta - get<3>(gMat[i][j]));
          double dphi = phi - get<4>(gMat[i][j]);
          if (dphi > pi) {dphi -= 2.0*pi;}
          if (dphi < -pi){dphi += 2.0*pi;}
          get<4>(gMat[i][j]) = dphi;//4/3
        } 
        sort(gMat[i].begin(), gMat[i].end(), comparelTuple);
      }
      fixDup(gMat);

      //Iterating over calos
      double ev = 0, ew = 0;
      vector<ushort> cal = {PCAL, ECIN, ECOUT};
      int e_s = ele[0]->cal(PCAL)->getSector();

      for (int k = 0; k < 3; k++){
        if (ev == 0){ev = ele[0]->cal(cal[k])->getLv();}
        if (ew == 0){ew = ele[0]->cal(cal[k])->getLw();}
      }

      //Rec digammma loop
      TString str = "";
      for (int i = 0; i < gs.size(); i++){
         int ga_s1;
         int g1_s[3];
         double g1v = 0, g1w = 0, g1b = gs[i]->par()->getBeta();
         for (int k = 0; k < 3; k++){
           if (g1v == 0){g1v = gs[i]->cal(cal[k])->getLv();}
           if (g1w == 0){g1w = gs[i]->cal(cal[k])->getLw();}    
           g1_s[k] = gs[i]->cal(cal[k])->getSector();
         }
             
        ga_s1 = findSector(g1_s[0], g1_s[1], g1_s[2]);

        for(int j = i+1; j < gs.size(); j++){

          //Fill ttrees if paricles pass cut
          if (goodEventRec(ele[0], gs[i], gs[j]) && q2Cut(p4_ele.Px(), p4_ele.Py(), p4_ele.Pz(), 10.6) ){
             int g2_s[3];
             double g2v = 0, g2w = 0, g2b = gs[j]->par()->getBeta();
             for (int k = 0; k < 3; k++){
               if (g2v == 0){g2v = gs[j]->cal(cal[k])->getLv();}
               if (g2w == 0){g2w = gs[j]->cal(cal[k])->getLw();}
                g2_s[k] = gs[j]->cal(cal[k])->getSector();
             }
          
             int ga_s2 = findSector(g2_s[0], g2_s[1], g2_s[2]); 

             // g1 and g2 generated matching
             // variables for matching 
             bool matchGammaGen = false;
             double delta_g1_gen = 1000, delta_g2_gen = 1000, delta_e_gen = 1000;
             int i_gen_i = 0, i_gen_j = 0;
             double g1x_gen=0, g1y_gen=0, g1z_gen=0;
             double g2x_gen=0, g2y_gen=0, g2z_gen=0;
            
             for (int i_gen = 0; i_gen < nRows; i_gen++){
              mcBank->setEntry(i_gen);
              
              //Parameters
              int pid         = mcBank->getPid();
              int index       = mcBank->getIndex();              
              //Save Generated pi0 and electron information


        //Save Generated pi0 and electron information
              if (pid == 11){
                TLorentzVector e1gen;
                
                e1gen.SetXYZM(mcBank->getPx(), mcBank->getPy(), mcBank->getPz(), m_e);
                double phi_diff_e = abs(p4_ele.Phi() - e1gen.Phi()) / (2 * 3.1415);
                double theta_diff_e = abs(p4_ele.Theta() - e1gen.Theta()) * 2 / 3.1415926;
                
                if (phi_diff_e + theta_diff_e < delta_e_gen){
                  delta_e_gen = phi_diff_e + theta_diff_e;
                  ex_gen = mcBank->getPx();
                  ey_gen = mcBank->getPy();
                  ez_gen = mcBank->getPz();
                  isThereElecGEN = 1;
                }
              }
                     
              if (pid == 22){
                double g_genX = mcBank->getPx();
                double g_genY = mcBank->getPy();
                double g_genZ = mcBank->getPz();

                TLorentzVector g1rec, g2rec, gGen;
                g1rec.SetXYZM(gs[i]->par()->getPx(), gs[i]->par()->getPy(), gs[i]->par()->getPz(), 0);
                g2rec.SetXYZM(gs[j]->par()->getPx(), gs[j]->par()->getPy(), gs[j]->par()->getPz(), 0);
                gGen.SetXYZM(mcBank->getPx(), mcBank->getPy(), mcBank->getPz(), 0);
  
                // angle distance in phi and theta between generated and rec gammas normalaized to the max diff range
                double phi_diff_g1 = abs(g1rec.Phi() - gGen.Phi()) / (2 * 3.1415926);
                double phi_diff_g2 = abs(g2rec.Phi() - gGen.Phi()) / (2 * 3.1415926);
                double theta_diff_g1 = abs(g1rec.Theta() - gGen.Theta()) * 2 / 3.1415926;
                double theta_diff_g2 = abs(g2rec.Theta() - gGen.Theta()) * 2 / 3.1415926;
  
                if (phi_diff_g1 + theta_diff_g1 < delta_g1_gen){
                  delta_g1_gen = phi_diff_g1 + theta_diff_g1;
                  i_gen_i = i_gen;
                  
                  g1x_gen = mcBank->getPx(); 
                  g1y_gen = mcBank->getPy(); 
                  g1z_gen = mcBank->getPz(); 
                }
                  
                if (phi_diff_g2 + theta_diff_g2 < delta_g2_gen){
                  delta_g2_gen = phi_diff_g2 + theta_diff_g2;
                  i_gen_j = i_gen;

                  g2x_gen = mcBank->getPx(); 
                  g2y_gen = mcBank->getPy(); 
                  g2z_gen = mcBank->getPz(); 
                }
               }
              }
              if ((i_gen_i != i_gen_j) && (delta_g1_gen < 0.15 && delta_g2_gen <0.15)){
                matchGammaGen = true;
              }
            
              // END: g1 and g2 generated matching
               
             add = {
                    // RECONSTRUCTED ELECS
                    p4_ele.Px(), p4_ele.Py(), p4_ele.Pz(), (double)e_s, ev, ew,

                    // RECONSTRUCTED GAMMAS
                    gs[i]->par()->getPx(), gs[i]->par()->getPy(), gs[i]->par()->getPz(), (double)ga_s1, g1v, g1w, (double)g1_s[0], g1b,
                    gs[j]->par()->getPx(), gs[j]->par()->getPy(), gs[j]->par()->getPz(), (double)ga_s2, g2v, g2w, (double)g2_s[0], g2b,

                    // GENERATED PARTICLES
                    get<0>(evec[0]).Px(), get<0>(evec[0]).Py(), get<0>(evec[0]).Pz(), 
                    (double)get<1>(evec[0]), getChisq(evec[0]), (double)get<6>(evec[0]), (double)get<7>(evec[0]),
                    get<0>(gMat[i][0]).Px(), get<0>(gMat[i][0]).Py(), get<0>(gMat[i][0]).Pz(), (double)get<1>(gMat[i][0]), 
                    (double)get<2>(gMat[i][0]), getChisq(gMat[i][0]), (double)get<6>(gMat[i][0]), (double)get<7>(gMat[i][0]), (double)get<8>(gMat[i][0]),                         
                    get<0>(gMat[j][0]).Px(), get<0>(gMat[j][0]).Py(), get<0>(gMat[j][0]).Pz(), (double)get<1>(gMat[j][0]), 
                    (double)get<2>(gMat[j][0]), getChisq(gMat[j][0]), (double)get<6>(gMat[j][0]), (double)get<7>(gMat[j][0]), (double)get<8>(gMat[j][0]),  
                    (double)ele[0]->getPid(), ele[0]->traj(DC,6)->getX(), ele[0]->traj(DC,6)->getY(), ele[0]->traj(DC,6)->getZ(),
                                      ele[0]->traj(DC,18)->getX(),ele[0]->traj(DC,18)->getY(),ele[0]->traj(DC,18)->getZ(),
                                      ele[0]->traj(DC,36)->getX(),ele[0]->traj(DC,36)->getY(),ele[0]->traj(DC,36)->getZ(),
                 
                    //Added by Valerii:
                    // for calor. Fid Cuts
                 	ele[0]->cal(PCAL)->getLu(),ele[0]->cal(PCAL)->getLv(),ele[0]->cal(PCAL)->getLw(),
                 	ele[0]->cal(ECIN)->getLu(),ele[0]->cal(ECIN)->getLv(),ele[0]->cal(ECIN)->getLw(),
                 	ele[0]->cal(ECOUT)->getLu(),ele[0]->cal(ECOUT)->getLv(),ele[0]->cal(ECOUT)->getLw(),
                 
                    gs[i]->cal(PCAL)->getLu(),gs[i]->cal(PCAL)->getLv(),gs[i]->cal(PCAL)->getLw(),
                 	gs[i]->cal(ECIN)->getLu(),gs[i]->cal(ECIN)->getLv(),gs[i]->cal(ECIN)->getLw(),
                 	gs[i]->cal(ECOUT)->getLu(),gs[i]->cal(ECOUT)->getLv(),gs[i]->cal(ECOUT)->getLw(),

                    gs[j]->cal(PCAL)->getLu(),gs[j]->cal(PCAL)->getLv(),gs[j]->cal(PCAL)->getLw(),
                 	gs[j]->cal(ECIN)->getLu(),gs[j]->cal(ECIN)->getLv(),gs[j]->cal(ECIN)->getLw(),
                 	gs[j]->cal(ECOUT)->getLu(),gs[j]->cal(ECOUT)->getLv(),gs[j]->cal(ECOUT)->getLw(),

                    // deposited Ener cal for SF cuts:
                    ele[0]->cal(PCAL)->getEnergy(),ele[0]->cal(ECIN)->getEnergy(),ele[0]->cal(ECOUT)->getEnergy(),
                    gs[i]->cal(PCAL)->getEnergy(),gs[i]->cal(ECIN)->getEnergy(),gs[i]->cal(ECOUT)->getEnergy(),
                    gs[j]->cal(PCAL)->getEnergy(),gs[j]->cal(ECIN)->getEnergy(),gs[j]->cal(ECOUT)->getEnergy(),

                    //Vz:
                    ele[0]->par()->getVz(),
                    (double)ele[0]->getSector(),
                    ele[0]->traj(DC,6)->getEdge(),ele[0]->traj(DC,18)->getEdge(),ele[0]->traj(DC,36)->getEdge(),
                   (double)matchGammaGen
                 };


            tt->Fill();
          }//Good event  
        }//j
      }//i
      countEvt++;
    }//events

    tt->Write();
    //tf->Print();
    //tf->Close();

  }//Files
}//gen_root_files_mcrec_match
