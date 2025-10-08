#include "gen_root_files.C"

void gen_root_files_dis(string inputPath, string hipoFile, string outputPath){

  HipoChain chain;
  chain.Add((inputPath + hipoFile).c_str());
  TString outputFile = Form("%s%s.root", outputPath.c_str(), hipoFile.c_str());
  outputFile.ReplaceAll("nSidis_","");
  outputFile.ReplaceAll(".hipo", "");

  TLorentzVector p4_ele(0, 0, 0, m_e);

  auto c12=chain.GetC12Reader(); //in case you want to configure, use this
  c12->applyQA("latest");//GETPASSSTRINGHERE="latest", "pass1, "pass2",...
  c12->db()->qadb_addQARequirement("MarginalOutlier");
  c12->db()->qadb_addQARequirement("TotalOutlier");
  c12->db()->qadb_addQARequirement("TerminalOutlier");
  c12->db()->qadb_addQARequirement("SectorLoss");
  c12->db()->qadb_addQARequirement("LowLiveTime");
  c12->addExactPid(11,1); //exactly 1 ele
  

    //Root file, tree, and variable initialization 
    auto tf = new TFile(outputFile.Data(), "RECREATE");
    auto tt = new TTree("h22", "h22");
    vector <string> bList = {"ex","ey","ez","esec","ev","ew","ePID", 
                             "eDC6x","eDC6y", "eDC6z","eDC18x", "eDC18y", "eDC18z", 
                             "eDC36x", "eDC36y", "eDC36z",
  
                              //Valerii:
                              //electron Cal
                              "e_pcal_Lu","e_pcal_Lv","e_pcal_Lw",
                              "e_ecin_Lu","e_ecin_Lv","e_ecin_Lw",
                              "e_ecout_Lu","e_ecout_Lv","e_ecout_Lw",
 
                              //deposited energy:
                              "e_pcalE","e_ecinE","e_ecoutE",
                              //electrons Vz:
                               "e_vz", "e_sec_DC",
                               "e_edge_R1", "e_edge_R2", "e_edge_R3"
  };
    
    vector <double> add(bList.size());
    for (int i = 0; i < bList.size(); i++){
      tt->Branch(bList[i].c_str(), &add[i], (bList[i] + "/D").c_str());
    }
    
    //Loop over files
    int countEvt=0;
    while(c12->next()==true){ 
      if(c12->getDetParticles().empty()){continue;}
      
      auto parts=c12->getDetParticles();
      auto ele=c12->getByID(11);

      if (!goodEventElectron(ele[0])){continue;}
      p4_ele.SetXYZM(ele[0]->par()->getPx(), ele[0]->par()->getPy(), ele[0]->par()->getPz(), m_e);
      
      auto qVec = (TLorentzVector){0, 0, sqrt(b_E*b_E - p4_ele.M2()), b_E} - p4_ele;
      TLorentzVector pro(0, 0, 0, m_p);
      if (qVec.E() / b_E >= 0.8){continue;}//y cut
      if ((pro + qVec).M() <= 2){continue;}//W cut
      
      int e_s = ele[0]->cal(PCAL)->getSector();

      //Iterating over calos
      vector<ushort> cal = {PCAL, ECIN, ECOUT};
      double ev = 0, ew = 0; 
      for (int k = 0; k < 3; k++){
        if (ev == 0){ev = ele[0]->cal(cal[k])->getLv();}
        if (ew == 0){ew = ele[0]->cal(cal[k])->getLw();}
      }
      add = {p4_ele.Px(), p4_ele.Py(), p4_ele.Pz(), (double)e_s, ev, ew,
            (double)ele[0]->getPid(), ele[0]->traj(DC,6)->getX(), ele[0]->traj(DC,6)->getY(), ele[0]->traj(DC,6)->getZ(),
            ele[0]->traj(DC,18)->getX(),ele[0]->traj(DC,18)->getY(),ele[0]->traj(DC,18)->getZ(),
            ele[0]->traj(DC,36)->getX(),ele[0]->traj(DC,36)->getY(),ele[0]->traj(DC,36)->getZ(),               
          
              //Added by Valerii:
                 // for calor. Fid Cuts
                 ele[0]->cal(PCAL)->getLu(),ele[0]->cal(PCAL)->getLv(),ele[0]->cal(PCAL)->getLw(),
                 ele[0]->cal(ECIN)->getLu(),ele[0]->cal(ECIN)->getLv(),ele[0]->cal(ECIN)->getLw(),
                 ele[0]->cal(ECOUT)->getLu(),ele[0]->cal(ECOUT)->getLv(),ele[0]->cal(ECOUT)->getLw(),

                  // deposited Ener cal for SF cuts:
                  ele[0]->cal(PCAL)->getEnergy(),ele[0]->cal(ECIN)->getEnergy(),ele[0]->cal(ECOUT)->getEnergy(),

                  //Vz:
                  ele[0]->par()->getVz(),
                  (double)ele[0]->getSector(),
                  ele[0]->traj(DC,6)->getEdge(),ele[0]->traj(DC,18)->getEdge(),ele[0]->traj(DC,36)->getEdge()
        };

      tt->Fill();
      countEvt++;
    }//Event
    tt->Write();
    //tf->Print();
    //tf->Close();
  //}//File loop
}//Function 
