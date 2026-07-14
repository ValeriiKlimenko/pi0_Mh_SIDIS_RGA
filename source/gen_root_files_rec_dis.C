#include "gen_root_files.C"

void gen_root_files_rec_dis(string inputPath, string hipoFile, string outputPath){
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


  //Gathering the file names from nSidis_files_mc_rec.txt
  HipoChain chain;
  chain.Add((inputPath + hipoFile).c_str());
  TString outputFile = Form("%s%s.root", outputPath.c_str(), hipoFile.c_str());
  outputFile.ReplaceAll(".hipo", "");

  //Returns pid
  auto getByFDID=[](int id, std::vector<region_part_ptr> parts)
       {return clas12::container_filter(parts, [id](region_part_ptr dr)
       {return dr->par()->getPid()==id;});
  };

  //create particles before looping to be more efficient
  auto db = TDatabasePDG::Instance();
  TLorentzVector p4_gamma1;
  TLorentzVector p4_gamma2;
  TLorentzVector p4_ele(0, 0, 0, m_e);
  TLorentzVector pro(0, 0, 0, m_p);
  vector<ushort> cal = {PCAL, ECIN, ECOUT};
  vector<double> res;

  //Loop over files
  for(int ifile=0;ifile<chain.GetNFiles();++ifile){
    clas12reader c12{chain.GetFileName(ifile).Data()};
    //c12.addAtLeastPid(11,0);    //exactly 1 ele

    //Root file, tree, and variable initialization 
    auto tf = new TFile(outputFile.Data(), "RECREATE");
    auto tt = new TTree("h22", "h22");
    
    vector <string> bList = {"ex","ey","ez","esec","ev","ew","ePID", "goodRecElec",
                             "ex_gen","ey_gen","ez_gen",

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

    //Loop over events in file
    int countEvt=0;
    while(c12.next()==true){ 
      if(c12.getDetParticles().empty()){continue;}
      double goodRecEle = 1;
      auto ele=c12.getByID(11);
      auto qVec = (TLorentzVector){0, 0, 0, 0};

      int e_s = 0;
      double ev = 0, ew = 0;

      //Rec electron selection
      if (ele.size() > 0){ 
        if (!goodEventElectron(ele[0])){goodRecEle = 0;}
        p4_ele.SetXYZM(ele[0]->par()->getPx(), ele[0]->par()->getPy(), ele[0]->par()->getPz(), m_e);
      
        qVec = (TLorentzVector){0, 0, sqrt(b_E*b_E - p4_ele.M2()), b_E} - p4_ele;

        if (qVec.E() / b_E >= 0.8){goodRecEle = 0;}//y cut
        if ((pro + qVec).M() <= 2){goodRecEle = 0;}//W cut
        e_s = ele[0]->cal(PCAL)->getSector();

        //Iterating over calos
        for (int k = 0; k < 3; k++){
          if (ev == 0){ev = ele[0]->cal(cal[k])->getLv();}
          if (ew == 0){ew = ele[0]->cal(cal[k])->getLw();}
        }
      }
      else{
        goodRecEle = 0;
      }

      //Generated electron and pi0 selection
      //Accessing mc banks
      auto mcBank = c12.mcparts();
      
      int nRows = mcBank->getRows();
      vector<TLorentzVector> evec;

      //Looping over particles in the event
      //Only intrested in electrons(pid=11) and neutral pions(pid==111)
      for (int i = 0; i < nRows; i++){
        mcBank->setEntry(i);
        
        //Parameters
        int pid         = mcBank->getPid();
        if (pid != 11){continue;}
        int dauIndex    = mcBank->getDaughter();
        TLorentzVector vec = {mcBank->getPx(), mcBank->getPy(), mcBank->getPz(), sqrt(pow(mcBank->getP(),2) + pow(mcBank->getMass(),2))};

        if (pid == 11 && dauIndex == 0){evec.push_back(vec);} 
      }//i
      //End generated banks

      //If there is at least 2 electrons(beam electron and scattered electron) and 1 pi0 we will fill the tree
      // it is from mc bank so the first line is the inc. elec.
      if (evec.size() > 1){
        countEvt++; 
        
          vector<double> res;
          // No good electron, so I put for momentum 0, 0, 0
          if (goodRecEle == 0){
            res = {0,0,0, (double)e_s, ev, ew, 0, goodRecEle,
                  // Generated momentum:
                  evec[1].Px(), evec[1].Py(), evec[1].Pz(),
                  
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 0,
                  0, 0, 
                  0, 0, 0,
              };
          }
          // there is a good rec electron:
          if (goodRecEle == 1){
            res = {p4_ele.Px(),p4_ele.Py(), p4_ele.Pz(), (double)e_s, ev, ew, (double)ele[0]->getPid(), goodRecEle,
                  // Generated elec P:
                  evec[1].Px(), evec[1].Py(), evec[1].Pz(), 
                
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
          }
          add = res;
          tt->Fill();
        
      }
      countEvt++;
   }//Good electron

   tt->Write();
   //tf->Print();
   //tf->Close();
  }//File loop
}//Function 