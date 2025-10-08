/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// DC fiducial cuts for CLAS12 (September 18 2020)
/// 
/// contact: sdiehl@jlab.org, Aron.Kripko@exp2.physik.uni-giessen.de
///
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Two versions (to be used for different cases):
///
/// DC_fiducial_cut_theta_phi(int j, int region)  --> based on local theta phi coordinate
/// --> Use this cut only for inbending hadrons! 
/// --> The cut parameters are not adjusted for any outbending hadron (Do not use it for outbending!)
///
/// DC_fiducial_cut_XY(int j, int region)  --> based on local x, y coordinate
/// --> USe this cut for inbending electrons and all outbending hadrons 
///
/// inbending electrons: DC_fiducial_cut_XY
/// inbending hadrons:  DC_fiducial_cut_theta_phi
/// all outbending: DC_fiducial_cut_XY
///

/// required avriables and functions:

// select the field setting in a global variable:

bool inbending = true;
bool outbending = false; 


// The variable arrays are defined as follows:

/// from REC::Partcile
/*Commented out by Marshall
part_pid[i] = vpart_pid->at(i);    

/// from REC::Traj

if(vTraj_pindex->at(i) >= 0 && vTraj_pindex->at(i) < BUFFER && vTraj_detID->at(i) == 6 && vTraj_layerID->at(i) == 6){    
  c1x[vTraj_pindex->at(i)] = vTraj_x->at(i);
  c1y[vTraj_pindex->at(i)] = vTraj_y->at(i);
  c1z[vTraj_pindex->at(i)] = vTraj_z->at(i);
}
if(vTraj_pindex->at(i) >= 0 && vTraj_pindex->at(i) < BUFFER && vTraj_detID->at(i) == 6 && vTraj_layerID->at(i) == 18){  
  c2x[vTraj_pindex->at(i)] = vTraj_x->at(i);
  c2y[vTraj_pindex->at(i)] = vTraj_y->at(i);
  c2z[vTraj_pindex->at(i)] = vTraj_z->at(i);
}
if(vTraj_pindex->at(i) >= 0 && vTraj_pindex->at(i) < BUFFER && vTraj_detID->at(i) == 6 && vTraj_layerID->at(i) == 36){  
  c3x[vTraj_pindex->at(i)] = vTraj_x->at(i);
  c3y[vTraj_pindex->at(i)] = vTraj_y->at(i);
  c3z[vTraj_pindex->at(i)] = vTraj_z->at(i);
}

part_DC_sector[i] = determineSector(i);
*/



int determineSector(vector<double> &c2){
  double phi = (180.0 / TMath::Pi() ) * atan2(c2[1] / sqrt(pow(c2[0], 2) + pow(c2[1], 2) + pow(c2[2], 2)), 
                                              c2[0] / sqrt(pow(c2[0], 2) + pow(c2[1], 2) + pow(c2[2], 2)));

  if     (phi < 30 && phi >= -30){   return 1;}
  else if(phi < 90 && phi >= 30){    return 2;}
  else if(phi < 150 && phi >= 90){   return 3;}
  else if(phi >= 150 || phi < -150){ return 4;}
  else if(phi < -90 && phi >= -150){ return 5;}
  else if(phi < -30 && phi >= -90){  return 6;}
  return 0;
}