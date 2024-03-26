#include<bits/stdc++.h>
#include<random>
#include<omp.h>
#include"auxFunctions.h"
#include"grids.h"

int main(){

  string pwd = filesystem::current_path();
  pwd += "/Grooves/";
  clearDirectory(pwd);
  
  // threads
  int thread_num;
  omp_set_num_threads(24);
  thread_num = omp_get_max_threads ();
  
  //SET DIMENSIONS OF THE SYSTEM
  double scale = 1;
  
  int dimX = 150*scale;
  int dimY = 150*scale;

  int dimZ = 100*scale;

  // 5 microns for the grooves -> 10 points
  // cell radius of 25 -> 50 points
  // nucleus radius of 10 -> 20 points >.<
  
  //MAKE GROOVES
  //geometry parameters
  int depth = 0;//10*scale;
  int spacing = 10*scale;
  int width = 10*scale;

  //global time
  double dt = 0.007; //small 0.007 original 0.05 
  
  //dynamic parameters
  double alpha = 0.001; //original 0.001
  double epsilon = 1; //original 0.001
  //construct
  cout<<"----CREATING GROOVES----"<<endl;  
  Groove mygrooves(dimX, dimY, dimZ, depth, spacing, width, epsilon, alpha);
  cout<<"----GROOVES DONE----"<<endl;  
  double time = 5;
  cout<<"----GROOVE STABILIZATION----"<<endl;
  mygrooves.stabilize(dt, time, 5);
  saveGridToVTI(pwd+"environment_0.vti",mygrooves.grid);
  cout<<"----GROOVE STABILIZATION----"<<endl;

  //MAKE CELL
  //geometry parameters
  depth = width; // for the plane case
  vector<double> radius = {15*scale,15*scale};
  double radius_nucleus = 9*scale;
  vector<double> center_coordinates = {(double)dimX/2,(double)dimY/2,(double)depth*2+(double)radius[1]+1};
  vector<double> nucleus_coordinates = {(double)dimX/2,(double)dimY/2,(double)depth*2+(double)radius[1]+1}; //{(double)dimX/2,(double)dimY/2,(double)depth*2+(double)radius_nucleus+10*scale};
  
  cout<<"----CREATING CELL----"<<endl;  
  Cell mycell (dimX, dimY, dimZ, center_coordinates, radius, nucleus_coordinates, radius_nucleus);
  saveGridToVTI(pwd+"cell_0.vti",mycell.grid);
  saveGridToVTI(pwd+"nucleus_0.vti",mycell.grid_nucleus);
  cout<<"----CELL DONE----"<<endl;  

  //evolve cel
  time = 50;

  mycell.epsilon = 1;         //surf->membrane thickness coefficient
  mycell.alpha = 0.001/12;    //volume coefficient //original 0.001
  mycell.gamma = 1;         //repulsion term //original 6/6
  mycell.eta = 2;          //adhesion coefficient //  2.6/6 will go down
  mycell.kappa = 1;        //surface tension fraction // este é o rácio
  double velocity = 0.;       //drag velocity
  double velocity_nuc = 0;
  
  cout<<"----CELL STABILIZATION----"<<endl;  
  evolve(mycell, mygrooves.grid, 0, 0, dt, 10, pwd);
  cout<<"----CELL STABILIZATION----"<<endl;

  saveGridToVTI(pwd+"cell_0.vti",mycell.grid);
  saveGridToVTI(pwd+"nucleus_0.vti",mycell.grid_nucleus);
  saveGridToVTI(pwd+"environment_0.vti",mygrooves.grid);

  //---------------------------------------------------------------//--------------------------
  double v0_cell = vol(mycell.grid,1e-6);
  double v0_nucleus = vol(mycell.grid_nucleus,1e-6);

  double blows = 0;
  if (v0_cell<1e-3) {
    blows = 1;
  }
  else{
    cout<<"before "<<vol(mycell.grid,1e-6)<<'\n';
    cout<<"----STARTING TIME EVOLUTION----"<<endl;
    evolve(mycell, mygrooves.grid, velocity, velocity_nuc,  dt, time, pwd, (int)((time/dt)/10));
    cout<<"----TIME EVOLUTION DONE----"<<endl;
  }
  
  saveGridToVTI(pwd+"cell_F.vti",mycell.grid);
  saveGridToVTI(pwd+"nucleus_F.vti",mycell.grid_nucleus);
  saveGridToVTI(pwd+"environment_F.vti",mygrooves.grid);
  
  cout<<"\n";

  //-------------------------------------------------------------------//----------------------------------
  cout<<"----POST-PROCESSING----\n";
  vector<double> volumes = DFS(mycell.grid, 2*depth, 0.4);
  vector<double> volumes_nuc = DFS(mycell.grid_nucleus, 2*depth, 0.4);
  double vf_cell = vol(mycell.grid,1e-6);
  double vf_nucleus = vol(mycell.grid_nucleus,1e-6);
  double vol_loss_cell = 1-vf_cell/v0_cell;
  double vol_loss_nucleus = 1-vf_nucleus/v0_nucleus;

  double blows_nucleus = 0;
  if (vf_cell<1e-3) blows = 1;
  if (vf_nucleus<1e-3) blows_nucleus = 1;
  
  cout<<"----POST-PROCESSING DONE----\n";
  cout<<'\n';
  
  cout << "Citoplasm:\n";
  cout << "Lowest Point: " << lowest(mycell.grid, depth) << "\n";
  cout << "Ratio of volume in Grooves: " << vol_in_grooves(mycell.grid, depth*2, 1e-3)/vol(mycell.grid, 1e-3) << "\n";
  cout<<"Volume Cell: "<<v0_cell<<" -> "<<vf_cell<<" corresponds to a loss of " << vol_loss_cell << "% in volume\n";
  cout<< "Groves penetrated: " << volumes.size() <<'\n';
  cout<< "Volumes in grooves: \n";
  for (auto x: volumes) cout<<x<<' ';
  cout<<'\n';
  cout<< "Got splited: " << isSplit(mycell.grid, 0.4) << '\n';
  cout<< "Blows up:" <<blows<<'\n';
  cout<<'\n';

  cout << "Nucleus:\n";
  cout << "Lowest Point: " << lowest(mycell.grid_nucleus, depth) << "\n";
  cout << "Ratio of volume in Grooves: " << vol_in_grooves(mycell.grid_nucleus, depth*2, 1e-3)/vol(mycell.grid_nucleus, 1e-3) << "\n";
  cout<<"Volume Nucleus: "<<v0_nucleus<<" -> "<<vf_nucleus<<" corresponds to a loss of " << 1-vf_nucleus/v0_nucleus << "% in volume\n";
  cout<< "Groves penetrated: " << volumes_nuc.size() <<'\n';
  cout<< "Volumes in grooves: \n";
  for (auto x: volumes_nuc) cout<<x<<' ';
  cout<<'\n';
  cout<< "Got splited: " << isSplit(mycell.grid_nucleus, 0.4) << '\n';
  cout<< "Blows up:" <<blows_nucleus<<'\n';
  cout<<'\n';
}
