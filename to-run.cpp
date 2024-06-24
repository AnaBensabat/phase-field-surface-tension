#include<bits/stdc++.h>
#include<random>
#include<omp.h>
#include"auxFunctions.h"
#include"grids.h"

void run_simulation(string pwd, double eta, double kappa){
  
  //SET DIMENSIONS OF THE SYSTEM
  double scale = 1;
  
  int dimX = 150*scale;
  int dimY = 150*scale;
  int dimZ = 70*scale;
  
  // 5 microns for the grooves -> 10 points
  // cell radius of 25 -> 50 points
  // nucleus radius of 10 -> 20 points 
  
  //MAKE GROOVES
  //geometry parameters
  int depth = 10*scale;
  int spacing = 10*scale;
  int width = 10*scale;
  
  //dynamic parameters
  double alpha = 0.001;
  double epsilon = 0.5;
  //construct
  cout<<"----CREATING GROOVES----"<<endl;  
  Groove mygrooves(dimX, dimY, dimZ, depth, spacing, width, epsilon, alpha);
  cout<<"----GROOVES DONE----"<<endl;  
  double dt = 0.07;
  double time = 5;
  cout<<"----GROOVE STABILIZATION----"<<endl;
  mygrooves.stabilize(dt, time, 5);
  cout<<"----GROOVE STABILIZATION----"<<endl;
  
  //MAKE CELL
  //geometry parameters
  vector<double> radius = {28*scale,15*scale};
  vector<double> radius_nucleus = {15*scale, 7*scale};
  depth = width;
  vector<double> center_coordinates = {(double)dimX/2,(double)dimY/2,(double)depth*2};
  vector<double> nucleus_coordinates = {(double)dimX/2-2,(double)dimY/2-2,(double)depth*2+1};
  
  cout<<"----CREATING CELL----"<<endl;
  Cell mycell (dimX, dimY, dimZ, center_coordinates, radius, nucleus_coordinates, radius_nucleus, false, true);
  saveGridToVTI(pwd+"cell_0.vti",mycell.grid);
  saveGridToVTI(pwd+"nucleus_0.vti",mycell.grid_nucleus);
  saveGridToVTI(pwd+"environment_0.vti",mycell.grid_nucleus);
  cout<<"----CELL DONE----"<<endl;  
  
  //cell parameters
  mycell.epsilon = 0.5;         //surf->membrane thickness coefficient
  mycell.alpha = 0.001/12;      //volume coefficient
  mycell.gamma = 1;             //repulsion term
  mycell.eta = eta;             //adhesion coefficient 
  mycell.kappa = kappa;         //surface tension fraction // este é o rácio
  mycell.chi = 100;             //non local term
  double velocity = 0.;         //drag velocity
  double velocity_nuc = 0;
  
  cout<<"----CELL STABILIZATION----"<<endl;
  time = 5;
  evolve(mycell, mygrooves.grid, 0, 0, dt, time, pwd);
  cout<<"----CELL STABILIZATION----"<<endl;
  
  string s1 = "cell_0_eta="+to_string(mycell.eta)+"_k="+to_string(mycell.kappa)+"_eps="+to_string(mycell.epsilon)+".vti";
  string s2 = "nucleus_0_eta="+to_string(mycell.eta)+"_k="+to_string(mycell.kappa)+"_eps="+to_string(mycell.epsilon)+".vti";
  string s3 = "environment_0_eta="+to_string(mycell.eta)+"_k="+to_string(mycell.kappa)+"_eps="+to_string(mycell.epsilon)+".vti";
  saveGridToVTI(pwd+s1,mycell.grid);
  saveGridToVTI(pwd+s2,mycell.grid_nucleus);
  saveGridToVTI(pwd+s3,mygrooves.grid);
  
  double v0_cell = vol(mycell.grid,1e-6);
  double v0_nucleus = vol(mycell.grid_nucleus,1e-6);
  
  double blows = 0;
  if (v0_cell<1e-3) {
    blows = 1;
    cout<<"----BLOWS UP----"<<endl;
  }
  else{
    time = 2000;

    cout<<"----STARTING TIME EVOLUTION----"<<endl;
    evolve(mycell, mygrooves.grid, velocity, velocity_nuc,  dt, time, pwd, 50000);
    
    s1 = "cell_F_eta="+to_string(mycell.eta)+"_k="+to_string(mycell.kappa)+".vti";
    s2 = "nucleus_F_eta="+to_string(mycell.eta)+"_k="+to_string(mycell.kappa)+".vti";
    s3 = "environment_F_eta="+to_string(mycell.eta)+"_k="+to_string(mycell.kappa)+".vti";
    saveGridToVTI(pwd+s1,mycell.grid);
    saveGridToVTI(pwd+s2,mycell.grid_nucleus);
    saveGridToVTI(pwd+s3,mygrooves.grid);
    cout<<"----TIME EVOLUTION DONE----"<<endl;
  }
  
  return ;
}

int main(){
  
  // threads
  int thread_num;
  omp_set_num_threads(32);
  thread_num = omp_get_max_threads ();

  string pwd = filesystem::current_path();
  filesystem::create_directory("output");
  pwd+="/output/";
  
  vector<double> etas = {0.05      , 0.21111111, 0.37222222, 0.53333333, 0.69444444, 0.85555556, 1.01666667, 1.17777778, 1.33888889, 1.5       };//{0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9};
  vector<double> kappas = {0.001, 0.0046, 0.0215, 0.1, 0.464, 2.15, 10, 46.1, 215.4, 1000};//{0.1, 0.5, 1, 1.5, 2, 2.5, 3};

  for(auto eta:etas){
    for(auto kappa:kappas){
      run_simulation(pwd, eta, kappa);
    }
  }
}
