// Example.C

#include "KF.hh"
//#include "KF.cc"

void Draw(map<int,TVectorD> &measurements,
	  KalmanFilter &myKF){
  
  // plot the state vectors 
  vector<TVectorD> states = myKF.GetStateVectors();
  int n = (int)states.size();
  double *x_state = new double[n];
  double *y_state = new double[n];
  double *z_state = new double[n];
  for(int i = 0; i < n; ++i){
    x_state[i] = states[i](0);
    y_state[i] = states[i](1);
    z_state[i] = states[i](2);
  }
  
  TGraph2D *g_state = new TGraph2D(n, x_state, y_state, z_state);
  delete [] x_state; delete [] y_state; delete [] z_state;
  g_state->SetMarkerStyle(8);
  g_state->SetMarkerColor(1);
  

  // plot the outlying vectors 
  vector<int> outliers = myKF.GetOutliers();  
  int m = (int)outliers.size();
  double *x_out = new double[m];
  double *y_out = new double[m];
  double *z_out = new double[m];
  for(int i = 0; i < m; ++i){
    x_out[i] = measurements[outliers[i]](0);
    y_out[i] = measurements[outliers[i]](1);
    z_out[i] = measurements[outliers[i]](2);
  }

  TGraph2D *g_out = new TGraph2D(m, x_out, y_out, z_out);
  delete [] x_out; delete [] y_out; delete [] z_out;
  g_out->SetMarkerStyle(4);
  g_out->SetMarkerColor(2);
  
  
  n = (int)measurements.size();
  double *x_m = new double[n];
  double *y_m = new double[n];
  double *z_m = new double[n];

  int counter = 0;
  map<int, TVectorD>::iterator it; // draw measurements
  for(it = measurements.begin(); it != measurements.end(); it++){
    
    TVectorD P = it->second;
    x_m[counter] = P(0);
    y_m[counter] = P(1);
    z_m[counter] = P(2);   
    counter++;
  }
  
   
  TGraph2D *g_m = new TGraph2D(n, x_m, y_m, z_m);
  delete [] x_m; delete [] y_m; delete [] z_m;
  g_m->SetMarkerStyle(8);
  g_m->SetMarkerColor(4);

 

  TCanvas* C = new TCanvas;  
  C->SetTitle("Example Candidate");
  g_m->SetTitle("Kalman Filter Example");
  g_m->Draw("P");
  g_state->Draw("SAME P");
  if(m != 0) g_out->Draw("SAME P");
}

void PlotResiduals(KalmanFilter &myKF){

  vector<TVectorD> residuals = myKF.GetResiduals();

  int n = (int)residuals.size();
  double *array = new double[n];
  double *x_array = new double[n];
  for(int i = 0; i < n; ++i){
    TVectorD res = residuals[i];
    double sum_of_res = pow(res(0), 2) + pow(res(1), 2) + pow(res(2), 2);
    array[i] = sqrt(sum_of_res);
    x_array[i] = i;
  }

  TGraph *g = new TGraph(n, x_array, array);
  delete [] array; delete [] x_array;
  g->SetMarkerStyle(8);
  g->SetMarkerColor(4);
  g->SetLineColor(4);

  TCanvas* C = new TCanvas;  
  C->SetTitle("Example Candidate");
  g->SetTitle("Residual Convergance");
  g->Draw("APL");
  
  g->GetXaxis()->SetTitle("Number of steps");
  g->GetYaxis()->SetTitle("Sum of sqaure residuals");
  C->Modified();
  C->Update();
  
}

void Example(double contamination_level = 10.0, double sigma = 0.5){

  /* 
     Generate a 3D linear trajectory containing 100 points
     randomly distributed around a line. It will also contain 
     some fraction of noise.
  */

  double r0_array[4] = {0,0,0,0};
  double v0_array[4] = {1,1,1,1};
  
  TVectorD r0; r0.Use(4, r0_array);
  TVectorD v; v.Use(4, v0_array);
  map<int, TVectorD> measurements;

  
  TRandom *rand = new TRandom();
  int number_of_contaminants = 0;
  for(int i = 0; i < 100 + contamination_level; i++){

    // create contaminants
    if(number_of_contaminants <= contamination_level){

      double prob = rand->Uniform(1);
      if(prob*100.0 < contamination_level/2.0){
	double p_array[4] = {rand->Uniform(100),
			     rand->Uniform(100),
			     rand->Uniform(100),
			     (double)i};
	
	TVectorD P; P.Use(4, p_array);
	measurements.insert(std::make_pair(i, P));
	number_of_contaminants++;
	continue;
      }
    }
    
    TVectorD P(4);
    P(0) = rand->Gaus(r0(0) + v(0)*i, sigma);
    P(1) = rand->Gaus(r0(1) + v(1)*i, sigma);
    P(2) = rand->Gaus(r0(2) + v(2)*i, sigma);
    P(3) = i;
    
    measurements.insert(std::make_pair(i, P));
  } delete rand;

  // Setup parameters
  float dt = 0.1;  // least significant time jump 
  TMatrixD F(6,6); // state dynamics (here we are using linear kinematics)
  TMatrixD Q(6,6); // stochastic processes
 
  for(int i = 0; i < 6; ++i){
    for(int j = 0; j < 6; ++j){
      Q(i,j) = 0;
   
      if(i == j) F(i,j) = 1;
      else if (j == i + 3) F(i,j) = dt;
      else F(i,j) = 0;
    }
  }
    
  TMatrixD H(3,6); // project state vector into measurement space
  TMatrixD R(3,3); // measurement covariance
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 6; ++j){
      if(i == j){
	if(i == 0 || i == 1 || i == 2) R(i,j) = pow(sigma,2);
	H(i,j) = 1;
      } else{
	if(j < 3) R(i,j) = 0;
	H(i,j) = 0;
      }
    }
  }

    
  TVectorD x0(6);  // initial state vector guess 
  x0(0) = 0;  // initial positon guess
  x0(1) = 0;
  x0(2) = 0;

  x0(3) = 1;  // initial state vector guess 
  x0(4) = 1;
  x0(5) = 1;


  TMatrixD P(6,6); // inital covariance matrix guess
  for(int i = 0; i < 6; ++i){
    if(i < 3) P(i,i) = pow(sigma,2);
    else P(i,i) = 1;
    
    for(int j = 0; j < 3; ++j){
      if(i != j){P(i,j) = 0; P(j,i) = 0;}
    }
  }

  double t0 = 0;
  KalmanFilter myKF(F,H,Q,R, dt,6);

  myKF.init(t0, x0,P);
  myKF.SetConfidance(0.01);
  myKF.RunFilter(measurements);
  Draw(measurements, myKF);
  


  // perform a filtering with a poor guess
  TRandom *rand2 = new TRandom();
  x0(0) = rand2->Uniform(100);  // initial positon guess
  x0(1) = rand2->Uniform(100);
  x0(2) = rand2->Uniform(100);

  x0(3) = rand2->Uniform(10);  // initial state vector guess 
  x0(4) = rand2->Uniform(10);
  x0(5) = rand2->Uniform(10);
  delete rand2;

  cout << "Randomly generated starting condition" << endl;
  x0.Print();
  
  KalmanFilter myMessyKF(F,H,Q,R, dt,6);

  myMessyKF.init(t0, x0,P);
  myMessyKF.SetConfidance(-1);
  myMessyKF.RunFilter(measurements);
  Draw(measurements, myMessyKF);
  PlotResiduals(myMessyKF);
}
