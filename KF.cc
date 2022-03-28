#include "KF.hh"

#include<Math/ProbFunc.h>
#include <TDecompChol.h>


//=====================================================
KalmanFilter::KalmanFilter(const TMatrixD &Fin, const TMatrixD &Hin,
			   const TMatrixD &Qin, const TMatrixD &Rin,
			   double Delta_t){

  n = Fin.GetNrows();
  m = Hin.GetNrows();

  (this->F).Use(n,n, Fin.GetMatrixArray());
  (this->H).Use(m,n, Hin.GetMatrixArray());
  (this->Q).Use(n,n, Qin.GetMatrixArray());
  (this->R).Use(m,m, Rin.GetMatrixArray());

  dt = Delta_t;
  initialize = false;
}


KalmanFilter::KalmanFilter(const TMatrixD &Fin, const TMatrixD &Hin,
			   const TMatrixD &Qin, const TMatrixD &Rin,
			   double Delta_t, int max){

  n = Fin.GetNrows();
  m = Hin.GetNrows();

  (this->F).Use(n,n, Fin.GetMatrixArray());
  (this->H).Use(m,n, Hin.GetMatrixArray());
  (this->Q).Use(n,n, Qin.GetMatrixArray());
  (this->R).Use(m,m, Rin.GetMatrixArray());
  
  dt = Delta_t;
  initialize = false;

  states.reserve(max + 1);
  covariance.reserve(max + 1);
  residuals.reserve(max);
}

void KalmanFilter::init(double t0, const TVectorD &x_in, const TMatrixD &P_in){

  (this->P0).Use(n,n, P_in.GetMatrixArray());
  (this->x0).Use(n, x_in.GetMatrixArray());
  t = t0;  t_last = t0;
  Chi2 = 0;
   
  states.push_back(x_in);
  covariance.push_back(P_in);
  initialize = true;
}



void KalmanFilter::Reset(){

  Chi2 = 0;
  states.clear();
  covariance.clear();
  residuals.clear();
  outliers.clear();

  initialize = false;
}



pair<TVectorD, TMatrixD> KalmanFilter::Predict(const TVectorD &xin, const TMatrixD &Pin){

  TMatrixD FT; // set the transpose using temp variable
  TMatrixD F_Temp(n,n, (this->F).GetMatrixArray());
  FT.Use(n,n, (F_Temp.T()).GetMatrixArray());
  
  TVectorD x_new = (this->F) * (xin);
  TMatrixD P_new = (this->F) * (Pin) * (FT) + (this->Q);
  t += dt;
 
  return std::make_pair(x_new, P_new); 
}


bool KalmanFilter::Check_Residual(const TVectorD &y, const TMatrixD &C){

  TMatrixD c_copy(m,m, C.GetMatrixArray());  
  
  TDecompChol cheol(c_copy);
  Bool_t ok;  cheol.Decompose();
  TVectorD y2 = cheol.Solve(y, ok);   
  Double_t Chi2_n = Dot(y,y2);



  // double sided chi2 test
  double upper_area = ROOT::Math::chisquared_cdf_c(Chi2_n, m-1);
  double lower_area = ROOT::Math::chisquared_cdf(Chi2_n, m-1);
  if((lower_area < alpha/2.0) || (upper_area < alpha/2.0)) return false;
  Chi2 += Chi2_n; return true; 
}

void KalmanFilter::SetConfidance(double new_alpha){
  if(new_alpha <= 0.0) check_outlier = false;
  else alpha = new_alpha;
}

tuple<bool, TVectorD, TMatrixD> KalmanFilter::Update(const TVectorD &z, const TVectorD &x, const TMatrixD &P){

  TMatrixD H_Temp(m,n, (this->H).GetMatrixArray());
  TMatrixD HT(n,m, (H_Temp.T()).GetMatrixArray());

  const TVectorD y_pre_fit =  z - (this->H)*x;  // measurment pre-fit residual
  TMatrixD S = (this->H)*P*HT + (this->R);      // pre-fit covariance

  
  // set identity matrix
  TMatrixD I(n,n); 
  for(int i = 0; i < n; ++i){
    for(int j = 0; j <= i; ++j){
      if(i == j) I(i,j) = 1;
      else I(i,j) = 0;
      I(i,j) = I(j,i); 
    }
  }

  const TMatrixD K = P*HT*S.Invert();           // Kalman gain
  const TVectorD x_new = x + K * y_pre_fit;     // a posteriori state 
  TMatrixD P_new = (I - K*this->H) * P;         // a posteriori covariance
  TVectorD y_post_fit = z - this->H * x_new;    // post-fit residual
  TMatrixD C = this->R - (this->H * P_new  * HT);

  
 
  if(check_outlier){
    if(!Check_Residual(y_post_fit, C)){// check if z it an outlier
      return std::make_tuple(false, states.back(), covariance.back());   
    }
  }
    
  residuals.push_back(y_post_fit);
  states.push_back(x_new);
  covariance.push_back(P_new);
  t_last = t;
  
  return std::make_tuple(true, x_new, P_new); 
} 


void KalmanFilter::RunFilter(map<int, TVectorD> &measurments){

  if(!initialize){
    cout << "Error! Request initialization! ..." << endl;
    return;
  }

  map<int, TVectorD>::iterator it; 
  pair<TVectorD, TMatrixD> XP = std::make_pair(this->x0, this->P0);
  for (it = measurments.begin(); it != measurments.end(); it++){
      
    // check for missing observations
    if(it->second(3) < t) goto state_update;
    while(fabs(it->second(3) - t) > dt){
      XP = Predict(XP.first, XP.second);
    }

    
  state_update:
    TVectorD z(3); // measured
    z(0) = it->second(0);  z(1) = it->second(1);  z(2) = it->second(2);
    
    tuple<bool, TVectorD, TMatrixD> Updated = Update(z, XP.first, XP.second);
    XP = std::make_pair(get<1>(Updated), get<2>(Updated));
    if(get<0>(Updated) == false){
      outliers.push_back(it->first);     
      t = t_last;
    }
  }

}


double KalmanFilter::GetChi2(){return Chi2;}
vector<TVectorD> KalmanFilter::GetStateVectors(){return states;}
vector<TVectorD> KalmanFilter::GetResiduals(){return residuals;}
vector<int> KalmanFilter::GetOutliers(){return outliers;}
