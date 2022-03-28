# 6D-Kalman-Filter-for-Particle-Tracking 
## Introduction 
This is a 6 dimensional discrete Kalman filter implemented for the purposes of tracking charged particles. The state vector keeps track of both the position and direction of the particle using time as the discrete quantity. 

## The Kalman Filter 
The Kalman filter is used to estimate the statistical behavior of a system in the presence of noise. This is done by updating a series of measurements taken with prior knowledge of the underlying physical state dynamics This iterative method approaches asymptotically approaches the desired behavior by applying corrasions to the *a priori* estimate. In addition, the Kalman filter simultaneously provides a quality score for removing false track candidates or "Ghost Tracks." Furthermore, the quality score introduces a robustness into the filtering routine against anomalous data and missing observations.  

If each individual measurement arrives discontinuously a discrete Kalman filter is used to estimate the dynamics. Here, the implementation of a 6-dimentional filter is described in detail. In general, the filter alternates between a **prediction** and **update** step to converge at a proper estimate of all state vectors.   

## Installation 
This program is dependent on the [ROOT](https://root.cern/) package v-1.6 developed by cern. With this dependency installed, the project can be cloned into a local folder.  

## Usage 
The filter requires 4 TMatrix objects as inputs. The martix Fin describes the state dynamics. Hin translates the *a postiriori* state estimate into measurement space so that the measurements and the states have the same dimensions. Qin describes any random processes within the state dynamics. Lastly, R describes the measurement covariance. Additionally, the width of the discrete layer Delta_t is required. This depends on the resolution of the measurements and can negatively affect the computation time.  

Initialization is ran by explicitly calling the init() function. This function requires an initial guess for time, position x0, and covariance P0.   

The confidence interval is used to determine outlying measurements using a double-sided chi-squared test. The function SetConfidance is used to set the alpha threshold. With the default set to 0.05, this subroutine may be ignored if alpha is set below 0. 

Lastly, a list of measurements is required to run the project using the function RunFilter(). This is an ordered set in the form of and std::map<int, TVectorD>. In this implementation the ordering key is determined by the final element of the vector and is used to keep track of measurements that fail the quality test. It is the user's responsibility to make sure that this ordering is performed properly. 

Since the program is constantly keeping track of the states and the covariance matrices, a function is provided to reset the object. This clears all underlying data structures from memory so that a new initial track guess can be provided. 

Additional functionality is provided by the functions GetStateVectors() and GetResiduals(). These functions return a vector of TVectorD objects for further analysis. The function GetOutliers() returns a vector of indices from the map of measurements. These measurements have failed the quality test and were not used to update the state dynamics. 

## Project Example 
The project example generates a series of points generated around the line originating at the origin in (1,1,1) direction. The level of contamination can be adjusted as well as the standard deviation of points around the line. In this example, the standard deviation was 0.5 and the model follows basic 3D linear kinematics such that the state vector is (x,y,z, vx,vy,vz). The contamination level was set to 10%. With an alpha threshold set to 0.01, the filter was successfully able to remove the abnormalities. The measurements are displayed in blue and the state vectors are black. The measurements that have failed the quality test are encircled in red.  

![image 1](/images/example_of_good_guess.png) 

Alternatively, if the quality test is ignored, the state vectors will converge onto the appropriate measurement. This is under the assumption that many measurements are provided. In this case, a starting guess of (15.43, 42.32, 79.05, 2.46, 3.98, 1.22) was randomly generated. Additionally, the sum of squared residuals is also presented for each iteration step.    

![image 2](/images/example_of_bad_guess.png) 

![image 3](/images/convergence.png) 

## Contributing 
I encourage any pull requests and experimentation; however, please open an issue thread in order to discuss possible modifications. 

## License and Copyright  
(C) Austin David Woolverton, The Institute of Quantum Computing at the University of Waterloo 
Licensed under the [MIT Licence](LICENCE.md).
