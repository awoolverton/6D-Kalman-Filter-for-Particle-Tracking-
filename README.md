# 6D-Kalman-Filter-for-Particle-Tracking

## Introduction
This is a 6 dimentional descrete Kalman filter implemented for the purposes of tracking charged particles. The state vector keeps track of both the position and direction of the particle using time as the descrete quantity.

## The Kalman Filter
The Kalman filter is used to estimate the statistical behavoir of a system in the presence of noise. This is done by updating a series of measurements taken with prior knowledge of the underlying physical state dynamics This iterative method approaches asymptotically approaches the desired behavior by applying corrations to the *a priori* estimate. In addition, the Kalman filter simultaneously provides a quality score for removing false track candidates or "Ghost Tracks." Furthermore, the quality score inroduces a robustness into the the filtering routine against anomalous data and missing observations. 

If each individual measurement arrives discontinuously a discrete Kalman filter is used to estimate the dynamics. Here, the implementation of a 6-dimentional filter is described in detail. In general, the filter alternates between a **prediction** and **update** step in order to converge at a proper estimate of all state vectors.  

## Installation
This program is dependant on the [ROOT](https://root.cern/) package v-1.6 developed by cern. With this dependancy properly installed, the project can be cloned into a local folder. 

## Usage
The filter requieres 4 TMatrix objects as inputs. The matex Fin describs the state dynamics. Hin translates the *a postiriori* state estimate into measurement space so that the measurements and the states have the same dimentions. Qin describs any random processes within the state dynamics. Lastly, R describes the measurement covariance. Additionaly, the width of the descrete layor Delta_t is required. This depends on the resolution of the measurements, and can negativly affect the computation time. 

Initialization is ran by explicitly calling the init() function. This function requires an initial guess for time, position x0, and covariance P0.  

The convidance interval is used to determine outlying measurements using a double-sided chi-squared test. The function SetConfidance is used to set the alpha theshold. With the default set to 0.05, this subroutine may be ignored if alpha is set below 0.

Lastly, a list of measurements is requred to run the project using the function RunFilter(). This is an ordered set in the form of and std::map<int, TVectorD>. In this implementation the ordering key is determined by the final element of the vector and is used to keep track of measurements that fail the quality test. It is the users responsability to make sure that this ordiering is performed properly.

Since the program is constantly keeping track of the states and the covaraince matricies, a fuction is provided to reset the object. This clears all underlying data structures from memory so that a new initial track guess can be provided.

Additional functionality is provided by the functions GetStateVectors() and GetResiduals(). These functions return a vector of TVectorD objects for further analysis. The function GetOutliers() returns a vector of indicies from the map of measurements. These measurements have failed the quality test and were not used to update the state dynamics.

## Project Example
The project example generates a series of points generated around the line originating at the origin in (1,1,1) direction. The level of contamination can be adjusted as well as the standard deviation of points aroind the line. In this example, the standard deviation was 0.5 and the model follows basic 3D linear kinematics such that the state vector is (x,y,z, vx,vy,vz). The contamination level was set to 10%. With an alph threshold set to 0.01, the filter was successfully able to remove the abnormalities. 

![image info](/images/example_of_good_guess.pdf)

## Contributing
I encurage any pull requests and experimentation; however, please open an issue thread in order to discuss possible modifications.

## Licence and Copywright 
(C) Austin David Woolverton, The Institute of Quantum Computing at the University of Waterloo
Licenced unter the [MIT Licence](LICENCE.md).
