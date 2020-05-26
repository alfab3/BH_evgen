#include <iostream>
#include "RandomReal.h"
#include <math.h>
#include "dSigmaBorn.h"

double E0, x, phi1, phi2, theta1, theta2, W, T, currVal, max0 = 0, max1 = 0, max2 = 0, max3 = 0, maxReturn;
double ElectronMass = 5.1099894761e-4, min_BeamEnergy = 8.1, max_BeamEnergy = 8.9, max_theta = 13.12 * M_PI/180, min_theta = 0.9 * M_PI/180;


double maxTest(){
    for(int n = 0; n <= 1000; n++){
		x = RandomReal(ElectronMass/min_BeamEnergy, 1 - ElectronMass/max_BeamEnergy);
		phi1 = RandomReal(0,2*M_PI);
		phi2 = RandomReal(0,2*M_PI);
		theta1 = 1/(RandomReal(1/max_theta,1/min_theta));
		theta2 = 1/(RandomReal(1/max_theta,1/min_theta));
		E0 = RandomReal(min_BeamEnergy, max_BeamEnergy);
        currVal = dSigmaBorn(x, theta1, theta2, phi1, phi2, E0);
        if(currVal > max0){
            max0 = currVal;
        }
    }
    for(int n = 0; n <= 10000; n++){
		x = RandomReal(ElectronMass/min_BeamEnergy, 1 - ElectronMass/max_BeamEnergy);
		phi1 = RandomReal(0,2*M_PI);
		phi2 = RandomReal(0,2*M_PI);
		theta1 = 1/(RandomReal(1/max_theta,1/min_theta));
		theta2 = 1/(RandomReal(1/max_theta,1/min_theta));
		E0 = RandomReal(min_BeamEnergy, max_BeamEnergy);
        currVal = dSigmaBorn(x, theta1, theta2, phi1, phi2, E0);
        if(currVal > max1){
            max1 = currVal;
        }
    }
    for(int n = 0; n <= 100000; n++){
		x = RandomReal(ElectronMass/min_BeamEnergy, 1 - ElectronMass/max_BeamEnergy);
		phi1 = RandomReal(0,2*M_PI);
		phi2 = RandomReal(0,2*M_PI);
		theta1 = 1/(RandomReal(1/max_theta,1/min_theta));
		theta2 = 1/(RandomReal(1/max_theta,1/min_theta));
		E0 = RandomReal(min_BeamEnergy, max_BeamEnergy);
        currVal = dSigmaBorn(x, theta1, theta2, phi1, phi2, E0);
        if(currVal > max2){
            max2 = currVal;
        }
    }
    for(int n = 0; n <= 1000000; n++){
		x = RandomReal(ElectronMass/min_BeamEnergy, 1 - ElectronMass/max_BeamEnergy);
		phi1 = RandomReal(0,2*M_PI);
		phi2 = RandomReal(0,2*M_PI);
		theta1 = 1/(RandomReal(1/max_theta,1/min_theta));
		theta2 = 1/(RandomReal(1/max_theta,1/min_theta));
		E0 = RandomReal(min_BeamEnergy, max_BeamEnergy);
        currVal = dSigmaBorn(x, theta1, theta2, phi1, phi2, E0);
        if(currVal > max3){
            max3 = currVal;
        }
    }
    if(max0 > max1 && max0 > max2 && max0 > max3){
        maxReturn = max0;
    }else if(max1 > max0 && max1 > max2 && max1 > max3){
        maxReturn = max1;
    }else if(max2 > max0 && max2 > max1 && max2 > max3){
        maxReturn = max2;
    }else if(max3 > max0 && max3 > max1 && max1 > max2){
        maxReturn = max1;
    }
    return(maxReturn);
}