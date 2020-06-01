#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "RandomReal.h"

struct density_rho0{
  double rho0, c_den[100], a_den[100]; //the overall normalization factor for nuclear charge densities
};

struct ZBQL0001{
  double ZBQLIX[43], B, C;
};

struct Brem_spect{
  double Eg[500], Br[500];
};


struct density_rho0 density_rho0;
struct ZBQL0001 zbql0001;
struct Brem_spect brem_spect;

double Brem(bool brem_init, bool cobrems, double E0, double Egamma);
double xsctn(double E0, int ztgt, double x, double theta1, double phi1, double theta2, double phi2, double pol, double m_part, bool nuc_FF);//units of nb/sr^2
double FF2(double q2, int ztgt, bool nuc_FF);
double analysis(double E0, double mtgt, double k1[3], double k2[3], double ktgt[3], double w_mumu, double t, double missing_mass, double m_part, bool pion_hypothesis);
void density_init(int ztgt);
double FF(double Q2, int ztgt);


//
//	To compile: "gfortran -ffixed-line-length-132 lepton_event_v17_4.f -o lepton_event_v17_4.exe”
//	To run: “./lepton_event_v17_4.exe”
//
//	version 6 of code weights the distribution of events by xs= acos( costheta ).  This is the most efficient version
//	of the code to run.  The nuclear and atomic form factors are included. 
//	version 8 allows for normal event weighting, or non-standard event weighting, theta=x**2, switchable with ‘standard_phase_space’. 
//		Non-standard is the most efficient. 
//	version 9: histograms J_T phi angle, uses simplified version of W_pol. 
//	version 10: histograms cross section arrays
//	version 11: histogram integrated cross section as a function of x and phi1. 
//	version 12: introduced logicals to do electron, muon, and muon with pion mass hypothesis histogramming. 
//	version 13: added some histogram options and controls
//	version 14: added some more controls, got rid of things in the code that weren’t useful 
//	version 15: brought over some features from the pi+pi- code, fixed some bugs.   Can control shape of the proton form factor. 
//	version 16: consistently use transverse momentum transfer squared everywhere, changed name to lepton_event
//	Uses rejection sampling: 
//	version 17: added feature to run one kinematic point
//	version 17_1: use weighting dcos theta/dx = (1-cos theta)**N , got rid of the option for bending the FF
//	version 17_2: use weighting theta = x**phase_space, with phase_space an integer >1 .   Set phase_space = 0 for dcos theta/dx =1. 
//		      makes a guess for the largest cross section, which seems to be correct
//	version 17_3: add option for log t plot, control proton rms radius independent of dipole FF parameter
//	version 17_4: read brem. file
//
//
double E0, theta_min, theta_max;
double pi, cos_max, cos_min, cross_sum, cross_max, cross_test, x_min, x_max, theta1, phi1, costheta1, theta2, phi2, costheta2, cross_section, total_xscn, k1[3], k2[3], m_e, m_part, m_muon, pol, x, q2, failure, w_mumu, xs1, xs2, xs_max, xs_min, Atgt[100], tgtlen[100], radlen[100], hbarc, W, jacobian, m_pi, phi_JT, delta_w, delta_x, delta_phi, mtgt, ktgt[3], Elab[2], missing_mass, t, delta_t, x_value, q2_T, W_pol, W_unpol, JS, JT[2], Rexp, xmax, theta1_max, theta2_max, phi1_max, phi2_max, temp, delta_log_t, frac_delta_t, delta_Egamma, data_array[200], error_array[200], Egamma_max, E_hi, E_lo, E_coherent, Egamma;
int iseed;
int i, itest[4], nevent, j, nfail, bad_max, i_array, j_array, ztgt, phase_space; 
bool hist_w, hist_x, hist_t, hist_phi_JT, hist_Egamma, output_event, hist_log_t, nuc_FF, muon, electron, pion_hypothesis, brem_init, cobrems;
//
//	Standard CPP configuration 
//	double ztgt = 82, E0 = 5.5, pol = 1.0, theta_min = 0.80,theta_max = 5.3; //Standard CPP configuration, min angle to TOF and max angle to MWPC
//
//	Standard GlueX configuration
//	data ztgt,E0,E_coherent,pol,theta_min,theta_max /1,11.0,8.7,1.0,0.90,13.12/	//standard GlueX config., min and max angles in deg. to TOF
int main(){
  std::cout << "i started the function\n";
    ztgt = 1, E0 = 11.0, E_coherent = 8.7, pol = 1.0, theta_min = 0.090, theta_max = 13.12;///standard GlueX config., min and max angles in deg. to TOF;
//	Set tagging interval
    E_hi = 8.8, E_lo = 8.6;
    itest[0] = 100000, itest[1] = 1000000, itest[2] = 10000000, itest[3] = 100000000, nevent = 10;
    m_e = 0.000511, m_muon = 0.105658, m_pi = 0.139570, hbarc = 0.197;
//	
//  Histogram parameters
    delta_w = 0.02, delta_t = 0.0002, delta_x = 0.02, delta_phi = 5.0, frac_delta_t = 0.2, delta_Egamma = 0.05; //units of GeV, GeV^2, //DIMENSIONless, degrees,
//											fractional bin width in t, GeV
//
//  Target information
    Atgt[0] = 1.0, tgtlen[0] = 0.0338, radlen[0] = 63.04;//tgtlen = # Rad lengths, radlen = rad length of material in g/cm^2;
//								this target is based on 30 cm LH2
    Atgt[5] = 12.0, tgtlen[5] = 0.05, radlen[5] = 42.70, density_rho0.c_den[5] = 2.45, density_rho0.a_den[5] = 0.524;//RL is in g/cm^2;
    Atgt[13] = 28.0, tgtlen[13] = 0.05, radlen[13] = 21.82, density_rho0.c_den[13] = 3.14, density_rho0.a_den[13] = 0.537;//units of c_den and a_den are fm;
    Atgt[19] = 40.0, tgtlen[19] = 0.05, radlen[19] = 16.12, density_rho0.c_den[19] = 3.51, density_rho0.a_den[19] = 0.563;//target length is in units of RL;
    Atgt[25] = 56.0, tgtlen[25] = 0.05, radlen[25] = 13.84, density_rho0.c_den[25] = 3.971,density_rho0.a_den[25] = 0.5935;//RL is g/cm^2;
    Atgt[49] = 116.0, tgtlen[49] = 0.05, radlen[49] = 8.82, density_rho0.c_den[49] = 5.416, density_rho0.a_den[49] = 0.552;
    Atgt[81] = 208.0, tgtlen[81] = 0.05, radlen[81] = 6.37;
    std::cout << "i finished array assignments\n";
//COMMON/density_rho0/rho0, c_den, a_den//the overall normalization factor for nuclear charge densities
//
    double zlo, zhi;
    zlo = 0.0, zhi = 1.0;
//
//       Initializes seed array etc. for random number generator.
//       The values below have themselves been generated using the
//       NAG generator.
//
    double ZBQLIX[43], B, C;
    ZBQLIX[0] = 8.001441, ZBQLIX[1] = 5.5321801, ZBQLIX[2] = 1.69570999, ZBQLIX[3] = 2.88589930, ZBQLIX[4] = 2.91581871, ZBQLIX[5] = 1.03842493, ZBQLIX[6] = 7.9952507,
    ZBQLIX[7] = 3.81202335, ZBQLIX[8] = 3.11575334, ZBQLIX[9] = 4.02878631, ZBQLIX[10] = 2.49757109, ZBQLIX[11] = 1.15192595, ZBQLIX[12] = 2.10629619, ZBQLIX[13] = 3.99952890,
    ZBQLIX[14] = 4.12280521, ZBQLIX[15] = 1.33873288, ZBQLIX[16] = 7.1345525, ZBQLIX[17] = 2.23467704, ZBQLIX[18] = 2.82934796, ZBQLIX[19] = 9.9756750, ZBQLIX[20] = 1.68564303,
    ZBQLIX[21] = 2.86817366, ZBQLIX[22] = 1.14310713, ZBQLIX[23] = 3.47045253, ZBQLIX[24] = 9.3762426, ZBQLIX[25] = 1.09670477, ZBQLIX[26] = 3.20029657, ZBQLIX[27] = 3.26369301,
    ZBQLIX[28] = 9.441177, ZBQLIX[29] = 3.53244738, ZBQLIX[30] = 2.44771580, ZBQLIX[31] = 1.59804337, ZBQLIX[32] = 2.07319904, ZBQLIX[33] = 3.37342907, ZBQLIX[34] = 3.75423178,
    ZBQLIX[35] = 7.0893571, ZBQLIX[36] = 4.26059785, ZBQLIX[37] = 3.95854390, ZBQLIX[38] = 2.0081010, ZBQLIX[39] = 5.9250059, ZBQLIX[40] = 1.62176640, ZBQLIX[41] = 3.20429173,
    ZBQLIX[42] = 2.63576576;
    std::cout << "i finished ZBQLIX\n";
    B = 4.294967291;
    C = 0.0;
//
    iseed = 0;//set equal to 0 for the // program to do a call to the system clock to initialize random number generator;
//			!for random number initialization
//
//  Initializations
//
    temp = -1;
    pi = acos(temp);
//
//Find the delta log t step
    delta_log_t = log10(1. + frac_delta_t);

//Set target mass in GeV
    mtgt = Atgt[ztgt]*0.931494;
    if (ztgt == 1) mtgt = 0.93828;
    i = 0;
    std::cout << "i got to the first while loop\n";
    while(i < 200){
        data_array[i] = 0;
        error_array[i] = 0;
        i++;
    }
    std::cout << "i finished the first while loop\n";
//Initialize Brem. distribution: select 1/Egamma or coherent Brems. file
    cobrems = true;//set true for scanfing coherent Brems. file, false for using a 1/Egamma distribution;
    if(cobrems == true){ // scanf coherent Brem. file
      brem_init = true;
      std::cout << "i am calling Brem function\n";
      temp = Brem(brem_init, cobrems, E0, Egamma);//scanf coherent brems file, then set brem_init = false;
    }
    int pause;
    std::cout << "this is the pause\n";
    //std::cin >> pause;
//Start logical assignments

//Only one of the histogram logicals can be true, the rest must be false.   The event output can be turned on independent of histograming.
    std::cout << "i got to the logicals\n";
    hist_w = false;
    hist_x = false;
    hist_t = false;
    hist_phi_JT = true;
    hist_log_t = false;
    hist_Egamma = false;
    output_event = true;
//Logical assignments
    phase_space = 4;//theta = x**phase_space, with int phase_space >1 . Note: phase_space = 4 seems to be fastest.;
//Set phase_space=0 for standard dcos theta/dx =1
    Rexp = float(phase_space);
    nuc_FF = true;
    muon = false;//this is how you change the particle type
//electron = .not.muon;
    if (muon == true){
        m_part = m_muon;
        pion_hypothesis = false; //set true if you want calculated invariant masses with pion assumption;
    }else{
        m_part = m_e;
        pion_hypothesis = false;//should always have this set false since we can distinguish e from mu;
    }
//
    if (nuc_FF) 
    {// setup the nuclear form factors
        density_rho0.c_den[ztgt] = density_rho0.c_den[ztgt]/hbarc;
        density_rho0.a_den[ztgt] = density_rho0.a_den[ztgt]/hbarc;
        density_init(ztgt); //for initializing the nuclear form factor(density_rho0.rho0 is density const.);
    }
//  
/*
    if((hist_w) || (hist_x) || (hist_t) || (hist_phi_JT) || (hist_log_t) || (hist_Egamma)){
        std::ofstream histFile;
        histFile.open("lepton_v17_4_hist.txt");
    }
    if (output_event){
        std::ofstream outFile;
        outFile.open("lepton_v17_4_event.txt");
    }/**/
//
//	End logical assignments
//
//***********************************
//	Evaluate the cross section at one kinematic point at coherent peak
    x = 0.5;
    theta1 = theta_min * (pi/180);
    theta2 = theta_min * (pi/180);
    phi1 = 90 * (pi/180);
    phi2 = 270 * (pi/180);
    cross_section = xsctn(E_coherent, ztgt, x, theta1, phi1, theta2, phi2, pol, m_part, nuc_FF);
    std::cout << E_coherent << " " << ztgt << " " << x << " " << theta1 << " " << phi1 << " " << theta2 << " " << phi2 << " " << pol << " " << m_part << " " << nuc_FF << " " << "\n";
    std::cout << " cross section nb/sr^2 = " << cross_section << "\n";
//*************************************
    theta_min = theta_min * pi/180;//switch to radians;
    theta_max = theta_max * pi/180;
    cos_max = cos(theta_min);
    cos_min = cos(theta_max);
// Limits on xs
    if (phase_space == 0){//dcos theta/dx = 1
        xs_max = cos_max;
        xs_min = cos_min;
    }else{
        xs_max = pow(theta_max, 1/Rexp);
        xs_min = pow(theta_min, 1/Rexp);
    }
//
//***************************************************************************
//
    j = 0;
    i = 0;
    std::cout << "i made it to the first loop\n";
    while(j < 3){//loop over 4 samplings of the phase space, each a factor of x10 larger, to see if the maximum
//cross section*Brem converges
        cross_max = 0;
        i = 0;
        while(i < itest[j]){//find maximum cross section in allowed phase space
            Egamma = RandomReal(E_lo, E_hi);//get tagged photon energy
            x = 0.5;//x_min + (x_max - x_min) RandomReal(zlo, zhi)//make a guess for the energy fraction
            phi1 = 90 * pi/180;//2.*pi RandomReal(zlo, zhi)//make a guess for phi1
            phi2 = 270 * pi/180;//2.*pi RandomReal(zlo, zhi)//make a guess for phi2
            xs2 = RandomReal(xs_min, xs_max);
            if (phase_space == 0) {// dcos theta/dx = 1
                theta1 = theta_min;//make a guess for theta1
                theta2 = acos(xs2);
                jacobian = 1;
            }else{
                theta1 = theta_min;//make a guess for theta1
                xs1 = pow(theta_min, (1/Rexp));
                theta2 = pow(xs2, phase_space);
                jacobian = (Rexp * pow(xs1, (phase_space - 1)) * sin(pow(xs1, phase_space))) * (Rexp * pow(xs2,(phase_space - 1)) * sin(pow(xs2,phase_space)));
            }
            cross_section = xsctn(Egamma, ztgt, x, theta1, phi1, theta2, phi2, pol, m_part, nuc_FF)*jacobian*Brem(brem_init, cobrems, E0, Egamma);
            if(cross_section > cross_max){
                cross_max = cross_section;
                Egamma_max = Egamma;
                xmax = x;
                theta1_max = theta1 * 180/pi;
                theta2_max = theta2 * 180/pi;
                phi1_max = phi1 * 180/pi;
                phi2_max = phi2 * 180/pi;
            
            }
            i++;
            std::cout << i << "\n";
        }
    
        std::cout << "test events " << itest[j] <<  " maximum xsctn*Brem " <<  cross_max << "\n";
        std::cout << "Egamma max " << Egamma_max << " x max " << xmax << " theta1 max " << theta1_max << " theta2 max " << theta2_max << " phi1 max " << phi1_max << " phi2 max " <<  phi2_max << "\n";
        j++;
        std::cout << j << "\n";
    }
    std::cout << "i have completed the maxtest loop\n";
//**********************************************************************************************************
// Loop over 4 samplings of the phase space at coherent peak, each a factor of x10 larger, to see if the integrated cross section converges
// Set limits on energy fraction
    j = 0;
    i = 0;
    x_max = (E_coherent - m_part)/E_coherent;
    x_min = m_part/E_coherent;
    std::cout << "i am starting the next loop\n";
    while(j < 3){//loop over 4 samplings of the phase space, each a factor of x10 larger, to see if the integrated
//cross section at the coherent peak converges
        cross_sum = 0;
        while(i < itest[j]){
            x = x_min + (x_max - x_min) * RandomReal(zlo, zhi);//energy fraction
            phi1 = 2*pi * RandomReal(zlo, zhi);
            phi2 = 2*pi * RandomReal(zlo, zhi);
            xs1 = RandomReal(xs_min, xs_max);
            xs2 = RandomReal(xs_min, xs_max);
            if (phase_space == 0){// dcos theta/dx = 1
                theta1 = acos(xs1);
                theta2 = acos(xs2);
                jacobian = 1;
            }else{
                theta1 = pow(xs1,phase_space);
                theta2 = pow(xs2,phase_space);
                jacobian = (Rexp * pow(xs1,(phase_space - 1)) * sin(pow(xs1,phase_space))) * (Rexp*pow(xs2,(phase_space-1))*sin(pow(xs2,phase_space)));
            }
            cross_section = xsctn(E_coherent, ztgt, x, theta1, phi1, theta2, phi2, pol, m_part, nuc_FF)*jacobian;
            cross_sum = cross_sum + cross_section;
            i++;
            std::cout << i << "\n";
        }
        total_xscn = cross_sum/float(itest[j])*pow((xs_max - xs_min), 2)*pow((2*pi), 2) * (x_max - x_min);
        std::cout << "test events " << itest[j] << " Egamma " <<  E_coherent << " total cross section nb " <<  total_xscn << "\n";
        j++;
    }
    std::cout << "i am about to start event generation\n";


//*************************************************************************
//
//	Start event generation
//
//Use the widest possible range in x by using the maximum accepted tagged photon energy, then test it
    x_max = (E_hi - m_part)/E_hi;//largest possible x
    x_min = m_part/E_hi;//smallest possible x
    //
    nfail = 0;
    bad_max = 0;
    i = 0;
    //
    while(i < nevent){
    g100:
        Egamma = RandomReal(E_lo, E_hi);//get tagged photon energy
        x = RandomReal(x_min, x_max);//energy fraction
  //	Test x to make sure it's within the allowed range for the photon energy Egamma
        if((x >= ((Egamma - m_part)/Egamma)) || (x <= (m_part/Egamma))) goto g100; // x is out of range, try again
        phi1 = 2 * pi * RandomReal(zlo, zhi);
        phi2 = 2 * pi * RandomReal(zlo, zhi);
        xs1 = RandomReal(xs_min, xs_max);
        xs2 = RandomReal(xs_min, xs_max);
        if(phase_space == 0){// dcos theta/dx = 1
            theta1 = acos(xs1);
            theta2 = acos(xs2);
            jacobian = 1;
        }else{
            theta1 = pow(xs1, phase_space);
            theta2 = pow(xs2, phase_space);
            jacobian = (Rexp * pow(xs1,(phase_space-1)) * sin(pow(xs1,phase_space))) * (Rexp*pow(xs2,(phase_space-1)) * sin(pow(xs2,phase_space)));
        }

        cross_section = xsctn(Egamma, ztgt, x, theta1, phi1, theta2, phi2, pol, m_part, nuc_FF)*jacobian*Brem(brem_init, cobrems, E0, Egamma);
        if (cross_section > cross_max){
            bad_max = bad_max + 1;//an occurrence of cross section larger than cross_max, not supposed to happen
            std::cout <<  "bad max cross section= " <<  cross_section << "\n";
        }
        cross_test = cross_max * RandomReal(zlo, zhi);
        if (cross_test > cross_section){//selection fails
            nfail = nfail + 1;
            goto g100;
        }
//	Event selection succeeds: 
//
        analysis(Egamma, mtgt, k1, k2, ktgt, w_mumu, t, missing_mass, m_part, pion_hypothesis);//analyze the event;
//							
//		Do the histogramming
//
        if(hist_w) i_array = int((w_mumu - 0.200)/delta_w);//w distribution;
        if(hist_x) i_array = int(x/delta_x);//x distribution
        if(hist_t) i_array = int(t/delta_t);//t distribution;
        if(hist_phi_JT) i_array = int(phi_JT*180/pi/delta_phi); //JT phi distribution in degrees;
        if(hist_log_t) i_array = int((log10(t) + 6)/delta_log_t);// 10^ - 6 GeV^2 is bin 0;
        if(hist_Egamma) i_array = int(Egamma/delta_Egamma);//photon energy distribution;
        if(i_array < 0) i_array = 0;
        if (i_array > 200) i_array = 200;
        data_array[i_array] = data_array[i_array] + 1;
//
//	3-momentum event output
        if(output_event){
            std::ofstream outputFile;
            outputFile.open("lepton_v17_4_event.txt");
            outputFile << Egamma << " " << k1[0] << " " << k1[1] << " " << k1[2] << " " << k2[0] << " " << k2[1] << " " << k2[2] << " " << ktgt[0] << " " << ktgt[1]<< " " << ktgt[2] << " " << "\n";
        }
	//g200:
// format(2x, f6.3, 1x, 9(f10.6, 1x))
//
        //if (i % 100 == 0) std::cout << ' event # ' << i;
        i++;
        std::cout << i << "\n";
    }


    std::cout << "i have finished the event generation\n";
//   Event generation ends
//	
    float failure = float(nfail)/float(nevent);
    std::cout <<  "Failures per event = " << failure << " Events with cross section exceeding max xsctn = " << bad_max << "\n";
    std::ofstream histFile;
    histFile.open("lepton_v17_4_hist.txt");
    int i = 0;
    while(i < 200){//printf out the arrays
        std::cout << "histogramming... ";
        if(hist_w){
        
            //std::ofstream histFile;
            //histFile.open("lepton_v17_4_hist.txt");
            x_value = float(i) * delta_w + 0.200;
            error_array[i] = sqrt(data_array[i]);
            histFile << x_value << " " << data_array[i] << " " << error_array[i] << "\n";
            //histFile.close();
        }else if (hist_x){
            //std::ofstream histFile;
            //histFile.open("lepton_v17_4_hist.txt");
            x_value = float(i) * delta_x;
            error_array[i] = sqrt(data_array[i]);
            histFile << x_value << " " << data_array[i] << " " << error_array[i] << "\n";
            //histFile.close();
        }else if(hist_t){
            //std::ofstream histFile;
            //histFile.open("lepton_v17_4_hist.txt");
            x_value = float(i)*delta_t;
            error_array[i] = sqrt(data_array[i]);
            histFile << x_value << " " << data_array[i] << " " << error_array[i] << "\n";
            //histFile.close();
        }else if(hist_phi_JT){
            //std::ofstream histFile;
            //histFile.open("lepton_v17_4_hist.txt");
            x_value = float(i)*delta_phi;
            if((x_value == 0.) || (x_value == 360.)) data_array[i] = 2.*data_array[i]; // this is a binning problem
            error_array[i] = sqrt(data_array[i]);
            histFile << x_value << " " << data_array[i] << " " << error_array[i] << "\n";
            //histFile.close();
        }else if(hist_log_t){
            //std::ofstream histFile;
            //histFile.open("lepton_v17_4_hist.txt");
            x_value = pow(10, (float(i) * delta_log_t-6));
            error_array[i] = sqrt(data_array[i])/(pow(10,(float(i+1)*delta_log_t-6)) - pow(10,(float(i-1) * delta_log_t - 6))) * 2;
            data_array[i] = data_array[i]/(pow(10,(float(i+1) * delta_log_t-6)) - pow(10,(float(i-1) * delta_log_t - 6))) * 2;
            histFile << x_value << " " << data_array[i] << " " << error_array[i] << "\n";
            //histFile.close();
        }else if(hist_Egamma){
            //std::ofstream histFile;
            //histFile.open("lepton_v17_4_hist.txt");
            x_value = float(i) * delta_Egamma;
            error_array[i] = sqrt(data_array[i]);
            histFile << x_value << " " << data_array[i] << " " << error_array[i] << "\n";
            //histFile.close();
        }
        i++;
    }
    histFile.close();


//

//
//exit(0);
}


//
//*******************************************************************
//
// --------------------------------------------
double Brem(bool brem_init, bool cobrems, double E0, double Egamma){
    double Eg[300], Br[300];
    int i, imax, ipoint;
    double bremOut;
    i = 0;
    FILE *CBD;
    if (brem_init == true){//fopen and scanf coherent Brem file
      CBD = fopen("CobremsDistribution.dat", "r");
      //std::cout << "i have opened the .dat file\n";
      while(i < 300){
	    fscanf(CBD,"%lf %lf", &brem_spect.Eg[i], &brem_spect.Br[i]);
	    //std::cout << i << "i have scanned .dat file\n";
	    i = i + 1;
      }
      //std::cout << "i have finished the scan\n";  
      imax = i - 1;
      fclose(CBD);
      //std::cout << "i have closed the file \n";
      brem_init = false;//done with initialization;
      //return bremOut;
    }
    if (!cobrems){
        bremOut = E0/Egamma;
        //return bremOut;
    }
    //
    if (cobrems){ //return coherent brems distribution
        ipoint = int((Egamma + .02)/.04);
        bremOut = brem_spect.Br[ipoint];
        //return bremOut;
    }
    //std::cout << "i am about to return\n";
    return bremOut;
}


//	
//******************************************************************
// The reference for this is Bakmaev et al., Physics Letters B 660 (2008) 494-500, eqn. 23
// Had to correct a mistake in Eqn. 23 of their paper.   The cos(phi1 + phi2) term should be 
// multiplied by 2.  You can see this by comparing Wp in Eqn. 22 with the vector current part of Eqn. 23
//
// --------------------------------------------
double xsctn(double E0, int ztgt, double x, double theta1, double phi1, double theta2, double phi2, double pol, double m_part, bool nuc_FF)//units of nb/sr^2;
{
//implicit none
  double Z, W_unpol, W_pol, q2_T;
    double alpha, hbarc;
    double xsctnOut;
    double pi, E1, E2, k1_mag, k2_mag, p1, p2, q2, c1, c2, JS, JT[2], FF_nuc, FF_TFM;
    int i;
    //
    alpha = 7.297352e-3, hbarc = 0.197;
    //
    pi = acos(-1);
    //
    Z = float(ztgt);
    E1 = E0 * x;
    E2 = E0 * (1 - x);
    k1_mag = sqrt(pow(E1, 2) - pow(m_part, 2));
    k2_mag = sqrt(pow(E2, 2) - pow(m_part, 2));
    //std::cout << "k1_mag " << k1_mag <<"\n";
    //std::cout << "k2_mag " << k2_mag <<"\n";
    k1[0] = k1_mag * sin(theta1) * cos(phi1);
    k1[1] = k1_mag * sin(theta1) * sin(phi1);
    k1[2] = k1_mag * cos(theta1);
    //std::cout << "k1[0] " << k1[0] <<"\n";
    //std::cout << "k1[1] " << k1[1] <<"\n";
    //std::cout << "k1[2] " << k1[2] <<"\n";
    k2[0] = k2_mag * sin(theta2) * cos(phi2);
    k2[1] = k2_mag * sin(theta2) * sin(phi2);
    k2[2] = k2_mag * cos(theta2);
    //std::cout << "k2[0] " << k2[0] <<"\n";
    //std::cout << "k2[1] " << k2[1] <<"\n";
    //std::cout << "k2[2] " << k2[2] <<"\n";
    //
    p1 = sqrt(pow(k1[0], 2) + pow(k1[1],2));//transverse momenta of muon #1, GeV
    //std::cout << "p1 " << p1 <<"\n";
    p2 = sqrt(pow(k2[0], 2) + pow(k2[1],2));//transverse momenta of muon #2, GeV
    //std::cout << "p2 " << p2 <<"\n";
    q2_T = pow((k1[0] + k2[0]),2) + pow((k1[1] + k2[1]),2); // this is transverse momentum transfer squared
    //std::cout << "q2_T " << q2_T <<"\n";
    q2 = q2_T + pow((E0 - k1[2] - k2[2]), 2); // this is 3 - momentum transfer squared
    //std::cout << "q2 " << q2 <<"\n";
    c1 = pow(p1, 2) + pow(m_part, 2);
    //std::cout << "c1 " << c1 <<"\n";
    c2 = pow(p2, 2) + pow(m_part, 2);
    //std::cout << "c2 " << c2 << "\n";
    JS = 1/c1 - 1/c2;//units of 1/GeV^2//scalar current, units of GeV^ - 2;
    //std::cout << "JS " << JS <<"\n";
    i = 0;
    while(i < 2){
      JT[i] = k1[i]/c1 + k2[i]/c2;//vector current, units of GeV^ - 1;
      //std::cout << "JT[i] " << JT[i] <<"\n";
      i++;
      
    }


    phi_JT = acos(JT[0]/sqrt(pow(JT[0],2) + pow(JT[1], 2)));//phi angle of JT wrt to x axis, radians
    if (JT[1] < 0) phi_JT = 2*pi - phi_JT;
    //std::cout << "phi_JT " << phi_JT <<"\n";
//
    W_unpol = pow(m_part,2) * pow(JS,2) + (pow(x,2) + pow((1 -x),2)) * (pow(JT[1],2) + pow(JT[2], 2));
    //std::cout << "W_unpol " << W_unpol <<"\n";
//     	W_pol = -2.*x*(1.-x)*((p2/c2)**2*cos(2.*phi2)+(p1/c1)**2*cos(2.*phi1)+2.*(p2/c2)*(p1/c1)*cos(phi1+phi2)) !the Bakmaev expression
//	xsctn=2.*alpha**3*Z**2*E0**4*x**2*(1.-x)**2/(pi**2*q2_T**2)*(W_unpol+pol*W_pol) ! note the absence of cos(2phi_JT) in this expression
//     &	*hbarc**2/100.*1.e9*FF2(q2_T,ztgt,nuc_FF) !units of nb/sr^2  The denominator uses the transverse 3-momentum transfer^2, 
    W_pol = -2 * x * (1 - x) * (pow(JT[1], 2) + pow(JT[2], 2));//this is my reduction of the Bakmaev equations;
    //std::cout << "W_pol " << W_pol <<"\n";
    xsctnOut = 2 * pow(alpha,3) * pow(Z, 2) * pow(E0, 4) * pow(x, 2) * pow((1 - x), 2)/(pow(pi,2) * pow(q2_T, 2)) * (W_unpol + pol * cos(2 * phi_JT) * W_pol) * pow(hbarc, 2)/100 * 1e9 * FF2(q2_T, ztgt, nuc_FF);
    //std::cout << xsctnOut << "\n";
//this contains the cos(2phi_JT) term*hbarc**2/100.*1.e9*FF2(q2_T, ztgt, nuc_FF) //units of nb/sr^2 The denominator uses the transverse 3 - momentum transfer^2
//
    return xsctnOut;
};


//
//****************************************************************************************
//
// --------------------------------------------
double FF2(double q2, int ztgt, bool nuc_FF)
{
    double hbarc, z, FF_nuc, FF_TFM, alpha[3], b[3], b0, m_e, c, outputFF2;
    int i;
    m_e = 0.511e-3, hbarc = 0.197;
    alpha[0] = 0.1; 
    alpha[1] = 0.55; 
    alpha[2] = 0.35;
    b[0] = 6.0; 
    b[1] = 1.2; 
    b[2] = 0.3;
    z = float(ztgt);
    c = m_e*pow(z,0.333);
    FF_nuc = FF(q2, ztgt);
    FF_TFM = 1;
    i = 0;
    while(i < 3){
      FF_TFM = FF_TFM - alpha[i]*q2/(q2 + pow((b[i]*c),2));
      i++;
    }
    if(nuc_FF){
	    outputFF2 = pow((FF_nuc - FF_TFM),2);
    }else{
	    outputFF2 = 1;
    }
    return outputFF2;
}
//
//****************************************************************************************
//
// --------------------------------------------
double analysis(double E0, double mtgt, double k1[3], double k2[3], double ktgt[3], double w_mumu, double t, double missing_mass, double m_part, bool pion_hypothesis)
{
    // implicit none
  double E1, E2, ks[3], m_x, analysisOut;
    double m_pi = 0.139570;
    pi = acos(-1);
    E1 = sqrt(pow(k1[0], 2) + pow(k1[1], 2) + pow(k1[2], 2) + pow(m_part, 2));//lepton energies
    E2 = sqrt(pow(k2[0], 2) + pow(k2[1], 2) + pow(k2[2], 2) + pow(m_part, 2));
    ks[0] = k1[0] + k2[0];//lepton summed momentum
    ks[1] = k1[1] + k2[1];
    ks[2] = k1[2] + k2[2];
    ktgt[0] = -ks[0];//target momentum;
    ktgt[1] = -ks[1];
    ktgt[2] = E0 - ks[2];
    missing_mass = sqrt(pow((E0 + mtgt - E1 - E2), 2) - pow(ktgt[0], 2) - pow(ktgt[1],2) - pow(ktgt[2], 2));
    t = pow(ks[0], 2) + pow(ks[1], 2) + pow((E0 - ks[2]), 2) - pow((E0 - E1 - E2), 2);//4 - momentum transfer squared to nucleus, this is positive;
    //
    //		mu mu invariant mass, possibly with pion hypothesis
    m_x = m_part;
    if (pion_hypothesis) m_x = m_pi;
    E1 = sqrt(pow(k1[0], 2) + pow(k1[1],2) + pow(k1[2], 2) + pow(m_x, 2));//need to put in the mass hypothesis
    E2 = sqrt(pow(k2[0], 2) + pow(k2[1], 2) + pow(k2[2], 2) + pow(m_x, 2));
    w_mumu = sqrt(pow(E1 + E2, 2) - pow(ks[0], 2) - pow(ks[1], 2) - pow(ks[2],2));
    return analysisOut;
}





//c******************************************************************
//c
// --------------------------------------------
void density_init(int ztgt)
{
    // implicit none
    double pi, c_den[100], a_den[100],rho0, w;
    int i;
    //COMMON/density_rho0/rho0, c_den, a_den
    //c
    pi = acos(-1);
    //c
    rho0 = 0;
    if((ztgt == 82) || (ztgt == 1)) return;
    //c	These equations have to do with Fermi distribution, reference? 
    w = 4 * pi * c_den[ztgt]/3 * (pow((pi * a_den[ztgt]), 2) + pow(c_den[ztgt], 2));
    i = 0;
    while(i < 10){
        w = w + 8 * pi * pow(a_den[ztgt], 3) * pow((-1), (i - 1)) * exp(-float(i)*c_den[ztgt]/a_den[ztgt])/pow(float(i), 3);
        i++;
    }
    density_rho0.rho0 = 1/w;
    return;
};


//
//******************************************************************

// --------------------------------------------
double FF(double Q2, int ztgt)
{
    // implicit none
    double q02, hbarc, Q, gamma, r[12], A[12], rho0, c_den[100], a_den[100], pi, norm, proton_rms;
    int i = 0;
    double returnFF;
    //c
    q02 = 0.71, proton_rms = 0.879;//proton dipole form factor parameter GeV^2, proton rms radius fm
    hbarc = 0.197;
    //COMMON/density_rho0/rho0, c_den, a_den
    //c
    r[0] = 0.1, A[0] = 0.003845;
    r[1] = 0.7, A[1] = 0.009724;
    r[2] = 1.6, A[2] = 0.033093;
    r[3] = 2.1, A[3] = 0.000120;
    r[4] = 2.7, A[4] = 0.083107;
    r[5] = 3.5, A[5] = 0.080869;
    r[6] = 4.2, A[6] = 0.139957;
    r[7] = 5.1, A[7] = 0.260892;
    r[8] = 6.0, A[8] = 0.336013;
    r[9] = 6.6, A[9] = 0.033637;
    r[10] = 7.6, A[10] = 0.018729;
    r[11] = 8.7, A[11] = 0.000020;
    gamma = 1.388;
    //
    //  Select the FF
    //
    pi = acos(-1);
    Q = sqrt(Q2);
    //
    if (ztgt == 1) 
    {//proton
        returnFF = 1/(pow((1. + Q2/q02),2)) + 2 * Q2/q02 - 1/6 * Q2 * pow(proton_rms, 2)/pow(hbarc, 2);
    }else if(ztgt == 82){//lead
        returnFF = 0;
	    do{
	    returnFF = returnFF + A[i] * (pow(gamma,2) * cos(Q * r[i]/hbarc) + 2 * r[i] * hbarc/Q * sin(Q * r[i]/hbarc))/(pow(gamma, 2) + 2 * pow(r[i], 2)) * exp(-Q2/4 * pow(gamma, 2)/pow(hbarc, 2));
	    i++;
        }while(i < 12);
    }else{ //for everything else use 2 - parameter fermi model, reference ?
      returnFF = 4 * pow(pi,2) * rho0 * pow(a_den[ztgt], 3)/(pow((Q * a_den[ztgt]), 2) * pow((sinh(pi * Q * a_den[ztgt])),2)) * (pi * Q * a_den[ztgt] * cosh(pi * Q * a_den[ztgt]) * sin(Q * c_den[ztgt]) - Q * c_den[ztgt] * cos(Q * c_den[ztgt]) * sinh(pi * Q * a_den[ztgt]));
      i = 0;
      do{
	    returnFF = returnFF + 8 * pi * rho0 * pow(a_den[ztgt], 3) * pow((-1), (i - 1)) * float(i) * exp(-float(i) * c_den[ztgt]/a_den[ztgt])/pow((pow(float(i), 2) + pow((Q * a_den[ztgt]),2)),2);
	    i++;
      }while(i < 10);
    }
return returnFF;
}
