//Date: 6/2/2020 9:08 PM
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "RandomReal.h"

struct density_rho0{
  double rho0, c_den[100], a_den[100]; //the overall normalization factor for nuclear charge densities
};


struct Brem_spect{
  double Eg[500], Br[500];
};


struct density_rho0 density_rho0;
struct Brem_spect brem_spect;

double Brem(bool brem_init, bool cobrems, double E0, double Egamma);
double xsctn(double E0, int ztgt, double x, double theta1, double phi1, double theta2, double phi2, double pol, double m_part, bool nuc_FF, double* phi_JT, double* k10, double* k11, double* k12, double* k20, double* k21, double* k22);//units of nb/sr^2
double FF2(double q2, int ztgt, bool nuc_FF);
double analysis(double E0, double mtgt, double missing_mass, double m_part, bool pion_hypothesis, double* w_mumu, double* t, double k10in, double k11in, double k12in, double k20in, double k21in, double k22in);
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
double mass_pi, delta_cos, delta_x_value, bin_width, x_value_hi, x_value_lo, total_xscn_old, delta_total_xscn, pi, cos_max, cos_min, cross_sum, cross_max, cross_test, x_min, x_max, theta1, phi1, costheta1, theta2, phi2, costheta2, cross_section, total_xscn, m_e, m_part, m_muon, pol, x, q2, failure, w_mumu, xs1, xs2, xs_max, xs_min, Atgt[100], tgtlen[100], radlen[100], hbarc, W, jacobian, phi_JT, delta_w, delta_x, delta_phi, mtgt, ktgt[3], Elab[2], missing_mass, t, delta_t, x_value, q2_T, W_pol, W_unpol, JS, JT[2], Rexp, xmax, theta1_max, theta2_max, phi1_max, phi2_max, temp, delta_log_t, frac_delta_t, delta_Egamma, data_array[200], error_array[200], Egamma_max, E_hi, E_lo, E_coherent, Egamma, data_array_w[200], data_array_x[200], data_array_t[200], data_array_phi_JT[200], data_array_nonlinear_t[200], data_array_Egamma[200], data_array_cos[200], cross_max_old, delta_cross_max, k10, k11, k12, k20, k21, k22, k10in, k11in, k12in, k20in, k21in, k22in;
int iseed;
int i, itest[4], nevent, j, nfail, bad_max, j_array, ztgt, phase_space; 
bool hist_w, hist_x, hist_t, hist_phi_JT, hist_Egamma, output_event, hist_nonlinear_t, muon, electron, pion_hypothesis, brem_init, cobrems, hist_cos, integral_xsctn, w_cut, nuc_FF, verbose_output;
int xi_array, ti_array, wi_array, phi_JTi_array, cos_maxi_array, nonlinear_ti_array, Egammai_array;
//	Standard CPP configuration 
//	double ztgt = 82, E0 = 5.5, pol = 1.0, theta_min = 0.80,theta_max = 5.3; //Standard CPP configuration, min angle to TOF and max angle to MWPC
//
//	Standard GlueX configuration
//	data ztgt,E0,E_coherent,pol,theta_min,theta_max /1,11.0,8.7,1.0,0.90,13.12/	//standard GlueX config., min and max angles in deg. to TOF
int main(){
    ztgt = 1, E0 = 11.0, E_coherent = 8.7, pol = 1.0, theta_min = 0.75, theta_max = 13.12;///standard GlueX config., min and max angles in deg. to TOF;
//	Set tagging interval
    double w_min = 0.25, w_max = 0.621;
    E_hi = 8.8, E_lo = 8.6;
    itest[0] = 100000, itest[1] = 1000000, itest[2] = 10000000, itest[3] = 100000000, nevent = 1000;
    m_e = 0.000511, m_muon = 0.105658, mass_pi = 0.139570, hbarc = 0.197;
//	
//  Histogram parameters
    delta_w = 0.02, delta_t = 0.0002, delta_x = 0.02, delta_phi = 5.0, frac_delta_t = 0.2, delta_Egamma = 0.05, delta_cos = 0.0001; //units of GeV, GeV^2, //DIMENSIONless, degrees,
//											fractional bin width in t, GeV
//
//  Target information
    Atgt[1] = 1.0, tgtlen[1] = 0.0338, radlen[1] = 63.04;//tgtlen = # Rad lengths, radlen = rad length of material in g/cm^2;
//								this target is based on 30 cm LH2
    Atgt[6] = 12.0, tgtlen[6] = 0.05, radlen[6] = 42.70, density_rho0.c_den[6] = 2.45, density_rho0.a_den[6] = 0.524;//RL is in g/cm^2;
    Atgt[14] = 28.0, tgtlen[14] = 0.05, radlen[14] = 21.82, density_rho0.c_den[14] = 3.14, density_rho0.a_den[14] = 0.537;//units of c_den and a_den are fm;
    Atgt[20] = 40.0, tgtlen[20] = 0.05, radlen[20] = 16.12, density_rho0.c_den[20] = 3.51, density_rho0.a_den[20] = 0.563;//target length is in units of RL;
    Atgt[26] = 56.0, tgtlen[26] = 0.05, radlen[26] = 13.84, density_rho0.c_den[26] = 3.971,density_rho0.a_den[26] = 0.5935;//RL is g/cm^2;
    Atgt[50] = 116.0, tgtlen[50] = 0.05, radlen[50] = 8.82, density_rho0.c_den[50] = 5.416, density_rho0.a_den[50] = 0.552;
    Atgt[82] = 208.0, tgtlen[82] = 0.05, radlen[82] = 6.37;
//COMMON/density_rho0/rho0, c_den, a_den//the overall normalization factor for nuclear charge densities
//
    double zlo, zhi;
    zlo = 0.0, zhi = 1.0;
//
//       Initializes seed array etc. for random number generator.
//       The values below have themselves been generated using the
//       NAG generator.
//
//  Initializations
//
//
//Find the delta log t step
    //delta_log_t = log10(1. + frac_delta_t);

//Set target mass in GeV
    mtgt = Atgt[ztgt]*0.931494;
    if (ztgt == 1) mtgt = 0.93828;
    i = 0;
    while(i < 200){
        data_array_w[i] = 0;
        data_array_x[i] = 0;
        data_array_t[i] = 0;
        data_array_phi_JT[i] = 0;
        data_array_nonlinear_t[i] = 0;
        data_array_Egamma[i] = 0;
        data_array_cos[i] = 0;
        error_array[i] = 0;
        i++;
    }
    //Initialize parameters for variable bin-width t distribution: bin_width_t = alpha*(t+beta)**gamma
    double bin_t_min = 0.0001;
    double alpha_t = 0.02;
    double gamma_t = 0.47;
    double beta_t = pow((bin_t_min/alpha_t),(1/gamma_t));

//Start logical assignments

//Can turn on any or all of these logicals.

    hist_w = true;
    hist_x = false;
    hist_t = false;
    hist_phi_JT = false;
    hist_nonlinear_t = false;
    hist_Egamma = false;
    hist_cos = false;
    output_event = true;     //for writing event file
    integral_xsctn = true;  //set true for outputting integrated cross sections
    w_cut = false;           //set true for applying w cuts to the data
    nuc_FF = true;           //set true for using a nuclear form factor, set false for FF=1;
    cobrems = true;          //set true for scanfing coherent Brems. file, false for using a 1/Egamma distribution;
    muon = false;            //set true for muons,false for electrons
    verbose_output = true;  //set true if you want to see test outputs printed to your screen
    electron = !muon;

//Initialize Brem. distribution: select 1/Egamma or coherent Brems. file
    
    if(cobrems == true){ // scanf coherent Brem. file
      brem_init = true;
      temp = Brem(brem_init, cobrems, E0, Egamma);//scanf coherent brems file, then set brem_init = false;
    }
//Logical assignments
    phase_space = 8;//theta = x**phase_space, with int phase_space >1 . Note: phase_space = 8 is the fastest for e+e-;
//Set phase_space=0 for standard dcos theta/dx =1
    Rexp = float(phase_space);

    if (muon == true){
        m_part = m_muon;
        pion_hypothesis = false; //set true if you want calculated invariant masses with pion assumption;
    }else{
        m_part = m_e;
        pion_hypothesis = false; //should always have this set false since we can distinguish e from mu;
    }
//
    if (nuc_FF){// setup the nuclear form factors
        density_rho0.c_den[ztgt] = density_rho0.c_den[ztgt]/hbarc;
        density_rho0.a_den[ztgt] = density_rho0.a_den[ztgt]/hbarc;
        density_init(ztgt); //for initializing the nuclear form factor(density_rho0.rho0 is density const., nothing to do with rho0 meson);
    }

//	Evaluate the cross section at one kinematic point at coherent peak
    x = 0.4;
    theta1 = theta_min * (M_PI/180);
    theta2 = theta_min * (M_PI/180);
    phi1 = 90 * (M_PI/180);
    phi2 = 270 * (M_PI/180);
    cross_section = xsctn(E_coherent, ztgt, x, theta1, phi1, theta2, phi2, pol, m_part, nuc_FF, &phi_JT, &k10, &k11, &k12, &k20, &k21, &k22);
    if(verbose_output == true){
        std::cout << E_coherent << " " << ztgt << " " << x << " " << theta1 << " " << phi1 << " " << theta2 << " " << phi2 << " " << pol << " " << m_part << " " << nuc_FF << " " << "\n";
        std::cout << " cross section nb/sr^2 = " << cross_section << "\n";
    }
//*************************************
    theta_min = theta_min * M_PI/180;//switch to radians;
    theta_max = theta_max * M_PI/180;
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
    while(j < 2){//loop over 4 samplings of the phase space, each a factor of x10 larger, to see if the maximum
//cross section*Brem converges
        cross_max_old = cross_max;
        i = 0;
        while(i < itest[j]){//find maximum cross section in allowed phase space
            Egamma = RandomReal(E_lo, E_hi);//get tagged photon energy
            x = 0.5;//x_min + (x_max - x_min) RandomReal(zlo, zhi)//make a guess for the energy fraction
            phi1 = 90 * M_PI/180;//2.*pi RandomReal(zlo, zhi)//make a guess for phi1
            phi2 = 270 * M_PI/180;//2.*pi RandomReal(zlo, zhi)//make a guess for phi2
            xs1 = RandomReal(xs_min, xs_max);
            xs2 = RandomReal(xs_min, xs_max);
            //std::cout << xs1 << " " << xs2 << "\n";
            if (phase_space == 0) {// dcos theta/dx = 1
                theta1 = acos(xs1);//make a guess for theta1
                theta2 = acos(xs2);
                jacobian = 1;
            }else{
                theta1 = pow(xs1, phase_space);//make a guess for theta1
                theta2 = pow(xs2, phase_space);
                jacobian = (Rexp * pow(xs1, (phase_space - 1)) * sin(pow(xs1, phase_space))) * (Rexp * pow(xs2,(phase_space - 1)) * sin(pow(xs2,phase_space)));
            }
            cross_section = xsctn(Egamma, ztgt, x, theta1, phi1, theta2, phi2, pol, m_part, nuc_FF, &phi_JT, &k10, &k11, &k12, &k20, &k21, &k22)*jacobian*Brem(brem_init, cobrems, E0, Egamma);
            if(cross_section > cross_max){
                cross_max = cross_section;
                Egamma_max = Egamma;
                xmax = x;
                theta1_max = theta1 * 180/M_PI;
                theta2_max = theta2 * 180/M_PI;
                phi1_max = phi1 * 180/M_PI;
                phi2_max = phi2 * 180/M_PI;
            
            }
            i++;
            std::cout << i << "\n";
        }
        delta_cross_max = (cross_max - cross_max_old)/cross_max*100;
        if(verbose_output == true){
            std::cout << "test events " << itest[j] <<  " maximum xsctn*Brem " <<  cross_max << "\n";
            std::cout << "Egamma max " << Egamma_max << " x max " << xmax << " theta1 max " << theta1_max << " theta2 max " << theta2_max << " phi1 max " << phi1_max << " phi2 max " <<  phi2_max << "\n";
        }
        j++;
        std::cout << j << "\n";
    }
//**********************************************************************************************************
// Loop over 4 samplings of the phase space at coherent peak, each a factor of x10 larger, to see if the integrated cross section converges
// Set limits on energy fraction
    if(integral_xsctn == true){
        j = 0;
        i = 0;
        x_max = (E_coherent - m_part)/E_coherent;
        x_min = m_part/E_coherent;
        total_xscn = 0;
        while(j < 2){//cross section at the coherent peak converges
            cross_sum = 0;
	    i = 0;
            while(i < itest[j]){
                x = x_min + (x_max - x_min) * RandomReal(zlo, zhi);//energy fraction
                phi1 = 2*M_PI * RandomReal(zlo, zhi);
                phi2 = 2*M_PI * RandomReal(zlo, zhi);
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
                cross_section = xsctn(E_coherent, ztgt, x, theta1, phi1, theta2, phi2, pol, m_part, nuc_FF, &phi_JT, &k10, &k11, &k12, &k20, &k21, &k22)*jacobian;
                k10in = k10;
                k11in = k11;
                k12in = k12;
                k20in = k20;
                k21in = k21;
                k22in = k22;
                analysis(Egamma, mtgt, missing_mass, m_e, pion_hypothesis , &w_mumu, &t, k10in, k11in, k12in, k20in, k21in, k22in);
                cross_sum = cross_sum + cross_section;
            i++;
            std::cout << i << "\n";
	    }
        total_xscn = cross_sum/float(itest[j])*pow((xs_max - xs_min), 2)*pow((2*M_PI), 2) * (x_max - x_min);
        delta_total_xscn = (total_xscn-total_xscn_old)/total_xscn*100;
        std::cout << "test events " << itest[j] << " Egamma " <<  E_coherent << " total cross section nb " <<  total_xscn << "\n";
        j++;
	std::cout << j << "\n";
        }
    }


//*************************************************************************
//
//	Start event generation
//
//Use the widest possible range in x by using the maximum accepted tagged photon energy, then test it
    x_max = (E_hi - m_part)/E_hi;//largest possible x
    x_min = m_part/E_hi;//smallest possible x
    nfail = 0;
    bad_max = 0;
    i = 0;
    while(i < nevent){
      std::cout << "nevent: " << i << "\n";

    g100:
        Egamma = RandomReal(E_lo, E_hi);//get tagged photon energy
        x = RandomReal(x_min, x_max);//energy fraction
  //	Test x to make sure it's within the allowed range for the photon energy Egamma
        if((x >= ((Egamma - m_part)/Egamma)) || (x <= (m_part/Egamma))) goto g100; // x is out of range, try again
        phi1 = 2 * M_PI * RandomReal(zlo, zhi);
        phi2 = 2 * M_PI * RandomReal(zlo, zhi);
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

        cross_section = xsctn(Egamma, ztgt, x, theta1, phi1, theta2, phi2, pol, m_part, nuc_FF, &phi_JT, &k10, &k11, &k12, &k20, &k21, &k22) * jacobian * Brem(brem_init, cobrems, E0, Egamma);
        if (cross_section > cross_max){
            bad_max = bad_max + 1;//an occurrence of cross section larger than cross_max, not supposed to happen
            std::cout <<  "bad max cross section= " <<  cross_section << "\n";
        }
        cross_test = cross_max * RandomReal(zlo, zhi);
        if (cross_test > cross_section){//selection fails
            nfail = nfail + 1;
            goto g100;
        }
        k10in = k10;
        k11in = k11;
        k12in = k12;
        k20in = k20;
        k21in = k21;
        k22in = k22;

        analysis(Egamma, mtgt, missing_mass, m_part, pion_hypothesis, &w_mumu, &t, k10in, k11in, k12in, k20in, k21in, k22in);//analyze the event;
        
        if(w_cut && (w_mumu < w_min || w_mumu > w_max)){
            goto g100;
        }

//  
//  ******************************Event selection succeeds:*********************************************
//		Do all the histogramming
//      
        wi_array = int(w_mumu/delta_w) + 1;//w distribution;
        //std::cout << "w_mumu " << w_mumu << "\n";
        //std::cout << "delta_w " << delta_w << "\n";
        std::cout << "wi_array " << wi_array << "\n";
        if(wi_array < 0){
            wi_array = 0;
        }
        if(wi_array > 200){
            wi_array = 200;
        }
        data_array_w[wi_array] = data_array[wi_array] + 1; 
        //data_array_w[i_array] = data_array_w[i_array] + 1;

        xi_array = x/delta_x; //x distribution
        if(xi_array < 0){
            xi_array = 0;
        }
        if(xi_array > 200){
            xi_array = 200;
        }
        //data_array_x[i_array] = data_array_x[i_array] + 1;
        data_array_x[i] = xi_array;

        ti_array = t/delta_t;//t distribution;
        if(ti_array < 0){
            ti_array = 0;
        }
        if(ti_array > 200){
            ti_array = 200;
        }
        //data_array_t[i_array] = data_array_t[i_array] + 1;
        data_array_t[i] = ti_array;

        phi_JTi_array = int(phi_JT*180/M_PI/delta_phi); //JT phi distribution in degrees;
        if(phi_JTi_array < 0){
            phi_JTi_array = 0;
        }
        if(phi_JTi_array > 200){
            phi_JTi_array = 200;
        }
        //data_array_phi_JT[i_array] = data_array_phi_JT[i_array] + 1;
        data_array_phi_JT[i] = phi_JTi_array;

        nonlinear_ti_array = (pow((t + beta_t), (1 + gamma_t)) - pow(beta_t, (1 - gamma_t)))/(alpha_t * (1 - gamma_t));  // variable t-bin width
	    if(nonlinear_ti_array < 0){
            nonlinear_ti_array = 0;
        }
        if(nonlinear_ti_array > 200){
            nonlinear_ti_array = 200;
        }	
	    //data_array_nonlinear_t[i_array] = data_array_nonlinear_t[i_array] + 1;
        data_array_nonlinear_t[i] = nonlinear_ti_array;

        Egammai_array = (Egamma - E_lo)/delta_Egamma;//photon energy distribution;
        if(Egammai_array < 0){
            Egammai_array = 0;
        }
        if(Egammai_array > 200){
            Egammai_array = 200;
        }
        //data_array_Egamma[i_array] = data_array_Egamma[i_array] + 1;
        data_array_Egamma[i] = Egammai_array;

        cos_maxi_array = (cos_max - cos(theta1))/delta_Egamma;
        if(cos_maxi_array < 0){
            cos_maxi_array = 0;
        }
        if(cos_maxi_array > 200){
            cos_maxi_array = 200;
        }
        //data_array_cos[i_array] = data_array_cos[i_array] + 1;
        data_array_cos[i] = cos_maxi_array;
//	3-momentum event output
        if(output_event){
            std::ofstream outputFile;
            outputFile.open("lepton_v17_4_event.txt");
            outputFile << Egamma << " " << k10 << " " << k11 << " " << k12 << " " << k20 << " " << k21 << " " << k22 << " " << "\n"; //<< ktgt[0] << " " << ktgt[1]<< " " << ktgt[2] << " " << "\n";
        }
        if(verbose_output && i % 50 == 0){
            std::cout << " event # " << i << "\n";
        }
        i++;
    }
//   Event generation ends
//	
    float failure = float(nfail)/float(nevent);
    if(verbose_output){
        std::cout << "Phase space parameter = " << phase_space <<  "Failures per event = " << failure << " Events with cross section exceeding max xsctn = " << bad_max << "\n";
    }
    int i = 0;
    if(hist_w){
        i = 0;
        std::ofstream w_histFile;
        w_histFile.open("lepton_w.txt");
        while(i < 200){
            x_value = w_min + (float(i) - 0.5) * delta_w;
            error_array[i] = sqrt(data_array_w[i]);
            //std::cout << "w " << data_array_w[i] << "\n";
            w_histFile << x_value << " " << data_array_w[i] << " " << error_array[i] << "\n";
            i++;
        }
        w_histFile.close();
    }
    if (hist_x){
        i = 0;
        std::ofstream x_histFile;
        x_histFile.open("lepton_x.txt");
        while(i < 200){
            x_value = (float(i)-0.5) * delta_x;
            error_array[i] = sqrt(data_array_x[i]);
            //std::cout << "x " << data_array_x[i] << "\n";
            x_histFile << x_value << " " << data_array_x[i] << " " << error_array[i] << "\n";
            i++;
        }
        x_histFile.close();
    }
    if(hist_t){
        i = 0;
        std::ofstream t_histFile;
        t_histFile.open("lepton_t.txt");
        while(i < 200){
            x_value = (float(i) - 0.5) * delta_t;
            error_array[i] = sqrt(data_array_t[i]);
            //std::cout << "t " << data_array_t[i] << "\n";
            t_histFile << x_value << " " << data_array_t[i] << " " << error_array[i] << "\n";
            i++;
        }
        t_histFile.close();
    }

    if(hist_phi_JT){
        i = 0;
        std::ofstream phi_JT_histFile;
        phi_JT_histFile.open("lepton_phi_JT.txt");
        while(i < 200){
            x_value = (float(i) - 0.5) * delta_phi;
            error_array[i] = sqrt(data_array_phi_JT[i]);
            //std::cout << "phi_JT " << data_array_phi_JT[i] << "\n";
            phi_JT_histFile << x_value << " " << data_array_phi_JT[i] << " " << error_array[i] << "\n";
            i++;
        }
        phi_JT_histFile.close();
    }
    if(hist_nonlinear_t == true){
        i = 0;
        std::ofstream nonlinear_t_histFile;
        nonlinear_t_histFile.open("lepton_nonlinear_t.txt");
        while(i < 200){
            x_value_lo = pow((alpha_t * (1 - gamma_t) * float(i-1) + pow(beta_t, 1 - gamma_t)), 1/(1-gamma_t)) - beta_t;
            x_value = pow((alpha_t * (1 - gamma_t) * float(i - 0.5) + pow(beta_t, 1 - gamma_t)), 1/(1-gamma_t)) - beta_t;
            x_value_hi = pow((alpha_t * (1 - gamma_t) * float(i) + pow(beta_t, 1 - gamma_t)), 1/(1-gamma_t)) - beta_t;
            bin_width = x_value_hi - x_value_lo;
            delta_x_value = bin_width/2;
            error_array[i] = sqrt(data_array_nonlinear_t[i])/bin_width;
            data_array_nonlinear_t[i] = data_array_nonlinear_t[i]/bin_width;
            //std::cout << "nonlin t " << data_array_nonlinear_t[i] << "\n";
            nonlinear_t_histFile << x_value << " " << delta_x_value << " "<< data_array_nonlinear_t[i] << " " << error_array[i] << "\n";
            i++;
        }
        nonlinear_t_histFile.close();
    }
    if(hist_Egamma == true){
        i = 0;
        std::ofstream Egamma_histFile;
        Egamma_histFile.open("lepton_Egamma.txt");
        while(i < 200){
            x_value = E_lo + (float(i) - 0.5) * delta_Egamma;
            error_array[i] = sqrt(data_array_Egamma[i]);
            //std::cout << "Egamma " << data_array_Egamma[i] << "\n";
            Egamma_histFile << x_value << " " << data_array_Egamma[i] << " " << error_array[i] << "\n";
            i++;
        }
        Egamma_histFile.close();
    }
    if(hist_cos == true){
        i = 0;
        std::ofstream cos_histFile;
        cos_histFile.open("lepton_cos.txt");
        while(i < 200){
            x_value = cos_max - (float(i) - 0.5) * delta_cos;
            error_array[i] = sqrt(data_array_cos[i]);
            //std::cout << "cos " << data_array_cos[i] << "\n";
            cos_histFile << x_value << " " << data_array_cos[i] << " " << error_array[i];
            i++;
        }
        cos_histFile.close();
    }
    return 0;
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
double xsctn(double E0, int ztgt, double x, double theta1, double phi1, double theta2, double phi2, double pol, double m_part, bool nuc_FF, double* phi_JT, double* k10, double* k11, double* k12, double* k20, double* k21, double* k22)//units of nb/sr^2;
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
    //
    Z = float(ztgt);
    E1 = E0 * x;
    E2 = E0 * (1 - x);
    k1_mag = sqrt(pow(E1, 2) - pow(m_part, 2));
    k2_mag = sqrt(pow(E2, 2) - pow(m_part, 2));
    //std::cout << "k1_mag " << k1_mag << " k2_mag " << k2_mag << "\n";
    *k10 = k1_mag * sin(theta1) * cos(phi1);
    *k11 = k1_mag * sin(theta1) * sin(phi1);
    *k12 = k1_mag * cos(theta1);
    //std::cout << "xsctn k1[0] " << *k10 << " xsctn k1[1] " << *k11 << " xsctn k1[2] " << *k12 <<"\n";
    *k20 = k2_mag * sin(theta2) * cos(phi2);
    *k21 = k2_mag * sin(theta2) * sin(phi2);
    *k22 = k2_mag * cos(theta2);
    //std::cout << "xsctn k2[0] " << *k20 << " xsctn k2[1] " << *k21 << " xsctn k2[2] " << *k22 <<"\n";

    p1 = sqrt(pow(*k11, 2) + pow(*k11,2));//transverse momenta of muon #1, GeV
    //std::cout << "p1 " << p1 <<"\n";
    p2 = sqrt(pow(*k20, 2) + pow(*k21,2));//transverse momenta of muon #2, GeV
    //std::cout << "p2 " << p2 <<"\n";
    q2_T = pow((*k10 + *k20),2) + pow((*k11 + *k21),2); // this is transverse momentum transfer squared
    //std::cout << "q2_T " << q2_T <<"\n";
    q2 = q2_T + pow((E0 - *k12 - *k22), 2); // this is 3 - momentum transfer squared
    //std::cout << "q2 " << q2 <<"\n";
    c1 = pow(p1, 2) + pow(m_part, 2);
    //std::cout << "c1 " << c1 <<"\n";
    c2 = pow(p2, 2) + pow(m_part, 2);
    //std::cout << "c2 " << c2 << "\n";
    JS = 1/c1 - 1/c2;//units of 1/GeV^2//scalar current, units of GeV^ - 2;
    //std::cout << "JS " << JS <<"\n";
    i = 0;
    JT[0] = *k10/c1 + *k20/c2;//vector current, units of GeV^ - 1;
    JT[1] = *k11/c1 + *k21/c2;
    JT[2] = *k12/c1 + *k22/c2;


    *phi_JT = acos(JT[0]/sqrt(pow(JT[0],2) + pow(JT[1], 2)));//phi angle of JT wrt to x axis, radians
    if (JT[1] < 0) *phi_JT = 2 * M_PI - *phi_JT;
    //std::cout << "phi_JT " << phi_JT <<"\n";
//
    W_unpol = pow(m_part,2) * pow(JS,2) + (pow(x,2) + pow((1 -x),2)) * (pow(JT[1],2) + pow(JT[2], 2));
    //std::cout << "W_unpol " << W_unpol <<"\n";
//     	W_pol = -2.*x*(1.-x)*((p2/c2)**2*cos(2.*phi2)+(p1/c1)**2*cos(2.*phi1)+2.*(p2/c2)*(p1/c1)*cos(phi1+phi2)) !the Bakmaev expression
//	xsctn=2.*alpha**3*Z**2*E0**4*x**2*(1.-x)**2/(pi**2*q2_T**2)*(W_unpol+pol*W_pol) ! note the absence of cos(2phi_JT) in this expression
//     &	*hbarc**2/100.*1.e9*FF2(q2_T,ztgt,nuc_FF) !units of nb/sr^2  The denominator uses the transverse 3-momentum transfer^2, 
    W_pol = -2 * x * (1 - x) * (pow(JT[1], 2) + pow(JT[2], 2));//this is my reduction of the Bakmaev equations;
    //std::cout << "W_pol " << W_pol <<"\n";
    xsctnOut = 2 * pow(alpha,3) * pow(Z, 2) * pow(E0, 4) * pow(x, 2) * pow((1 - x), 2)/(pow(M_PI,2) * pow(q2_T, 2)) * (W_unpol + pol * cos(2 * *phi_JT) * W_pol) * pow(hbarc, 2)/100 * 1e9 * FF2(q2_T, ztgt, nuc_FF);
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
double analysis(double E0, double mtgt, double missing_mass, double m_part, bool pion_hypothesis, double* w_mumu, double* t, double k10in, double k11in, double k12in, double k20in, double k21in, double k22in)
{
    // implicit none
    double E1, E2, ks[3], m_x, analysisOut, ktgt[3];
    double mass_pi = 0.139570;
    //std::cout << "analysis k1[0] " << k10in <<"\n";
    //std::cout << "analysis k1[1] " << k11in <<"\n";
    //std::cout << "analysis k1[2] " << k12in <<"\n";
    //std::cout << "analysis k2[0] " << k20in <<"\n";
    //std::cout << "analysis k2[1] " << k21in <<"\n";
    //std::cout << "analysis k2[2] " << k22in <<"\n";
    E1 = sqrt(pow(k10in, 2) + pow(k11in, 2) + pow(k12in, 2) + pow(m_part, 2));//lepton energies
    E2 = sqrt(pow(k20in, 2) + pow(k21in, 2) + pow(k22in, 2) + pow(m_part, 2));
    ks[0] = k10in + k20in;//lepton summed momentum
    ks[1] = k11in + k21in;
    ks[2] = k12in + k22in;
    ktgt[0] = -ks[0];//target momentum;
    ktgt[1] = -ks[1];
    ktgt[2] = E0 - ks[2];
    missing_mass = sqrt(pow((E0 + mtgt - E1 - E2), 2) - pow(ktgt[0], 2) - pow(ktgt[1],2) - pow(ktgt[2], 2));
    *t = pow(ks[0], 2) + pow(ks[1], 2) + pow((E0 - ks[2]), 2) - pow((E0 - E1 - E2), 2);//4 - momentum transfer squared to nucleus, this is positive;
    //
    //		mu mu invariant mass, possibly with pion hypothesis
    m_x = m_part;
    if (pion_hypothesis){
        m_x = mass_pi;
    }
    E1 = sqrt(pow(k10in, 2) + pow(k11in,2) + pow(k12in, 2) + pow(m_x, 2));//need to put in the mass hypothesis
    E2 = sqrt(pow(k20in, 2) + pow(k21in, 2) + pow(k22in, 2) + pow(m_x, 2));
    //std::cout << "E1 " << E1 << " E2 " << E2 <<  " ks[0] " << ks[0] << " ks[1] " << ks[1] << " ks[2] "<< ks[2] << "\n";
    *w_mumu = sqrt(pow((E1 + E2), 2) - pow(ks[0], 2) - pow(ks[1], 2) - pow(ks[2],2));
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
    //c
    rho0 = 0;
    if((ztgt == 82) || (ztgt == 1)) return;
    //c	These equations have to do with Fermi distribution, reference? 
    w = 4 * M_PI * c_den[ztgt]/3 * (pow((M_PI * a_den[ztgt]), 2) + pow(c_den[ztgt], 2));
    i = 0;
    while(i < 10){
        w = w + 8 * M_PI * pow(a_den[ztgt], 3) * pow((-1), (i - 1)) * exp(-float(i)*c_den[ztgt]/a_den[ztgt])/pow(float(i), 3);
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
    double q02, hbarc, Q, gamma, r[12], A[12], rho0, c_den[100], a_den[100], norm, proton_rms;
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
      returnFF = 4 * pow(M_PI,2) * rho0 * pow(a_den[ztgt], 3)/(pow((Q * a_den[ztgt]), 2) * pow((sinh(M_PI * Q * a_den[ztgt])),2)) * (M_PI * Q * a_den[ztgt] * cosh(M_PI * Q * a_den[ztgt]) * sin(Q * c_den[ztgt]) - Q * c_den[ztgt] * cos(Q * c_den[ztgt]) * sinh(M_PI * Q * a_den[ztgt]));
      i = 0;
      do{
	    returnFF = returnFF + 8 * M_PI * rho0 * pow(a_den[ztgt], 3) * pow((-1), (i - 1)) * float(i) * exp(-float(i) * c_den[ztgt]/a_den[ztgt])/pow((pow(float(i), 2) + pow((Q * a_den[ztgt]),2)),2);
	    i++;
      }while(i < 10);
    }
return returnFF;
}
