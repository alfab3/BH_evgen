//Date: 6/10/2020 12:06 PM
//	To compile: "g++ lepton_event_v17_5.cpp -lgsl -lgslcblas”
//	To run: “./a.out”
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <sys/time.h>


struct density_rho0{
  double c_den[100], a_den[100]; //the overall normalization factor for nuclear charge densities
};


struct Brem_spect{
  double Eg[500], Br[500];
};

std::vector<int> E_lookup = {2,6,10,14,18,22,26,30,34,38,42,46,50,54,58,62,66,70,74,78,82,86,90,94,98,102,106,110,114,118,122,126,130,134,138,142,146,150,154,158,162,166,170,174,178,182,186,190,194,198,202,206,210,214,218,222,226,230,234,238,242,246,250,254,258,262,266,270,274,278,282,286,290,294,298,302,306,310,314,318,322,326,330,334,338,342,346,350,354,358,362,366,370,374,378,382,386,390,394,398,402,406,410,414,418,422,426,430,434,438,442,446,450,454,458,462,466,470,474,478,482,486,490,494,498,502,506,510,514,518,522,526,530,534,538,542,546,550,554,558,562,566,570,574,578,582,586,590,594,598,602,606,610,614,618,622,626,630,634,638,642,646,650,654,658,662,666,670,674,678,682,686,690,694,698,702,706,710,714,718,722,726,730,734,738,742,746,750,754,758,762,766,770,774,778,782,786,790,794,798,802,806,810,814,818,822,826,830,834,838,842,846,850,854,858,862,866,870,874,878,882,886,890,894,898,902,906,910,914,918,922,926,930,934,938,942,946,950,954,958,962,966,970,974,978,982,986,990,994,998,1002,1006,1010,1014,1018,1022,1026,1030,1034,1038,1042,1046,1050,1054,1058,1062,1066,1070,1074,1078,1082,1086,1090,1094,1098,1102,1106,1110,1114,1118,1122,1126,1130,1134,1138,1142,1146,1150,1154,1158,1162,1166,1170,1174,1178,1182,1186,1190,1194,1198};
std::vector<int> EnergyCounts = {0,0,0,0,0,0,0,0,0,0,0,0,811,1699,1770,1852,1917,1777,1809,1818,1777,1687,1623,1559,1440,1390,1375,1272,1168,1072,1154,1016,969,889,926,907,844,861,771,765,747,731,707,648,673,700,702,639,627,559,577,523,622,542,550,530,516,519,515,484,474,460,431,480,433,423,420,413,411,408,394,395,405,378,353,369,383,338,337,350,376,348,361,313,332,304,314,274,308,261,292,305,290,277,297,272,282,259,264,274,233,220,243,259,271,258,212,234,243,223,217,245,208,194,242,203,185,227,178,229,200,190,194,218,184,190,224,193,176,168,194,173,203,191,179,201,163,185,181,172,189,168,158,166,188,179,168,151,142,155,161,170,157,149,151,155,167,151,185,158,168,160,145,161,159,155,155,168,152,152,170,148,148,180,166,188,165,181,167,189,168,178,175,180,179,196,183,201,203,211,222,223,239,224,258,244,269,273,290,284,284,339,298,334,360,347,360,383,378,432,459,466,472,447,512,508,566,558,610,585,445,96,120,117,120,103,124,114,115,119,120,118,135,109,144,139,155,141,158,144,133,159,154,168,171,179,194,163,175,201,202,104,109,119,123,105,115,131,108,102,127,121,120,109,105,92,107,109,94,82,84,79,89,106,86,72,109,104,94,115,142,138,155,148,109,108,113,95,105,99,0,0,0,0,0,0,0,0,0,0};

struct density_rho0 density_rho0;
struct Brem_spect brem_spect;

#define _ISEED_ 0

#define _OVERFLOW_PROTECTION_ 1


double Brem(bool brem_init, bool cobrems, double E0, double Egamma);
double xsctn(double E0, int ztgt, double x, double theta1, double phi1, double theta2, double phi2, double pol, double m_part, bool nuc_FF, double* phi_JT, double* k10, double* k11, double* k12, double* k20, double* k21, double* k22);//units of nb/sr^2
double FF2(double q2, int ztgt, bool nuc_FF);
double analysis(double E0, double mtgt, double missing_mass, double m_part, bool pion_hypothesis, double* w_mumu, double* t, double k10in, double k11in, double k12in, double k20in, double k21in, double k22in);
double density_init(int ztgt);
double FF(double Q2, int ztgt);
void    ZBQLINI(int seed, double ZBQLIX[43+1], double *bb);
double  ZBQLUAB(double x1, double x2);
double  ZBQLU01();



double E0, theta_min, theta_max;
double raddeg, degrad, bremMain, zlo, zhi, w_min, w_max, mass_pi, delta_cos, delta_x_value, bin_width, x_value_hi, x_value_lo, total_xscn_old, delta_total_xscn, pi, cos_max, cos_min, cross_sum, cross_max, cross_test, x_min, x_max, theta1, phi1, costheta1, theta2, phi2, costheta2, cross_section, total_xscn, m_e, m_part, m_muon, pol, x, q2, failure, w_mumu, xs1, xs2, xs_max, xs_min, Atgt[100], tgtlen[100], radlen[100], hbarc, W, jacobian, phi_JT, delta_w, delta_x, delta_phi, mtgt, ktgt[3], Elab[2], missing_mass, t, delta_t, x_value, q2_T, W_pol, W_unpol, JS, JT[2], Rexp, xmax, theta1_max, theta2_max, phi1_max, phi2_max, temp, delta_log_t, frac_delta_t, delta_Egamma, data_array[200], error_array[200], Egamma_max, E_hi, E_lo, E_coherent, Egamma, data_array_w[200], data_array_x[200], data_array_t[200], data_array_phi_JT[200], data_array_nonlinear_t[200], data_array_Egamma[200], data_array_cos[200], data_array_log_t[200], cross_max_old, delta_cross_max, k10, k11, k12, k20, k21, k22, k10in, k11in, k12in, k20in, k21in, k22in;
int iseed, i_array, imax;
int i, itest[4], nevent, j, nfail, bad_max, j_array, ztgt, phase_space, ipoint; 
bool hist_w, hist_x, hist_t, hist_phi_JT, hist_Egamma, output_event, hist_nonlinear_t, hist_log_t, muon, electron, pion_hypothesis, brem_init, cobrems, hist_cos, integral_xsctn, w_cut, nuc_FF, verbose_output;
int main(){
    //Standard CPP configuration 
    //ztgt = 82, E0 = 5.5, pol = 1.0, theta_min = 0.80,theta_max = 5.3; //Standard CPP configuration, min angle to TOF and max angle to MWPC
    //Standard GlueX config.
    ztgt = 1, E0 = 11.0, E_coherent = 8.7, pol = 1.0, theta_min = 0.75, theta_max = 13.12;//Standard GlueX config., min and max angles in deg. to TOF;
    //Set invariant mass cut
    w_min = 0.25, w_max = 0.621;
    //Set tagging interval
    E_hi = 8.8, E_lo = 8.6;
    itest[0] = 100000, itest[1] = 1000000, itest[2] = 10000000, itest[3] = 100000000, nevent = 1000000;
    m_e = 0.000511, m_muon = 0.105658, mass_pi = 0.139570, hbarc = 0.197;

    //Histogram parameters
    delta_w = 0.02, delta_t = 0.0002, delta_x = 0.02, delta_phi = 5.0, frac_delta_t = 0.2, delta_Egamma = 0.05, delta_cos = 0.0001; //units of GeV, GeV^2, //DIMENSIONless, degrees,

    //Target information
    Atgt[1] = 1.0,//tgtlen = # Rad lengths, radlen = rad length of material in g/cm^2;
    Atgt[6] = 12.0, density_rho0.c_den[6] = 2.45, density_rho0.a_den[6] = 0.524;//RL is in g/cm^2;
    Atgt[14] = 28.0, density_rho0.c_den[14] = 3.14, density_rho0.a_den[14] = 0.537;//units of c_den and a_den are fm;
    Atgt[20] = 40.0, density_rho0.c_den[20] = 3.51, density_rho0.a_den[20] = 0.563;//target length is in units of RL;
    Atgt[26] = 56.0, density_rho0.c_den[26] = 3.971,density_rho0.a_den[26] = 0.5935;//RL is g/cm^2;
    Atgt[50] = 116.0, density_rho0.c_den[50] = 5.416, density_rho0.a_den[50] = 0.552;
    Atgt[82] = 208.0;

    zlo = 0.0, zhi = 1.0;
//
//  Initializations
//
//Find the delta log t step
    delta_log_t = log10(1. + frac_delta_t);

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
    degrad = M_PI/180.0;
    raddeg = 180.0/M_PI;

//Start logical assignments

//Can turn on any or all of these logicals.
    hist_w = true;
    hist_x = true;
    hist_t = true;
    hist_phi_JT = true;
    hist_nonlinear_t = false; //not currently operational
    hist_log_t = true;
    hist_Egamma = false; 
    hist_cos = false;
    output_event = true;     //for writing event file
    integral_xsctn = true;  //set true for outputting integrated cross sections
    w_cut = false;           //set true for applying w cuts to the data
    nuc_FF = true;           //set true for using a nuclear form factor, set false for FF=1;
    cobrems = false;          //set true for scanfing coherent Brems. file, false for using a 1/Egamma distribution;
    muon = false;            //set true for muons,false for electrons
    verbose_output = true;  //set true if you want to see test outputs printed to your screen
    electron = !muon;

//Initialize Brem. distribution: select 1/Egamma or coherent Brems. file
    
    if(cobrems == true){ // scanf coherent Brem. file
      brem_init = true;
      temp = Brem(brem_init, cobrems, E0, Egamma);//scanf coherent brems file, then set brem_init = false;
    }

    FILE *CBD;
    if (brem_init == true){//fopen and scanf coherent Brem file
        CBD = fopen("CobremsDistribution.dat", "r");

        while(i < 300){
	        fscanf(CBD,"%lf %lf", &brem_spect.Eg[i], &brem_spect.Br[i]);
	        i++;
        }
        imax = i - 1;
        fclose(CBD);
        brem_init = false;//done with initialization;
    }


    phase_space = 8;//theta = x**phase_space, with int phase_space >1 . Note: phase_space = 8 is the fastest for e+e-;
    //Set phase_space=0 for standard dcos theta/dx =1
    Rexp = phase_space;

    if(muon){
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
    theta1 = theta_min * (degrad);
    theta2 = theta_min * (degrad);
    phi1 = 90 * (degrad);
    phi2 = 270 * (degrad);
    cross_section = xsctn(E_coherent, ztgt, x, theta1, phi1, theta2, phi2, pol, m_part, nuc_FF, &phi_JT, &k10, &k11, &k12, &k20, &k21, &k22);
    if(verbose_output == true){
        std::cout << E_coherent << " " << ztgt << " " << x << " " << theta1 << " " << phi1 << " " << theta2 << " " << phi2 << " " << pol << " " << m_part << " " << nuc_FF << " " << "\n";
        std::cout << " cross section nb/sr^2 = " << cross_section << "\n";
    }
//*************************************
    theta_min = theta_min * degrad;//switch to radians;
    theta_max = theta_max * degrad;
    cos_max = cos(theta_min);
    cos_min = cos(theta_max);

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
    g25:
        i = 0;
        while(i < itest[j]){//find maximum cross section in allowed phase space
            Egamma = ZBQLUAB(E_lo, E_hi);//get tagged photon energy
            x = 0.5;//x_min + (x_max - x_min) ZBQLUAB(zlo, zhi)//make a guess for the energy fraction
            phi1 = 90 * degrad;//2.*pi ZBQLUAB(zlo, zhi)//make a guess for phi1
            phi2 = 270 * degrad;//2.*pi ZBQLUAB(zlo, zhi)//make a guess for phi2
            xs1 = ZBQLUAB(xs_min, xs_max);
            xs2 = ZBQLUAB(xs_min, xs_max);
            if (phase_space == 0) {// dcos theta/dx = 1
                theta1 = acos(xs1);//make a guess for theta1
                theta2 = acos(xs2);
                jacobian = 1;
            }else{
                theta1 = pow(xs1, phase_space);//make a guess for theta1
                theta2 = pow(xs2, phase_space);
                jacobian = Rexp * pow(xs1, phase_space - 1) * sin(pow(xs1, phase_space)) * Rexp * pow(xs2,phase_space - 1) * sin(pow(xs2,phase_space));
            }
            cross_section = xsctn(Egamma, ztgt, x, theta1, phi1, theta2, phi2, pol, m_part, nuc_FF, &phi_JT, &k10, &k11, &k12, &k20, &k21, &k22) * jacobian * Brem(brem_init, cobrems, E0, Egamma);

            //analysis(Egamma, mtgt, missing_mass, m_part, pion_hypothesis, &w_mumu, &t, k10in, k11in, k12in, k20in, k21in, k22in);
            //if(w_cut && (w_mumu < w_min || w_mumu > w_max)){
            //    goto g25;
            //}

            if(cross_section > cross_max){
                cross_max = cross_section;
                Egamma_max = Egamma;
                xmax = x;
                theta1_max = theta1 * raddeg;
                theta2_max = theta2 * raddeg;
                phi1_max = phi1 * raddeg;
                phi2_max = phi2 * raddeg;

            }
            i++;
        }
        if(verbose_output == true){
            std::cout << "test events " << itest[j] <<  " maximum xsctn*Brem " <<  cross_max << "\n";
            std::cout << "Egamma max " << Egamma_max << " x max " << xmax << " theta1 max " << theta1_max << " theta2 max " << theta2_max << " phi1 max " << phi1_max << " phi2 max " <<  phi2_max << "\n";
        }
        j++;
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
                x = x_min + (x_max - x_min) * ZBQLUAB(zlo, zhi);//energy fraction
                phi1 = 2 * M_PI * ZBQLUAB(zlo, zhi);
                phi2 = 2 * M_PI * ZBQLUAB(zlo, zhi);
                xs1 = ZBQLUAB(xs_min, xs_max);
                xs2 = ZBQLUAB(xs_min, xs_max);
                if (phase_space == 0){// dcos theta/dx = 1
                    theta1 = acos(xs1);
                    theta2 = acos(xs2);
                    jacobian = 1;
                }else{
                    theta1 = pow(xs1,phase_space);
                    theta2 = pow(xs2,phase_space);
                    jacobian = Rexp * pow(xs1,phase_space - 1) * sin(pow(xs1,phase_space)) * Rexp * pow(xs2,phase_space-1)*sin(pow(xs2,phase_space));
                }
                cross_section = xsctn(E_coherent, ztgt, x, theta1, phi1, theta2, phi2, pol, m_part, nuc_FF, &phi_JT, &k10, &k11, &k12, &k20, &k21, &k22)*jacobian;
                //k10in = k10;
                //k11in = k11;
                //k12in = k12;
                //k20in = k20;
                //k21in = k21;
                //k22in = k22;
                //analysis(Egamma, mtgt, missing_mass, m_e, pion_hypothesis , &w_mumu, &t, k10in, k11in, k12in, k20in, k21in, k22in);
                cross_sum = cross_sum + cross_section;
            i++;
	    }
        total_xscn = cross_sum/double(itest[j])*pow((xs_max - xs_min), 2)*pow((2*M_PI), 2) * (x_max - x_min);
        //delta_total_xscn = (total_xscn-total_xscn_old)/total_xscn*100;
        std::cout << "test events " << itest[j] << " Egamma " <<  E_coherent << " total cross section nb " <<  total_xscn << "\n";
        j++;
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
    int j;
    double rrz;
    std::ofstream outputFile;
    outputFile.open("lepton_event.txt");
    while(i < nevent){

    g100:

        Egamma = ZBQLUAB(E_lo, E_hi);//get tagged photon energy
        x = ZBQLUAB(x_min, x_max);//energy fraction
        //Test x to make sure it's within the allowed range for the photon energy Egamma
        if(x >= (Egamma - m_part)/Egamma || x <= m_part/Egamma) continue;

        phi1 = 2 * M_PI * ZBQLUAB(zlo, zhi);
        phi2 = 2 * M_PI * ZBQLUAB(zlo, zhi);
        xs1 = ZBQLUAB(xs_min, xs_max);
        xs2 = ZBQLUAB(xs_min, xs_max);

        if(phase_space == 0){// dcos theta/dx = 1
            theta1 = acos(xs1);
            theta2 = acos(xs2);
            jacobian = 1;
        }else{
            theta1 = pow(xs1, phase_space);
            theta2 = pow(xs2, phase_space);
            jacobian = Rexp * pow(xs1,phase_space-1) * sin(pow(xs1,phase_space)) * Rexp*pow(xs2,phase_space-1) * sin(pow(xs2,phase_space));
        }

        cross_section = xsctn(Egamma, ztgt, x, theta1, phi1, theta2, phi2, pol, m_part, nuc_FF, &phi_JT, &k10, &k11, &k12, &k20, &k21, &k22) * jacobian * Brem(brem_init, cobrems, E0, Egamma);

        if (cross_section > cross_max){
            bad_max = bad_max + 1;//an occurrence of cross section larger than cross_max, not supposed to happen
            std::cout <<  "bad max cross section= " <<  cross_section << "\n";
        }
        rrz = ZBQLUAB(zlo,zhi);
        cross_test = cross_max * rrz;
        
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
        //std::cout << w_mumu << "\n";
        //if(w_cut && (w_mumu < w_min || w_mumu > w_max)){
        //    goto g100;
        //}

//  
//  ******************************Event selection succeeds:*********************************************
//		Do all the histogramming
//      
        if(hist_w) i_array = int(w_mumu/delta_w) + 1;//w distribution;
        if(i_array < 0) i_array = 0;
        if(i_array > 200) i_array = 200;
        data_array_w[i_array] = data_array_w[i_array] + 1;

        if(hist_x) i_array = int(x/delta_x); //x distribution
        if(i_array < 0) i_array = 0;
        if(i_array > 200) i_array = 200;
        data_array_x[i_array] = data_array_x[i_array] + 1;

        if(hist_t) i_array = int(t/delta_t);//t distribution;
        if(i_array < 0) i_array = 0;
        if(i_array > 200) i_array = 200;
        data_array_t[i_array] = data_array_t[i_array] + 1;

        if(hist_phi_JT) i_array = int(phi_JT*180/M_PI/delta_phi); //JT phi distribution in degrees;
        if(i_array < 0) i_array = 0;
        if(i_array > 200) i_array = 200;
        data_array_phi_JT[i_array] = data_array_phi_JT[i_array] + 1;

        /*
        if(hist_nonlinear_t) i_array = int((pow((t + beta_t), (1 + gamma_t)) - pow(beta_t, (1 - gamma_t)))/(alpha_t * (1 - gamma_t)));  // variable t-bin width
        std::cout << "nonlint " << i_array << "\n";
        if(i_array < 0) i_array = 0;
        if(i_array > 200) i_array = 200;
	    data_array_nonlinear_t[i_array] = data_array_nonlinear_t[i_array] + 1;
        /**/

        if(hist_log_t) i_array = int((log10(t) + 6.0)/delta_log_t) + 1;
        if(i_array < 0) i_array = 0;
        if(i_array > 200) i_array = 200;
        data_array_log_t[i_array] = data_array_log_t[i_array] + 1;


        if(hist_Egamma) i_array = int((Egamma/delta_Egamma)) + 1;//photon energy distribution;
        if(i_array < 0) i_array = 0;
        if(i_array > 200) i_array = 200;
        std::cout << "Egamma " << i_array << "\n";
        data_array_Egamma[i_array] = data_array_Egamma[i_array] + 1;

        if(hist_cos) i_array = int((cos_max - cos(theta1))/delta_Egamma);
        if(i_array < 0) i_array = 0;
        if(i_array > 200) i_array = 200;
        std::cout << "cosmax " << i_array << "\n";
        data_array_cos[i_array] = data_array_cos[i_array] + 1;


        if(output_event){
            outputFile << Egamma << " " << k10 << " " << k11 << " " << k12 << " " << k20 << " " << k21 << " " << k22 << " " << " " << ktgt[0] << " " << ktgt[1]<< " " << ktgt[2] << " " << "\n";
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
    if(hist_w){
        i = 0;
        std::ofstream w_histFile;
        w_histFile.open("lepton_w.txt");
        while(i < 200){
            x_value = w_min + (float(i) - 0.5) * delta_w;
            error_array[i] = sqrt(data_array_w[i]);
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
            phi_JT_histFile << x_value << " " << data_array_phi_JT[i] << " " << error_array[i] << "\n";
            i++;
        }
        phi_JT_histFile.close();
    }

    if(hist_log_t){
        i = 0;
        std::ofstream log_t_histFile;
        log_t_histFile.open("lepton_log_t.txt");
        while(i < 200){
            x_value = pow(10.0, float(i)*delta_log_t - 6.0);
            error_array[i] = sqrt(data_array_log_t[i])/(pow(10.0, (i+1) * delta_log_t - 6.0) - pow(10.0, (i-1) * delta_log_t - 6.0)) * 2;
            data_array_log_t[i] = data_array_log_t[i]/(pow(10.0, (i+1) * delta_log_t - 6.0) - pow(10.0, (i-1) * delta_log_t - 6.0)) * 2;
            log_t_histFile << x_value << " " << data_array_log_t[i] << " " << error_array[i] << "\n";
            i++;
        }
    }

    /*
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
            nonlinear_t_histFile << x_value << " " << delta_x_value << " " << data_array_nonlinear_t[i] << " " << error_array[i] << "\n";
            i++;
        }
        nonlinear_t_histFile.close();
    }
    /**/

    if(hist_Egamma == true){
        i = 0;
        std::ofstream Egamma_histFile;
        Egamma_histFile.open("lepton_Egamma.txt");
        while(i < 200){
            x_value = E_lo + (float(i) - 0.5) * delta_Egamma;
            error_array[i] = sqrt(data_array_Egamma[i]);
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
            cos_histFile << x_value << " " << data_array_cos[i] << " " << error_array[i] << "\n";
            i++;
        }
        cos_histFile.close();
    }
    return 0;
}


double Brem(bool brem_init, bool cobrems, double E0, double Egamma){
    int i;
    double bremOut;
    int Ecompare = std::ceil(100*Egamma);
    /*
    if (brem_init && !cobrems){
        bremMain = E0/Egamma;
    }
    if (brem_init && cobrems){ //return coherent brems distribution
        ipoint = int((Egamma + .02)/.04);
        bremMain = brem_spect.Br[ipoint];
    }
    /**/
    i = 0;
    while(E_lookup[i] < Ecompare){
        i++;
    }
    bremOut = EnergyCounts[i];
    return bremOut;
}

//******************************************************************
// The reference for this is Bakmaev et al., Physics Letters B 660 (2008) 494-500, eqn. 23
// Had to correct a mistake in Eqn. 23 of their paper.   The cos(phi1 + phi2) term should be 
// multiplied by 2.  You can see this by comparing Wp in Eqn. 22 with the vector current part of Eqn. 23
//
// --------------------------------------------
double xsctn(double E0, int ztgt, double x, double theta1, double phi1, double theta2, double phi2, double pol, double m_part, bool nuc_FF, double* phi_JT, double* k10, double* k11, double* k12, double* k20, double* k21, double* k22)//units of nb/sr^2;
{
    double Z, W_unpol, W_pol, q2_T;
    double alpha, hbarc;
    double xsctnOut, xsctn_point;
    double E1, E2, k1_mag, k2_mag, p1, p2, q2, c1, c2, JS, JT[2], FF_nuc, FF_TFM, FF2G;
    int i;
    alpha = 7.297352e-3, hbarc = 0.197;
    Z = float(ztgt);
    E1 = E0 * x;
    E2 = E0 * (1 - x);
    k1_mag = sqrt(pow(E1, 2) - pow(m_part, 2));
    k2_mag = sqrt(pow(E2, 2) - pow(m_part, 2));
    *k10 = k1_mag * sin(theta1) * cos(phi1);
    *k11 = k1_mag * sin(theta1) * sin(phi1);
    *k12 = k1_mag * cos(theta1);
    *k20 = k2_mag * sin(theta2) * cos(phi2);
    *k21 = k2_mag * sin(theta2) * sin(phi2);
    *k22 = k2_mag * cos(theta2);

    p1 = sqrt(pow(*k10, 2) + pow(*k11,2));//transverse momenta of muon #1, GeV
    p2 = sqrt(pow(*k20, 2) + pow(*k21,2));//transverse momenta of muon #2, GeV
    q2_T = pow((*k10 + *k20),2) + pow((*k11 + *k21),2); // this is transverse momentum transfer squared
    q2 = q2_T + pow((E0 - *k12 - *k22), 2); // this is 3 - momentum transfer squared
    c1 = pow(p1, 2) + pow(m_part, 2);
    c2 = pow(p2, 2) + pow(m_part, 2);
    JS = 1./c1 - 1./c2;//units of 1/GeV^2//scalar current, units of GeV^ - 2;
    JT[0] = *k10/c1 + *k20/c2;//vector current, units of GeV^ - 1;
    JT[1] = *k11/c1 + *k21/c2;
    //JT[2] = *k12/c1 + *k22/c2;


    *phi_JT = acos(JT[0]/sqrt(pow(JT[0],2) + pow(JT[1], 2)));//phi angle of JT wrt to x axis, radians
    if (JT[1] < 0) *phi_JT = 2 * M_PI - *phi_JT;

    W_unpol = pow(m_part,2) * pow(JS,2) + (pow(x,2) + pow((1 -x),2)) * (pow(JT[0],2) + pow(JT[1], 2));

//  W_pol = -2.*x*(1.-x)*((p2/c2)**2*cos(2.*phi2)+(p1/c1)**2*cos(2.*phi1)+2.*(p2/c2)*(p1/c1)*cos(phi1+phi2)) !the Bakmaev expression
//	xsctn=2.*alpha**3*Z**2*E0**4*x**2*(1.-x)**2/(pi**2*q2_T**2)*(W_unpol+pol*W_pol) ! note the absence of cos(2phi_JT) in this expression
//     &	*hbarc**2/100.*1.e9*FF2(q2_T,ztgt,nuc_FF) !units of nb/sr^2  The denominator uses the transverse 3-momentum transfer^2, 
    W_pol = -2 * x * (1 - x) * (pow(JT[0], 2) + pow(JT[1], 2));//this is my reduction of the Bakmaev equations;
    xsctn_point = 2 * pow(alpha,3) * pow(Z, 2) * pow(E0, 4) * pow(x, 2) * pow((1 - x), 2)/(pow(M_PI,2) * pow(q2_T, 2)) * (W_unpol + pol * cos(2 * *phi_JT) * W_pol) * pow(hbarc, 2)/100 * 1e9;
    FF2G = FF2(q2_T,ztgt,nuc_FF);
    xsctnOut = xsctn_point * FF2G;
    //std::cout << "xsctn " << xsctnOut << "\n";
//this contains the cos(2phi_JT) term*hbarc**2/100.*1.e9*FF2(q2_T, ztgt, nuc_FF) //units of nb/sr^2 The denominator uses the transverse 3 - momentum transfer^2
    return xsctnOut;
}


double FF2(double q2, int ztgt, bool nuc_FF)
{
    float hbarc, z, FF_nuc, FF_TFM, alpha[3], b[3], b0, m_e, c, outputFF2;
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
    if(nuc_FF){
        FF_nuc = FF(q2, ztgt);
    }else{
        FF_nuc = 1;
    }
    i = 0;
    FF_TFM = 1;
    while(i < 3){
      FF_TFM = FF_TFM - alpha[i]*q2/(q2 + pow((b[i]*c),2));
      i++;
    }
    outputFF2 = pow((FF_nuc - FF_TFM), 2);
    return outputFF2;
}

double analysis(double E0, double mtgt, double missing_mass, double m_part, bool pion_hypothesis, double* w_mumu, double* t, double k10in, double k11in, double k12in, double k20in, double k21in, double k22in)
{
    double E1, E2, ks[3], m_x, analysisOut, ktgt[3];
    double mass_pi = 0.139570;
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
    *w_mumu = sqrt(pow((E1 + E2), 2) - pow(ks[0], 2) - pow(ks[1], 2) - pow(ks[2],2));
}


double density_init(int ztgt)
{
    double pi, c_den[100], a_den[100],rho0, w;
    int i;
    rho0 = 0;
    if((ztgt == 82) || (ztgt == 1)){
        return rho0;
    }
    //	These equations have to do with Fermi distribution, reference? 
    w = 4 * M_PI * c_den[ztgt]/3 * (pow((M_PI * a_den[ztgt]), 2) + pow(c_den[ztgt], 2));
    i = 0;
    while(i < 10){
        w = w + 8 * M_PI * pow(a_den[ztgt], 3) * pow((-1), (i - 1)) * exp(-float(i)*c_den[ztgt]/a_den[ztgt])/pow(float(i), 3);
        i++;
    }
    rho0 = 1/w;
    return rho0;
}



double FF(double Q2, int ztgt)
{
    double q02, hbarc, Q, gamma, r[12], A[12], rho0, c_den[100], a_den[100], norm, proton_rms;
    int i = 0;
    double returnFF, returnFF_1, returnFF_2, returnFF_3;
    q02 = 0.71, proton_rms = 0.879;//proton dipole form factor parameter GeV^2, proton rms radius fm
    hbarc = 0.197;
    rho0 = density_init(ztgt);
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
        const double q02 = 0.71, proton_rms = 0.879;
        double rq = Q2/q02, pm = proton_rms * proton_rms/hbarc/hbarc;
        returnFF = 1./(1.+rq)/(1.+ rq) + 2. *rq - 1./6. * Q2 * pm;
    
    }
    else if(ztgt == 82){//lead
        returnFF = 0;
	    while(i < 12){
	    returnFF = returnFF + A[i] * (pow(gamma,2) * cos(Q * r[i]/hbarc) + 2 * r[i] * hbarc/Q * sin(Q * r[i]/hbarc))/(pow(gamma, 2) + 2 * pow(r[i], 2)) * exp(-Q2/4 * pow(gamma, 2)/pow(hbarc, 2));
	    i++;
        }
    }else{ //for everything else use 2 - parameter fermi model, reference ?
      returnFF = 4 * pow(M_PI,2) * rho0 * pow(density_rho0.a_den[ztgt], 3)/(pow((Q * density_rho0.a_den[ztgt]), 2) * pow((sinh(M_PI * Q * density_rho0.a_den[ztgt])),2)) * (M_PI * Q * density_rho0.a_den[ztgt] * cosh(M_PI * Q * density_rho0.a_den[ztgt]) * sin(Q * density_rho0.c_den[ztgt]) - Q * density_rho0.c_den[ztgt] * cos(Q * density_rho0.c_den[ztgt]) * sinh(M_PI * Q * density_rho0.a_den[ztgt]));
      i = 0;
      while(i < 10){
	    returnFF = returnFF + 8 * M_PI * rho0 * pow(density_rho0.a_den[ztgt], 3) * pow((-1), (i - 1)) * float(i) * exp(-float(i) * density_rho0.c_den[ztgt]/density_rho0.a_den[ztgt])/pow((pow(float(i), 2) + pow((Q * density_rho0.a_den[ztgt]),2)),2);
	    i++;
      }
    }
    return returnFF;
}


void ZBQLINI(int seed, double ZBQLIX[43+1], double *bb) {

  static int init = 0;
  if(init) {
    if(init==1) printf("***WARNING**** You have called routine ZBQLINI more than once. Ignoring any subsequent calls.\n");
    ++init;
    return;
  } else  {init = 1;}

  double B = 4.294967291e9;

  if(!seed) {
    struct timespec tt;
    clock_gettime(CLOCK_REALTIME,&tt);
    ZBQLIX[1] = fmod(double(tt.tv_nsec)*4.e-9*B,B);
  } else {
    ZBQLIX[1] = fmod(double(seed),B);
  }

  for(int i = 1; i < 43; ++i) ZBQLIX[i+1] = fmod(ZBQLIX[i]*30269.,B);

  *bb = B;
  return;
}

//
//  Returns a uniform random number between 0 & 1, using
//  a Marsaglia-Zaman type subtract-with-borrow generator
//

double  ZBQLU01() {

  static int init = 0;
  static double  ZBQLIX[43+1], B;

  if(!init) {
    double  zz[43+1], bb;
    int iseed = _ISEED_; ZBQLINI(iseed,zz,&bb);
    B = bb;
    memcpy(ZBQLIX,zz,sizeof(zz));
    init = 1;
  }

  static int curpos = 1, id22 = 22, id43 = 43, C = 0.;
  double x, B2 = B, BINV = 1./B;

  while(1) {
    x = ZBQLIX[id22] - ZBQLIX[id43] - C;
    if(x<0.) {x += B; C = 1.;} else {C = 0.;}
    ZBQLIX[id43] = x;
    --curpos; --id22; --id43;
    if(!curpos) {
      curpos = 43;
    } else {
      if(!id22) {id22 = 43;} else {if(!id43) id43 = 43;}
    }
    if(x<BINV) {B2 *= B;} else {break;}
  }

   return x/B2;
}

//
//  Returns a random number uniformly distributed on (x1,x2)
//  Even if x1 > x2, this will work as x2-x1 will then be -ve
//
double ZBQLUAB(double x1, double x2) {

  if(x1==x2)
    printf("****WARNING**** (function ZBQLUAB) Upper and lower limits on uniform distribution are identical\n");

#if _OVERFLOW_PROTECTION_
  double z = -1.;
  const double eps = 2.5e-6;
  while(z<eps||z>1.-eps) z = ZBQLU01();
  z = x1+(x2-x1)*z;
#else
  double z = x1+(x2-x1)*ZBQLU01();
#endif

  return z;
}
