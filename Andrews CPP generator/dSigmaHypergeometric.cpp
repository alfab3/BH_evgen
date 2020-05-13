//g++ -std=c++17 -lgsl `root-config --cflags` `root-config --libs` -lflint -lmpfr -lgmp -larb dSigmaHypergeometric.cpp RandomReal.cpp EnergyLookup.cpp -o Hyper_e5
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include <math.h>
#include "RandomReal.h"
#include "arb.h"
#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include <iostream>
#include <ctime>
#include "EnergyLookup.h"


const int N = 1e5; 
//1e4 events took 2 minutes and 10 seconds for hypergeometric. About a factor of 4 faster than muldrows code. 
//1e5 events took 17.4005 minutes.

const double PI = 3.14159265;
//Experimental Parameters
int Z = 1; // Atomic Number
double max_theta = 13.12 * PI/180.; //It's the perpendicular "corner" from beam axis of FCAL
double min_theta = 0.1 * 0.9 * PI/180.; //The 0.1 is to make it 10% of the minimum TOF angle. To really throw a bunch of these events down the beam pipe. 
double min_BeamEnergy  = 8.12;
double max_BeamEnergy  = 8.88;
const double xi3 = 1.0; 

//physics constants
const double fs_alpha = 0.0072973525664; //fine structure constant
const double ElectronMass = 5.109989461e-4; //Mass of electron
//User Defined Constants
const double nu = Z * fs_alpha;
const double alpha1 = 0.1;
const double alpha2 = 0.55;
const double alpha3 = 0.35;
const double b1 = 6.0;
const double b2 = 1.2;
const double b3 = 0.3;
const double mu1 = (ElectronMass * pow(Z,1/3)) * b1;
const double mu2 = (ElectronMass * pow(Z,1/3)) * b2;
const double mu3 = (ElectronMass * pow(Z,1/3)) * b3;
const double E0 = 8.7; //8.78 Energy of Incoming Photon (Beam Energy)
//cross section is probably 1 over gev squared in Bakmaev. To turn that into a cross section you need to multiply
//by two powers of hbar c
const double hbarc_sq = 3.89379e5; // (hbar * c)^2   in units of GeV*nm

/*
Bakmaev gives the cross section in inverse GeV^2. Need to convert it to an area (nanobarns.)

hbar * c = 197.3269804 MeV * fm 
1 barn = 1 bn = 100 fm^2   -> 10^(-2) bn = 1 fm^2
1 nanobarn = 10^(-9) bn    -> 10^9 nbn = 1 bn

10^7 nbn = 1 fm^2
1 MeV^2 = 10^(-3) GeV * 10^(-3) GeV = 10^(-6) GeV^2

1 MeV^2 fm^2 = 10^(-6) GeV^2 * 10^7 nbn = 10 GeV^2 nbn

(hbar * c)^2 = 38937.93 MeV^2 * fm^2 = 389379.3 GeV^2 nbn = 3.893793 * 10^(5) nbn GeV^2



*/




using namespace ROOT::Math;
using namespace std;

//FUNCTION For building 3 Vectors
XYZVector ThreeVec(double Theta, double Phi)
{
	XYZVector k(  ( sin(Theta) * cos(Phi)),  sin(Theta)*sin(Phi), cos(Theta) ) ; //three momentum for lepton
	return k;
}

//FUNCTION For the Born Cross Section
double dSigmaBorn(double x, double theta1, double theta2, double phi1, double phi2, double E0){
 

	XYZVector k1 = ThreeVec(theta1, phi1) * sqrt( pow(x,2) * pow(E0,2) - pow(ElectronMass,2)); 
	XYZVector k2 = ThreeVec(theta2, phi2) * sqrt( pow(1 - x,2) * pow(E0,2) - pow(ElectronMass,2));
	XYVector p1(k1.X(), k1.Y()); //transverse Momentum of Lepton 1
    XYVector p2(k2.X(), k2.Y()); //transverse Momentum of Lepton 2 
    XYVector q = p1 + p2; //the sum of lepton 1's transverse momentum and lepton 2's transverse momentum

	// c1 and c2 are defined in Bakmaev (4)
	double c1 = pow(k1.X(), 2) + pow(k1.Y(),2) + pow(ElectronMass,2);  
	double c2 = pow(k2.X(), 2) + pow(k2.Y(),2) + pow(ElectronMass,2);

	double z = 1 - pow(ElectronMass,2) * q.Mag2() / (c1 * c2); //defined on page 497 of Bakmaev

	double AtomicFF = 1 - q.Mag2() * (  alpha1/( pow(mu1,2) + q.Mag2() ) 
							       +    alpha2/( pow(mu2,2) + q.Mag2() )
							       +    alpha3/( pow(mu3,2) + q.Mag2() ) 
							          ); //What Muldrow calls Fatomic
	double FullOnFF = abs( pow( 1 + q.Mag2() /pow(0.71,2) , -2 ) - AtomicFF ); //Full on FF in units of GeV

 	double S1 = (1/c1) - (1/c2);
 	XYVector T1 = (p1/c1) + (p2/c2);

 	double dSigmaBorn = 2 * pow(fs_alpha,3) * pow(Z,2) * pow(E0,4) * pow(x,2) * pow(1 - x, 2)/(pow(PI,2)*pow(q.Mag2(),2))
 					    * (  pow(ElectronMass,2)*pow(S1,2) + ( pow(x,2) + pow(1 - x, 2) ) * T1.Mag2()
 					       - 2*x*(1-x)*xi3*(  
 					       					(p2.Mag2()/pow(c2,2))*cos(2*phi2)  + (p1.Mag2()/pow(c1,2))*cos(2*phi1) 
 					       				     + 2 * cos(phi1 + phi2) * sqrt(p2.Mag2() * p1.Mag2())/(c1*c2)  
 					       				   )	
 					      )
 					    * FullOnFF
 					    * abs(pow(theta1,2) * sin(theta1)) * abs(pow(theta2,2) * sin(theta2));;

  return dSigmaBorn;
}


//FUNCTION For the hypergeometric differential cross section
double dSigmaHypergeometric(double x, double theta1, double theta2, double phi1, double phi2, double E0){

	XYZVector k1 = ThreeVec(theta1, phi1) * sqrt( pow(x,2) * pow(E0,2) - pow(ElectronMass,2)); 
	XYZVector k2 = ThreeVec(theta2, phi2) * sqrt( pow(1 - x,2) * pow(E0,2) - pow(ElectronMass,2));
	XYVector p1(k1.X(), k1.Y()); //transverse Momentum of Lepton 1
    XYVector p2(k2.X(), k2.Y()); //transverse Momentum of Lepton 2 
    XYVector q = p1 + p2; //the sum of lepton 1's transverse momentum and lepton 2's transverse momentum

	// c1 and c2 are defined in Bakmaev (4)
	double c1 = pow(k1.X(), 2) + pow(k1.Y(),2) + pow(ElectronMass,2);  
	double c2 = pow(k2.X(), 2) + pow(k2.Y(),2) + pow(ElectronMass,2);
	double z = 1 - pow(ElectronMass,2) * q.Mag2() / (c1 * c2); //defined on page 497 of Bakmaev
	
	double AtomicFF = 1 - q.Mag2() * (  alpha1/( pow(mu1,2) + q.Mag2() ) 
							       +    alpha2/( pow(mu2,2) + q.Mag2() )
							       +    alpha3/( pow(mu3,2) + q.Mag2() ) 
							          ); //What Muldrow calls Fatomic

	double FullOnFF = abs( pow( 1 + q.Mag2() /pow(0.71,2) , -2 ) - AtomicFF ); //Full on FF in units of GeV
	//This is what Muldrow calls F
	


		/************* acb functions--gamma and hypergeometrics-- ***********/


		acb_t hres1; acb_t ha1; acb_t hb1; acb_t hc1; acb_t hd1; //Hypergeometric variable initialization I
	    acb_t hres2; acb_t hb2; acb_t hc2; //Hypergeometric variable initialization II
	    acb_t factor;
	    // We will need to create new acb_t variables for our previously generated p1 and p2 vectors.
		acb_t p1x; acb_t p1y; 
		acb_t p2x; acb_t p2y;
		acb_t Px; acb_t Py; acb_t Pxconj; acb_t Pyconj;
		acb_t cx_fac1; acb_t cx_fac2; acb_t Complex2For_cx_fac2; acb_t cx_fac2conj;
		acb_t MassSq_c1; acb_t MassSq_c2;
		acb_t WunpTerm1res; acb_t WunpTerm2res; acb_t WunpTerm2resConj;
				//Introducing the terms for Wpol
		acb_t WpolTerm1; acb_t WpolTerm2; acb_t WpolTerm3;
		acb_t cosine2phi1; acb_t cosine2phi2; acb_t cosinephi1_plus_phi2;
		acb_t hres2conj; acb_t hres2squared;


		acb_init(hres1); acb_init(ha1); acb_init(hb1); acb_init(hc1); acb_init(hd1);
		acb_init(hres2); acb_init(hb2); acb_init(hc2); 
		acb_init(factor);


		acb_set_d_d(ha1,0, nu); acb_set_d_d(hb1, 0, -nu); acb_set_d(hc1,1); acb_set_d(hd1, z); //Setting values for Gauss HYP1
		acb_set_d_d(hb2,1,-nu); acb_set_d(hc2,2);  //Setting values for Gauss HYP2
		acb_set_d_d(factor,1,-nu);
	


		acb_hypgeom_2f1(hres1, ha1, hb1, hc1, hd1, ACB_HYPGEOM_2F1_ABC, 30); //Computation of the first hypergeometric
		acb_hypgeom_2f1(hres2, ha1, hb2, hc2, hd1, ACB_HYPGEOM_2F1_ABC, 30); //Computation of the second hypergeometric. Results are stored in hres1, hres2 respectively.
		/*
		cout << "prior hres2 = ";
		acb_printn(hres2,30,0);
		cout <<"\n";
		*/
		acb_mul(hres2, hres2, factor, 30);	//Bakmaev has a factor of (1 - i) in front of the second hypergeometric, so let's add it:

		/*
		cout << "hypergeometric results" << "\n";
		acb_printn(hres1,30,0);
		cout << "\n";
		acb_printn(hres2,30,0);
		cout << "\n";
		*/


		//We are now ready to the unpolarized and polarized weighting factors, Wunp, Wpol. 

		// For Wunp we want: 
		// Wunp = m^2 (x^2 + (1 - x)^2)  |F2 p1/c1 + (2*F1 - F2) p2/c2 |^2       +       |F2 - F1 + (2*F1 - F2)m^2/c2 - F2 m^2/c1|^2
		// 						call what's directly above WunpTerm1res								call this one WunpTerm2res
	


		acb_init(p1x); acb_init(p1y); acb_init(p2x); acb_init(p2y); //the components of the transverse vectors
		acb_init(Px); acb_init(Py); acb_init(Pxconj); acb_init(Pyconj);  
		acb_init(cx_fac1); acb_init(cx_fac2); acb_init(Complex2For_cx_fac2); acb_init(cx_fac2conj);
		acb_init(WunpTerm1res); acb_init(WunpTerm2res); acb_init(WunpTerm2resConj);
		acb_init(MassSq_c1); acb_init(MassSq_c2);

		acb_set_d(p1x, p1.X()/c1); acb_set_d(p1y, p1.Y()/c1); 
		acb_set_d(p2x, p2.X()/c2); acb_set_d(p2y, p2.Y()/c2);
		acb_set_d(Complex2For_cx_fac2,2);


	    acb_mul(cx_fac2, Complex2For_cx_fac2, hres1, 30); acb_sub(cx_fac2, cx_fac2, hres2, 30); //This is 2F1 - F2. Can call it again in the second term.
/*	   
	    cout << "2F_1 - F_2 = ";
	    acb_printn(cx_fac2, 30, 0);
	    cout << "\n";
*/
		acb_mul(p1x, p1x, hres2, 20); acb_mul(p1y, p1y, hres2, 20); //This is F2 p1/c1's two components
		acb_mul(p2x, p2x, cx_fac2, 20); acb_mul(p2y, p2y, cx_fac2, 20); //This is (2F1 - F2)p2/c2's two components
		acb_add(Px, p1x, p2x, 20); acb_add(Py, p1y, p2y, 20); //Px is [F2 p1/c1 + (2*F1 - F2) p2/c2].X(), Py is [F2 p1/c1 + (2*F1 - F2) p2/c2].Y()
		acb_conj(Pxconj, Px); acb_conj(Pyconj, Py); //Sets Pxconj to be the complex conjugate of Px. Does the same for Py.

		acb_mul(WunpTerm1res, Px, Pxconj, 20); acb_addmul(WunpTerm1res, Py, Pyconj, 20); //First computes Px times Px*, then computes Py times Py* and adds it to the result of Px times Px*

	    // Now for the |F2 - F1 + (2*F1 - F2)m^2/c2 - F2 m^2/c1|^2
		acb_set_d(MassSq_c1, pow(ElectronMass,2)/c1 ); acb_mul(MassSq_c1, MassSq_c1, hres2,20); //F2 m^2/c1
		acb_set_d(MassSq_c2, pow(ElectronMass,2)/c2 ); acb_mul(MassSq_c2, MassSq_c2, cx_fac2,20); // (2*F1 - F2)m^2/c2



	    acb_sub(WunpTerm2res, hres2, hres1, 20); //F2 - F1
	    acb_add(WunpTerm2res, WunpTerm2res, MassSq_c2,20); //F2 - F1 +    (2*F1 - F2)m^2/c2
	    acb_sub(WunpTerm2res, WunpTerm2res, MassSq_c1, 20); //F2 - F1 + (2*F1 - F2)m^2/c2  -      F2 m^2/c1
	    acb_conj(WunpTerm2resConj, WunpTerm2res); acb_mul(WunpTerm2res, WunpTerm2res, WunpTerm2resConj, 20);



		double Wunp = pow(ElectronMass, 2) * ( pow(x, 2) + pow(1 - x, 2) ) * 
			(arf_get_d(arb_midref(acb_realref(WunpTerm1res)), ARF_RND_NEAR) + arf_get_d(arb_midref(acb_realref(WunpTerm2res)), ARF_RND_NEAR));
		

		//Now we're ready to do Wpol. 
		//W_pol = -2x(1-x)m^2[  |2*F1 - F2|^2 * (p_2/c_2)^2 cos(2*phi_2)   		term 1
		//					  + |F2|^2 * (p_1/c_1)^2 cos(2*phi_1) 				term 2
		//					  + 2 Re(conj(F2)(2*F1 - F2)) |p_2||p_1|/(c_2 c_1) cos(phi_1 + phi_2) 	term 3



		acb_init(WpolTerm1); acb_init(WpolTerm2); acb_init(WpolTerm3);
		acb_init(cosine2phi1); acb_init(cosine2phi2); acb_init(cosinephi1_plus_phi2);
		acb_init(hres2conj); acb_init(hres2squared);

		//term 1
		acb_set_d(WpolTerm1, (pow(p2.X(),2) + pow(p2.Y(),2))/pow(c2,2)); //(p2/c2)^2
		acb_set_d(cosine2phi2, cos(2*phi2)); //cos(2*phi_2)
		acb_mul(WpolTerm1, WpolTerm1, cosine2phi2,20); //(p2/c2)^2 * cos(2*phi_2)
		acb_conj(cx_fac2conj, cx_fac2);
		acb_mul(WpolTerm1, WpolTerm1, cx_fac2conj, 20); //I think it's okay to do this. I think this falls under associativity? 
		acb_mul(WpolTerm1, WpolTerm1, cx_fac2, 20); //|2*F1 - F2|^2 * (p_2/c_2)^2 cos(2*phi_2) finished
		

		/*
		cout << "WpolTerm1 = ";
		acb_printn(WpolTerm1, 30,0);
		cout << "\n";
		//term 1 matches muldrow
		*/

		//term 2
		acb_set_d(WpolTerm2, (pow(p1.X(), 2) + pow(p1.Y(),2) )/pow(c1,2) ); //(p1/c1)^2
		acb_set_d(cosine2phi1, cos(2*phi1)); //cos(2*phi_1)
		acb_mul(WpolTerm2, WpolTerm2,cosine2phi1,20); //(p_1/c_1)^2 cos(2*phi_1)
		acb_conj(hres2conj, hres2);  
		acb_mul(hres2squared, hres2conj, hres2,20);
		acb_mul(WpolTerm2, WpolTerm2, hres2squared, 20); //|F2|^2 * (p_1/c_1)^2 cos(2*phi_1)  finished

		/*
		cout << "WpolTerm2 = ";
		acb_printn(WpolTerm2, 30, 0);
		cout << "\n";
		*/


		//term 3
		acb_set_d(WpolTerm3, cos(phi1 + phi2) * sqrt(p1.Mag2())*sqrt(p2.Mag2())/( c1 * c2 ) ); // |p_2||p_1|/(c_2 c_1) cos(phi_1 + phi_2) 
	 	acb_mul(WpolTerm3, WpolTerm3, cx_fac2, 20); //(2*F1 - F2)) |p_2||p_1|/(c_2 c_1) cos(phi_1 + phi_2)
	 	acb_mul(WpolTerm3, WpolTerm3, hres2conj, 20); //conj(F2)(2*F1 - F2)) |p_2||p_1|/(c_2 c_1) cos(phi_1 + phi_2) 

	 	double Wpol = -2*x*(1 - x)*pow(ElectronMass,2)* (
	 				  + (arf_get_d(arb_midref(acb_realref(WpolTerm1)), ARF_RND_NEAR))
	 				  + (arf_get_d(arb_midref(acb_realref(WpolTerm2)), ARF_RND_NEAR))
	 				  + 2 * (arf_get_d(arb_midref(acb_realref(WpolTerm3)), ARF_RND_NEAR)) );
	 											



		//acb variables that are from Wunp and Wpol
		acb_clear(hres1); acb_clear(ha1); acb_clear(hb1); acb_clear(hc1); acb_clear(hd1); 
		acb_clear(hres2); acb_clear(hres2conj); acb_clear(hres2squared); 
		acb_clear(hb2); acb_clear(hc2); acb_clear(factor); 
		acb_clear(p1x); acb_clear(p1y); acb_clear(p2x); acb_clear(p2y); 
		acb_clear(Px); acb_clear(Pxconj); acb_clear(Pyconj); acb_clear(Py);
		acb_clear(cx_fac1); acb_clear(cx_fac2); acb_clear(Complex2For_cx_fac2);
		acb_clear(WunpTerm1res); acb_clear(WunpTerm2res); acb_clear(WunpTerm2resConj);
		acb_clear(MassSq_c1); acb_clear(MassSq_c2);

		acb_clear(WpolTerm1); acb_clear(WpolTerm2); acb_clear(WpolTerm3);
		acb_clear(cosine2phi1); acb_clear(cosine2phi2); acb_clear(cosinephi1_plus_phi2);

		//Gamma Function calculation
		//Need to compute |Gamma(1 - i nu)|^4; 	let z = 1 - i nu
		acb_t AbsGamma4; acb_t GammaConj; acb_t arg_of_gamma;
		acb_init(AbsGamma4); acb_init(GammaConj); acb_init(arg_of_gamma);
		acb_set_d_d(arg_of_gamma, 1, - nu);
		acb_gamma(AbsGamma4, arg_of_gamma, 20); //Just one power of it
		acb_conj(GammaConj, AbsGamma4);
		acb_mul(AbsGamma4, AbsGamma4, GammaConj, 20); //|Gamma(z)|^2  
		acb_mul(AbsGamma4, AbsGamma4, AbsGamma4, 20); //|Gamma(z)|^4

		double AbsGamma4th =  arf_get_d(arb_midref(acb_realref(AbsGamma4)), ARF_RND_NEAR);

		acb_clear(AbsGamma4); acb_clear(GammaConj); acb_clear(arg_of_gamma);
	/******** Done with ARB ************/
	//And now, the moment you've been waiting for: the differential cross section!
	double result = (2*fs_alpha* pow(nu,2) * pow(E0,4) /(pow(PI,2) * pow(ElectronMass,2)) )
								  *(pow(x,2) * pow(1-x,2)/(pow(q.Mag2(),2)))
								  * AbsGamma4th
								  * (Wunp  + xi3 * Wpol)
								  * pow(FullOnFF,2)
	  							  * abs(pow(theta1,2) * sin(theta1)) * abs(pow(theta2,2) * sin(theta2))
								  * hbarc_sq ;

/* diagnostics
cout << "k1 = {" << k1.X() << "," << k1.Y() << "," << k1.Z() << "} \n";
cout << "k2 = {" << k2.X() << "," << k2.Y() << "," << k2.Z() << "} \n";
cout << "c1 = " << c1 << "\n"; 
cout << "c2 = " << c2 << "\n";
cout << "z = " << z << "\n";
cout << "AtomicFF = " << AtomicFF << "\n";
cout << "FullOnFF = " << FullOnFF << " \n ";
cout << "Wpol = " << Wpol << "\n";
cout << "Wunp = " << Wunp << "\n";
*/

	return result;
 }


//Metropolis Method
 

int main(){
	std::clock_t start;
	double duration;
	start = std::clock();


   //create a Tree file tree1.root                                                                                                                                                                   

   //create the file, the Tree and a few branches                                                                                                                                                    
   TFile f("HyperTreeBrem.root","recreate");
   TTree t1("HyperTree","Hypergeometric Tree");
   Double_t x_t, theta1_t, theta2_t, phi1_t, phi2_t, E0_t;
   Int_t ev_t, Attempts_t, HyperGeometricCrossSection_t;
   Double_t CrossSectionRatio_t;
   t1.Branch("x",&x_t,"x/D");
   t1.Branch("theta1",&theta1_t,"theta1/D");
   t1.Branch("theta2",&theta2_t,"theta2/D");
   t1.Branch("phi1",&phi1_t,"phi1/D");
   t1.Branch("phi2",&phi2_t,"phi2/D");
   t1.Branch("E0",&E0_t,"E0/D");
   t1.Branch("CrossSectionRatio",&CrossSectionRatio_t, "CrossSectionRatio/D");
   t1.Branch("HyperGeometricCrossSection",&HyperGeometricCrossSection_t, "HyperGeometricCrossSection_t/D");
   t1.Branch("Attempts",&Attempts_t,"Attempts/I");
   t1.Branch("ev",&ev_t,"ev/I");

XYVector E0(RandomReal(min_BeamEnergy, max_BeamEnergy), 
			RandomReal(min_BeamEnergy, max_BeamEnergy)); //Energy of Incoming Photon (Beam Energy)


XYVector x(RandomReal(ElectronMass/min_BeamEnergy, 1 - (ElectronMass/max_BeamEnergy) ), 
	       RandomReal(ElectronMass/min_BeamEnergy, 1 - (ElectronMass/max_BeamEnergy) ));

//Make sure our two initial points in x and E0 are good:
while(x.X() < (ElectronMass/E0.X()) || x.X() > (1 - (ElectronMass/E0.X()) ))
{
 x.SetX(RandomReal(ElectronMass/min_BeamEnergy, 1 - (ElectronMass/max_BeamEnergy) ));
}

while(x.Y() < (ElectronMass/E0.Y()) || x.Y() > (1 - (ElectronMass/E0.Y()) ))
{
 x.SetY(RandomReal(ElectronMass/min_BeamEnergy, 1 - (ElectronMass/max_BeamEnergy) ));
}


XYVector theta1( 1/(RandomReal( 1/max_theta, 1/min_theta) ),
				 1/(RandomReal( 1/max_theta, 1/min_theta) ) );
XYVector theta2( 1/(RandomReal( 1/max_theta, 1/min_theta) ),
				 1/(RandomReal( 1/max_theta, 1/min_theta) ) ); 
XYVector phi1(RandomReal(0, 2*PI),RandomReal(0, 2*PI));
XYVector phi2(RandomReal(0, 2*PI),RandomReal(0, 2*PI));


int i = 1;
/*cout << "first x = " << x.X() << "     Previous Generated x = " << x.Y() << "     iteration = " << i << "\n" ;
cout << "first theta1 = " << theta1.X() << "     Previous theta1 = " << theta1.Y() << "\n";
cout << "first theta2 = " << theta2.X() << "     Previous theta2 = " << theta2.Y() << "\n";
cout << "first phi1 = " << phi1.X() << "     Previous phi1 = " << phi1.Y() << "\n";
cout << "first phi2 = " << phi2.X() << "     Previous phi2 = " << phi2.Y() << "\n";
cout << "break condition number " << breakCondition << "\n";*/
//cout << "initial Hypergeometric result " << dSigmaHypergeometric(x.X(), theta1.X(), theta2.X(), phi1.X(), phi2.X(), E0) / 
	 // dSigmaHypergeometric(x.Y(), theta1.Y(), theta2.Y(), phi1.Y(), phi2.Y(), E0) << "\n";
//while(dSigmaHypergeometric(x.X(), theta1.X(), theta2.X(), phi1.X(), phi2.X(), E0) / 
//	  dSigmaHypergeometric(x.Y(), theta1.Y(), theta2.Y(), phi1.Y(), phi2.Y(), E0)    <= breakCondition  ){
//cout << "Hypergeometric result " << dSigmaHypergeometric(x.X(), theta1.X(), theta2.X(), phi1.X(), phi2.X(), E0) / 
	  //dSigmaHypergeometric(x.Y(), theta1.Y(), theta2.Y(), phi1.Y(), phi2.Y(), E0) << "\n";


for(int n = 1; n <= N; n++){
//cout << phi2.Y() << "\n";	
//double breakCondition = RandomReal(0,1);
i = 1;
	while(dSigmaHypergeometric(x.X(), theta1.X(), theta2.X(), phi1.X(), phi2.X(), E0.X()) * EnergyLookup(E0.X())  / 
			  (dSigmaHypergeometric(x.Y(), theta1.Y(), theta2.Y(), phi1.Y(), phi2.Y(), E0.Y()) *EnergyLookup(E0.Y()) * RandomReal(0,1)) < 1 ) {
	//while( (dSigmaBorn(x.X(), theta1.X(), theta2.X(), phi1.X(), phi2.X(), E0.X()) * EnergyLookup(E0.X()) ) / 
	//	 (dSigmaBorn(x.Y(), theta1.Y(), theta2.Y(), phi1.Y(), phi2.Y(), E0.Y()) * EnergyLookup(E0.Y()) * RandomReal(0,1)) < 1 ) {	
		/*  	
		cout << "Numerator = " << dSigmaHypergeometric(x.X(), theta1.X(), theta2.X(), phi1.X(), phi2.X(), E0) << "\n";
		cout << "Denominator = " << dSigmaHypergeometric(x.Y(), theta1.Y(), theta2.Y(), phi1.Y(), phi2.Y(), E0) << "\n";
		cout << "Ratio = " << dSigmaHypergeometric(x.X(), theta1.X(), theta2.X(), phi1.X(), phi2.X(), E0) / 
	                          dSigmaHypergeometric(x.Y(), theta1.Y(), theta2.Y(), phi1.Y(), phi2.Y(), E0) << "\n";
	      */                    
	/*while(dSigmaBorn(x.X(), theta1.X(), theta2.X(), phi1.X(), phi2.X(), E0)  / 
		  dSigmaBorn(x.Y(), theta1.Y(), theta2.Y(), phi1.Y(), phi2.Y(), E0)    <= RandomReal(0,1) ){
		cout << "Numerator = " << dSigmaBorn(x.X(), theta1.X(), theta2.X(), phi1.X(), phi2.X(), E0) << "\n";
		cout << "Denominator = " << dSigmaBorn(x.Y(), theta1.Y(), theta2.Y(), phi1.Y(), phi2.Y(), E0) << "\n";
		cout << "Ratio = " << dSigmaBorn(x.X(), theta1.X(), theta2.X(), phi1.X(), phi2.X(), E0) / 
	                          dSigmaBorn(x.Y(), theta1.Y(), theta2.Y(), phi1.Y(), phi2.Y(), E0) << "\n";
*/

	//cout << "current x = " << x.X() << "     Previous Generated x = " << x.Y() << "     iteration = " << i << "\n" ;
	//cout << "current theta1 = " << theta1.X() << "     Previous theta1 = " << theta1.Y() << "\n";
	//cout << "current theta2 = " << theta2.X() << "     Previous theta2 = " << theta2.Y() << "\n";
	//cout << "current phi1 = " << phi1.X() << "     Previous phi1 = " << phi1.Y() << "\n";
	//cout << "current phi2 = " << phi2.X() << "     Previous phi2 = " << phi2.Y() << "\n";

	//x.SetY(x.X()); //move to the previous value. 
	//cout << "Numerator = " << dSigmaHypergeometric(x.X(), theta1.X(), theta2.X(), phi1.X(), phi2.X(), E0) << "\n";
	//cout << "Denominator = " << dSigmaHypergeometric(x.Y(), theta1.Y(), theta2.Y(), phi1.Y(), phi2.Y(), E0) << "\n";
	//cout << "Ratio = " << dSigmaHypergeometric(x.X(), theta1.X(), theta2.X(), phi1.X(), phi2.X(), E0) / 
	//	  dSigmaHypergeometric(x.Y(), theta1.Y(), theta2.Y(), phi1.Y(), phi2.Y(), E0) << "\n";


	E0.SetX(RandomReal(min_BeamEnergy, max_BeamEnergy));

	x.SetX(RandomReal(ElectronMass/min_BeamEnergy, 1 - (ElectronMass/max_BeamEnergy) )); //make a new value for x
		while(x.X() < (ElectronMass/E0.X()) || x.X() > (1 - (ElectronMass/E0.X()) ))
		{
 			x.SetX(RandomReal(ElectronMass/min_BeamEnergy, 1 - (ElectronMass/max_BeamEnergy) ));
		}
	//theta1.SetY(theta1.X());
	theta1.SetX(1/(RandomReal( 1/max_theta, 1/min_theta) ));

	//theta2.SetY(theta2.X());
	theta2.SetX(1/(RandomReal( 1/max_theta, 1/min_theta) ));

	//phi1.SetY(phi1.X());
	phi1.SetX(RandomReal(0, 2*PI));

	//phi2.SetY(phi2.X());
	phi2.SetX(RandomReal(0, 2*PI));

	i++; 
	}
//cout << phi2.X() << "\n";	
 x_t = x.X();
 theta1_t = theta1.X(); theta2_t = theta2.X();
 phi1_t = phi1.X(); phi2_t = phi2.X();
 E0_t = E0.X(); Attempts_t = i;
 ev_t = n; 

//CrossSectionRatio_t = dSigmaHypergeometric(x.X(), theta1.X(), theta2.X(), phi1.X(), phi2.X(), E0) / 
//		  dSigmaHypergeometric(x.Y(), theta1.Y(), theta2.Y(), phi1.Y(), phi2.Y(), E0); 

//HyperGeometricCrossSection_t = dSigmaHypergeometric(x.X(), theta1.X(), theta2.X(), phi1.X(), phi2.X(), E0);

 t1.Fill();

//cout << "\n"; 
//cout << "\n"; cout << "\n"; cout << "\n"; cout << "\n"; cout << "\n"; cout << "\n"; cout << "\n"; 
//cout << "done!" << "\n";
//cout << i << "\n";
//cout << dSigmaHypergeometric(x.X(), theta1.X(), theta2.X(), phi1.X(), phi2.X(), E0) / 
//		  dSigmaHypergeometric(x.Y(), theta1.Y(), theta2.Y(), phi1.Y(), phi2.Y(), E0) << "\n";

E0.SetX(RandomReal(min_BeamEnergy, max_BeamEnergy));


x.SetY(x.X()); //move to the previous value. 
x.SetX(RandomReal(ElectronMass/min_BeamEnergy, 1 - (ElectronMass/max_BeamEnergy) )); //make a new value for x
		while(x.X() < (ElectronMass/E0.X()) || x.X() > (1 - (ElectronMass/E0.X()) ))
		{
 			x.SetX(RandomReal(ElectronMass/min_BeamEnergy, 1 - (ElectronMass/max_BeamEnergy) ));
		}

theta1.SetY(theta1.X());
theta1.SetX(1/(RandomReal( 1/max_theta, 1/min_theta) ));

theta2.SetY(theta2.X());
theta2.SetX(1/(RandomReal( 1/max_theta, 1/min_theta) ));

phi1.SetY(phi1.X());
phi1.SetX(RandomReal(0, 2*PI));


phi2.SetY(phi2.X());
phi2.SetX(RandomReal(0, 2*PI));


}

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout<<"printf: "<< duration <<'\n';

t1.Write();



return 0;
}

//t1->Draw("CrossSectionRatio", "CrossSectionRatio < 10")
