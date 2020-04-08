//g++ -std=c++17 -lgsl `root-config --cflags` `root-config --libs` -lflint -lmpfr -lgmp -larb dSigmaBorn.cpp RandomReal.cpp EnergyLookup.cpp -o Born_e7
#include <iostream>
#include <ctime>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"

#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include <math.h>

#include "arb.h"
#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

#include "EnergyLookup.h"
#include "RandomReal.h"


const int N = 1e5; 
//1e6 took 6.12 minutes


const double PI = 3.14159265;
//Experimental Parameters
int Z = 1; // Atomic Number
double max_theta = 13.12 * PI/180.; //It's the perpendicular "corner" from beam axis of FCAL
double min_theta = 0.9 * PI/180.;
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
	//double FullOnFF = abs( pow( 1 + q.Mag2() /pow(0.71,2) , -2 ) - AtomicFF ); //Full on FF in units of GeV
	//double FullOnFF = abs( pow( 1 + q.Mag2() /pow(0.71,2) , -2 ) - AtomicFF );
	double FullOnFF = abs( pow( 1 + (q.Mag2() /0.71) , -2) - AtomicFF );
 	double S1 = (1/c1) - (1/c2);
 	XYVector T1 = (p1/c1) + (p2/c2);

 	double dSigmaBorn = 2 * pow(fs_alpha,3) * pow(Z,2) * pow(E0,4) * pow(x,2) * pow(1 - x, 2)/(pow(PI,2)*pow(q.Mag2(),2))
 					    * (  pow(ElectronMass,2)*pow(S1,2) + ( pow(x,2) + pow(1 - x, 2) ) * T1.Mag2()
 					       - 2*x*(1-x)*xi3*(  
 					       					(p2.Mag2()/pow(c2,2))*cos(2*phi2)  + (p1.Mag2()/pow(c1,2))*cos(2*phi1) 
 					       				     + 2 * cos(phi1 + phi2) * sqrt(p2.Mag2() * p1.Mag2())/(c1*c2)  
 					       				   )	
 					      )
 					    * pow(FullOnFF,2)
 					    * abs(pow(theta1,2) * sin(theta1)) * abs(pow(theta2,2) * sin(theta2));

  return dSigmaBorn;
}



//Metropolis Method
 

int main(){
	std::clock_t start;
	double duration;
	start = std::clock();


   //create a Tree file tree1.root                                                                                                                                                                   

   //create the file, the Tree and a few branches                                                                                                                                                    
   TFile f("BornTreeBremNoNFF.root","recreate");
   TTree t1("BornTreeNoNFF","Born Tree no NFF");
   Double_t x_t, theta1_t, theta2_t, phi1_t, phi2_t, E0_t, t_t;
   Int_t ev_t, Attempts_t, HyperGeometricCrossSection_t;

   TH1D* Hist_t = new TH1D("t", ";Momentum Transfer Squared", 234, 0, .14 );

   XYZVector k1, k2;
   TLorentzVector p1, p2, BeamP4;

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
   t1.Branch("t", &t_t, "MomentumTransferSq/D");

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
	//while(dSigmaHypergeometric(x.X(), theta1.X(), theta2.X(), phi1.X(), phi2.X(), E0) / 
	//	  (dSigmaHypergeometric(x.Y(), theta1.Y(), theta2.Y(), phi1.Y(), phi2.Y(), E0) * RandomReal(0,1)) < 1 ) {
	while( (dSigmaBorn(x.X(), theta1.X(), theta2.X(), phi1.X(), phi2.X(), E0.X()) * EnergyLookup(E0.X()) ) / 
		 (dSigmaBorn(x.Y(), theta1.Y(), theta2.Y(), phi1.Y(), phi2.Y(), E0.Y()) * EnergyLookup(E0.Y()) ) < RandomReal(0,1) ) {	
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
 k1 = ThreeVec(theta1.X(), phi1.X()) * sqrt( pow(x.X(),2) * pow(E0.X(),2) - pow(ElectronMass,2));
 k2 = ThreeVec(theta2.X(), phi2.X()) * sqrt( pow(1 - x.X(),2) * pow(E0.X(),2) - pow(ElectronMass,2));
 p1.SetXYZM(k1.X(),k1.Y(),k1.Z(),ElectronMass); //v=(x,y,z,e=Sqrt(x*x+y*y+z*z+m*m))
 p2.SetXYZM(k2.X(),k2.Y(),k2.Z(),ElectronMass);
 BeamP4.SetXYZT(0,0,E0.X(),E0.X());
 t_t = -(BeamP4 - p1 - p2).Mag2();
 Hist_t->Fill(t_t);

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
