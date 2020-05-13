#include "EnergyLookup.h"
#include <cmath>
#include <vector>

std::vector<int> E_lookup = {810, 814, 818, 822, 826, 830, 834, 838, 842, 846, 850, 854, 858, 862, 866, 870, 874, 878, 882, 886, 890};

std::vector<int> EnergyCounts = {298, 334, 360, 347, 360, 383, 378, 432, 459, 466, 472, 447, 512, 508, 566, 558, 610, 585, 445, 96, 120};

int EnergyLookup(double E0){
	int Ecompare = std::floor(100*E0);
	if(Ecompare < 810 || Ecompare > 890){
		return 0;
	}
	int i = 0;
	while(E_lookup[i] < Ecompare){
		i++;
	}
	return EnergyCounts[i-1];
}


