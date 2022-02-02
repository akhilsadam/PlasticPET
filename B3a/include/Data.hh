#ifndef DATA_TABLES
#define DATA_TABLES

#include <vector>
using namespace std;
// SiPM QE (50um HAMAMATSU S14160-3050HS):
	// https://www.hamamatsu.com/us/en/product/optical-sensors/mppc/mppc_mppc-array/index.html
	const vector<double> sipm_QE_lambda{274.7,279.18,284.55,287.15,289.58,294.72,300,319.65,331.32,341.18,361.86,374.38,406.66,445.33,474.2,498.65,543.01,595.5,685.93,798.8,898.97};
	const vector<double> sipm_QE{2.614,3.307,4.618,8.011,17.034,25.671,29.835,35.54,37.004,38.699,39.931,43.091,48.332,50.41,49.173,45.083,37.674,29.725,18.378,9.805,4.396};

#endif 
