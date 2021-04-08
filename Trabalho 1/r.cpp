#include <iostream>
#include <iomanip>

constexpr double kappa {0.6};
constexpr int rho {600};
constexpr int cp {1200};
constexpr double h {15.0};
constexpr int g {100000};
constexpr double T0 {20.0};
constexpr double TL {20.0};
constexpr auto alpha = kappa/(rho*cp);


int main (int argc, char* argv[]){

	double dt {4};
	double r = 0.0;
	std::cout << std::setprecision(8) << "\tk\t dt \t dx \t r_max\n" ;
	for (double dx = 0.001; dx <= 0.03; dx +=0.001){
		r = (alpha*dt)/2*(dx*dx);
		std::cout << std::setw(10) << std::setprecision(8)<< kappa << '\t' << dt << '\t' << dx << '\t' << r << std::endl;
	}
}
