#include <iostream>
#include <cmath>
#include <array>
#include <string>   //stof

int main (int argc, char* argv[]){

	double ti = 0.05;
	if (argc > 1)
		ti = std::stof(argv[1]);
	std::cout << "argc = " << argc << std::endl;

	std::cout << "ti = " << ti << std::endl;

	double dt {0.1};
	double dx {0.1};
	double r_max {0.5};
	const double k {0.025};

	for (double r = 0.1; r <= r_max; r +=ti){
		dt = (r * dx * dx) / k;
		std::cout << "r = " << r << " , dt = " << dt << std::endl;
	}
}