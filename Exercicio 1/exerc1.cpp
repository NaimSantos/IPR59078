#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "gsiedel.h"

double f(double x);

constexpr double dx = {0.1};
constexpr double dt = {0.05};

int main (int argc, char* argv[]){
	auto x = 0.0;
	for (int i = 0; i < 100; i++)
		std:: cout << f(x + 0.01*i) << ' ' ;
}

double f(double x){
	return (x <= 0.5) ? (x) : (1 - x);
}
