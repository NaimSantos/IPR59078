#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>	//std::fill()

/*
	TDMA: TriDiagonal Matrix Algorithm ("Thomas")
	a = diagonal inferior
	b = diagonal principal
	c = diagonal superior
	d = termos independentes
*/

void tdma_solver(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, std::vector<double>& d){
	auto n = static_cast<int>(d.size()-1);

	auto c_temp = c;
	
	c_temp[0] /= b[0];
	d[0] /= b[0];

	for (int i = 1; i < n; i++){
		c_temp[i] = (c_temp[i] ) / (b[i] - a[i]*c_temp[i-1]);
		d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c_temp[i-1]);
	}

	d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c_temp[n-1]);

	for (int i = n; i-- > 0;){
		d[i] = d[i] - (c_temp[i]*d[i+1]);
	}
}
/*
int main(int argc, char* arg[]){
	std::vector<double> a1 = {0, 1, 2, 3};	//diagonal inferior
	std::vector<double> b1 = {2, 3, 5, 8};	//diagonal principal
	std::vector<double> c1 = {1, 2, 1, 0};	//diagonal superior

	std::vector<double> d1 = {7, 19, 31, 52};

	tdma_solver(a1, b1, c1, d1);

	std::cout << "T = ";
	for (auto &e: d1){
		std::cout << e << ' ';
	}
}
*/