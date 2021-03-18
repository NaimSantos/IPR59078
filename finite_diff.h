/*
	12 implementações de Métodos de Diferenças Finitas para avaliar derivadas
*/
#pragma once

#include <functional>

//Diferenças finitas avançada para primeira derivada.
double fdiff_a_1st(double x, double h, std::function<double (double)> f){
	double res = (f(x+h) - f(x)) / h;
	return res;
}
double fdiff_a_1st_v2(double x, double h, std::function<double (double)> f){
	double res = (-f(x+2*h) + 4*f(x+h) - 3*f(x)) / (2*h);
	return res;
}

//Diferenças finitas avançada para segunda derivada.
double fdiff_a_2nd(double x, double h, std::function<double (double)> f){
	double res = (f(x+2*h) - 2*f(x+h) + f(x)) / (h*h);
	return res;
}
double fdiff_a_2nd_v2(double x, double h, std::function<double (double)> f){
	double res = (-f(x+3*h) - 4*f(x+2*h) - 5*f(x+h) + 2*f(x)) / (h*h);
	return res;
}

//Diferenças finitas recuada para primeira derivada.
double fdiff_r_1st(double x, double h, std::function<double (double)> f){
	double res = (f(x) - f(x-h)) / h;
	return res;
}
double fdiff_r_1st_v2(double x, double h, std::function<double (double)> f){
	double res = (3*f(x) - 4*f(x-h) + f(x-2*h)) / (2*h);
	return res;
}

//Diferenças finitas recuada para segunda derivada:
double fdiff_r_2nd(double x, double h, std::function<double (double)> f){
	double res = (f(x) - 2*f(x-h) + f(x-2*h)) / (h*h);
	return res;
}
double fdiff_r_2nd_v2(double x, double h, std::function<double (double)> f){
	double res = (2*f(x) - 5*f(x-h) + 4*f(x-2*h) + f(x-3*h)) / (h*h);
	return res;
}

//Diferenças finitas centrada  para primeira derivada:
double fdiff_c_1st(double x, double h, std::function<double (double)> f){
	double res = (f(x+h) - f(x-h)) / (2*h);
	return res;
}
double fdiff_c_1st_v2(double x, double h, std::function<double (double)> f){
	double res = (-f(x+2*h) + 8*f(x+h) - 8*f(x-h) + f(x-2*h)) / (12*h);
	return res;
}

//Diferenças finitas centrada  para segunda derivada:
double fdiff_c_2nd(double x, double h, std::function<double (double)> f){
	double res = (f(x+h) - 2*f(x) + f(x-h)) / (h*h);
	return res;
}
double fdiff_c_2nd_v2(double x, double h, std::function<double (double)> f){
	double res = (-f(x+2*h) + 16*f(x+h) - 30*f(x) + 16*f(x-h) - f(x-2*h)) / (12*(h*h));
	return res;
}