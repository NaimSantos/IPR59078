/*
	12 implementações de Métodos de Diferenças Finitas para avaliar derivadas
*/
#pragma once

#include <functional>

//Diferenças finitas avançada para a primeira derivada:
double fdiff_a_1st(double x, double h, std::function<double (double)> f){
	double res = (f(x+h) - f(x)) / h;
	return res;
}
//Diferenças finitas avançada para a segunda derivada:
double fdiff_a_2nd(double x, double h, std::function<double (double)> f){
	double res = (f(x+2*h) - 2*f(x+h) + f(x)) / (h*h);
	return res;
}
//Diferenças finitas recuada para a primeira derivada:
double fdiff_r_1st(double x, double h, std::function<double (double)> f){
	double res = (f(x) - f(x-h)) / h;
	return res;
}
//Diferenças finitas recuada para a segunda derivada:
double fdiff_r_2nd(double x, double h, std::function<double (double)> f){
	double res = (f(x) - 2*f(x-h) + f(x-2*h)) / (h*h);
	return res;
}
//Diferenças finitas centrada para a primeira derivada:
double fdiff_c_1st(double x, double h, std::function<double (double)> f){
	double res = (f(x+h) - f(x-h)) / (2*h);
	return res;
}
//Diferenças finitas centrada para a segunda derivada:
double fdiff_c_2nd(double x, double h, std::function<double (double)> f){
	double res = (f(x+h) - 2*f(x) + f(x-h)) / (h*h);
	return res;
}
