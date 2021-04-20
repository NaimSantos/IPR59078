#pragma once

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>         //std::fabs
#include <algorithm>     //std::max
#include <functional>

constexpr auto NPI {4*std::atan(1)};

// Registro do tempo de execução de uma parte do código usando o escopo de objetos
class CustomTimer
{
	public:
		CustomTimer(){
			startpoint = std::chrono::steady_clock::now();
		}
		~CustomTimer(){
			auto endpoint = std::chrono::steady_clock::now();
			auto start = std::chrono::time_point_cast<std::chrono::milliseconds>(startpoint).time_since_epoch().count();
			auto end = std::chrono::time_point_cast<std::chrono::milliseconds>(endpoint).time_since_epoch().count();

			auto total = end - start;
			double s = total*0.001;
			std::cout << "Tempo decorrido: " << s << " segundos (" << total << " milisegundos)" << std::endl;
		}

	private:
		std::chrono::time_point<std::chrono::steady_clock> startpoint;
};


// Implementação para o método de Gauss-Siedel
constexpr double eps = 0.000001;
constexpr unsigned int MAX_ITER {50};

void GS_Solver(const std::vector<std::vector<double>>& A,  std::vector<double>& B){

	auto n = static_cast<int>(B.size());
	auto C = B;
	std::vector<double> X (n, 0.0);
	size_t counter {1};
	bool teste {false};
	double erro_max {1};

	while((erro_max >= eps) && counter<MAX_ITER){
		for (int i = 0; i < n; i++){
			auto x_old = C[i];
			C[i] = (B[i] / A[i][i]);
			for (int j = 0; j < n; j++){
				if (j==i)
					continue;
				C[i] = C[i] - ((A[i][j] / A[i][i]) * X[j]);
			}
			auto x_new = C[i];
			auto erro_new = std::fabs((x_old - x_new) / x_new);
			(i == 0) ? (erro_max = erro_new) : (erro_max = std::max(erro_new, erro_max));
			X[i] = C[i];
		}
		counter++;
	}
	B = C;
}

// Integração numérica pela regra 1/3 de Simpson
double int_simpson(double a, double b, double h, std::function<double (double)> f){	
	const auto n = static_cast<int>( std::floor((std::fabs(b - a)) / h));
	//n - 1 pontos internos
	double sum_odds {0.0};
	for (int i = 1; i < n; i += 2){
		sum_odds += f(a + i*h);
	}
	double sum_evens {0.0};
	for (int i = 2; i < n; i += 2){
		sum_evens += f(a + i*h);
	}
	return (f(a) + f(b) + 2*sum_evens + 4*sum_odds) * (h/3);
}

//A derivada é avaliada por diferencas finitas centrada. Erro: O(h^4)
double f_diff(double x, double h, std::function<double (double)> f){
	double res = (-f(x+2*h) + 8*f(x+h) - 8*f(x-h) + f(x-2*h)) / (12*h);
	return res;
}

//Raiz numérica via Newton Raphson
double solveNewton(double x, std::function <double (double)> fnwtn){
	auto h = fnwtn(x) / f_diff(x, 0.001, fnwtn);
	unsigned int i {0};
	while (std::fabs(h) >= 0.00001 && i<100){
		h = fnwtn(x) / f_diff(x, 0.001, fnwtn);
		x = x - h;
		i++;
	}

	return x;
}