#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

//#include "gsiedel.hpp"
//#include "tdma.hpp"

double f(double x);
//void linspacefill(std::vector<double>& Vec, const int Num, const double xi = 0.0, const double xf = 1.0);
void explictsolver(std::vector<double>& Temp, const double r = 0.5);
void printvec(const std::vector<double>& Vec);


//Dados do domínio:
constexpr double dx {0.1}; //refinamento da malha espacial
constexpr double dt {1};	//passo de tempo
constexpr double L {1.0};	//comprimento total do domínio
constexpr double ti {0.0};	//tempo inicial da simulação
constexpr double tf {100.0};	//tempo final de simulação
constexpr auto N  = L/dx + 1; //número de nós da malha (intervalos + 1)
constexpr double I = (tf - ti) / dt; //número de passos de tempo

//Dados do problema:
constexpr double kappa {0.6};
constexpr double rho {600};
constexpr double cp {1200};
constexpr double h {15.0};
constexpr auto g {100000};
constexpr double T_zero {20.0};

int main (int argc, char* argv[]){

	const auto alfa = kappa / (rho * cp);
	const auto r = alfa * dt / (dx * dx);
	std::cout << "dt = " << dt << ", dx = " << dx  << ", alfa = " << alfa << ", CFL (r) = " << r << std::endl;


	std::vector<double> Temperature(N, 0.0);

	//Preencher os temperaturas iniciais a partir da função fornecida no problema:
	auto n = Temperature.size();
	for (int i = 0; i < n; i++){
		Temperature[i] = f(i*dx);
	}

	printvec(Temperature);
	for (int i = 1; i < 100; i++){
		explictsolver(Temperature, 0.25);
		printvec(Temperature);
	}

}
void explictsolver(std::vector<double>& Temp, const double r){
	for (int i = 1; i < N - 1; i++ ){
		Temp[i] = Temp[i] + r*(Temp[i-1] - 2.0*Temp[i] + Temp[i+1] );
	}
}

double f(double x){
	return (x <= 0.5) ? (x) : (1 - x);
}
void printvec(const std::vector<double>& Vec){
	for (auto& x : Vec)
		std::cout << std::setw(10) << std::setprecision(5) << x << ' ';
	std::cout << std::endl;
}



/*
void linspacefill(std::vector<double>& Vec, const int Num, const double xi, const double xf){
	auto h = (xf - xi) / Num;
	auto n = Vec.size();
	for (int i = 0; i < n; i++){
		Vec[i] = xi + i*h;
	}
}
//Distribuições espaciais e temporais:
//linspacefill(X, N, 0.0, L);
//linspacefill(t, I, ti, tf);
*/
