#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

using std::vector;

#include "tdma.hpp"

double f(double x);
void explictsolver(vector<double>& Temp, const double r = 0.5, const int tempo = 10);
void implicitsolver(vector<double>& U, const vector<double>& a, const vector<double>& b, const vector<double>& c, const int tempo = 10);
void printvec(const vector<double>& Vec);
void reset_first_and_last(vector<double>& Vec);

//Dados do domínio:
constexpr double dx {0.1};              //refinamento da malha espacial
constexpr double dt {1};                //passo de tempo
constexpr double L {1.0};               //comprimento total do domínio
constexpr double ti {0.0};              //tempo inicial da simulação
constexpr double tf {100.0};            //tempo final de simulação
constexpr auto N  = L/dx + 1;           //número de nós da malha (intervalos + 1)
constexpr double I = (tf - ti) / dt;    //número de passos de tempo

//Dados do problema:
constexpr double kappa {0.6};
constexpr double rho {600};
constexpr double cp {1200};
constexpr double h {15.0};
constexpr auto g {100000};
constexpr double T_zero {20.0};

int main (int argc, char* argv[]){

	const auto alfa = kappa / (rho * cp);
	const auto r = (alfa * dt) / (dx * dx);
	//std::cout << "dt = " << dt << ", dx = " << dx  << ", alfa = " << alfa << ", CFL (r) = " << r << std::endl;
	const auto rtest = 0.5;

	// -----------------Solução 1 ---------------------//
	vector<double> Temperature(N, 0.0);

	//Preencher os temperaturas iniciais a partir da função fornecida no problema:
	auto n = Temperature.size();
	for (int i = 0; i < n; i++){
		Temperature[i] = f(i*dx);
	}
	auto T2 = Temperature; //cópia do vetor para o esquema implícito

	//Invocar o solver explícito:
	std::cout << "\nSolucao explicita:\n";
	explictsolver(Temperature, rtest, 10);


	// -----------------Solução 2 ---------------------//
	vector<double> Dsup (N, (-rtest));          //Diagonal superior, inicializada em -r
	vector<double> Dmain (N, (1 + 2*rtest));    //Diagonal principal, inicializada em 1 + 2r
	vector<double> Dinf (N, (-rtest));          //Diagonal inferior, inicializada em -r


	//Corrigimos os termos cuja inicialização é diferente:
	Dmain[0] = Dmain[n-1] = 1.0;
	Dsup[0] = Dinf[n-1] = 0.0;

	//Invoca o solver implicito:
	std::cout << "\nSolucao implicita:\n";
	implicitsolver(T2, Dsup, Dmain, Dinf, 10);
}

double f(double x){
	return (x <= 0.5) ? (x) : (1 - x);
}

void explictsolver(vector<double>& Temp, const double r, const int tempo){
	printvec(Temp);
	for (int iter = 1; iter <= tempo; iter++){
		for (int i = 1; i < N - 1; i++ ){
			Temp[i] = Temp[i] + r*(Temp[i-1] - 2.0*Temp[i] + Temp[i+1] );
		}
		printvec(Temp);
	}
}

void implicitsolver(vector<double>& U, const vector<double>& a, const vector<double>& b, const vector<double>& c, const int tempo){
	printvec(U);
	for (int iter = 1; iter <= tempo; iter++){
		tdma_solver(a, b, c, U);
		reset_first_and_last(U);
		printvec(U);
	}
}

void printvec(const vector<double>& Vec){
	for (auto& x : Vec)
		std::cout << std::setw(10) << std::setprecision(5) << x << ' ';
	std::cout << std::endl;
}

void reset_first_and_last(vector<double>& V){
	V[0] = V[V.size()-1] = 0.0;
}
