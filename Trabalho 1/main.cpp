#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

#include "gsiedel.hpp"

using std::vector;

void printvec(const vector<double>& Vec);
void implicit(vector<vector<double>>& A, vector<double>& B, const double r);
void resumesaving(std::fstream& printer, const vector<double>& B, const int iter);
void printmatriz(const vector<vector<double>>& A);

// Variáveis do domínio da simulação:
constexpr double L {0.03};                            // comprimento total da placa
constexpr int N {100};                                 // número de nós da malha
constexpr double ti {0.0};                            // tempo inicial da simulação
constexpr double tf {500.0};                          // tempo final da simulação
constexpr double dt {0.1};                            // passo de tempo
constexpr auto dx { L / (N - 1)};                      // comprimento do intervalo
constexpr auto nsteps = static_cast<int>((tf-ti)/dt);  // número de passos de tempo

// Dados do problema:
constexpr double kappa {0.6};
constexpr int rho {600};
constexpr int cp {1200};
constexpr double h {15.0};
constexpr int g {100000};
constexpr double T0 {20.0};
constexpr double TL {20.0};
constexpr auto alpha = kappa/(rho*cp);
constexpr auto r1 = (alpha*dt)/(dx*dx);               // coeficiente r do método implícito
constexpr auto r2 = (alpha*dt)/(2*dx*dx);             // coeficiente r do Crank-Nicolson
constexpr auto gamma = 3.0 + ((2*h*dx)/kappa);
constexpr auto lambda = (g*dt)/(rho*cp);


int main (int argc, char* argv[]){
	
	vector<vector<double>> A (N, std::vector<double>(N, 0.0));
	vector<double> B (N, T0);


	implicit(A, B, r1);
}

void implicit(vector<vector<double>>& A, vector<double>& B, const double r){
	
	std::fstream printer {"DadosDeTemperatura.dat", std::ios::app};
	printer << "Perfil de Temperatura via solver implicito.\n";
	printer << "i t x T\n";
	
	int step = 0;
	double t = ti;
	printer << step << ' ' << t << ' ';
		for (double x = 0.0; x <= L; x=x+dx)
			printer << x << ' ';
	printer << T0 << '\n';
	
	//Preenchemos a matriz A (B já foi preenchido):
	A[0][0] = -3.0;
	A[0][1] = 4.0;
	A[0][2] = -1.0;
	int k = 0;
	for (int i = 1; i < N-1; i++){
			A[i][k] = - r;
			A[i][k+1] = 1 + 2*r;
			A[i][k+2] = -r;
			k++;
		}
    A[N-1][N-3] = 1.0;
	A[N-1][N-2] = -4.0;
	A[N-1][N-1] = gamma;

	//printmatriz(A);
	//Os passos iterativos do método:
	for (step = 1; step < nsteps; step++){
		// Corrige B:
		B[0] = 0.0;
        for (int i = 1; i < N-1; i++){
            B[i] = B[i] + lambda;
		}
		B[N-1] = (2*h*dx*T0)/kappa;
		//Resolve o sistema:
		GS_Solver(A, B);
		//resumesaving(printer, B, step);
	}
	resumesaving(printer, B, step);
	std::cout << "\nSolucao obtida: " << std::endl;
	printvec(B);
}

void printvec(const vector<double>& Vec){
	for (auto& x : Vec)
		std::cout << std::setw(10) << std::setprecision(8) << x << ' ';
	std::cout << std::endl;
}
void resumesaving(std::fstream& printer, const vector<double>& B, const int iter){
	printer << iter << ' ' << iter*dt << ' ';
	for (int i = 0; i < N; i++){
		printer << B[i] << ' ';
	}
	printer << '\n';
}
void printmatriz(const vector<vector<double>>& A){
	for (int k=0; k < N; k++){
		for (int w = 0; w < N; w++)
			std::cout <<  std::setw(10) << std::setprecision(5)  << A[k][w] << ' ';
		std::cout << std::endl;
	}
}
/*
	Teste:

	vector<vector<double>> A { {2, 1, 1},
							   {3, 5, 2},
							   {2, 1, 4} };
	vector<double> B {5, 15, 8};
	std::cout << "Antes do Gauss Siedel" << std::endl;
	printvec(B);
	
	GS_Solver(A, B);
	std::cout << "\nDepois do Gauss Siedel" << std::endl;
	printvec(B);
*/