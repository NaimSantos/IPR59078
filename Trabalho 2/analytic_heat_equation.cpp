#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

#include "utilities.hpp"	//Gauss-Siedel e registro de tempo

using std::vector;

void nicolson_diff(vector<vector<double>>& A, vector<double>& B, const double r);
void resumesaving(std::fstream& printer, const vector<double>& B, const int iter);
void linspace(vector<double>& Vec, const int Num, const double xf = 1.0, const double xi = 0.0);
void printvec(const vector<double>& Vec);
void printmatriz(const vector<vector<double>>& A);

// Variáveis do domínio da simulação:
constexpr double L {0.03};                             // comprimento total da placa
constexpr int N {17};                                  // número de nós da malha
constexpr double ti {0.0};                             // tempo inicial da simulação
constexpr double tf {500.0};                          // tempo final da simulação
constexpr auto dx  = L/(N-1);                          // comprimento do intervalo
constexpr int nsteps {2048*4};                          // número de passos de tempo
constexpr auto dt = (tf-ti)/nsteps;                    // passo de tempo
//constexpr auto dt {0.0025};                            // passo de tempo
//constexpr auto nsteps = static_cast<int>((tf-ti)/dt);  // número de passos de tempo

// Dados do problema:
constexpr double kappa {0.6};
constexpr int rho {600};
constexpr int cp {1200};
constexpr double h {15.0};
constexpr int g {100000};
constexpr double T0 {20.0};
constexpr double TL {20.0};
constexpr auto alpha = kappa/(rho*cp);
//constexpr auto r1 = (alpha*dt)/(dx*dx);               // coeficiente r do método implícito
constexpr auto r2 = (alpha*dt)/(2*dx*dx);             // coeficiente r do Crank-Nicolson
constexpr auto eta = 3.0 + ((2*h*dx)/kappa);
constexpr auto sigma = (g*dt)/(rho*cp);


int main (int argc, char* argv[]){

	vector<vector<double>> A (N, std::vector<double>(N, 0.0));
	vector<double> B (N, T0);

	std::cout << "Solucao da Equacao Difusivo-Advectiva" << std::endl;
	std::cout << "tf = " << tf << " s, nsteps = " << nsteps << ", dt = " << dt << ", N = " << N << ", dx = " << dx << std::endl;

	nicolson_diff(A, B, r2);

	// Salvar em um arquivo a comparação entre todos os métodos:
	vector<double> X(N, 0.0);
	linspace(X, N, L, 0.0);
	std::fstream allprint {"dados.dat", std::ios::out|std::ios::trunc};
	allprint << "tf = " << tf << " s, nsteps = " << nsteps << ", dt = " << dt << ", N = " << N << ", dx = " << dx << std::endl;
	allprint << "X CN_DIF\n";
	for (int i = 0; i < N; i++){
		allprint << X[i] << ' ' << B[i] << '\n';
	}
}

void nicolson_diff(vector<vector<double>>& A, vector<double>& B, const double r){

	std::fstream printer {"Temperatura_Nicolson_Diff.dat", std::ios::out|std::ios::trunc};
	printer << "Perfil de Temperatura via Crank-Nicolson, contorno via diferencas finitas.\n";
	printer << "i t x T\n";

	int step = 0;
	printer << step << ' ' << ti << ' ';
	for (double x = 0.0; x <= L; x = x+dx)
		printer << x << ' ';
	printer << T0 << '\n';

	// Preenchemos a matriz A:
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
	A[N-1][N-1] = eta;

	CustomTimer time_cn_dif;
	// Processo Iterativo do Esquema Crank-Nicolson:
	for (step = 1; step < nsteps; step++){
		//Corrige B, incializado em T0
		for (int i = 1; i < N-1; i++){
			B[i] = r*B[i-1] + (1-2*r)*B[i] + r*B[i+1] + sigma;
		}
		B[0] = 0.0;
		B[N-1] = (2*h*dx*T0)/kappa;

		// Resolve o sistema:
		GS_Solver(A, B);
	}
	resumesaving(printer, B, step);
	std::cout << "\nCrank-Nicolson com contorno de diferencas finitas e r = " << r << ": " << std::endl;
	//printvec(B);
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
			std::cout << std::setw(10) << std::setprecision(7) << A[k][w] << ' ';
		std::cout << std::endl;
	}
}

void linspace(vector<double>& Vec, const int Num, const double xf, const double xi){
	auto h = (xf - xi) / (Num-1);
	auto n = static_cast<int>(Vec.size());
	for (int i = 0; i < n; i++){
		Vec[i] = xi + i*h;
	}
}
