#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

#include "gsiedel.hpp"

using std::vector;

void implicit_diff(vector<vector<double>>& A, vector<double>& B, const double r);
void implicit_fic(vector<vector<double>>& A, vector<double>& B, const double r);
void nicolson_diff(vector<vector<double>>& A, vector<double>& B, const double r);
void nicolson_fic(vector<vector<double>>& A, vector<double>& B, const double r);
void resumesaving(std::fstream& printer, const vector<double>& B, const int iter);
void linspacefill(vector<double>& Vec, const int Num, const double xi = 0.0, const double xf = 1.0);
void printvec(const vector<double>& Vec);
void printmatriz(const vector<vector<double>>& A);

// Variáveis do domínio da simulação:
constexpr double L {0.03};                            // comprimento total da placa
constexpr int N {10};                                 // número de nós da malha
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

	vector<vector<double>> A1 (N, std::vector<double>(N, 0.0));
	vector<double> B1 (N, T0);
	vector<vector<double>> A2 (N, std::vector<double>(N, 0.0));
	vector<double> B2 (N, T0);
	vector<vector<double>> A3 (N, std::vector<double>(N, 0.0));
	vector<double> B3 (N, T0);
	vector<vector<double>> A4 (N, std::vector<double>(N, 0.0));
	vector<double> B4 (N, T0);
	std::cout << "Solucao da Equacao Difusivo-Advectiva" << std::endl;
	std::cout << "tf = " << tf << " s, nsteps = " << nsteps << ", dt = " << dt << ", N = " <<  N << ", dx = " << dx << std::endl;

	implicit_diff(A1, B1, r1);
	nicolson_diff(A2, B2, r2);
	implicit_fic(A3, B3, r1);
	nicolson_fic(A4, B4, r2);
}

void implicit_diff(vector<vector<double>>& A, vector<double>& B, const double r){

	std::fstream printer {"Temperatura_Implicit_Diff.dat", std::ios::app};
	printer << "Perfil de Temperatura via solver implicito, contorno via diferencas finitas.\n";
	printer << "i t x T\n";

	int step = 0;
	printer << step << ' ' << ti << ' ';
	for (double x = 0.0; x <= L; x = x+dx)
		printer << x << ' ';
	printer << T0 << '\n';

	// Preenchemos a matriz A (B já foi preenchido com os valores de T0):
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

	// Os passos iterativos do método:
	for (step = 1; step < nsteps; step++){
		// Corrige B:
		B[0] = 0.0;
		for (int i = 1; i < N-1; i++){
			B[i] = B[i] + lambda;
		}
		B[N-1] = (2*h*dx*T0)/kappa;
		// Resolve o sistema:
		GS_Solver(A, B);
	}
	resumesaving(printer, B, step);
	std::cout << "\nEsquema implicito, com contorno de diferencas finitas e  r = " << r << ": " << std::endl;
	printvec(B);
}

void implicit_fic(vector<vector<double>>& A, vector<double>& B, const double r){

	std::fstream printer {"Temperatura_Implicit_Fic.dat", std::ios::app};
	printer << "Perfil de Temperatura via solver implicito, contorno via nos fictios.\n";
	printer << "i t x T\n";

	int step = 0;
	printer << step << ' ' << ti << ' ';
	for (double x = 0.0; x <= L; x = x+dx)
		printer << x << ' ';
	printer << T0 << '\n';

	// Preenchemos a matriz A:
	A[0][0] = 1 + 2*r;
	A[0][1] = -2*r;
	int k = 0;
	for (int i = 1; i < N-1; i++){
			A[i][k] = - r;
			A[i][k+1] = 1 + 2*r;
			A[i][k+2] = -r;
			k++;
		}
	A[N-1][N-2] = -2*r;
	A[N-1][N-1] = 1 + 2*r + (2*r*dx*h)/kappa;

	// Os passos iterativos do método:
	for (step = 1; step < nsteps; step++){
		// Corrige B:
		for (int i = 0; i < N-1; i++){
			B[i] = B[i] + lambda;
		}
		B[N-1] = B[N-1] + lambda + (2*r*h*dx*T0)/kappa;
		// Resolve o sistema:
		GS_Solver(A, B);
	}
	resumesaving(printer, B, step);
	std::cout << "\nEsquema implicito, com contorno via nos ficticios e  r = " << r << ": " << std::endl;
	printvec(B);
}

void nicolson_diff(vector<vector<double>>& A, vector<double>& B, const double r){

	std::fstream printer {"Temperatura_Nicolson_Diff.dat", std::ios::app};
	printer << "Perfil de Temperatura via Crank-Nicolson, contorno via diferencas finitas.\n";
	printer << "i t x T\n";

	int step = 0;
	printer << step << ' ' << ti << ' ';
	for (double x = 0.0; x <= L; x = x+dx)
		printer << x << ' ';
	printer << T0 << '\n';

	// Preenchemos a matriz A (B já foi preenchido):
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

	// Os passos iterativos do método:
	for (step = 1; step < nsteps; step++){
		// Corrige os termos internos de B:
		for (int i = 1; i < N-1; i++){
			B[i] = r*B[i-1] + (1-2*r)*B[i] + r*B[i+1] + lambda;
		}
		// Cprrige os termos no contorno:
		B[0] = 0.0;
		B[N-1] = (2*h*dx*T0)/kappa;
		// Resolve o sistema:
		GS_Solver(A, B);
	}
	resumesaving(printer, B, step);
	std::cout << "\nCrank-Nicolson com contorno de diferencas finitas e r = " << r << ": " << std::endl;
	printvec(B);
}

void nicolson_fic(vector<vector<double>>& A, vector<double>& B, const double r){

	std::fstream printer {"Temperatura_Nicolson_Fic.dat", std::ios::app};
	printer << "Perfil de Temperatura via Crank-Nicolson, contorno via nos fictios.\n";
	printer << "i t x T\n";

	int step = 0;
	printer << step << ' ' << ti << ' ';
	for (double x = 0.0; x <= L; x = x+dx)
		printer << x << ' ';
	printer << T0 << '\n';

	// Preenchemos a matriz A :
	A[0][0] = 1 + 2*r;
	A[0][1] = -2*r;
	int k = 0;
	for (int i = 1; i < N-1; i++){
			A[i][k] = - r;
			A[i][k+1] = 1 + 2*r;
			A[i][k+2] = -r;
			k++;
	}
	A[N-1][N-2] = -2*r;
	A[N-1][N-1] = 1 + 2*r + (2*r*dx*h)/kappa;

	// Os passos iterativos do método:
	for (step = 1; step < nsteps; step++){
		// Corrige os termos internos de B:
		for (int i = 1; i < N-1; i++){
			B[i] = r*B[i-1] + (1-2*r)*B[i] + r*B[i+1] + lambda;
		}
		// Corrige os termos no contorno de B:
		B[0] = (1 - 2*r)*B[0] + (2*r)*B[1] + lambda;
		B[N-1] = (2*r)*B[N-2] + (1 - 2*r - 2*r*dx*h/kappa)*B[N-1] + (4*r*h*dx*T0)/kappa + lambda;
		// Resolve o sistema:
		GS_Solver(A, B);
	}
	resumesaving(printer, B, step);
	std::cout << "\nCrank-Nicolson com contorno via nos ficticios e r = " << r << ": " << std::endl;
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
			std::cout << std::setw(10) << std::setprecision(5) << A[k][w] << ' ';
		std::cout << std::endl;
	}
}

void linspacefill(vector<double>& Vec, const int Num, const double xi, const double xf){
	auto h = (xf - xi) / Num;
	auto n = Vec.size();
	for (int i = 0; i < n; i++){
		Vec[i] = xi + i*h;
	}
}