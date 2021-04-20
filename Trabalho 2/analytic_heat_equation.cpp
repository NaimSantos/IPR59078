#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

#include "utilities.hpp"	//Gauss-Siedel, registro de tempo e integração numérica

using std::vector;

void analytic_solver(vector<vector<double>>& A, const vector<double>& B);
void nicolson_diff(vector<vector<double>>& A, vector<double>& B, const double r);
void resumesaving(std::fstream& printer, const vector<double>& B, const int iter);
void linspace(vector<double>& Vec, const int Num, const double xf = 1.0, const double xi = 0.0);
void printvec(const vector<double>& Vec);
void printmatriz(const vector<vector<double>>& A);
double f1(double x);
double f2(double x);
double Txt(double x, double t);
double int_trapz(double a, double b, const double i);
double f_p1(const double x, const int i);
double f_p2(const double x, const int i);

// Variáveis do domínio da simulação:
constexpr double L {1};                                 // comprimento total da placa
constexpr int N {17};                                      // número de nós da malha
constexpr double ti {0.0};                                 // tempo inicial da simulação
constexpr double tf {500.0};                               // tempo final da simulação
constexpr auto dx = L/(N-1);                               // comprimento do intervalo
constexpr auto dt {0.005};                                 // passo de tempo
constexpr auto nsteps = static_cast<int>((tf-ti)/dt);      // número de passos de tempo
constexpr auto dt_f {0.1};                                 // passo de tempo na série de Fourier
constexpr auto nsteps_f = static_cast<int>((tf-ti)/dt_f);  // número de passos de tempo usando a série de Fourier
constexpr int NF {5};                                      // número de elementos na série de Fourier
constexpr auto dx_f = L/(NF-1);                            // comprimento do intervalo na série de Fourier
constexpr auto N_i {L/2};                                  // coeficiente na série

// Dados do problema:
constexpr double kappa {0.6};
constexpr int rho {600};
constexpr int cp {1200};
constexpr double h {15.0};
constexpr int g {100000};
constexpr double T0 {0.0};
constexpr double TL {0.0};
constexpr auto alpha = kappa/(rho*cp);
constexpr auto r2 = (alpha*dt)/(2*dx*dx);             // coeficiente r do Crank-Nicolson
constexpr auto eta = 3.0 + ((2*h*dx)/kappa);
constexpr auto sigma = (g*dt)/(rho*cp);


int main (int argc, char* argv[]){
	
	//solução Analítica
	vector<vector<double>> C (nsteps_f, std::vector<double>(NF, 0.0));
	vector<double> D(NF, 0.0);
	linspace(D, NF, L, 0.0);

	analytic_solver(C, D);
	//Crank-Nicolson:
	/*
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
	*/
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

	CustomTimer TimerNC;
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

void analytic_solver(vector<vector<double>>& T, const vector<double>& X){
	//A is a 2D-vector: rows for time, collumns for position
	//B is a 1D-vector with the positions
	
	for (int k = 0; k < NF; k++)
		T[0][k] = f1(X[k]);
	
	std::fstream printer {"Temperatura_Analitica.dat", std::ios::out|std::ios::trunc};
	printer << "Perfil de Temperatura via solucao analitica\n";
	printer << "t";
	for (int k = 0; k < NF; k++)
		printer << " X" << k;
	printer <<"\n ";
	for (int k = 0; k < NF; k++)
		printer << ' ' << X[k] ;
	printer <<"\n ";
	for (int k = 0; k < NF; k++)
		printer << ' ' <<T[0][k];
	printer <<"\n ";
	
	
	for (int i = 1; i < N; i++){
		for (int j = 0; j < NF; j++)
			T[i][j]=Txt(j*dx_f, i*dt);
	}
	/*
	for (int i = 1; i <= N; i++){
		auto intg = int_trapz(0, L, 0.00001, f);
		auto total {0.0};
		total = std::exp(-alpha*(mu_i*mu_i) * t) * std::sin(mu_i) * intg / N_i;
	}
	*/
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

double f1(double x){
	return (x <= 0.5*L) ? (x) : (L - x);
}
double f2(double x){
	auto res = std::sin(x);
	return (x <= 0.5) ? (x * res) : ((1 - x)*res);
}
//mu_i = i * pi / L
double Txt(double x, double t){
	double res {0.0}; double res1{0.0}; double res2{0.0}; double res3{0.0};

	for (int i = 1; i <= 5; i++){
		res1 = exp(- alpha * t * (std::pow(i*NPI/L, 2)));
		res2 = std::sin(x*i*NPI/L);
		res3 = int_trapz(0, L, i);
	}
	return res;
}
double fxsenx(double x, double i){
	if (x <= 0.5*L) 
		return x*std::sin(x*i*NPI/L);
	else
		return (L - x)*std::sin(x*i*NPI/L);
		
}
// Integração numérica pela regra do trapézio
double int_trapz(double a, double b, const double i){
	const double h = 0.000001;           // passo
	const auto L2 = (b - a) /2;	        // metade do intervalo (L/2)
	
	// Primeira parte da integral, de 0 a L/2:
	double res1 {0.0};
	const auto n = static_cast<int>( std::floor((std::fabs(L2 - a)) / h));
	for (int k = 0; k < n - 1; k++){
		res1 += f_p1(a + k*h, i);
	}
	res1 += (f_p1(a, i) + f_p1(L2, i) ) / 2;
	res1 *= h;

	// Segunda parte da integral, de L/2 a L:
	double res2 {0.0};
	const auto m = static_cast<int>( std::floor((std::fabs(b - L2)) / h));	
	for (int k = 0; k < m - 1; k++){
		res2 += f_p2(L2 + k*h, i);
	}
	res2 += (f_p2(L2, i) + f_p2(b, i) ) / 2;
	res2 *= h;
	return res1 + res2;
}

double f_p1(const double x, const int i){
	return x*std::sin(i*NPI*x/L);
}
double f_p2(const double x, const int i){
	return (L - x)*std::sin(i*NPI*x/L);
}