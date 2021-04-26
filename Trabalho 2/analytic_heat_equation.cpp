#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

using std::vector;

void analytic_solver(vector<vector<double>>& A, const vector<double>& B);
void linspace(vector<double>& Vec, const int Num, const double xf = 1.0, const double xi = 0.0);
void printparameters();
void savedata(const vector<vector<double>>& T, const vector<double>& X);
double f1(double x);
double f2(double x);
double FourierAdjust(double x, double t);
double int_trapz(double a, double b, const double i);
double f_p1(const double x, const int i);
double f_p2(const double x, const int i);

// Estimativa para pi:
constexpr auto NPI {4*std::atan(1)};

// Variáveis do domínio da simulação:
constexpr double L {1.0};                                  // comprimento total da placa
constexpr double ti {0.0};                                 // tempo inicial da simulação
constexpr double tf {50.0};                               // tempo final da simulação
constexpr auto dt {0.1};                                   // passo de tempo na série de Fourier
constexpr auto nsteps = static_cast<int>((tf-ti)/dt) + 1;      // número de passos de tempo usando a série de Fourier
constexpr int N {5};                                       // número de elementos na série de Fourier
constexpr auto dx {0.1};                                  // comprimento do intervalo na série de Fourier
constexpr auto N_i = L/2;                                  // coeficiente na série
constexpr auto Npoints = L/dx + 1;

// Dados do problema:
constexpr double kappa {0.6};
constexpr int rho {600};
constexpr int cp {1200};
constexpr double h {15.0};
constexpr int g {100000};
constexpr double T0 {0.0};
constexpr double TL {0.0};
constexpr auto alpha = kappa/(rho*cp);


int main (int argc, char* argv[]){

	printparameters();
	// Solução Analítica
	vector<vector<double>> T (nsteps, std::vector<double>(Npoints, 0.0));
	vector<double> X(Npoints, 0.0);
	linspace(X, Npoints, L, 0.0);       // escreve em D uma distribuição dos N pontos de avaliação

	analytic_solver(T, X);
}

void analytic_solver(vector<vector<double>>& T, const vector<double>& X){
	//A is a 2D-vector: rows for time, collumns for position
	//B is a 1D-vector with the positions

	for (int i = 0; i < nsteps; i++){
		for (int j = 0; j < Npoints; j++){
			//std::cout << "i = " << i << ", j = " << j <<std::endl;
			T[i][j]=FourierAdjust(j*dx, i*dt);
		}
	}

	// Salvar os resultados em um arquivo:
	savedata(T, X);
}

double f1(double x){
	return (x <= 0.5*L) ? (x) : (L - x);
}

double f2(double x){
	auto res = std::sin(x);
	return (x <= 0.5) ? (x * res) : ((1 - x)*res);
}

// mu_i = i * pi / L
double FourierAdjust(double x, double t){
	double res {0.0};
	double res1{0.0};
	double res2{0.0};
	double res3{0.0};

	for (int i = 1; i <= N; i++){
		res1 = std::exp(- alpha * t * (std::pow(i*NPI/L, 2)));
		res2 = std::sin(x*i*NPI/L);
		res3 = int_trapz(0, L, i);
		res += (1/N_i)*(res1 * res2 * res3);
	}
	return res;
}


// Integração numérica pela regra do trapézio
double int_trapz(double a, double b, const double i){
	//std::cout << "Integral evaluation. a=" << a << ", b=" << b << ", i=" << i << std::endl;
	const double h = 0.0001;           // passo
	const auto L2 = (b - a)/2;          // metade do intervalo (L/2)

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
	return 10*x*std::sin(i*NPI*x/L);
}

double f_p2(const double x, const int i){
	return 10*(L - x)*std::sin(i*NPI*x/L);
}

void linspace(vector<double>& Vec, const int Num, const double xf, const double xi){
	auto h = (xf - xi) / (Num-1);
	auto n = static_cast<int>(Vec.size());
	for (int i = 0; i < n; i++){
		Vec[i] = xi + i*h;
	}
}

void printparameters(){	
	std::cout << "L = " << L << "\tdx = " << dx << "\tNpoints = " << Npoints
				<< "\ntempo total = " << tf << "\tdt = " << dt	<< "\tnsteps = "
				<< nsteps << "\nElementos na serie (N) = " << N	<<  std::endl;
}				

void savedata(const vector<vector<double>>& T, const vector<double>& X){
	std::cout << "SaveData function called" << std::endl;
	// Salvar os resultados em um arquivo:
	std::fstream printer {"Temperatura_Analitica.dat", std::ios::out|std::ios::trunc};
	printer << "Perfil de Temperatura via solucao analitica\n";
	printer << "t";
	for (int k = 0; k < Npoints; k++)
		printer << " X" << k;
	printer <<"\n-";
	for (int k = 0; k < Npoints; k++)
		printer << ' ' << X[k] ;
	double t = 0.0;
	for (int i = 0; i < nsteps; i++){
		printer <<"\n" << t;
		for (int j = 0; j < Npoints; j++){
			printer <<' ' << T[i][j];
		}
		t+=dt;
		
	}
	
}