#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <functional>

using std::vector;

void analytic_solver(vector<vector<double>>& A, const vector<double>& B);
double FourierAdjust(double x, double t, const vector<double>& EVal);
double int_trapz(double a, double b, const double i);
double f_init(const double x, const double i, bool onintegral = false);
double Tfx(const double x);
double invers_N(const double Bn);
double eigenvalue(const double Bn);
double solveNewton(double x, std::function<double (double)> f);
double f_diff(double x, double h, std::function<double (double)> f);
void linspace(vector<double>& Vec, const int Num, const double xf = 1.0, const double xi = 0.0);
void printparameters();
void savedata(const vector<vector<double>>& T, const vector<double>& X);


constexpr auto NPI = 4*std::atan(1);

// Variáveis do domínio da simulação:
constexpr double L {0.03};                                   // comprimento total da placa
constexpr double ti {0.0};                                   // tempo inicial da simulação
constexpr double tf {50.0};                                  // tempo final da simulação
constexpr auto dt {0.1};                                     // passo de tempo na série de Fourier
constexpr auto nsteps = static_cast<int>((tf-ti)/dt) + 1;    // número de passos de tempo usando a série de Fourier
constexpr auto dx {0.0025};                                  // comprimento do intervalo na série de Fourier
constexpr auto Npoints = L/dx + 1;

// Dados do problema:
constexpr double kappa {0.6};
constexpr int rho {600};
constexpr int cp {1200};
constexpr double h {15.0};
constexpr int g {100000};
constexpr double T_inf {20.0};
constexpr auto alpha = kappa/(rho*cp);
constexpr auto H2 = h/kappa;

int main (int argc, char* argv[]){
	printparameters();

	// Solução Analítica
	vector<vector<double>> T (nsteps, std::vector<double>(Npoints, 0.0));
	vector<double> X(Npoints, 0.0);
	linspace(X, Npoints, L, 0.0);       // escreve em X uma distribuição dos Npoints pontos de avaliação

	analytic_solver(T, X);
}

void analytic_solver(vector<vector<double>>& T, const vector<double>& X){
	// A is a 2D-vector: rows for time, collumns for position
	// B is a 1D-vector with the positions

	// Calcular os autovalores:
	int nroots = 10;                        // número de auto vetores (truncamento da série de Fourier)
	vector<double> EVal (nroots, 0.0);      // array onde nroots autovalores serão armazenados.

	EVal[0] = solveNewton(20, eigenvalue);
	for (int i = 1; i < nroots; i++){
		EVal[i] =  solveNewton(i*100, eigenvalue);
	}
	//for (const auto& x : EVal) std::cout << std::setprecision(6) << x << ' ';

	for (int i = 0; i < nsteps; i++){
		for (int j = 0; j < Npoints; j++){
			T[i][j] = FourierAdjust(j*dx, i*dt, EVal);
			// A Solução T(x,t) é a soma das soluções (filtrada + filtro): T(x,t) = T*(x,t) + Tf(x)
			T[i][j] = T[i][j] + Tfx(j*dx);
		}
	}
	// Salvar os resultados em um arquivo:
	savedata(T, X);
}

double FourierAdjust(double x, double t, const vector<double>& EVal){
	double res {0.0}, res1{0.0}, res2{0.0}, res3{0.0}, res4 {0.0};
	auto m = EVal.size(); // número de autovalores usados;

	for (int j = 0; j < m; j++){
		auto i = EVal[j];
		res1 = invers_N(i);
		res2 = std::exp(- alpha * std::pow(i*NPI/L, 2) * t);
		res3 = std::cos(x*i*NPI/L);
		res4 = int_trapz(0, L, i);

		res += (res1 * res2 * res3 * res4);
	}
	return res;
}

// Integração numérica pela regra do trapézio
double int_trapz(double a, double b, const double i){
	double res {0.0};
	const double h = 0.0001;                // passo da integração numérica
	const auto n = static_cast<int>( std::floor((std::fabs(b - a)) / h));

	for (int k = 0; k < n - 1; k++){
		res += f_init(a + k*h, i, true);
	}
	res += (f_init(a, i, true) + f_init(b, i, true) ) / 2;
	res *= h;
	return res;
}

// Função inicial do problema:
double f_init(const double x, double i, bool onintegral){
	double res = T_inf - ( (-g*x*x)/(2*kappa) + T_inf + (g*L*L)/(2*kappa) + g*L/h);
	if (onintegral)
		return res*std::cos(x*i*NPI/L);
	else
		return res;
}

double Tfx(const double x){
	return (-g*x*x)/(2*kappa) + T_inf + (g*L*L)/(2*kappa) + g*L/h;
}

double invers_N(const double Bn){
	return 2 * (Bn*Bn + H2*H2) / (L* (Bn*Bn + H2*H2) + H2);
}

double eigenvalue(const double Bn){
	return (Bn*std::tan(Bn*L) - H2);
}

// Cálculo de raiz via Newton-Rhaphson
double solveNewton(double x, std::function <double (double)> f){
	auto h = f(x) / f_diff(x, 0.001, f);
	unsigned int i {0};
	while (std::fabs(h) >= 0.00001 && i<100){
		h = f(x) / f_diff(x, 0.001, f);
		x = x - h;
		i++;
	}
	return x;
}

// Diferenças finitas centras para a derivada no Newton-Rhaphson
double f_diff(double x, double h, std::function<double (double)> f){
	double res = (-f(x+2*h) + 8*f(x+h) - 8*f(x-h) + f(x-2*h)) / (12*h);
	return res;
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
		<< "\ntempo total = " << tf << "\tdt = " << dt << "\tnsteps = "
		<< nsteps << std::endl;
}

void savedata(const vector<vector<double>>& T, const vector<double>& X){
	std::cout << "\nSaveData function called" << std::endl;
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