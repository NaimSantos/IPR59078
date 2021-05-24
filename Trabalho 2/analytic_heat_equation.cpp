#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <functional>
#include <string>

using std::vector;
using std::string;

void analytic_solver(vector<vector<double>>& A, const vector<double>& B, bool fullmatrix);
void eigenvalueTreatment(vector<double>& EVal, const int nroots);
double FourierAdjust(double x, double t, const vector<double>& EVal);
double int_trapz(double a, double b, const double n);
double Tx_xt_0(const double x);
double Tfx(const double x);
double f_intg(const double x, const double value);
double invers_N(const double Bn);
double eigen_value_function(const double Bn);
double solveNewton(double x, std::function<double (double)> f);
double f_diff(double x, double h, std::function<double (double)> f);
void linspace(vector<double>& Vec, const int Num, const double xf = 1.0, const double xi = 0.0);
void printparameters();
void savedata(const vector<vector<double>>& T, const vector<double>& X, const string& filename);
void savedataL2(const vector<double>& T, const string&);
void savetoplot(const vector<vector<double>>& T, const vector<double>& X, const vector<int>& Tpoints, const string& filename);
void analytic_solverL2(vector<double>& T);

// Variáveis do domínio da simulação:
constexpr double L {0.03};                                   // comprimento total da placa
constexpr double ti {0.0};                                   // tempo inicial da simulação
constexpr double tf {10000.0};                                // tempo final da simulação
constexpr auto dt {10};                                      // passo de tempo na série de Fourier
constexpr auto nsteps = static_cast<int>((tf-ti)/dt) + 1;    // número de passos de tempo usando a série de Fourier
constexpr auto dx {0.0001};                                   // comprimento do intervalo na série de Fourier
constexpr auto Npoints = L/dx + 1;                           // número de pontos avaliados para x 

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
	vector<double> X(Npoints, 0.0);
	linspace(X, Npoints, L, 0.0);            // escreve em X uma distribuição dos Npoints pontos de avaliação

	// Chamada com todos os ts (para gerar a distribuição):
	vector<vector<double>> T (nsteps, std::vector<double>(Npoints, 0.0)); // linhas para tempo, colunas para posições
	//analytic_solver(T, X, true);
	//const string nome {"dados_totais.dat"};
	savedata(T, X, nome);
	
	// Chamada para valores específicos:
	// Um vetor para os tempos: 0, 10, 50, 100, 500, 750, 1000 e 20000
	vector<vector<double>> T2 (8, std::vector<double>(Npoints, 0.0));
	analytic_solver(T2, X, false);
	
	//Chamada para l/2
	vector<double> TL2 (nsteps, 0);
	analytic_solverL2(TL2);


	
}
void analytic_solverL2(vector<double>& T){
	// Calcular os autovalores:
	int nroots = 5;                                 // número de autovalores (truncamento da série de Fourier)
	vector<double> EVal (nroots, 0.0);               // array onde nroots autovalores serão armazenados.
	eigenvalueTreatment(EVal, EVal.size());          // calcular os autovalores

	// Solução em L/2
	for (int i = 0; i < nsteps; i++){
		T[i] = FourierAdjust(L/2, i*dt, EVal);
		// A Solução T(x,t) é a soma das soluções (filtrada + filtro): T(x,t) = T*(x,t) + Tf(x)
		T[i] = T[i] + Tfx(L/2);
	}
	const string nome_f {"dados_L2.dat"};
	savedataL2(T, nome_f);
}

void analytic_solver(vector<vector<double>>& T, const vector<double>& X, bool fullmatrix){
	// Calcular os autovalores:
	int nroots = 10;                                 // número de autovalores (truncamento da série de Fourier)
	vector<double> EVal (nroots, 0.0);               // array onde nroots autovalores serão armazenados.
	eigenvalueTreatment(EVal, EVal.size());          // calcular os autovalores


	// A Solução T(x,t) é a soma das soluções (filtrada + filtro): T(x,t) = T*(x,t) + Tf(x)
	static int count = 0;
	count++;
	std::cout << "\nAnalytic Solver called " << count << std::endl;
	if (fullmatrix){
		// Solução em todos os tempos
		for (int i = 0; i < nsteps; i++){
			for (int j = 0; j < Npoints; j++){
				T[i][j] = FourierAdjust(j*dx, i*dt, EVal);
				// A Solução T(x,t) é a soma das soluções (filtrada + filtro): T(x,t) = T*(x,t) + Tf(x)
				T[i][j] = T[i][j] + Tfx(j*dx);
			}
		}
	}
	else {
		// Avaliar a solução via ajuste de Fourier nos tempos: 0, 10, 50, 100, 500, 750 e 1000:
		std::cout << "\nSolucao para um unico tempo" << std::endl; 
		// Solução apenas para tempos requeridos:
		//vector<int> Tpoints {0, 100};
		vector<int> Tpoints {0, 10, 50, 100, 500, 750, 1000, 20000};
		auto tam = Tpoints.size();
		
		for (int i = 0; i < tam; i++){
			auto t = Tpoints[i];
			//std::cout << "Tpoints[i] = " << t << std::endl;
			for (int j = 0; j < Npoints; j++){
				T[i][j] = FourierAdjust(j*dx, t, EVal);
				T[i][j] = T[i][j] + Tfx(j*dx);
			}
		}
		string nome {"dados_parciais_v2.dat"};
		savetoplot(T, X, Tpoints, nome);
	}	
}
void eigenvalueTreatment(vector<double>& EVal, const int nroots){
	EVal[0] = solveNewton(20, eigen_value_function);
	for (int i = 1; i < nroots; i++){
		EVal[i] =  solveNewton(i*100, eigen_value_function);
	}
	std::cout << "\n" << EVal.size() << " autovalores obtidos: " ;
	for (const auto& x : EVal) std::cout << std::setprecision(8) << x << ' ';
}

double FourierAdjust(double x, double t, const vector<double>& EVal){
	double res {0.0}, res1{0.0}, res2{0.0}, res3{0.0}, res4 {0.0};
	auto total = EVal.size();                        // número de autovalores usados

	for (int n = 0; n < total; n++){
		auto eigval = EVal[n];
		res1 = invers_N(eigval);
		res2 = std::exp(- alpha * eigval * eigval * t);
		res3 = std::cos(x * eigval);
		res4 = int_trapz(0, L, eigval);
		res += (res1 * res2 * res3 * res4);
	}
	//std::cout <<total << " autovalores, x = " << x << ", t = " << t << ", res = "<< res << std::endl;
	
	return res;
}

// Integração numérica pela regra do trapézio
double int_trapz(double a, double b, const double value){
	double res {0.0};
	const double h = 0.0001;                         // passo da integração numérica
	const auto m = static_cast<int>( std::floor((std::fabs(b - a)) / h));

	for (int k = 0; k < m - 1; k++){
		res += f_intg(a + k*h, value);
	}
	res += (f_intg(a, value) + f_intg(b, value)) / 2;
	res *= h;
	return res;
}

// Função inicial do problema:
double Tx_xt_0(const double x){
	return (g*x*x)/(2*kappa) - (g*L*L)/(2*kappa) - g*L/h;
}

// Função na integral dos coeficientes:
double f_intg(const double x, const double value){
	return Tx_xt_0(x)*std::cos(x*value);
}

// Solução filtro:
double Tfx(const double x){
	return T_inf -(g*x*x)/(2*kappa) + (g*L*L)/(2*kappa) + (g*L)/h;
}

// Coeficiente 1/N:
double invers_N(const double Bn){
	return 2 * ((Bn*Bn + H2*H2) / (L*(Bn*Bn + H2*H2) + H2));
}

// Função cujas raízes são autovalores:
double eigen_value_function(const double Bn){
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
	std::cout << "L = " << L << " dx = " << dx << " Npoints = " << Npoints;
}

void savedata(const vector<vector<double>>& T, const vector<double>& X, const string& filename){
	std::cout << "\nSaveData function called" << std::endl;
	// Salvar os resultados em um arquivo:
	std::fstream printer {filename, std::ios::out|std::ios::trunc};
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
void savedataL2(const vector<double>& T, const string& filename){
	std::cout << "\nFunction SaveDataL2 called" << std::endl;
	// Salvar os resultados em um arquivo:
	std::fstream printer {filename, std::ios::out|std::ios::trunc};
	printer << "Temperatura em x = L/2\n";
	printer << "t T";
	double t = 0.0;
	for (int i = 0; i < nsteps; i++){
		printer <<"\n" << t <<' ' << T[i];
		t+=dt;
	}
}
void savetoplot(const vector<vector<double>>& T, const vector<double>& X, const vector<int>& Tpoints, const string& filename){

	std::cout << "\nSaving data to plot..." << std::endl;
	std::fstream printer {filename, std::ios::out|std::ios::trunc};

	printer << "Perfil de Temperatura via solucao analitica\n";
	printer << "t = ";
	auto amount = Tpoints.size();
	for (int w = 0; w < amount; w++){
		printer << " " << Tpoints[w];
	}
	printer << "\nPoint Temperature\n";
	
	auto m = (T[0]).size();	
	auto n = T.size();

	for (int i = 0; i < m; i++){
		printer << X[i] << ' ';
		for (int j = 0; j < n; j++){
			printer << T[j][i] << ' ';
		}
		printer << "\n";
	}

}
