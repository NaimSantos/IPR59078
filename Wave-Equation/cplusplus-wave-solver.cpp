#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <functional>
#include <string>

using std::vector;

void print_parameters();
double int_trapz(double a, double b, const double n);
void fill_eingenvalues(vector<double>& Values, const int total);
double fourier_adjust(double x, double t, const vector<double>& Values);
double target_function(const double x, const int n);
double coef_an(const int n);
double int_trapz(double a, double b, const int n);
void linspace(vector<double>& Vec, const int Num, const double xf = 1.0, const double xi = 0.0);
void save_data(const vector<vector<double>>& U, const vector<double>& X, const std::string& filename);
void save_plot(const vector<vector<double>>& U, const vector<double>& X, const std::string& filename);

constexpr double beta {0.5};
constexpr double L {1.0};
constexpr double ti {0.0};                                  // tempo inicial da simulação
constexpr double tf {5.0};                                  // tempo final da simulação
constexpr auto dt {1};                                      // passo temporal
constexpr auto nsteps = static_cast<int>((tf-ti)/dt) + 1;   // número de passos de tempo
constexpr auto dx {0.01};                                   // passo espacial
constexpr auto npoints = L/dx + 1;                          // número de pontos avaliados para x
constexpr int N {5};                                        // número de termos no somatório da Série de Fourier
constexpr auto NPI = 4*std::atan(1);

int main (int argc, char* argv[]){

	print_parameters();
	vector<vector<double>> U_xt (nsteps, std::vector<double>(npoints, 0.0)); // vetor para armazenamento dos resultados
	vector<double> Values (N, 0.0);                                          // vetor com os auto valores
	vector<double> X (npoints, 0.0);                                         // vetor com os pontos de avaliação de em X (para plotar);

	linspace(X, npoints, L);
	fill_eingenvalues(Values, N);

	// Solução em todos os tempos
	for (int i = 0; i < nsteps; i++){
		std::cout << "Current time step  = " << i << std::endl;
		for (int j = 0; j < npoints; j++){
			U_xt[i][j] = fourier_adjust(j*dx, i*dt, Values);
			//std::cout << U_xt[i][j] << " " ;
		}
	}
	std::string file_name {"data_full.dat"};
	std::string file_name2 {"dados_N5.dat"};
	//save_data(U_xt, X, file_name);
	save_plot(U_xt, X, file_name2);
}

void fill_eingenvalues(vector<double>& Values, const int total){
	std::cout << "Entered the EigenValues" << std::endl;
	for (int i {0}; i < total; i++){
		Values[i] = std::pow((i+1)*NPI, 2);
	}
	for (auto& e : Values)
		std::cout << e << " ";
}

double fourier_adjust(double x, double t, const vector<double>& Values){
	double res {0.0}, res_0 {0.0}, res_1 {0.0}, res_2 {0.0}, res_3 {0.0}, res_4 {0.0}, res_5{0.0}, lambda_n {0.0}, w_n {0.0};
	auto ntotal = Values.size();

	//std::cout << "Total of Eigenvalues: " << ntotal << std::endl;
	res_0 = std::exp(-beta*t);
	for (int n = 0; n < ntotal; n++){
		lambda_n = Values[n];
		w_n = std::sqrt(lambda_n - beta*beta);
		res_1 = coef_an(n);
		res_2 = std::cos(w_n*t);
		res_3 = (beta/w_n)*res_1;
		res_4 = std::sin(w_n*t);
		res_5 = std::sin(n*NPI*x);
		res += (res_1*res_2 + res_3*res_4)*res_5;
		//std::cout << w_n << " ";
	}
	//std::cout << std::endl;
	res *= res_0;
	return res;
}
double target_function(const double x, const int n){
	return std::sin(n*NPI*x)*x*(1 - x*x);
}
double coef_an(const int n){
	return 2*int_trapz(0, L, n);
}
// Integração numérica pela regra do trapézio
double int_trapz(double a, double b, const int n){
	double res {0.0};
	const double h = 0.00005;                                             // passo da integração numérica
	const auto m = static_cast<int>( std::floor((std::fabs(b - a)) / h));

	for (int k = 0; k < m - 1; k++){
		res += target_function(a + k*h, n);
	}
	res += (target_function(a, n) + target_function(b, n)) / 2;
	res *= h;
	return res;
}

// Torna um vector linearmente espaçado
void linspace(vector<double>& Vec, const int Num, const double xf, const double xi){
	auto h = (xf - xi) / (Num-1);
	auto n = static_cast<int>(Vec.size());
	for (int i = 0; i < n; i++){
		Vec[i] = xi + i*h;
	}
}
void save_data(const vector<vector<double>>& T, const vector<double>& X, const std::string& filename){
	std::fstream printer {filename, std::ios::out|std::ios::trunc};

	auto m = (T[0]).size();
	auto n = T.size();
	
	printer << "Solucao analitica da equacao da onda\n";
	printer << "t";
	for (int k=0; k < m; k++)
		printer << " X" << k;

	for (int i = 0; i < m; i++){
		printer << "\n" << dt*i << ' ';
		for (int j = 0; j < n; j++){
			printer << std::setw(12)<< T[j][i] << ' ';
		}
	}
}

void save_plot(const vector<vector<double>>& T, const vector<double>& X, const std::string& filename){
	std::fstream printer {filename, std::ios::out|std::ios::trunc};

	auto m = (T[0]).size();
	auto n = T.size();

	printer << "Solucao analitica da equacao da onda\n";
	printer << "X";
	for (int k=0; k < n; k++)
		printer << std::setw(12) << " T=" << k*dt;

	for (int i = 0; i < m; i++){
		printer << "\n" << X[i] << ' ';
		for (int j = 0; j < n; j++){
			printer <<  std::setw(12) << T[j][i] << ' ';
		}
	}
	
}

void print_parameters(){
	std::cout << "\nDAMPED WAVE-EQUATION SOLVER" << std::endl;
	std::cout << "\nL = " << L << " m";
	std::cout << "\ndx = " << dx << " m";
	std::cout << "\nt = " << tf << " s";
	std::cout << "\ndt = " << dt << " s";
	std::cout << "\nnx_points = " << npoints;
	std::cout << "\nn_steps = " << nsteps << std::endl;
}
