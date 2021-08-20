#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <functional>
#include <string>

using std::vector;
using std::string;

double int_trapz(std::function<double (double, int)> f, double a, double b, int t);
double int_trapz(std::function<double (double, int)> f, double a, double b, int index, bool teste);
double fi_int(double x, int index);
double psin(double x, int index);
double f(double t, int index);
double beta(const double x);
double bissection(double a, double b);
double integrate(std::function<double (double, int)>f, double a, double b, int index, int int_step);
double psi(double x, int index);
double psi2(double x, int index);
double Tf(double x);
double run_t(const double t);
void savedata(const vector<double>& T, const vector<double>& X, const string& filename);
void linspace(vector<double>& Vec, const int Num, const double xf = 1.0, const double xi = 0.0);
template <typename T> void print_vector(const vector<T>& Vec);

double k {0.6};
double rho {600.0};
double cp {1200.0};
double T0 {20.0};
double Tinf {20.0};
const double h {15.0};
double L {0.03};
const double g {100000.0};
auto H = h/k;
auto P = g/k;
auto F = T0;
auto alfa0 {0.0};
auto beta0 {1.0};
auto phi0 {0.0};
auto alfal = h/k;
auto betal = 1.0;
auto phil = h*Tinf/k;
double kx = 1.0;
double t_end = 500;
double t_init = 0;
const int N = 3;
auto alpha = k/(rho*cp);
auto w = 1.0/alpha;
const int nodes_x = 31;
const int nodes_t = 31;
auto step_x = L/(nodes_x - 1);
auto step_t = (t_end - t_init)/(nodes_t - 1);


vector<double> X (nodes_x, 0.0);
vector<double> t (nodes_t, 0.0);
vector<vector<double>> T (nodes_x, std::vector<double>(nodes_t, 0.0));
vector<double> T_f (nodes_x, 0.0);
vector<double> EVal (N, 0.0);
vector<double> N_i (N, 0.0);
vector<double> fi (N, 0.0);
vector<double> giaa (N, 0.0);
vector<double> gia (N, 0.0);
vector<double> gi (N, 0.0);
vector<double> A_int (N, 0.0);
vector<double> Tf_vec (nodes_x, 0.0);
vector<double> T_t (nodes_t+1, 0.0);

int main (int argc, char* argv[]){

	linspace(X, nodes_x, L);
	linspace(t, nodes_t, t_end);
	print_vector(X);
	print_vector(t);

	int count = 0;
	double eig = 0;
	int i = 0;

	std::cout << "Ponto 1" << std::endl;
	while(eig < N){
		if(((beta(i) * beta(i+1)) < 0) && (count==0)){
			EVal[eig] = bissection(i, i+20);
			count++;
			eig++;
		}
		else if (((beta(i) * beta(i+1)) < 0) and (count==1)){
			count = 0;
		}
		i++;
	}

	std::cout << "Ponto 2" << std::endl;
	for (int i = 0; i < N; i++){
		N_i[i] = w*((L* (EVal[i]*EVal[i] + H*H) + H)/ (2.0 * (EVal[i]*EVal[i] + H*H)));
	}

	std::cout << "Ponto 3" << std::endl;
	for (int i = 0; i < N; i++){
		auto I = int_trapz(fi_int, 0, L, i);
		fi[i] = I;
	}

	std::cout << "Ponto 4" << std::endl;
	//g i barra asterisco asterisco (t)
	for (int i = 0; i < N; i++){
		auto I = int_trapz(psin, 0, L, i);
		giaa[i] = P * I;
	}
	for (auto x : giaa)
		std::cout << x << " ";
	std::cout << "Ponto 5" << std::endl;
	double derivate_step = 0.00000000000001;
	//g i barra asterisco (t)
	for (int i = 0; i < N; i++){
		gia[i] = phil * (kx * (psin(L+derivate_step, i) - psin(L, i))/derivate_step - psin(L, i)) / (alfal + betal);
		gia[i] = gia[i] + phi0 * (kx * (psin(L+derivate_step, i) - psin(L, i))/derivate_step - psin(L, i)) / (alfa0 + beta0);
	}

	std::cout << "Ponto 6" << std::endl;
	for (int i = 0; i < N; i++){
		gi[i] = giaa[i] - gia[i];
	}

	std::cout << "Ponto 7" << std::endl;
	for (int i = 0; i < N; i++){
		std::cout << "N= " << i << std::endl;
		auto I = int_trapz(f, 0, t_end, i, true);
		A_int[i] = std::exp(-alpha*EVal[i]*EVal[i]*t_end)*(fi[i] + gi[i] * I);
	}

	std::cout << "Ponto 8" << std::endl;
	double sum {0.0};
	for (int i = 0; i < nodes_x; i++){
		for (int i = 0; i < N; i++){
			sum = sum + psin(X[i], k) * A_int[k];
		}
		T[i][1] = sum;
		sum = 0;
	}

	std::cout << "Ponto 9" << std::endl;
	for (int i = 0; i < nodes_x; i++){
		T_f[i] = T[i][1] + Tf(X[i]);
	}

	std::cout << "Ponto 10" << std::endl;
	//Armazena em cada Tf_vec 0.0, para manter consistência
	for (int i = 0; i < nodes_x; i++){
		Tf_vec[i] = Tf(X[i]);
	}
	std::cout << "Ponto 11" << std::endl;
	for (int i = 0; i <= (nodes_t); i++){
		T_t[i] = run_t(t[i]);
	}
	string arquivo {"Dados_Transformada.dat"};
	savedata(T_f, X, arquivo);
	
	return 0;
}

double int_trapz(std::function<double (double, int)> f, double a, double b, int t){
	double res {0.0};
	const double h {0.00001};
	const auto n = static_cast<int>( std::floor((std::fabs(b - a)) / h));
	for (int i = 0; i < n - 1; i++){
		res += f(a+i*h, t);
	}
	res += (f(a, t) + f(b, t) ) / 2;
	res *= h;
	return res;
}
double int_trapz(std::function<double (double, int)> f, double a, double b, int index, bool teste){
	static int count = 1;
	std::cout << "Integral via trapezio V2 invocada " << count << " vezes" << std::endl;
	double res {0.0};
	const auto n = 512;
	static long long iter = 1;
	for (int i = 0; i < n - 1; i++){
		res += f(a+i*h, index);
		iter++;
	}
	std::cout << "Exited the for loop" << std::endl;
	res += (f(a, index) + f(b, index) ) / 2;
	res *= h;
	count++;
	return res;
	
}
double fi_int(double x, int index){
	return w * psin(x, index) * (F - Tf(x));
}
//phi tilde i (x)
double psin(double x, int index){
	return psi(x, index) / std::sqrt(N_i[index]);
}
double f(double t, const int index){
	return w*std::exp(t*alpha*EVal[index]*EVal[index]);
}
//phi(x)
double beta(const double x){
	return x * std::tan(x*L) - H;
}
double bissection(double a, double b){
	double c = (a + b) / 2;
	while(std::fabs(beta(c)) > 10e-8)
	{
		c = (a + b)/2.0;
		if(beta(a) * beta(c) < 0)
			b = c;
		else
			a = c;
	}
	return c;
}
double integrate(std::function<double (double, int)>f, double a, double b, int index, int int_step){
	auto size = (b - a)/int_step;
	double sum {0.0};
	for (int i = 0; i < int_step-1; i++)
	{
		sum = sum + (f(i*size, index) + f((i+1)*size, index))/2;
	}
	sum = sum * size;
	return sum;
}
double psi(double x, int index){
	return std::cos(EVal[index]*x);
}
double psi2(double x, int index){
	return std::pow(std::cos(EVal[index]*x), 2);
}
double Tf(double x){
	return 0.0;
}
double run_t(const double t){
	double res {0.0};
	for (int i = 0; i < N; i++){
		res = integrate(f, 0.0, t, i, 1024*4);
		A_int[i] = std::exp(-alpha*EVal[i]*EVal[i]*t)*(fi[i] + gi[i] * res);
	}
	double sum {0.0};
	for (int i = 0; i < N; i++){
		sum = sum + psin(X[nodes_x/2], k) * A_int[k];
	}
	return sum + Tf_vec[1];
}
void savedata(const vector<double>& T, const vector<double>& X, const string& filename){
	std::cout << "\nSaveData function called" << std::endl;
	// Salvar os resultados em um arquivo:
	std::fstream printer {filename, std::ios::out|std::ios::trunc};
	printer << "Perfil de Temperatura via solucao analitica\n";
	
	auto xtotal = X.size();
	auto ttotal = T.size();
	
	std::cout << "xtotal, ttotal = " << xtotal << ", " << ttotal <<std::endl;
	for (int k = 0; k < xtotal; k++)
		printer << " X_" << k;
	printer <<"\n-";
	for (int k = 0; k < xtotal; k++)
		printer << ' ' << X[k] ;
	printer <<"\n-";
	for (int j = 0; j < ttotal; j++)
		printer <<' ' << T[j];
}
void linspace(vector<double>& Vec, const int Num, const double xf, const double xi){
	auto h = (xf - xi) / (Num-1);
	auto n = static_cast<int>(Vec.size());
	for (int i = 0; i < n; i++){
		Vec[i] = xi + i*h;
	}
}
template <typename T>
void print_vector(const vector<T>& Vec){
	static int ammount = 1;
	std::cout << "\nExibindo vetor " << ammount << " : " << std::endl;
	for (const auto& x: Vec){
		std::cout << x << " ";
	}
	ammount++;
	std::cout << std::endl;
}