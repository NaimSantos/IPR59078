#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <functional>
#include <string>

using std::vector;
using std::string;


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
const int N = 5;
auto alpha = k/(rho*cp);
auto w = 1.0/alpha;
const int nodes_x = 31;
const int nodes_t = 31;
auto dx = L/(nodes_x - 1);
auto dt = (t_end - t_init)/(nodes_t - 1);

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

	std::cout << "Ponto 1" << std::endl;
	// Calcula os autovalores:
	fill_bis();
	
	
}

void fill_bis(){
	int eig =0;
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
}
double beta(const double x){
	return x*std::tan(x*L) - H;
}
double Xi(double x, int i){
	return std::cos(EVal[i]*x);
}
double Ni(double x, int i){
	return w*((L* (EVal[i]*EVal[i] + H*H) + H)/ (2.0 * (EVal[i]*EVal[i] + H*H)));
}