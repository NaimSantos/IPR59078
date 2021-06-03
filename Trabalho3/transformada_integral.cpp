#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <functional>
#include <string>
#include <cassert>

using std::vector;
using std::string;

constexpr auto k {0.6};
constexpr auto rho {600};
constexpr auto cp {1200};
constexpr auto T0 {20.0};
constexpr auto Tinf {20.0};
constexpr auto h {15.0};
constexpr auto l {0.03};
constexpr auto g {100000};

constexpr auto H = h/k;
constexpr auto P {0.0};
constexpr auto F = T0;

constexpr auto alfa0 {0.0};
constexpr auto beta0 {1.0};
constexpr auto phi0 {0.0};
constexpr auto alfal {h};
constexpr auto betal {k};
constexpr auto phil {0.0};

constexpr auto t_end {11000};
constexpr auto t_init {0};
constexpr auto N {3};

constexpr auto alpha = k/(rho * cp);
constexpr auto w = 1.0/alpha;

constexpr int nodes_x {31};
constexpr int nodes_t {31};

int main (int argc, char* argv[]){
	auto step_x = l / (nodes_x - 1);
	auto step_t = (t_end - t_init)/ (nodes_t - 1);
	
	//x = np.arange(0, l + step_x, step_x)
	//t = np.arange(t_init, t_end + step_t, step_t)
	
	vector<double> x (nodes_x, 0.0);
	vector<double> t (nodes_t, 0.0);
	
	vector<vector<double>> T (nodes_x, std::vector<double>(nodes_t, 0.0));
	vector<double> T_f (nodes_x, 0.0);

	int count = 0;
	int eig = 0;
	int i = 0;

	vector<double> D (N, 0.0);
	while(eig < N){
		if((beta(i)*beta(i+1) < 0) && count == 0){
			D[eig] = bissection(i, i+30);
			count = count + 1;
			eig = eig + 1;
		}
		else if(beta(i) * beta(i+1) < 0 and count == 1){
			count = 0;
			i = i + 1;
		}
	}

	vector<double> ni (N, 0.0);	
	for (int i = 0; i < N; i++){
	//ni[i] = (1./alpha) * integrate(psi2, 0, l, i, 1024*4)
		ni[i] = w * ((l * (D[i]*D[i] + H*H) + H)/ (2.0 * (D[i]*D[i] + H*H)));
	}

	vector<double> Tf_vec (nodes_x, 0.0);
	for (int i=0; i < nodes_x; i++){
		Tf_vec[i] = Tf(x[i]);
	}

	vector<double> fi (N, 0.0);
	for (int i=0; i < N; i++){
		//fi[i] = integrate(fi_int, 0, l, i, 1024*4)
		I = quad(fi_int, 0, l, args=(i));
		fi[i] = I[0];
	}

	vector<double> giaa (N, 0.0);
	for (int i=0; i < N; i++){
		giaa[i] = P * integrate(psin, 0, l, i, 1024*4);
	}

	vector<double> gia (N, 0.0);
	auto derivate_step = 0.0000000001;
	
	
	for (int i=0; i < N; i++){
		gia[i] = phil * (k * (psin(l+derivate_step, i) - psin(l, i))/derivate_step - psin(l, i)) / (alfal + betal);
		gia[i] = gia[i] + phi0 * (k * (psin(l+derivate_step, i) - psin(l, i))/derivate_step - psin(l, i)) / (alfa0 + beta0);
	}
	
	vector<double> gi (N, 0.0);
	for (int i=0; i < N; i++){
		gi[i] = giaa[i] - gia[i];
	}
	vector<double> T_t (nodes_t + 1, 0.0);
	vector<double> A_int (nodes_t + 1, 0.0);

	for (int i=0; i < N; i++){
		T_t[i] = run_t(t[i], A_int);
	}
	/*
	A_int = np.zeros(N)
	for (int i=0; i < N; i++){
		A_int[i] = std::exp(-alpha * D[i]*D[i] * t_end) * (fi[i] + gi[i] * integrate(f, 0.0, t_end, i, 1024*4));
	}
	*/
	
	for (int i=0; i < nodes_x; i++){
		double sum = 0.0;
		for (int k=0; k < N; k++){
			sum = sum + psin(x[i], k) * A_int[k];
		}
		T[i][1] = sum;
	}

	
	for (int i=0; i < nodes_x; i++){
		T_f[i] = T[i][1] + Tf_vec[i];
	}

}

double beta(double x){
	return x * std::tan(x * l) - H;
}
double bissection(double a, double b){
	double c = (a + b)/2;
	while(std::fabs(beta(c)) > 0.00000001){
		c = (a + b)/2.0;
		if(beta(a) * beta(c) < 0)
			b = c;
		else
			a = c;
	}
	return c;
}
double integrate(std::function<double>(double) f, double a, double b, int index, double int_step){
	auto size = (b - a)/int_step;
	double sum = 0.0;
	for (int i = 0; i < int_step-1; i++){
		sum = sum + (f(i * size, index) + f((i+1)*size, index))/2;
	}
	sum = sum * size;
	return sum;
}
double psi(double x, int index, vector<double>& D){
	return std::cos(D[index] * x);
}
double psi2(double x, int index, vector<double>& D){
	return std::pow(std::cos(D[index]*x), 2);
}
double psin(double x, int index, vector<double>& D, const vector<double>& ni){
	return psi(x, index, D) / std::sqrt(ni[index]);
}
double Tf(double x){
	//return Tinf/k + (g*l/h) + g * (l**2 - x**2)/(2 * k)
	return (2.0*g*k*l + g*h*(std::pow(l, 2.0)) + 2.0*h*k*Tinf - g*h*(std::pow(x, 2.0))/(2.0*h*k));
}

double fi_int(double x, int index, vector<double>& D, const vector<double>& ni){
	return (w) * psin(x, index, D, ni) * (F - Tf(x));
}
double f(double t, int index, vector<double>& D){
	return std::exp(t * alpha * D[index]*D[index]);		
}
double run_t(double t, vector<double>& A_int){
	for (int i=0; i < N; i++){
		A_int[i] = std::exp(-alpha * D[i]*D[i] * t) * (fi[i] + gi[i] * integrate(f, 0., t, i, 1024*4));
	}


	sum = 0
	for k in range(N):
		sum = sum + psin(x[nodes_x/2], k) * A_int[k]
		#print(A_int[k])

	return sum + Tf_vec[1]
}

/*
double integrate(std::function<double>(double)f, a, b, index, int_step){
	auto size = (b - a)/int_step;
	double sum = 0;
	for i in range(int_step - 1):
		sum = sum + (f(i * size, index) + f((i+1)*size, index))/2

	sum = sum * size;
	return sum;
}
*/