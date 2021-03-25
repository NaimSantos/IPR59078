#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>

/*
Método de Gauss-Siedel para resolver um sistema de equações lineares
Parâmetros:
	A = Matriz dos coeficientes
	m,n = dimensões de A
	X = Matriz das incógnitas
	B = Matriz dos termos independentes
	n_eq = dimensão de B e X (numero de equações)
	eps = define o erro permitido, usado como critério de parada
No final, X é sobreescrito.
*/

constexpr double eps = 0.00001;

void GS_Solver(double** A, const int m, const int n, double* B, const unsigned int n_eq, const float eps, double* X){

	double* Y = new double [n_eq];	//matriz auxiliar
	double* E = new double [n_eq];	//necessária para estimar o erro de uma iteração a outra
	for (int i = 0; i < n_eq; i++)
		E[i] = X[i];
	
	unsigned int counter = 1; //Contar iterações apenas pro caso da tolerancia nao ser atingida. Se a matriz é diagonal dominante, a convergência é garantida
	bool teste = false;

	while(!teste && counter<20){
		teste = true;
		std::cout << "Iteracao " << std::setprecision(10) << counter << '\n';
		for (int i = 0; i < m; i++){
			Y[i] = (B[i] / A[i][i]);
			for (int j = 0; j < n; j++){
				if (j==i)
					continue;
				Y[i] = Y[i] - ((A[i][j] / A[i][i]) * X[j]);
				X[i] = Y[i]; //escreve em X a estimativa encontrada
			}
			auto res = std::fabs(((X[i] - E[i]) / X[i])) <= eps;
			teste = teste & res;
			std::cout<< "x" << i + 1 << " = " << Y[i] << '\n';
			E[i] = X[i];
		}
		counter++;
		std::cout << '\n';
	}
	delete[] E;
	delete[] Y;
}
