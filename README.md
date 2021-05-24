# MNED-doutorado


Um repositório com minhas implementações de exercícios e trabalhos da disciplina Métodos Numéricos para Equações Diferenciais, ministrada pelo professor Diego Knupp, do Doutorado em Modelagem Computacional do IPRJ - UERJ.

## Trabalho 1

* Solução numérica da equação do calor (difusão-advecção).
	* Com geração de calor interna (fonte).
	* Condições de contorno: segundo tipo (Neumann) no lado esquerdo, terceiro tipo (Robin) no lado direito.
* Discretização do domínio por diferenças finitas.
* Discretização dos contornos via diferenças finitas e via nós-fictícios.
* Esquema implícito simples e esquema da Cranck-Nicolson.


## Trabalho 2:

* Solução analítica da equação do calor (difusão-advecção)
	* Com geração de calor interna (fonte).
	* Condições de contorno: segundo tipo (Neumann) no lado esquerdo, terceiro tipo (Robin) no lado direito.
* Formulação via separação de variáveis
	* Tratamento via solução filtro e solução filtrada.
* Solução via série de Fourier
	* Integral no cálculo dos coeficientes da série via integração numérica por regra do trapézio.
	* Cálculo de autovalores na equação transcendental via Newton-Rhapson.
	
## Wave Equation

* Solução analítica da equação da onda
	* d2u/dx2 = d2u/dt2 + 2 beta * du/dt
	* u(0,t) = 0
	* u(1,t) = 0
	* u(x,0) = x(1-x^2)
	* du/dt = 0 em t=0
	* 0 <= x <= 1
	* 0 < beta < 1

* Solução via série de Fourier
	* Integral no cálculo dos coeficientes da série via integração numérica por regra do trapézio.
	* Cálculo de autovalores na equação transcendental via Newton-Rhapson.
