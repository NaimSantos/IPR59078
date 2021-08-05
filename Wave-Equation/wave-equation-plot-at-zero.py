import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math


# Variáveis do domínio da simulação e do problema:
beta = 0.9                        # coeficiente de amortecimento
L = 1.0                           # comprimento total da corda
ti = 0.0                          # tempo inicial da simulação
tf = int(1.0)                     # tempo final da simulação
N = 3                             # número de elementos na série de Fourier
dx = 0.05                         # intervalo em x
dt = 0.1                          # passo de tempo
npoints  = int(L/dx + 1)          # número de pontos de avaliação de x
nsteps = int((tf - ti) / dt)      # número de passos de tempo

# Arrays utilizados
T = np.zeros((nsteps, npoints))   # nsteps instantes de tempo por npoints pontos em x
X = np.linspace(0.0, L, npoints)  # posições em x
ts = np.linspace(ti, tf, nsteps)  # tempos
V = np.zeros(N)                   # autovalores



XValues = np.linspace(0.0, L, npoints*2)
YValues = np.linspace(0.0, L, npoints*2)

def polinomio(x):
    return x *(1 - x**2)

def fixXYyoplot():
    p = 0
    while p < (npoints*2):
        YValues[p] = polinomio(XValues[p])
        p += 1

fixXYyoplot()
def plotfxy(eixo_x, eixo_y):
    plt.plot(eixo_x, eixo_y, label = "Série (N=3)", color = 'b')
    plt.plot(XValues, YValues, label = "u(x,t) ", color = 'r')
    plt.legend()
    plt.grid("on")
    plt.xlabel("x (m)", fontsize = 13)
    plt.ylabel("u (m)", fontsize = 13)
    plt.savefig('resultado_N3.png')
    plt.show()

def fill_eigen_values(V) :
    print("\nPreenchendo auto valores...")
    n = len(V)
    print("\nNúmero de autovalores usados: ", n)
    i = 0
    while i < n :
        V[i]  = ((i+1)*(math.pi))**2
        i += 1

def function_target(x, n) :
    res = (x - x**3) * math.sin(math.pi*n*x)
    return res

def int_trapz(a, b, n) :
    res  = 0.0
    h = 0.0001
    m = int((b - a)/h)

    k = 0
    while k < (m-1) :
        res += function_target(a + k*h, n)
        k += 1

    res += (function_target(a, n) + function_target(b, n)) / 2
    res *= h
    return res

def coef_an(n) :
    return 2*int_trapz(0, L, n)

def solver() :
    print("\nInicializando o solver...")
    i = 0
    while i < (nsteps) :
        print("Instante de tempo = ", i*dt)
        j = 0
        while j < (npoints) :
            T[i][j] = fourier_adjust(j*dx, i*dt)
            j += 1
        i += 1

def fourier_adjust(x, t) :
    res0, res1, res2, res3, res4, res5 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    w_n, lambda_n = 0.0, 0.0
    res = 0.0
    n = 0
    res0 = math.exp(-beta * t)
    while n < len(V) :
        lambda_n = V[n]
        w_n = math.sqrt(lambda_n - beta*beta)
        res1 = coef_an(n)
        res2 = math.cos(w_n*t)
        res3 = (beta/w_n)*res1
        res4 = math.sin(w_n*t)
        res5 = math.sin(math.pi * n * x)
        res = res + (((res1*res2) + (res3*res4))*res5)
        n += 1
    res = res* res0
    return res
    

fill_eigen_values(V)
solver()


plotfxy(X, T[0, :])
