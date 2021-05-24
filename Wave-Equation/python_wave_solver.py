import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math


# Variáveis do domínio da simulação e do problema:
beta = 0.5                       # coeficiente de amortecimento
L = 1.0                          # comprimento total da corda
ti = 0.0                         # tempo inicial da simulação
tf = int(50.0)                   # tempo final da simulação
N = 5                            # número de elementos na série de Fourier
dx = 0.05                        # intervalo em x
dt = 1                           # passo de tempo
npoints  = int(L/dx + 1)         # número de pontos de avaliação de x
nsteps = int((tf - ti) / dt)     # número de passos de tempo

# Arrays utilizados
T = np.zeros((nsteps, npoints))   # N elementos, 100 tempos
X = np.linspace(0.0, L, npoints)  # posições em x
ts = np.linspace(ti, tf, nsteps)  # tempos
V = np.zeros(N)                   # autovalores


def plotfxy(eixo_x, eixo_y):
    plt.plot(eixo_x, eixo_y, "r")
    plt.title("Corda em repouso teste")
    plt.xlabel("x (m)", fontsize = 13)
    plt.ylabel("u (m)", fontsize = 13)
    plt.savefig('resultado_t_0.png')
    #plt.show()

def fill_eigen_values(V) :
    n = len(V)
    print("Numero de autovalores usados: ", n)
    i = 0
    while i < n :
        V[i]  = ((i+1)*(math.pi))**2
        i += 1

def function_target(x, n) :
    res = (x - x**3) * math.sin(math.pi*n*x)
    return res

def int_trapz(a, b, n) :
    res  = 0.0
    h = 0.00005
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
    i = 0
    while i < (nsteps) :
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

fill_eigen_values(V);
print(V)
solver();

plotfxy(X, T[0, :])

fig = plt.figure()
plt.title("Solucao da equacao da onda")
plt.xlabel("Comprimento", fontsize = 11)
plt.ylabel("Deslocamento", fontsize = 11)
ax = plt.axes(xlim=(0, 1), ylim=(0, 0.45))
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
line, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

# animation function.  This is called sequentially
def animate(i):
    x = X
    y = T[i]
    time_text.set_text('t = % .01f s' % ts[i])
    line.set_data(x, y)
    return line, time_text

anim = animation.FuncAnimation(fig, animate, init_func=init,frames=nsteps, interval=dt, blit=True)

anim.save('temperatura.gif')
plt.show()
