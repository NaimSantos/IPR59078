import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math


# Variáveis do domínio da simulação e do problema:
beta = 0.1                        # coeficiente de amortecimento
L = 6                           # comprimento total da corda
ti = 0.0                          # tempo inicial da simulação
tf = int(24.0)                    # tempo final da simulação
N = 100                             # número de elementos na série de Fourier
dx = 0.2                         # intervalo em x
dt = 0.25                          # passo de tempo
npoints  = int(L/dx + 1)          # número de pontos de avaliação de x
nsteps = int((tf - ti) / dt)      # número de passos de tempo

# Arrays utilizados
T = np.zeros((nsteps, npoints))   # nsteps instantes de tempo por npoints pontos em x
X = np.linspace(0.0, L, npoints)  # posições em x
ts = np.linspace(ti, tf, nsteps)  # tempos
V = np.zeros(N)                   # autovalores


def plotfxy(eixo_x, eixo_y):
    plt.plot(eixo_x, eixo_y, "r")
    plt.title("Configuração da corda em t = 0 s")
    plt.xlabel("x (m)", fontsize = 13)
    plt.ylabel("u (m)", fontsize = 13)
    plt.savefig('resultado_t_0.png')
    plt.show()

def fill_eigen_values(V) :
    n = len(V)
    print("Número de autovalores usados: ", n)
    i = 0
    while i < n :
        V[i]  = ((i+1)*(math.pi))**2
        i += 1

def function_target(x, n) :
    if (x < 2) :
        res =  x/2 * math.sin(math.pi*n*x/L)
    elif (x >= 2 and x < 4) :
        res =  math.sin(math.pi*n*x/L)
    else:
        res = (- (x-6) / 2) * math.sin(math.pi*n*x/L)
    return res

def int_trapz(a, b, n) :
    res  = 0.0
    h = 0.001
    m = int((b - a)/h)

    k = 0
    while k < (m-1) :
        res += function_target(a + k*h, n)
        k += 1

    res += (function_target(a, n) + function_target(b, n)) / 2
    res *= h
    return res

def coef_an(n) :
    return (2/L)*(int_trapz(0, 2, n) + int_trapz(2, 4, n) + int_trapz(4, 6, n))

def solver() :
    i = 0
    while i < (nsteps) :
        j = 0
        while j < (npoints) :
            T[i][j] = fourier_adjust(j*dx, i*dt)
            j += 1
        i += 1

def fourier_adjust(x, t) :
    print("Fourier Adjust Called... x, t", x, t)
    res1, res2, res3, = 0.0, 0.0, 0.0,
    res = 0.0
    n = 1
    while n < N :
        res1 = coef_an(n)
        res2 = math.cos(n*math.pi*t/L)
        res3 = math.sin(n*math.pi*x/L)
        res = res + (res1*res2*res3)
        n += 1
    return res

#fill_eigen_values(V);
#print(V)
solver();
print(T[0, :])
print(T[1, :])
plotfxy(X, T[0, :])

fig = plt.figure()
plt.title("Solução da equação da onda (N=100)")
plt.xlabel("Comprimento (m)", fontsize = 11)
plt.ylabel("Deslocamento (m)", fontsize = 11)
ax = plt.axes(xlim=(0, 6), ylim=(-1.5, 1.5))
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
line, = ax.plot([], [], lw=2)

# função de inicialização, chamada a cada frame
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

# função para a animação
def animate(i):
    x = X
    y = T[i]
    time_text.set_text('t = % .01f s' % ts[i])
    line.set_data(x, y)
    return line, time_text

anim = animation.FuncAnimation(fig, animate, init_func=init,frames=nsteps, interval=100, blit=True)

anim.save('resultado_beta_01.gif')
plt.show()
