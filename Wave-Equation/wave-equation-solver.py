import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math


# Variáveis do domínio da simulação e do problema:
beta = 0.1                        # coeficiente de amortecimento
L = 1.0                           # comprimento total da corda
ti = 0.0                          # tempo inicial da simulação
tf = int(10.0)                     # tempo final da simulação
N = 5                             # número de elementos na série de Fourier
dx = 0.05                         # intervalo em x
dt = 0.1                          # passo de tempo
npoints  = int(L/dx + 1)          # número de pontos de avaliação de x
nsteps = int((tf - ti) / dt)      # número de passos de tempo

# Arrays utilizados
T = np.zeros((nsteps, npoints))   # nsteps instantes de tempo por npoints pontos em x
X = np.linspace(0.0, L, npoints)  # posições em x
ts = np.linspace(ti, tf, nsteps)  # tempos
V = np.zeros(N)                   # autovalores
X025 = np.zeros(nsteps)
X05 = np.zeros(nsteps)
X075 = np.zeros(nsteps)

def plotfxy(eixo_x, eixo_y):
    plt.plot(eixo_x, eixo_y, "r")
    plt.title("Configuração da corda em t = 0 s")
    plt.xlabel("x (m)", fontsize = 13)
    plt.ylabel("u (m)", fontsize = 13)
    plt.savefig('resultado_t_0.png')
    plt.show()

def plot_over_time(eixo_x, eixo_y1, eixo_y2, eixo_y3):
    plt.plot(eixo_x, eixo_y1, label = " x = L/4", color = "red")
    plt.plot(eixo_x, eixo_y2, label = " x = L/2", color = "blue")
    plt.plot(eixo_x, eixo_y3, label = " x = 3L/4", color = "darkgreen")
    plt.title("Evolução no tempo (N=5, β=0.1)")
    plt.legend()
    plt.grid("on")
    plt.xlabel("t (s)", fontsize = 13)
    plt.ylabel("u (m)", fontsize = 13)
    plt.savefig('Plot_Over_Time_beta01.png')
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
        print("\nInstante de tempo = ", i*dt)
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

def fillXmed() :
    x_025 = math.floor(npoints/4) - 1
    x_05 = math.floor(npoints/2)
    x_075 = 3*math.floor(npoints/4) - 1
    print(" x 0,25 = ", x_025)
    print(" x 0,5 = ", x_05)
    print(" x 0,75 = ", x_075)
    k = 0
    while k < nsteps :
        X025[k] = T[k][x_025]
        X05[k]  = T[k][x_05]
        X075[k] = T[k][x_075]
        k += 1

fill_eigen_values(V)
solver()

fillXmed()
plot_over_time(ts, X025, X05, X075)

ymax = np.amax(T)
ymin = np.amin(T)
xmax = np.amax(X)
xmin = np.amin(X)

print("\nx max = ", xmax)
print("x min = ", xmin)
print("y max = ", ymax)
print("y min = ", ymin)
#plotfxy(X, T[0, :])

fig = plt.figure()
plt.title("Solução da equação da onda (N=5, β=0.1)")
plt.xlabel("x (m)", fontsize = 11)
plt.ylabel("u (m)", fontsize = 11)
ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))
time_text = ax.text(0.8, 0.95, '', transform=ax.transAxes)
ax.grid("on")
ax.axis("off")
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
