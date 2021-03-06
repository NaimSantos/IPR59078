import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani


# Variáveis do domínio da simulação:
r = 0.5             # condição CFL
dx = 0.05            # intervalo em x
dt = 1               # passo de tempo
L = 1.0              # comprimento total da placa
ti = 0.0             # tempo inicial da simulação
tf = int(100.0)      #tempo final da simulação
N  = int(L/dx + 1)   # número de nós da malha (intervalos + 1)
I = (tf - ti) / dt   # número de passos de tempo

# Dados do problema:
kappa = 0.6
rho = 600
cp = 1200
h = 15.0
g = 100000
T_zero = 20.0


def f(x):
    if x <= 0.5:
        return (x)
    else:
        return (1 - x)


def fillbounds(T, f, N):
    j = 0
    while j < N:
        T[0][j]  = f(j*dx)
        j += 1

def plotfxy(eixo_x, eixo_y):
    plt.plot(eixo_x, eixo_y, "r") # red diamonds
    plt.title("Perfil de temperatura da placa unidimensional")
    plt.xlabel("Comprimento", fontsize = 13)
    plt.ylabel("Temperatura", fontsize = 13)
    plt.show()


# Solver explicito para a temperatura:
def explictsolver(Arr, cfl, N, tf):
    i = 1 # i = tempo = linhas
    while i < tf :
        j = 1 #comecamos no segundo elemento
        while j < (N - 1):
            Arr[i][j] = Arr[i-1][j] + cfl*(Arr[i-1][j-1] - 2.0*Arr[i-1][j] + Arr[i-1][j+1] )
            j += 1
        i += 1

T = np.zeros((tf, N)) #Array para temperaturas, N elementos, 100 tempos
X = np.linspace(0.0, L, N)
ts = np.arange(0, 100, 1)


explictsolver(T, 0.5, N, tf)
# fillbounds(T, f, N)
# print(T)

x_ax = X 
y_ax = T[0, :]

fig = plt.figure()
plt.title("Perfil de temperatura da placa unidimensional")
plt.xlabel("Comprimento", fontsize = 11)
plt.ylabel("Temperatura", fontsize = 11)
ax = plt.axes(xlim = (0, 1.0), ylim=(0, 0.5))
line, = ax.plot([], [], lw=2)
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

def animate(i):
    x = X
    y = T[i, :]
    time_text.set_text('t = % .01f s' % ts[i])
    line.set_data(x, y)
    return line, time_text


anim = ani.FuncAnimation(fig, animate, init_func=init, frames=100, interval=1, blit=True)
anim.save('temperatura.gif')
plt.show()


res = np.linalg.solve()