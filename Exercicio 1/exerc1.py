import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani


# Variáveis do domínio da simulação:
r = 0.5 # condição CFL
dx = 0.1 # intervalo em x
dt = 1 # passo de tempo
L = 1.0 # comprimento total da placa
ti = 0.0 # tempo inicial da simulação
tf = int(20.0) #tempo final da simulação
N  = int(L/dx + 1) # número de nós da malha (intervalos + 1)
I = (tf - ti) / dt # número de passos de tempo


# Dados do problema:
kappa = 0.6
rho = 600
cp = 1200
h = 15.0
g = 100000
T_zero = 20.0


# Solver explicito para a temperatura:
def explictsolver(Arr, cfl, N, tf):
    i = 1 # i = tempo = linhas
    while i < tf :
        j = 1 #comecamos no segundo elemento
        while j < (N - 1):
            Arr[i][j] = Arr[i-1][j] + cfl*(Arr[i-1][j-1] - 2.0*Arr[i-1][j] + Arr[i-1][j+1] )
            j += 1
        i += 1

def f(x):
    if x <= 0.5:
        return (x)
    else:
        return (1 - x)


T = np.zeros((tf, N)) #Array para temperaturas, N elementos, 100 tempos
X = np.linspace(0.0, L, N)
ts = np.arange(0, 100, 1)
print(X)
print(ts)

j = 0 
while j < N:
    T[0][j]  = f(j*dx)
    j += 1


explictsolver(T, 0.5, N, tf)
# print(T)

x_ax = X 
y_ax = T[0,  :]
# plt.plot(x_ax, y_ax, "r") # red diamonds
# plt.title("Perfil de temperatura da placa unidimensional")
# plt.xlabel("Comprimento", fontsize = 13)
# plt.ylabel("Temperatura", fontsize = 13)

# plt.show()

fig = plt.figure()
plt.title("Perfil de temperatura da placa unidimensional")
plt.xlabel("Comprimento", fontsize = 11)
plt.ylabel("Temperatura", fontsize = 11)
ax = plt.axes(xlim = (0, 1.0), ylim=(0, 0.5))
line, = ax.plot([], [], lw=2)

def init():
    line.set_data([], [])
    return line,

def animate(i):
    x = X
    y = T[i, :]
    line.set_data(x, y)
    return line,


#animate(X, T)
anim = ani.FuncAnimation(fig, animate, init_func=init, frames=20, interval=1, blit=True)
anim.save('temperatura.gif')
plt.show()