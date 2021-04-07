import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani


# Variáveis do domínio da simulação:
L = 0.03                  # comprimento total da placa
N  = 100                   # número de nós da malha
ti = 0.0                  # tempo inicial da simulação
tf = 500.0                # tempo final da simulação
dx = L/(N-1)              # comprimento do intervalo
dt = 0.01                  # passo de tempo
nsteps = int((tf-ti)/dt)  # número de passos de tempo

# Dados do problema:
kappa = 0.6
rho = 600
cp = 1200
h = 15.0
g = 100000
T0 = 20.0
TL = 20.0
r1 = (kappa*dt)/(rho*cp*dx*dx)         # coeficiente r do método implícito
r2 = (kappa*dt)/(rho*cp*2*dx*dx)       # coeficiente r do Crank-Nicolson
eta = 3.0 + ((2*h*dx)/kappa)
mu = (g*dt)/(rho*cp)

print("r implicito: ", r1)
print("r Nicolson: ", r2)
def plotfxy(eixo_x, y1, y2):
    plt.plot(eixo_x, y1, 'r', label= 'Implícito')
    plt.plot(eixo_x, y2, 'blue', label= 'Crank-Nicolson')
    plt.title("Perfil de temperatura da placa unidimensional (t=500 s)")
    plt.xlabel("Comprimento (m)", fontsize = 11)
    plt.ylabel("Temperatura (° C)", fontsize = 11)
    plt.legend(loc='upper center', fontsize=9)
    plt.grid(True, 'major', 'both')
    plt.savefig('Grafico1.png')
    plt.show()
def init():
    line.set_data([], [])
    time_text.set_text('t = 0.0 s')
    return line, time_text
def animate(i):
    x = X
    y = T[i, :]
    ymin=np.min(y)
    ymax=np.max(y)
    time_text.set_text('t = % .1f s' % ts[i])
    ax.set_ylim(ymin, ymax)
    line.set_data(x, y)
    return line, time_text

# Matrizes para armazenamento dos resultados
T1 = np.zeros((nsteps, N))   # Array para temperaturas, com N elementos por linha em nsteps linhas
T2 = np.zeros((nsteps, N))  # Copia para o Cranck Nicolson
X = np.linspace(0.0, L, N)  # Vetor das posições linearmente espaçado

# Método totalmente implícito:
def implictsolver(A, B, T):
    T[0] = B.reshape(1, N)
    t = 1
    while t < nsteps:
        B[0][0] = 0.0
        B[N-1][0] = (2*h*dx*T0)/kappa
        i = 1
        while i < N-1 :
            B[i][0] = B[i][0] + mu
            i = i + 1
        B = np.linalg.solve(A, B)
        T[t] = B.reshape(1, N)
        t = t + 1
def solveimplicitly(r, T):
    # Preenchimento da matriz de termos independentes:
    B = np.full((N, 1), T0)
    # Preenchimento da matriz de coeficientes:
    A = np.zeros((N,N))
    A[0][0] = -3.0
    A[0][1] = 4.0
    A[0][2] = -1.0
    A[N-1][N-3] = 1.0
    A[N-1][N-2] = -4.0
    A[N-1][N-1] = eta
    i = 1
    j = 0
    while i < N-1 :
        A[i][j] = - r
        A[i][j+1] = 1 + 2*r
        A[i][j+2] = -r
        j = j+1
        i = i+1
    implictsolver(A, B, T)

# Crank-Nicolson:
def nicolsonsolver(A, B, T, r):
    T[0] = B.reshape(1, N)
    t = 1
    while t < nsteps:
        i = 1
        # Corrigimos os Bs internos primeiro:
        while i < N-1 :
            B[i][0] = r*B[i-1][0] + (1-2*r)*B[i][0] + r*B[i+1][0] + mu
            i = i + 1
        # Corrigimos os Bs do contorno:  
        B[0][0] = 0.0
        B[N-1][0] = (2*h*dx*T0)/kappa  
        B = np.linalg.solve(A, B)
        T[t] = B.reshape(1, N)
        t = t + 1
def solvebyNicolson(r, T):
    # Preenchimento da matriz de termos independentes:
    B = np.full((N, 1), T0)
    # Preenchimento da matriz de coeficientes:
    A = np.zeros((N,N))
    A[0][0] = -3.0
    A[0][1] = 4.0
    A[0][2] = -1.0
    A[N-1][N-3] = 1.0
    A[N-1][N-2] = -4.0
    A[N-1][N-1] = eta
    i = 1
    j = 0
    while i < N-1 :
        A[i][j] = -r
        A[i][j+1] = 1 + 2*r
        A[i][j+2] = -r
        j = j+1
        i = i+1
    nicolsonsolver(A, B, T, r)



solveimplicitly(r1, T1)
print("T[x=0] =", T1[nsteps-1][0])
print("T[x=0.03] =", T1[nsteps-1][N-1])

solvebyNicolson(r2, T2)
print("T[x=0] =", T2[nsteps-1][0])
print("T[x=0.03] =", T2[nsteps-1][N-1])

plotfxy(X, T1[nsteps-1], T2[nsteps-1])

# Animação :
# x_ax = X
# y_ax = T[0, :]
# ts = np.linspace(0, tf, nsteps)

# fig = plt.figure()
# plt.title("Perfil de temperatura da placa unidimensional")
# plt.xlabel("Comprimento", fontsize = 10)
# plt.ylabel("Temperatura", fontsize = 10)
# ax = plt.axes(xlim = (0.0, L), ylim=(0, 21), autoscale_on=False)
# line, = ax.plot([], [], lw=2,)
# time_text = ax.text(0.9, 0.9, 'tempo', transform=ax.transAxes)

# anim = ani.FuncAnimation(fig, animate, init_func=init, frames=nsteps,interval=1, blit=False)
# anim.save('temperatura.gif')
# plt.show()





