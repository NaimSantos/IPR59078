import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani


# Variáveis do domínio da simulação:
L = 0.03                  # comprimento total da placa
N  = 25                   # número de nós da malha (intervalos + 1)
ti = 0.0                  # tempo inicial da simulação
tf = 500.0                # tempo final da simulação
dx = L / (N - 1)          # comprimento do intervalo
dt = 0.1                  # passo de tempo
nsteps = int((tf-ti)/dt)  # número de passos de tempo

# Dados do problema:
kappa = 0.6
rho = 600
cp = 1200
h = 15.0
g = 100000
T0 = 20.0
TL = 20.0
alpha = kappa/(rho*cp)
r1 = (alpha*dt)/(dx*dx)         # coeficiente r do método implícito
r2 = (alpha*dt) / (2*dx*dx)     # coeficiente r do Crank-Nicolson
gamma = 3.0 + ((2*h*dx)/kappa)
llambda = (g*dt)/(rho*cp)

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
    plt.plot(eixo_x, eixo_y, "r")
    plt.title("Perfil de temperatura da placa unidimensional")
    plt.xlabel("Comprimento", fontsize = 12)
    plt.ylabel("Temperatura", fontsize = 12)
    plt.show()

# Solver explícito para a temperatura:
def explictsolver(Arr, cfl, N, tf):
    i = 1
    while i < tf :
        j = 1 # começamos no segundo elemento
        while j < (N - 1):
            Arr[i][j] = Arr[i-1][j] + cfl*(Arr[i-1][j-1] - 2.0*Arr[i-1][j] + Arr[i-1][j+1] )
            j += 1
        i += 1

T = np.zeros((nsteps, N))   # Array para temperaturas, com N elementos por linha em nsteps linhas
T2 = np.zeros((nsteps, N))  # copia para o Cranck Nicolson
X = np.linspace(0.0, L, N)  # Vetor das posições linearmente espaçado

# Solver implícito:
def implictsolver(A, B, C):
    T[0] = B.reshape(1, N)
    t = 1
    while t < nsteps:
        B[0][0] = 0.0
        B[N-1][0] = (2*h*dx*T0)/kappa
        i = 1
        while i < N-1 :
            B[i][0] = B[i][0] + llambda
            i = i + 1
        B = np.linalg.solve(A, B)
        T[t] = B.reshape(1, N)
        t = t + 1

def solveimplicitly(r):
    # Preenchimento da matriz de termos independentes:
    B = np.full((N, 1), T0)
    C = np.full((N, 1), T0)
    # Preenchimento da matriz de coeficientes:
    A = np.zeros((N,N))
    A[0][0] = -3.0
    A[0][1] = 4.0
    A[0][2] = -1.0
    A[N-1][N-3] = 1.0
    A[N-1][N-2] = -4.0
    A[N-1][N-1] = gamma
    i = 1 # linha
    j = 0 # coluna
    while i < N-1 :
        A[i][j] = - r
        A[i][j+1] = 1 + 2*r
        A[i][j+2] = -r
        j = j+1
        i = i+1
    #print("A: ", A)
    #print("B antes: ", B)
    implictsolver(A, B, C)



solveimplicitly(r1)
plotfxy(X, T[nsteps-1])

# explictsolver(T, 0.5, N, tf)
# fillbounds(T, f, N)
# print(T)


x_ax = X 
y_ax = T[0, :]
ts = np.arange(0, nsteps, 1)

fig = plt.figure()
plt.title("Perfil de temperatura da placa unidimensional")
plt.xlabel("Comprimento", fontsize = 11)
plt.ylabel("Temperatura", fontsize = 11)
ax = plt.axes(xlim = (0, 0.03), ylim=(20, 100))
line, = ax.plot([], [], lw=2)
time_text = ax.text(0.01,01.0, '', transform=ax.transAxes)

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

#anim = ani.FuncAnimation(fig, animate, init_func=init, frames=100, #interval=1, blit=True)
#anim.save('temperatura.gif')
#plt.show()
