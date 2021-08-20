import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import scipy.integrate as scpy

h = 15.
L = 0.03
g = 5000.0
k = 0.6
rho = 600.0
cp = 1200.0
T0 = 20.0
Tinf = 20.0
H = h/k
P = g/k
F = T0
alfa0 = 0.0
beta0 = 1.0
phi0 = 0.0
alfal = h/k
betal = 1.0
phil = h*Tinf/k
kx = 1
t_end = 11000
t_init = 0
N = 3
alpha = k/(rho*cp)
w = 1.0/alpha
nodes_x = 100
nodes_t = 10
dx = L / (nodes_x - 1)
dt = (t_end-t_init) / (nodes_t-1)

x = np.arange(0, L + dx, dx)
t = np.arange(t_init, t_end + dt, dt)
T = np.zeros((nodes_x, nodes_t))
T_f = np.zeros(nodes_x)

def beta(x):
	return x * np.tan(x * L) - H

def rootsbybissect(a, b):
	c = (a + b)/2
	while(np.abs(beta(c)) > 10**(-6)):
		c = (a + b)/2.0
		if(beta(a) * beta(c) < 0):
			b = c
		else:
			a = c
	return c

ammount = 0
eigenvalue = 0
i = 0
Roots = np.zeros(N)
while(eigenvalue < N):
	if(beta(i) * beta(i+1) < 0 and ammount == 0):
		Roots[eigenvalue] = rootsbybissect(i, i+20)
		ammount = ammount + 1
		eigenvalue = eigenvalue + 1
	elif(((beta(i)*beta(i+1)) < 0) and (ammount==1)):
		ammount = 0
	i = i + 1

def integrate(f, a, b, indice, int_step):
	h = (b-a)/int_step
	sum = 0.0
	for i in range(int_step - 1):
		sum = sum + (f(i * h, indice) + f((i+1)*h, indice))/2

	sum = sum*h
	return sum

def psi(x, indice):
	return np.cos(Roots[indice]*x)
def psi_quad(x, indice):
	return np.cos(Roots[indice]*x)**2
def Tf(x):
	return 0

ni = np.zeros(N)
for i in range(N):
	ni[i] = w * ((L * (Roots[i]**2 + H**2) + H)/ (2.0 * (Roots[i]**2 + H**2)))

def psin(x, indice):
	return psi(x,indice)/np.sqrt(ni[indice])
def fi_int(x, indice):
	return (w) * psin(x, indice) * (F - Tf(x))

fi = np.zeros(N)
for i in range(N):
	res = scpy.quad(fi_int, 0, L, args=(i))
	fi[i] = res[0]

giaa = np.zeros(N)
for i in range(N):
	res = scpy.quad(psin, 0, L, args=(i))
	giaa[i] = res[0] * P

gia = np.zeros(N)

limite_step = 0.00000000000001
for i in range(N):
	# contribuição do contorno x=L:
	gia[i] = phil * (kx * (psin(L+limite_step,i) - psin(L, i))/limite_step - psin(L, i)) / (alfal + betal)
	# contrinuição em L + contribuição do contorno x=0:
	gia[i] = gia[i] + phi0 * (kx * (psin(L+limite_step,i) - psin(L, i))/limite_step - psin(L, i)) / (alfa0 + beta0)

gi = np.zeros(N)
for i in range(N):
	gi[i] = giaa[i] - gia[i]
def f(t, indice):
	return np.exp(t * alpha * Roots[indice]**2)

A_int = np.zeros(N)
for i in range(N):
	I = scpy.quad(f, 0, t_end, args=(i))
	A_int[i] = np.exp(-alpha * Roots[i]**2. * t_end) * (fi[i] + gi[i] * I[0])

sum = 0
for i in range(nodes_x):
	for k in range(N):
		sum = sum + psin(x[i], k) * A_int[k]
	T[i][1] = sum
	sum = 0

for i in range(nodes_x):
	T_f[i] = T[i][1] + Tf(x[i])
	print(T_f[i])

Tf_vec = np.zeros(nodes_x)
for i in range(nodes_x):
	Tf_vec[i] = Tf(x[i])

T_t = np.zeros(nodes_t + 1)
	
def evaluator(t):
	sum = 0
	for i in range(N):
		A_int[i] = np.exp(-alpha * Roots[i]**2. * t) * (fi[i] + gi[i] * integrate(f, 0., t, i, 1024*4))

	for k in range(N):
		sum = sum + psin(x[0], k) * A_int[k]

	return sum + Tf_vec[1]

#for i in range(nodes_t + 1):
#	T_t[i] = evaluator(t[i])


def plotfxy(eixo_x, eixo_y):
    plt.plot(eixo_x, eixo_y, 'b', linewidth=2)
    plt.xlabel("Comprimento (m)", fontsize = 11)
    plt.ylabel("Temperatura (° C)", fontsize = 11)
    plt.legend(loc='upper center', fontsize=9)
    plt.grid(True, 'major', 'both')
    plt.savefig('Grafico1.png')
    plt.show()
    
plotfxy(x, T_f)

