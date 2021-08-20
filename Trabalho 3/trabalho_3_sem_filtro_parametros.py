# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from scipy.integrate import quad

k = 0.6
rho = 600.
cp = 1200.
T0 = 20.
Tinf = 20.
h = 15.
l = 0.03
g = 100000.

H = h/k
P = g/k
F = T0

alfa0 = 0.
beta0 = 1.
phi0 = 0.

alfal = h/k
betal = 1.
phil = h * Tinf/k

kx = 1

t_end = 11000
t_init = 0

N = 3

alpha = k/(rho * cp)
w = 1./alpha

nodes_x = 31
nodes_t = 31

step_x = l / (nodes_x - 1)
step_t = (t_end - t_init) / (nodes_t - 1)

x = np.arange(0, l + step_x, step_x)
t = np.arange(t_init, t_end + step_t, step_t)

#x = [l/4., l/2., 3.*l/4.]

T = np.zeros((nodes_x, nodes_t))
T_f = np.zeros(nodes_x)

def beta(x):
	return x * np.tan(x * l) - H

def bissection(a, b):
	c = (a + b)/2.
	while(np.abs(beta(c)) > 10**(-8)):
		c = (a + b)/2.0
		if(beta(a) * beta(c) < 0):
			b = c 
		else:
			a = c
	return c

count = 0
eig = 0
i = 0

D = np.zeros(N)
while(eig < N):
	if(beta(i) * beta(i+1) < 0 and count == 0):
		D[eig] = bissection(i, i+20)
		count = count + 1
		eig = eig + 1

	elif(beta(i) * beta(i+1) < 0 and count == 1):
		count = 0
	i = i + 1

def integrate(f, a, b, index, int_step):
	size = (b - a)/int_step
	sum = 0
	for i in range(int_step - 1):
		sum = sum + (f(i * size, index) + f((i+1)*size, index))/2

	sum = sum * size
	return sum

def psi(x, index):
	return np.cos(D[index] * x)

def psi2(x, index):
	return np.cos(D[index] * x)**2

def Tf(x):
	return 0

ni = np.zeros(N)

for i in range(N):
	ni[i] = w * ((l * (D[i]**2 + H**2) + H)/ (2. * (D[i]**2 + H**2)))
	#ni[i] = (w) * integrate(psi2, 0, l, i, 1024*4)
	#print(ni[i])

def psin(x, index):
	return psi(x,index)/np.sqrt(ni[index])

def fi_int(x, index):
	return (w) * psin(x, index) * (F - Tf(x))    

fi = np.zeros(N)
for i in range(N):
	#fi[i] = 1./(alpha) * T0  * integrate(psin, 0, l, i, 1024*4)
	I = quad(fi_int, 0, l, args=(i))
	fi[i] = I[0]
	#fi[i] = integrate(fi_int, 0, l, i, 1024*64)
	#print(fi[i])

giaa = np.zeros(N)
for i in range(N):
	I = quad(psin, 0, l, args=(i))
	giaa[i] = P * I[0]
	#print(giaa[i])

gia = np.zeros(N)

derivate_step = 0.00000000000001
for i in range(N):
	gia[i] = phil * (kx * (psin(l+derivate_step,i) - psin(l, i))/derivate_step - psin(l, i)) / (alfal + betal)
	gia[i] = gia[i] + phi0 * (kx * (psin(l+derivate_step,i) - psin(l, i))/derivate_step - psin(l, i)) / (alfa0 + beta0)
	#print(gia[i])

gi = np.zeros(N)
for i in range(N):
	gi[i] = giaa[i] - gia[i]
	#print(gi[i])

def f(t, index):
	return np.exp(t * alpha * D[index]**2)

A_int = np.zeros(N)
for i in range(N):
	I = quad(f, 0, t_end, args=(i))
	A_int[i] = np.exp(-alpha * D[i]**2. * t_end) * (fi[i] + gi[i] * I[0])
	#A_int[i] = np.exp(-alpha * D[i]**2. * t_end) * (fi[i] + gi[i] * integrate(f, 0., t_end, i, 1024*64))
	#A_int[i] = np.exp(-alpha * D[i]**2 * t_end) * (fi[i] + gi[i] * (np.exp(alpha * D[i]**2 * t_end)/(alpha * D[i]**2) - 1/(alpha * D[i]**2)))
	#print(A_int[i])

sum = 0
for i in range(nodes_x):
	for k in range(N):
		#sum = sum + psin(x[i], k) * A_int[k]
		sum = sum + psin(x[i], k) * A_int[k]
	T[i][1] = sum
	sum = 0


print("\n")
for i in range(nodes_x):
	T_f[i] = T[i][1] + Tf(x[i])
	print(T_f[i])

Tf_vec = np.zeros(nodes_x)
for i in range(nodes_x):
	Tf_vec[i] = Tf(x[i])

T_t = np.zeros(nodes_t + 1)
	
def run_t(t):	
	for i in range(N):
		A_int[i] = np.exp(-alpha * D[i]**2. * t) * (fi[i] + gi[i] * integrate(f, 0., t, i, 1024*4))
		#A_int[i] = np.exp(-alpha * D[i]**2 * t_end) * (fi[i] + gi[i] * (np.exp(alpha * D[i]**2 * t_end)/(alpha * D[i]**2) - 1/(alpha * D[i]**2)))
		#print((np.exp(alpha * D[i]**2 * t_end)/(alpha * D[i]**2) - 1/(alpha * D[i]**2)))

	sum = 0
	for k in range(N):
		sum = sum + psin(x[nodes_x/2], k) * A_int[k]
		#print(A_int[k])

	return sum + Tf_vec[1]

for i in range(nodes_t + 1):
	T_t[i] = run_t(t[i])
	

#print("%d & %.5f & %.5f & %.5f \\\\" % (N, T_t[0], T_t[1], T_t[2]))


#################################################################

chart = plt.figure()    
plt.plot(x, T_f,'r', linewidth=1)
plt.ylabel(r'Temperatura $(^oC)$')
plt.grid(True)

chart.savefig("/Users/Libotte/Desktop/chart_citt.pdf", bbox_inches='tight')

plt.show()


#################################################################
'''
chart = plt.figure()    
plt.plot(t, T_t,'r', linewidth=1)
plt.xlabel(r'Tempo $(s)$')
plt.ylabel(r'Temperatura $(^oC)$')
plt.grid(True)

chart.savefig("/Users/Libotte/Desktop/chart_citt.pdf", bbox_inches='tight')

plt.show()
'''