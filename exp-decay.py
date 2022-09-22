import random as rn
import numpy as np
import matplotlib.pyplot as plt
import plot
from scipy import integrate

#SETUP
parameters=[0]*5
with open('setup.txt', 'r') as f:
	i = 0
	for line in f:
		try:
			parameters[i] = line.split("\t")[1]
			i += 1
		except:
			print('invalid setup input')

N = [int(nucleos) for nucleos in parameters[0].split()]
tau= float(parameters[1])
t_min = float(parameters[2])
t_max = float(parameters[3])*tau
dt = tau/int(parameters[4])
n = int(np.floor((t_max-t_min)/dt))
t_max = n*dt
print('setup parameters:')
print('N:\t', N)
print('tau:\t', tau)
print('t min:\t', t_min)
print('t max:\t', t_max)
print('dt:\t', round(dt,3), '\n')

#función a integrar
def prob(t):
	return np.exp(-t/tau)/tau

#probabilidad acumulativa
p_acum = [0]*n
for i in range(n):
	a = t_min + i*dt
	b = t_min + (i+1)*dt
	if i == 0:
		p_acum[0] = integrate.quad(prob, a, b)[0]
	else:
		p_acum[i] = p_acum[i-1] + integrate.quad(prob, a, b)[0]

#normalización
p_acum = [prob / p_acum[-1] for prob in p_acum]

#guardar p-acum
with open('datos/exp-acum.txt', 'w') as f:
    for i in range(n):
        print(p_acum[i], file = f)

#resultados
print('resultados:')
print('exp-acum final:\t',p_acum[-1])

#histograma
histograma = np.array([0]*n)
rn.seed(10)

for nucleos in N:
	print('N:', nucleos)
	for nucleo in range(nucleos):
		r = rn.random()
		if r <= p_acum[0]:
			histograma[0] += 1
			continue
		for Dt in range(n):
			try:
				if p_acum[Dt] < r  and r <= p_acum[Dt+1]:
					bingo = Dt + 1
					break
			except:
				continue
		histograma[bingo] += 1

	total = 0
	for freq in histograma: total += freq
	print('Núcleos decaidos exp:\t', total)

	#----Guardar resultados---
	with open('datos/exp-hist-'+str(nucleos)+'.txt', 'w') as f:
	    for i in range(n) :
	        print(histograma[i], file = f)
	histograma = histograma*0
