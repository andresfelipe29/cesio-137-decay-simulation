import random as rn
import numpy as np
import matplotlib.pyplot as plt
import plot

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
p=dt/tau
print('setup parameters:')
print('N:\t', N)
print('tau:\t', tau)
print('t min:\t', t_min)
print('t max:\t', t_max)
print('dt:\t', round(dt,3))
print('p:\t', p, '\n')

#---------Acumulativa-----
x_a=[0,1]
y_a=[p,1]

#---------Conteo----------
rn.seed(20)
for nucleos in N:
    t = 0
    k = 1
    decaimientos = 0
    histograma = [0]*n
    total = 0
    print('N:',nucleos)

    while t<t_max:

        for nucleo in range(nucleos):
            r = rn.random()
            if r <= p:
                decaimientos += 1

        nucleos -= decaimientos
        total += decaimientos
        histograma[k-1] = decaimientos
        if nucleos == 0: break

        t=t_min+k*dt
        k+=1
        decaimientos=0

    print('Núcleos decaidos bin:\t', total)

    #----Guardar resultados---
    f_name='datos/bin-hist-'+str(nucleos+total)+'.txt'
    with open(f_name, 'w') as f:
        for i in range(n):
            print(histograma[i], file = f)
