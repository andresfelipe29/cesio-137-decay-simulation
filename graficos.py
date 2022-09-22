import matplotlib.pyplot as plt
import numpy as np
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

#densidad de probabilidad P(t)=exp(-t/tau)/tau
def Acumulativa(t): #integral de P(t) hasta un t'
    return -np.exp(-t/tau) + 1

def Actividad(t, n): #numero de atomos decaidos por unidad de tiempo
    return (n/tau)*np.exp(-t/tau)

def Restantes(t, n): #numero de atomos restantes en la muestra
    return n*np.exp(-t/tau)

#acumulativa exponencial simulada
exp_acum = [0]*n
with open('datos/exp-acum.txt', 'r') as f:
    i = 0
    for line in f:
        exp_acum[i] = float(line)
        i += 1

x = np.arange(t_min, t_max, step=dt)
plot.set_plot(ticks=[t_min,t_max], ticks_step=30, labels=['t [años]','F(t)', 'Acumulativa exponencial'])
plt.step(x, exp_acum, where='mid', label='F(t)') #simulada
plt.plot(x, Acumulativa(x), color='red') #teorica
x_norm = [t_min, t_max]
y_norm = [1 , 1]
plt.plot(x_norm, y_norm, alpha=0.6) #y=1
plt.legend()
plt.savefig('graficos/exp-acum.pdf')

#acumulativa binomial
x_a=[0,1]
y_a=[1-p,1]
plot.set_plot(ticks=[0,1], ticks_step=1, labels=['Variable aleatoria','F(t)','Acumulativa binomial'])
plt.step(x_a, y_a, where='mid', label='F(t)')

plt.legend()
plt.savefig('graficos/bin-acum.pdf')

#resultados
print('resultados:')
print('exp-acum final:\t', exp_acum[-1])

for nucleos in N:
    #plot Actividad exponencial
    print('Nucleos totales:\t',nucleos)
    exp_hist = [0]*n
    with open('datos/exp-hist-'+str(nucleos)+'.txt', 'r') as f:
        i = 0
        for line in f:
            exp_hist[i] = int(line)
            i += 1

    total = 0
    for freq in exp_hist: total += freq
    print('Núcleos decaidos exp:\t', total)

    plot.set_plot(ticks=[t_min,t_max], ticks_step=30, labels=['t [años]','A(t)','Actividad para N='+str(nucleos)])
    plt.yscale('log')
    plt.step(x, exp_hist, where='mid', label='A(t)')
    plt.plot(x, Actividad(x, nucleos), color='red', label='Teórico')
    plt.legend()
    plt.savefig('graficos/exp-'+str(nucleos)+'.pdf')

    #plot Actividad binomial
    bin_hist = [0]*n
    with open('datos/bin-hist-'+str(nucleos)+'.txt', 'r') as f:
        i = 0
        for line in f:
            bin_hist[i] = int(line)
            i += 1

    total = 0
    for freq in bin_hist: total += freq
    print('Núcleos decaidos bin:\t', total)    

    plot.set_plot(ticks=[t_min,t_max], ticks_step=30, labels=['t [años]','A(t)','Actividad para N='+str(nucleos)])
    plt.yscale('log')
    plt.step(x, bin_hist, where='mid', label='A(t) Binomial')
    plt.step(x, exp_hist, where='mid', label='A(t) Exponencial')
    plt.plot(x, Actividad(x, nucleos), color='red', label='Teórico')
    plt.legend()
    plt.savefig('graficos/bin-'+str(nucleos)+'.pdf')

    #plot nucleos restantes binomial y exponencial
    bin_rest = [0]*n #nucleos restantes en simulacion binomial
    for i in range(n):
        if i==0:
            bin_rest[i]=nucleos-bin_hist[i]
        else:
            bin_rest[i]=bin_rest[i-1]-bin_hist[i]

    exp_rest=[0]*n #nucleos restantes en simulacion exponencial
    for i in range(n):
        if i==0:
            exp_rest[i]=nucleos-exp_hist[i]
        else:
            exp_rest[i]=exp_rest[i-1]-exp_hist[i]

    plot.set_plot(ticks=[t_min,t_max], ticks_step=30, labels=['t [años]','N(t)','Nucleos restantes para N='+str(nucleos)])
    plt.yscale('log')
    plt.step(x, bin_rest, where='mid', label='N(t)')
    plt.step(x, exp_rest, where='mid', label='N(t) Exponencial')
    plt.plot(x, Restantes(x, nucleos), color='red', label='Teórico')
    plt.legend()
    plt.savefig('graficos/bin-res-'+str(nucleos)+'.pdf')

    print('Núcleos restantes exp:\t', exp_rest[-1])
    print('Núcleos restantes bin:\t', bin_rest[-1], '\n')