"""
Makes plots of maximal energy of pion depending on momentum of proton.
"""
import matplotlib.pyplot as plt
import math
import scipy
import numpy
data = open('test.f14', 'r')
k = 0
energies = []
number_of_events = 10000
for line in data:
    line = line.split(' ') 
    temp_l = []
    k += 1
    for j in line:#delete '0' and '\n' elements from sublists of data_l
        if len(j) != 0 and j != '\n':
            temp_l.append(j)
    line = temp_l
    if k == 20:
        if line[0] == 'UQMD':
            k = 0
            continue
        if line[9] == '1':
            energy = float(line[4])
            energies.append(energy)
        k = 19
delta_e = 0.5
number_of_elements = int(max(energies)//delta_e + 1)
energy_interval = []
number_of_particles = []
for i in range(number_of_elements):
    number_of_particles.append(0)
    energy_interval.append(delta_e*(i+1))
for energy in energies:
    n = int(energy//delta_e)
    number_of_particles[n] += 1
x = []
y = []
for i in range(len(number_of_particles)):
    if energy_interval[i] > 0.938:
        x.append(energy_interval[i])
        y.append(math.log(number_of_particles[i]/4/math.pi/(math.sqrt(energy_interval[i]**2 - 0.938**2))/delta_e)/number_of_events)
f = open('data1.txt', 'w')
for i in x[int(len(x)/5):int(len(x)/1.25)]:
    f.write(str(i))
    f.write('\n')
f.write('\n')
for i in y[int(len(x)/5):int(len(x)/1.25)]:
    f.write(str(i))
    f.write('\n')
f.close()
plt.figure(1)
plt.scatter(x, y)
def f(t):
    return -(4.4973011*10**(-4))*t+0.0012144
t = numpy.linspace(0.9, 2, 2)
p = f(t)
plt.plot(t, p, label = 'y = -(4.4973011*10^(-4))*x+0.0012144 - approximation1')
def f(t):
    return -(1.9934093*10**(-5))*t+(3.7990867*10**(-4))
t = numpy.linspace(2, 5, 2)
p = f(t)
plt.plot(t, p, label = 'y = -(1.9934093*10^(-5))*x+(3.7990867*10^(-4)) - approximation2')
plt.xlabel('Energy, GeV')
plt.legend()
#plt.axis([0.8, 6.2, 0, 0.0008])
plt.show()

"""
plt.figure(2)
plt.scatter(x[0:int(len(x)/2.5)], y[0:int(len(x)/2.5)])
plt.xlabel('Energy, GeV')
plt.show()
plt.figure(3)
plt.scatter(x[int(len(x)/2.5):int(len(x)/1.25)], y[int(len(x)/2.5):int(len(x)/1.25)])
plt.xlabel('Energy, GeV')
plt.show()
"""


