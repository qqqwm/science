import matplotlib.pyplot as plt
import math
from matplotlib.ticker import AutoMinorLocator
delta_e = 0.01
number_of_events = 500000000
nucleons = []
nucleon_energies = open("proton_energies_final_3.txt", 'r')
def plot(array, label = 'None', m = 0.940):
    "draws a plot"
    number_of_elements = int(1.5126698//delta_e + 1)
    energy_interval = []
    number_of_particles = []
    for i in range(number_of_elements):
        number_of_particles.append(0)
        energy_interval.append(delta_e*(i+1))
    for energy in nucleon_energies:
        if float(energy) > 0.094:
            n = int(float(energy)//delta_e)
            number_of_particles[n] += 1
    x = []
    y = []
    for i in range(len(number_of_particles)):
        if energy_interval[i] > m:
            if number_of_particles[i] != 0:
                x.append(energy_interval[i])
                y.append(number_of_particles[i]/(math.sqrt(energy_interval[i]**2 - m**2))/delta_e/number_of_events)
    plt.step(x, y, label = label)
    plt.xlabel('$E$, GeV', size = 16)
    plt.ylabel('$(pN_{ev})^{-1}dN/dE$, c/GeV$^2$', size = 16)
#for i in nucleon_energies:
#    nucleons.append(float(i))
plt.figure(1)
plot(nucleons, 'pi')
'''
plot(pions_plus)
plot(pions_minus)
plot(pions_0)
plt.yscale('log')
plt.figure(2)
plot(resonances, 'resonances')
plot(strings, 'strings')
'''
plt.axes().get_yaxis().set_minor_locator(AutoMinorLocator())
plt.axes().get_xaxis().set_minor_locator(AutoMinorLocator())
plt.tick_params(which = 'major', direction='in', bottom ='on', top ='on', right ='on', left ='on', length=6)
plt.tick_params(which = 'minor', direction='in', bottom ='on', top ='on', right ='on', left ='on', length=3)
plt.tick_params(axis = 'both', which = 'major', labelsize = 12 )
plt.yscale('log')
#plt.legend(prop = {'size' : 12})
plt.show()
