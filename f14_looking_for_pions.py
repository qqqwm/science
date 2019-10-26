import matplotlib.pyplot as plt
import math
data = open('pions.f14', 'r')
k = 0

heading_id = 'UQMD'
pions = []
pions_plus = []
pions_minus = []
pions_0 = []
strings = []
resonances = []
number_of_events = 100000
for line in data:
    line = line.split(' ') 
    temp_l = []
    k += 1
    for j in line:#delete '0' and '\n' elements from sublists of data_l
        if len(j) != 0 and '\n' not in j:
            temp_l.append(j)
        elif '\n' in j:
            temp_l.append(j[0:len(j)-1])
    line = temp_l
    if k == 20:
        if line[0] == heading_id:
            k = 0
            continue
        if line[9] == '101' and (float(line[5])**2 + float(line[6])**2)**(1/2)/float(line[7])\
        <= 0.105104235266 and float(line[7]) < 0:
            pions.append(float(line[4]))
            if line[10] == '2':
                pions_plus.append(float(line[4]))
            elif line[10] == '-2':
                pions_minus.append(float(line[4]))
            elif line[10] == '0':
                pions_0.append(float(line[4]))
            if line[14] == '20':
                resonances.append(float(line[4]))
            elif line[14] in ['15', '23', '24', '27', '28']:
                strings.append(float(line[4]))
        k = 19


print(pions)

delta_e = 0.01


def plot(array, label = 'None', m = 0.140):
    "draws a plot"
    number_of_elements = int(max(array)//delta_e + 1)
    energy_interval = []
    number_of_particles = []
    for i in range(number_of_elements):
        number_of_particles.append(0)
        energy_interval.append(delta_e*(i+1))
    for energy in array:
        n = int(energy//delta_e)
        number_of_particles[n] += 1
    x = []
    y = []
    for i in range(len(number_of_particles)):
        if energy_interval[i] > m:
            if number_of_particles[i] != 0:
                x.append(energy_interval[i])
                y.append(number_of_particles[i]/(math.sqrt(energy_interval[i]**2 - m**2))/delta_e/number_of_events)
    plt.step(x, y, label = label)


plt.figure(1)
plot(pions, 'pi')
plot(pions_plus)
plot(pions_minus)
plot(pions_0)
plt.yscale('log')
plt.figure(2)
plot(pions, 'pi')
plot(resonances, 'resonances')
plot(strings, 'strings')
plt.yscale('log')
plt.legend()
plt.show()
