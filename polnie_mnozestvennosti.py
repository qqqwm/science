#import matplotlib.pyplot as plt
#import math
#import scipy
#import numpy
data = open('test31100.f14', 'r')
number_of_events = 10000
number_of_inelastic_events = 0
protons = 0
antiprotons = 0
pi_minus = 0
pi_plus = 0
k_minus = 0
k_plus = 0
lambd = 0
identifiers = [[]]#list of sublists of tuples, each tuple contains particle_ID and iso3, each sublist contains tuples that correspond to all particles in one event
indentifiers_inelastic = []#similarly to identifiers, elastic collisions excluded
string_counter = 0
list_counter = 0
for line in data:
    line = line.split(' ') 
    temp_l = []
    string_counter += 1
    for j in line:#delete '0' and '\n' elements from sublists of data_l
        if len(j) != 0 and j != '\n':
            temp_l.append(j)
    line = temp_l
    if string_counter == 20:
        if line[0] == 'UQMD':
            string_counter = 1
            list_counter += 1
            identifiers.append([])
            continue
        identifiers[list_counter].append((float(line[9]), float(line[10])))
        string_counter = 19
for index in identifiers:
    if len(index) != 2:
         indentifiers_inelastic.append(index)
         number_of_inelastic_events += 1
for i in indentifiers_inelastic:#number of different particles counting
    for j in i:
        if j[0] == 1 and j[1] == 1:
            protons += 1
        elif j[0] == -1 and j[1] == -1:
            antiprotons += 1
        elif j[0] == 106 and j[1] == 1:
            k_plus += 1
        elif j[0] == -106 and j[1] == -1:
            k_minus += 1
        elif j[0] == 101:
            if j[1] == -2:
                pi_minus += 1
            elif j[1] == 2:
                pi_plus += 1
        elif j[0] == 27:
            lambd += 1
print('<p>', protons/number_of_inelastic_events)
print('<anti_p>', antiprotons/number_of_inelastic_events)
print('<pi^+>', pi_plus/number_of_inelastic_events)
print('<pi^->', pi_minus/number_of_inelastic_events)
print('<k^+>', k_plus/number_of_inelastic_events)
print('<k^->', k_minus/number_of_inelastic_events)
print('<lambda>', lambd/number_of_inelastic_events)
print('sum', protons + antiprotons + pi_plus + pi_minus + k_plus + k_minus)

