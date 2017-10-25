import math
data = open('testprotons.f14', 'r')
k = 0
number = 1
numbers = open("numbers.txt", 'a')
pi_energies = open("proton_energies_1.txt", 'a')
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
        if line[0] == 'UQMD':
            k = 0
            number += 1
            continue
        if line[9] == '1' and (float(line[5])**2 + float(line[6])**2)**(1/2)/float(line[7])\
        <= 0.105104235266 and float(line[7]) < 0:
            pi_energies.write(str(float(line[4])))
            pi_energies.write('\n')
        k = 19
pi_energies.close()
numbers.write(str(number))
numbers.close()
data.close()
