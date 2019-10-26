import math
data = open('pions.f14', 'r')
k = 0

heading_id = 'UQMD'
'''
pions = []
pions_plus = []
pions_minus = []
pions_0 = []
strings = []
resonances = []
'''
pi_energies = open("pion_energies.txt", 'a')
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
            pi_energies.write(str(float(line[4])))
            pi_energies.write('\n')
            '''
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
                strings.append(float(line[4]))'''
        k = 19
pi_energies.close()
