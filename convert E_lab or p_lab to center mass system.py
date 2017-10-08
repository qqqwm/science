m = 0.938 #mass of proton
"""x = input("Press 'pE' to convert p_lab to E_cm, \n\
press 'pp' to convert p_lab to p_cm, \n\
press 'EE' to convert E_lab to E_cm, \n\
press 'Ep' to convert E_lab to p_cm. \n")
if x == 'pE' or x == 'pp':
    p_lab = input('Enter p_lab in GeV: \n')
    if x == 'pp':
        pass
    else:
        E_cm = 2*m/(1 - p_lab**2/(m + (p_lab**2 + m**2)**(1/2)))**(1/2)
        print(E_cm)
if x == 'EE' or x == 'Ep':
    E_lab = input('Enter E_lab in GeV: \n')"""
p_lab = (20, 31, 40, 80, 158)
E_cm = []
for p in p_lab:
    E = 2*m/(1 - p**2/(m + (p**2 + m**2)**(1/2))**2)**(1/2)
    E_cm.append(E)
print(p_lab)
print(E_cm)
