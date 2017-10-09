"""
Converts energies in laboratory system to momenta in center of mass system.
"""
m = 0.938 #mass of proton
p_lab = (20, 31, 40, 80, 158)
E_cm = []
for p in p_lab:
    E = 2*m/(1 - p**2/(m + (p**2 + m**2)**(1/2))**2)**(1/2)
    E_cm.append(E)
print(p_lab)
print(E_cm)
