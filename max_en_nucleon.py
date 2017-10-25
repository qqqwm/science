nucleon_energies = open("proton_energies_final_3.txt", 'r')
max_en = 0
for energy in nucleon_energies:
    if float(energy) > max_en:
        max_en = float(energy)
print(max_en)
