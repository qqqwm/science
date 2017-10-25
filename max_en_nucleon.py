nucleon_energies = open("proton_energies_final_3.txt", 'r')
maximal_energy = 0
for energy in nucleon_energies:
    if float(energy) > maximal_energy:
        max_en = float(energy)
print(max_en)
nucleon_energies.close()
