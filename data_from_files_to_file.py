st = "testdata\proton_energies.txt"
nucleon_energies = open(st, 'r')
energies = open("proton_energies_2.txt", 'a')
for line in nucleon_energies:
    energies.write(line)
nucleon_energies.close()
energies.close()
