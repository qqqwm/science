energies = open("proton_energies_final.txt", 'a')
numbers = open("numbers_final.txt", 'a')
for i in range (1, 26):
    name = "urqmd_" + str(i) + "/proton_energies_1.txt"
    nucleon_energies = open(name, 'r')
    name2 = "urqmd_" + str(i) + "/numbers.txt"
    number = open(name2, 'r')
    for line in nucleon_energies:
        energies.write(line)
nucleon_energies.close()
energies.close()
