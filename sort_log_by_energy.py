import numpy as np
import os
import re
import matplotlib.pyplot as plt

def log2xyz(gaussian_file):
    # Open the Gaussian output file
    with open(gaussian_file, 'r') as f:
        lines = f.readlines()

    # Extract the final energy from the file
        #energy_pattern = re.compile("Sum of electronic and zero-point Energies=")
        #energy_pattern = re.compile("Wavefunction amplitudes converged")
        energy_pattern = re.compile("SCF Done:")
        for line in lines:
            energy_match = energy_pattern.search(line)
            if energy_match:
                #energy = float(line.split()[6])
                energy = float(line.split()[4])
    return energy

    
E = []
files = os.listdir("b3lyp_aug_cc_pvtz")
#print(files)
for i in files:
    if i != ".DS_Store":
        E.append((i, log2xyz(f'b3lyp_aug_cc_pvtz/{i}')))
sorted_array = sorted(E, key=lambda x: x[1])
print("")
print(sorted_array)
print("")
E_diff = []
for j in range(len(sorted_array)):
    E_diff.append(round((sorted_array[0][1]-sorted_array[j][1])*-27, 3))
print(E_diff)
plt.plot(E_diff, "*")
plt.ylabel("dE in eV")
plt.show()

