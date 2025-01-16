import re
import os

def log2xyz(gaussian_file, nam):
    name = nam
    atom_number_to_symbol = {1: "H", 8: "O", 47: "Ag"}
    # Open the Gaussian output file
    with open(gaussian_file, 'r') as f:
        lines = f.readlines()

    # Extract the final energy from the file
    #energy_pattern = re.compile(r'SCF Done:  E\([A-Z]\) =')
    energy_pattern = re.compile("SCF Done:")
    for line in lines:
        energy_match = energy_pattern.search(line)
        if energy_match:
            energy = float(line.split()[4])
            break


    # Extract the final coordinates from the file
    coordinates_pattern = re.compile(r'Standard orientation:')
    Natoms_pattern = re.compile(r'NAtoms=')
    Natoms = 0
    for i, line in enumerate(lines):
        if Natoms_pattern.search(line):
                Natoms = int(line.split()[1])
        if coordinates_pattern.search(line):
            coordinates_start = i + 5
    coordinates = []
    atoms = []
    for line in lines[coordinates_start:coordinates_start+Natoms]:
        coordinates.append([line.split()[3], line.split()[4], line.split()[5]])
        atoms.append(int(line.split()[1]))
    
    # Write the final energy and coordinates to an XYZ file
    j = 0
    with open(f'strucs/{name}.xyz', 'w') as f:
        f.write(str(Natoms) + '\n')
        try:
            f.write('Energy = ' + str(energy) + '\n')
        except:
            f.write(f'Energy = {str(energy)}\n')
        for coord in coordinates:
            symbol = atom_number_to_symbol.get(atoms[j], "X")
            f.write(symbol + " " + coord[0] + ' ' + coord[1] + ' ' + coord[2] + '\n')
            j += 1


files = os.listdir("5_strucs_gatt")
for i in files:
    if i != ".DS_Store":
        print(i)
        log2xyz(f"5_strucs_gatt/{i}", i)

        