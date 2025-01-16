import re
import os

def log2xyz(gaussian_file,o, nname):
    n = nname
    atom_number_to_symbol = {1: "H", 8: "O", 47: "Ag"}
    # Open the Gaussian output file
    with open(gaussian_file, 'r') as f:
        lines = f.readlines()

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
    with open(f"out/{n}.cc.com", 'w') as f:
        f.write("%nproc=20\n")
        f.write("%Mem=190GB\n")
        f.write(f"%chk={n}.cc.chk\n")
        f.write("#CCSD(T)/genecp\n")
        f.write("\n")
        f.write(f"CCSD(T) with aug-cc-pvtz-pp for {n}.com\n")
        f.write("\n")
        f.write("1 1\n")
        for coord in coordinates:
            symbol = atom_number_to_symbol.get(atoms[j], "X")
            f.write(symbol + " " + coord[0] + ' ' + coord[1] + ' ' + coord[2] + '\n')
            j += 1
        f.write(" ")

    with open("aug_cc_pVQZ.txt", "r") as source_file, open(f"out/{n}.cc.com", "a") as destination_file:
        # read the data from the source file
        data = source_file.read()
        # write (append) the data to the destination file
        destination_file.write(data)



files = os.listdir("6_to_calc")
l = 0
for i in files:
    if i != ".DS_Store" and i != "4.rrk":
        print(i)
        log2xyz(f"6_to_calc/{i}", l, i)
        l += 1