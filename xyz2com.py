import re
import os

def log2xyz(gaussian_file,o):
    n = o
    with open(gaussian_file, 'r') as f:
        content = f.read()
    with open(f"out/{n}.5.blyp.com", 'w') as f:
        f.write("%nproc=4\n")
        f.write(f"%chk={n}.chk\n")
        f.write("#BLYP/def2SVP scf=yqc DensityFit freq opt\n")
        f.write("\n")
        f.write(f"BLYP optimization and freq for {n}\n")
        f.write("\n")
        f.write("1 1\n")
        f.write(content)
        f.write(" ")
        

files = os.listdir("isomers")
l = 0
for i in files:
    if i != ".DS_Store":
        log2xyz(f"isomers/{i}", l)
        l += 1