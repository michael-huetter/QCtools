#this script will look at all .log files in the dir and look for negative freq (Gaussian)
import os

def log2xyz(gaussian_file):
    # Open the Gaussian output file
    with open(gaussian_file, 'r') as f:
        lines = f.readlines()

    freq = []
    for line in lines:
        content = line.split()
        if len(content) > 0:
            if content[0] == "Frequencies":
                for k in range(len(content)-2):
                    freq.append(float(content[k+2]))
    return freq
              
def has_negative(arr):
    return any(num < 0 for num in arr)

   
f = []
files = os.listdir("mycalc")
#print(files)
for i in files:
    if i != ".DS_Store":
        f.append((i, log2xyz(f'mycalc/{i}')))
sorted_array = sorted(f, key=lambda x: x[1])
print(sorted_array)
print("")

for i in range(len(sorted_array)):
    to_test = sorted_array[i][1]
    if has_negative(to_test) == True:
        print(f"Negative frequancys found in: {sorted_array[i][0]}")
        print(to_test)
print("")