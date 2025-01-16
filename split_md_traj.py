#split a MD traj into three parts, the initail frame is always the same for correct alignment of the trajectories,
#the rest is split in the ratio given in lines 17-19
#20 a.u of time to ps -> *0.0004837768508

input_file = "co4_md/CO4W1.B3LYP.augdz.30kJ.movie.xyz"  # replace with your input file name
output_file1 = "run1.xyz"
output_file2 = "R3_30_1"
output_file3 = "R3_30_2"

# Read the input file
with open(input_file, "r") as f:
    lines = f.readlines()

n_atoms = int(lines[0].strip())
lines_per_frame = n_atoms + 2
total_frames = len(lines) // lines_per_frame

#total_time = total_frames*0.0004837768508
#fraction_to_remove = 0.2/total_time
#if fraction_to_remove >= 1:
#    exit("time to remove is greater than time in trajectory")

# Calculate the number of frames for each output file
frames_file1 = int(total_frames * 0.02)
frames_file2 = int(total_frames * 0.49)
frames_file3 = total_frames - frames_file1 - frames_file2

print("")
print(f"Times in frame in ps:  1: {(frames_file1*0.0004837768508)}, 2: {(frames_file2*0.0004837768508)}, 3: {(frames_file3*0.0004837768508)}")
print("")

initial_frame = lines[:lines_per_frame]

# Write the output files
with open(output_file1, "w") as f:
    f.writelines(initial_frame + lines[:lines_per_frame * frames_file1])

with open(output_file2, "w") as f:
    f.writelines(initial_frame + lines[lines_per_frame * frames_file1:lines_per_frame * (frames_file1 + frames_file2)])

with open(output_file3, "w") as f:
    f.writelines(initial_frame + lines[lines_per_frame * (frames_file1 + frames_file2):])