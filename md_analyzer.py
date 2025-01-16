#MD trajectory analyzer (from xyz file). For alighnment all atoms labeled as tc are taken into account. Do this with the rename_md_traj.py script. Voulumes is calculated with the convex hull method. All the functions are named with water hydrogen etc. becaus i wrote the script initially for them, but works for any atom typ.

#Changes might be needed in lines marked with <-

import numpy as np
#import os
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import MDAnalysis as mda
from MDAnalysis.analysis import align
from rename_md_traj import rename_atoms_Ag, rename_atoms_co
#from scipy.spatial import distance_matrix
#import alphashape
from scipy.spatial.qhull import ConvexHull
#from scipy.spatial import Voronoi
#import seaborn as sns
from sklearn.neighbors import KernelDensity
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys
#from scipy.integrate import nquad
#from scipy.integrate import quad_vec
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from scipy.spatial import cKDTree


#------------------------------------------------------------
#some stuff i aleso tried but dont need atm


def bounding_box(points, coverage):
    sorted_points = [np.sort(points[:, i]) for i in range(3)]

    num_points = len(points)
    num_points_to_cover = int(coverage * num_points)

    min_lengths = [np.inf] * 3
    min_indices = [0] * 3
    max_indices = [0] * 3

    for i in range(3):
        for j in range(num_points - num_points_to_cover + 1):
            length = sorted_points[i][j + num_points_to_cover - 1] - sorted_points[i][j]
            if length < min_lengths[i]:
                min_lengths[i] = length
                min_indices[i] = j
                max_indices[i] = j + num_points_to_cover - 1

    min_coords = [sorted_points[i][min_indices[i]] for i in range(3)]
    max_coords = [sorted_points[i][max_indices[i]] for i in range(3)]

    x_length, y_length, z_length = [max_coords[i] - min_coords[i] for i in range(3)]

    return min_coords, max_coords, x_length, y_length, z_length

def visualize_2d(trajectory):
    all_coords = np.vstack([frame[:, 1:3] for frame in trajectory])

    # Calculate the convex hull using the combined points
    convex_hull = ConvexHull(all_coords)

    fig, ax = plt.subplots()

    # Plot the convex hull
    convex_hull_plot_2d(convex_hull, ax=ax)

    # Plot the atoms
    for frame in trajectory:
        for atom in frame:
            element, x, y, z = atom[0], atom[1], atom[2], atom[3]
            ax.scatter(x, y, label=element, zorder=10)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.title('2D Projection of Molecules with Convex Hull Contour')

    plt.show()    

#------------------------------------------------------------




def read_xyz(filename):
    trajectory = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        num_atoms = int(lines[0].strip())
        frame = []
        for i, line in enumerate(lines[2:]):
            if not line.startswith(("Ag", "O", "H", "C")):  #<-
                continue

            frame.append([line.strip().split()[0]] + list(map(float, line.strip().split()[1:])))
            if (i + 1) % num_atoms == 0:
                trajectory.append(frame)
                frame = []
    return np.array(trajectory)


def get_water_molecules_indexes(trajectory):
    #only consider atoms labeled as tc
    water_molecules_indexes = []
    for i, atom in enumerate(trajectory[0]):
        if atom[0] == "tc":
            water_molecules_indexes.append(i)
    return water_molecules_indexes


def get_water_coords(trajectory, water_index):
    water_coords = []
    for frame in trajectory:
        water_coords.append(frame[water_index][1:])
    return np.array(water_coords)


def downsample_data(x, y, z, normalized_density, fraction):
    num_points = len(x)
    num_sampled_points = int(fraction * num_points)
    sampled_indices = np.random.choice(num_points, num_sampled_points, replace=False)

    x_sampled = x[sampled_indices]
    y_sampled = y[sampled_indices]
    z_sampled = z[sampled_indices]
    normalized_density_sampled = normalized_density[sampled_indices]
    #print(f"sampled_indices dtype: {sampled_indices.dtype}, shape: {sampled_indices.shape}")

    return x_sampled, y_sampled, z_sampled, normalized_density_sampled


def downsample_dense_points(chunk_indices, coordinates, downsampling_distance):
    chunk_coordinates = coordinates[chunk_indices]
    chunk_tree = cKDTree(chunk_coordinates)
    downsampled_indices = chunk_tree.query_pairs(downsampling_distance)
    indices_to_remove = np.unique([pair[1] for pair in downsampled_indices])
    #print(f"indices_to_remove dtype: {indices_to_remove.dtype}, shape: {indices_to_remove.shape}")

    return np.delete(chunk_indices, indices_to_remove)


def chunk_indices(indices, num_chunks):
    chunk_size = len(indices) // num_chunks
    for i in range(0, len(indices), chunk_size):
        yield indices[i:i + chunk_size]

def compute_density(chunk_xyz, kde):
    return np.exp(kde.score_samples(chunk_xyz))

def visualize(trajectory, water_molecules_coords, water_oxy_molecules_coords, hydrogen_molecules_coords, create_fig, name_fig, n,k):
    
    fig = plt.figure(); ax = fig.add_subplot(121, projection='3d')

    # Plot the initial molecule structure (frame[0])
    for atom in trajectory[0]:
        element, x, y, z = atom[0], atom[1], atom[2], atom[3]
        if element == "O" or element == "tc" or element == "Ag" or element == "C":   #<-
            color = "black"
            ax.scatter(x, y, z, color=color, marker='o', s=60)
            #ax.text(x, y, z, str(oxygen_counter), fontsize=12, color='k')
            #oxygen_counter += 1

    # save the trajectory of the oxygen atoms
    x = []; y = []; z = []
    for molecule_coords, color in zip([water_oxy_molecules_coords], ["purple"]):  #<-
        for coords in molecule_coords:
            x.extend(coords[:, 0]); y.extend(coords[:, 1]); z.extend(coords[:, 2])

    #Compute the point density of the oxygen atom trajectory using a kernel density estimation
    xyz = np.vstack([x, y, z]).T
    kde = KernelDensity(kernel='gaussian', bandwidth=0.5).fit(xyz)   # <---------- at the moment the proble super slow
    density = np.exp(kde.score_samples(xyz))  
 
    # Normalize the density values to the range [0, 1]
    normalized_density = (density - density.min()) / (density.max() - density.min())    
    
    x = np.array(x); y = np.array(y); z = np.array(z)
    #Downsample the data to reduce the number of points plotted
    density_threshold = 0.8
    dense_points_indices = np.where(normalized_density > density_threshold)[0]
    downsampling_distance = 0.004

    num_workers = multiprocessing.cpu_count()   
    print(f"Used CPUs: {num_workers}")
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        downsampled_chunks = list(executor.map(downsample_dense_points, chunk_indices(dense_points_indices, num_workers), [xyz] * num_workers, [downsampling_distance] * num_workers))
    downsampled_dense_points_indices = np.concatenate(downsampled_chunks)
    
    mask = np.ones_like(x, dtype=bool)
    mask[dense_points_indices] = False; mask[downsampled_dense_points_indices] = True
    x = x[mask]; y = y[mask]; z = z[mask]    
    time_traj = (len(trajectory)*0.0004837768508)
    normalized_density = normalized_density[mask]

    #plot the oxygen atom trajectory with a color map representing the point density
    sc = ax.scatter(x, y, z, c=normalized_density, cmap='viridis', s=50)

    # Add color bar to show the density scale
    axins = inset_axes(ax,
                    width="5%",  # width of the colorbar
                    height="50%",  # height of the colorbar
                    loc='center right',  # location of the colorbar
                    #bbox_to_anchor=(1.1, 0., 1, 1),  # specify the anchor point
                    #bbox_to_anchor=(1, 0., 0.8, 0.8),  # specify the anchor point
                    bbox_transform=ax.transAxes,  # set the coordinate system
                    borderpad=-3,  # padding between the colorbar and the plot
                    )


    cbar = plt.colorbar(sc, cax=axins)
    #cbar = fig.colorbar(sc, cax=axins)
    cbar.set_label('normalized point density')
    
    #set a initaiton view point
    ax.view_init(elev=90, azim=0)

    ax.set_xlabel("x in Å")
    ax.set_ylabel('y in Å')
    #ax.set_zlabel('z in Å')
    ax.zaxis.set_ticklabels([])
    ax.zaxis.set_ticks([])

    # Plot histograms of distances moved by oxygen atoms relative to the initial frame
    #ax_hist = fig.add_subplot(122)
    #for i, water_coords in enumerate(water_oxy_molecules_coords):
    #    initial_coords = water_coords[0]
    #    distances = np.sqrt(np.sum((water_coords - initial_coords)**2, axis=1))
    #    print(f"mean displacement: {np.mean(distances):.2f} Å")
    #    ax_hist.hist(distances, bins=30, color="purple", alpha=0.8, label=f"O displacement")

    #ax_hist.set_xlabel("Displacement in Å")
    #ax_hist.set_ylabel("counts")
    #ax_hist.legend()

    #change figure size
    ax.grid(False)
    fig.set_size_inches(13, 5) 
      
    if create_fig == True:
        fig.savefig(f"{name_fig}.pdf")
    #plt.savefig("plot.pdf")
    plt.show()



def random_points_on_sphere_surface(center, radius, num_points):
    points = np.random.normal(size=(num_points, 3))
    points /= np.linalg.norm(points, axis=1)[:, np.newaxis]
    points *= radius
    points += center
    return points

def get_volume(coords):
    #compute the smallest convex polyhedron that encompasses all the points - atomes are modeld using the van der Waals radius
    if len(coords) < 4:  # At least 4 points are required to define a volume
        return 0

    hull = ConvexHull(coords)
    return hull.volume


def expanded_volume(trajectory, vdw_radius, num_points=100):

    # Generate random points on the surface of spheres around each atom
    all_points = []
    for point in trajectory:
        x, y, z = point[0], point[1], point[2]
        center = np.array([x, y, z])
        radius = vdw_radius
        points = random_points_on_sphere_surface(center, radius, num_points)
        all_points.extend(points)

    all_points = np.array(all_points)

    # Calculate the convex hull using the combined points
    convex_hull = ConvexHull(all_points)

    # Compute the volume of the convex hull
    volume = convex_hull.volume
    return volume, all_points

def get_alpha_shape_volume(coords, alpha):
    #calculate valume using a alpha shape while controling the alpha value wiht the average vdw radius
    if len(coords) < 4:  # At least 4 points are required to define a volume
        return 0

    # Calculate the alpha shape (concave hull)
    alpha_shape = alphashape.alphashape(coords, alpha)

    # Calculate the volume of the alpha shape
    return alpha_shape.volume

def voronoi_cell_volume(vor, idx):
    cell_vertices = vor.vertices[vor.regions[vor.point_region[idx]]]
    if not cell_vertices.size or len(cell_vertices) < 4:  # Empty, infinite cell, or not enough vertices
        return np.nan  # Return a special value to indicate an invalid volume
    return ConvexHull(cell_vertices).volume


def vol_KDA(water_oxy_molecules_coords):

    # Data preparation
    x = []
    y = []
    z = []
    for molecule_coords, color in zip([water_oxy_molecules_coords], ["purple"]):
        for coords in molecule_coords:
            x.extend(coords[:, 0])
            y.extend(coords[:, 1])
            z.extend(coords[:, 2])

    # Compute the point density of the oxygen atom trajectory using a kernel density estimation
    xyz = np.vstack([x, y, z]).T
    kde = KernelDensity(kernel='gaussian', bandwidth=0.5).fit(xyz)

    # Set a density threshold
    threshold = 0.05  # Adjust this value according to your application

    # Define the function to integrate
    def integrand(x, y, z):
        point_density = np.exp(kde.score_samples(np.column_stack((x, y, z))))
        return np.where(point_density >= threshold, 1, 0)

    # Define the integration limits for x, y, and z
    limits = ((min(x), max(x)), (min(y), max(y)), (min(z), max(z)))

    # Monte Carlo integration parameters
    n_samples = 10000  # Adjust this value based on the desired accuracy and computational resources

    # Generate random samples for the Monte Carlo integration
    random_samples = [np.random.rand(n_samples) * (high - low) + low for low, high in limits]

    total_volume = np.prod([high - low for low, high in limits])
    # Compute the occupied volume
    integral_values = integrand(*random_samples)
    occupied_volume = np.prod([high - low for low, high in limits]) * np.mean(integral_values)
    print("Occupied Volume MC-INT:", occupied_volume / total_volume)



def main(input_file, create_fig, name_fig, vis, n, k):


    output_file = "test.renamed.xyz"
    to_rename = [5, 4, 2, 3] #indexes corrisponding to the atom ordering like in chemcraft (starting from 1) CURRANTLY NOT USED
    n_atoms = 8 #number of atoms in the molecule

    rename_atoms_co(input_file, output_file, to_rename, 8)
    #rename_atoms_Ag(input_file, output_file, to_rename, 13)

    trajectory_filename = output_file
    original_trajectory_filename = input_file


    # Load the trajectory using MDAnalysis
    u = mda.Universe(trajectory_filename, format="XYZ", topology_format="XYZ")

    # Align the trajectory to the first frame using the water oxygen atoms
    aligner = align.AlignTraj(u, u, select="name tc", in_memory=True)
    aligner.run()

    # Convert the aligned trajectory back to a list
    trajectory = []
    for ts in u.trajectory:
        frame = []
        for atom, position in zip(u.atoms, ts.positions):
            element = atom.name
            frame.append([element] + list(position))
        trajectory.append(frame)

    print(f"Time in trajectory: {len(trajectory)*0.0004837768508}")

    #coordinates of tc atoms
    water_molecules_indexes = get_water_molecules_indexes(trajectory)
    water_molecules_coords = [get_water_coords(trajectory, index) for index in water_molecules_indexes]

    # Get the water oxygen molecules of the water
    water_oxy_molecules_indexes = [i for i, atom in enumerate(trajectory[0]) if atom[0] == "O"]
    water_oxy_molecules_coords = [get_water_coords(trajectory, index) for index in water_oxy_molecules_indexes]
    # Get hydrogen positions
    hydrogen_molecules_indexes = [i for i, atom in enumerate(trajectory[0]) if atom[0] == "H"]
    hydrogen_molecules_coords = [get_water_coords(trajectory, index) for index in hydrogen_molecules_indexes]
    #get all coords
    all_indexes = [i for i, atom in enumerate(trajectory[0]) if atom[0] == "O" or atom[0] == "tc"]
    all_coords = [get_water_coords(trajectory, index) for index in all_indexes]

    for i, coords in enumerate(water_oxy_molecules_coords):
        volume = get_volume(coords)
        print(f"Volume of oxygen (convex hall) {i + 1}: {volume:.3f} A^3")

    for i, coords in enumerate(water_oxy_molecules_coords):
        volume, all_points = expanded_volume(coords, 1.52)  #1.52 is the van der Waals radius of oxygen
        print(f"Expanded volume: \033[1m{volume:.3f} A^3\033[0m")

    
    #vol_KDA(water_oxy_molecules_coords)

        

    #get alpha volume
    #for i, coords in enumerate(water_oxy_molecules_coords):
    #    volume = get_alpha_shape_volume(coords, alpha=1)
    #    print(f"Volume of oxygen (alpha) {i + 1}: {volume:.3f} A^3")


    
    if vis == True:
        visualize(trajectory, water_molecules_coords, water_oxy_molecules_coords, hydrogen_molecules_coords, create_fig, name_fig, n, k)
    
    return volume



if __name__ == "__main__":

    #rename trajectory file to fix the atoms for rotational correction of the trajectory
    #input_file = ["m_1", "m_2"]
    #input_file = ["co4md/R3_40_1", "co4md/R3_40_2"]
    input_file = ["calc/40_tot.xyz"]
    for i in range(1, len(sys.argv)-1):
        input_file.append(str(sys.argv[i]))
    vis = True
    create_fig = True
    name_fig = str(sys.argv[-1])


    exp_vol = []

    n = len(input_file)
    k = 0
    for i in input_file:
        volume = main(i, create_fig, name_fig, vis, n, k)
        exp_vol.append(volume)
        try:
            print(f"Err in %: {abs((1 - (exp_vol[0] / exp_vol[1])) * 100)}")
        except:
            pass
        k = k+1

