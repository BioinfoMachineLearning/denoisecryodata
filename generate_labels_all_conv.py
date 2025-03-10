"""
@author: nabin

Generates labels for density maps for different tasks:
1. Regression (no part of map (zero part), type1, type 2) => type1 is for atoms found in PDB, type2 is for neighbors of atoms. 
values are copied from simulated map
2. Classification (no part of map (zero part), type 1, type 2 ) => (0, 1, 2)
3. Classification ( no part of map (zero part), Ca, CB, C, O, N ) => (0, 1, 2, 3, 4, 5)

"""

import math
import numpy as np
import mrcfile
import os
import subprocess
import warnings
from Bio import PDB
from copy import deepcopy
import sys
import shutil

from scipy.signal import convolve

chimera_path = '/home/ngzvh/bin/ChimeraX'

# Suppress all warnings
warnings.filterwarnings("ignore")



def get_index(cord, origin, voxel):
    return math.ceil(math.floor(cord - origin) / voxel)


def get_xyz(idx, voxel, origin):
    return (idx * voxel) - origin



def get_neighbors(indices, radius=1):
    i, j, k = indices
    x_range = np.arange(i - radius, i + radius + 1)
    y_range = np.arange(j - radius, j + radius + 1)
    z_range = np.arange(k - radius, k + radius + 1)

    neighbors = np.array(np.meshgrid(x_range, y_range, z_range)).T.reshape(-1, 3)
    neighbors = neighbors[~np.all(neighbors == [i, j, k], axis=1)]

    return neighbors.tolist()



def get_neighbors_1(indices):
    neighbors = list()
    
    radius = 12
    
    i, j, k = indices 
    for x in [i-radius, i, i+radius]:
        for y in [j-radius, j, j+radius]:
            for z in [k-radius, k, k+radius]:
                if (x, y, z) != (i, j, k):  
                    neighbors.append((x, y, z))
    return neighbors




def generate_label_type1(regression_label_path, classification_label_path, experimental_map, experimental_map_cord_idx, simulated_map_cord_idx ):

    common_cords = experimental_map_cord_idx.keys() & simulated_map_cord_idx.keys()
    common_cords_vals = {key: simulated_map_cord_idx[key] for key in common_cords}

    org_map = mrcfile.open(experimental_map, mode='r')
    regression_data = np.zeros(org_map.data.shape, dtype=np.int16)
    regression_data = regression_data.astype('float32')

    classification_data = np.zeros(org_map.data.shape, dtype=np.int16)
    
    x_origin = org_map.header.origin['x']
    y_origin = org_map.header.origin['y']
    z_origin = org_map.header.origin['z']
    x_voxel = org_map.voxel_size['x']
    y_voxel = org_map.voxel_size['y']
    z_voxel = org_map.voxel_size['z']

    for key, value in common_cords_vals.items():
        x, y, z = key
        iz = int(get_index(z, z_origin, z_voxel))
        jy = int(get_index(y, y_origin, y_voxel))
        kx = int(get_index(x, x_origin, x_voxel))
        if value != 0:
            regression_data[iz,jy,kx] = value
            classification_data[iz,jy,kx] = 1


    with mrcfile.new(regression_label_path, overwrite=True) as mrc:
        mrc.set_data(regression_data)
        mrc.voxel_size = org_map.voxel_size
        mrc.nstart = org_map.nstart
        mrc.header.origin = org_map.header.origin
        mrc.close()



    with mrcfile.new(classification_label_path, overwrite=True) as mrc:
        mrc.set_data(classification_data)
        mrc.voxel_size = org_map.voxel_size
        mrc.nstart = org_map.nstart
        mrc.header.origin = org_map.header.origin
        mrc.close()

    org_map.close()


def fill_zeros_near_nonzeros(reg, classification, fill_arr, radius=6):
    type2_count = 0
    kernel = np.ones((2*radius + 1, 2*radius + 1, 2*radius + 1), dtype=int)
    conv_result = convolve((reg != 0).astype(int), kernel, mode='same', method='direct')

    # Identify zero voxels near non-zero voxels using the convolution result
    zero_indices = np.argwhere((reg == 0) & (conv_result > 0))
    # print(zero_indices)

    # Fill zero voxels with corresponding values from sim map
    for i, j, k in zero_indices:
        try:
            # if fill_arr[i, j, k] != 0:
            reg[i, j, k] = fill_arr[i, j, k]
            classification[i, j, k] = 2
            type2_count += 1
        except IndexError:
            pass
    print("Type 2 Count:", type2_count)
    return reg, classification




def generate_label_type2(regression_label_path, classification_label_path, simulated_map, experimental_map_cord_idx, simulated_map_cord_idx):
    # common_cords_exp_sim = experimental_map_cord_idx.keys() & simulated_map_cord_idx.keys()
    # common_cords_vals_exp = {key: experimental_map_cord_idx[key] for key in common_cords_exp_sim}
    # labeled_voxels = np.array(list(common_cords_vals_exp.values()))
    
    sim_map = mrcfile.open(simulated_map, mode='r')
    sim_map_data = deepcopy(sim_map.data)
    

    regression_map = mrcfile.open(regression_label_path, mode='r+')
    classification_map = mrcfile.open(classification_label_path, mode='r+')
    regression_map_data = deepcopy(regression_map.data)
    classification_map_data = deepcopy(classification_map.data)

    # print(labeled_voxels)
    print("NUM NON ZERO IN SIM MAP:",len(np.nonzero(sim_map_data)[0]))
    print("NUM NON ZERO IN REGRESSION MAP TYPE 1: ", len(np.nonzero(regression_map_data)[0]))



    # Fill zero voxels near non-zero voxels in regression_map.data with corresponding values from sim_map_data
    regression_map_data1, classification_map_data1 = fill_zeros_near_nonzeros(regression_map_data, classification_map_data, sim_map_data, radius=6)

    regression_map.set_data(regression_map_data1) 
    classification_map.set_data(classification_map_data1)

    print("NUM NON ZERO IN REGRESSION MAP TYPE 1 and TYPE 2: ", len(np.nonzero(regression_map.data)[0]))

    # Close the maps
    sim_map.close()
    regression_map.close()
    classification_map.close()


def generate_label_classification_types(experimental_map, pdb_structure, classification_types_label_path):
    # read the map and structure files
    error_list = set()
    org_map = mrcfile.open(experimental_map, mode='r')
    data = np.zeros(org_map.data.shape, dtype=np.int16)
    data = data.astype('float32')
    x_origin = org_map.header.origin['x']
    y_origin = org_map.header.origin['y']
    z_origin = org_map.header.origin['z']
    x_voxel = org_map.voxel_size['x']
    y_voxel = org_map.voxel_size['y']
    z_voxel = org_map.voxel_size['z']
    parser = PDB.PDBParser()
    struct = parser.get_structure("CA", pdb_structure)
    for model in struct:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    x, y, z = atom.get_coord()
                    iz = int(get_index(z, z_origin, z_voxel))
                    jy = int(get_index(y, y_origin, y_voxel))
                    kx = int(get_index(x, x_origin, x_voxel))
                    if atom.get_name() == "CA":
                        try:
                            data[iz, jy, kx] = 1
                        except IndexError as error:
                            error_list.add([iz, jy, kx])
                    elif atom.get_name() == "CB":
                        try:
                            data[iz, jy, kx] = 2
                        except IndexError as error:
                            error_list.add([iz, jy, kx])
                    elif atom.get_name() == "C":
                        try:
                            data[iz, jy, kx] = 3
                        except IndexError as error:
                            error_list.add([iz, jy, kx])
                    elif atom.get_name() == "O":
                        try:
                            data[iz, jy, kx] = 4
                        except IndexError as error:
                            error_list.add([iz, jy, kx])
                    elif atom.get_name() == "N":
                        try:
                            data[iz, jy, kx] = 5
                        except IndexError as error:
                            error_list.add([iz, jy, kx])
                    else:
                        pass

    
    with mrcfile.new(classification_types_label_path, overwrite=True) as mrc:
        mrc.set_data(data)
        mrc.voxel_size = x_voxel
        mrc.header.origin = org_map.header.origin
        mrc.nstart = org_map.nstart
        mrc.close()
    if len(error_list) !=0:
        print('Index Errors occurred at the following indices:')
        print(error_list)


def get_pdb_indices(density_map, pdb_structure):
    cord_idx = dict()
    open_density_map = mrcfile.open(density_map, mode='r')
    x_origin = open_density_map.header.origin['x']
    y_origin = open_density_map.header.origin['y']
    z_origin = open_density_map.header.origin['z']
    x_voxel = open_density_map.voxel_size['x']
    y_voxel = open_density_map.voxel_size['y']
    z_voxel = open_density_map.voxel_size['z']

    parser = PDB.PDBParser()
    struct = parser.get_structure("PROTEIN", pdb_structure)
    for model in struct:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    x, y, z = atom.get_coord()
                    iz = int(get_index(z, z_origin, z_voxel))
                    jy = int(get_index(y, y_origin, y_voxel))
                    kx = int(get_index(x, x_origin, x_voxel))
                    cord_idx[(x,y,z)] = (iz, jy, kx)
    return cord_idx

def get_simulated_map_values(density_map, pdb_structure):
    index_error = 0
    count = 0
    cord_vals = dict()
    open_density_map = mrcfile.open(density_map, mode='r')
    density_map_data = deepcopy(open_density_map.data)

    x_origin = open_density_map.header.origin['x']
    y_origin = open_density_map.header.origin['y']
    z_origin = open_density_map.header.origin['z']
    x_voxel = open_density_map.voxel_size['x']
    y_voxel = open_density_map.voxel_size['y']
    z_voxel = open_density_map.voxel_size['z']

    parser = PDB.PDBParser()
    struct = parser.get_structure("CA", pdb_structure)
    for model in struct:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    x, y, z = atom.get_coord()
                    count += 1
                    iz = int(get_index(z, z_origin, z_voxel))
                    jy = int(get_index(y, y_origin, y_voxel))
                    kx = int(get_index(x, x_origin, x_voxel))
                    try:
                        cord_vals[(x,y,z)] = density_map_data[iz, jy, kx]
                    except IndexError:
                        index_error += 1

    print(f"Total Index Error : {index_error} out of {count}")
    open_density_map.close()
    return cord_vals

def run_chimera(experimental_map, simulated_map, path_save_name, label_save_path, m):
    sim_map = mrcfile.open(simulated_map, mode='r+')
    exp_map = mrcfile.open(experimental_map, mode='r')
    
    sim_map.header.origin = exp_map.header.origin
    sim_map.nstart = exp_map.nstart
    sim_map.voxel_size = exp_map.voxel_size
    # sim_map.flush()
    sim_map.close()

    
    chimera_resample_script = f'{label_save_path}/{m}_resample.cxc'
    chimera_scripts = open(f'{chimera_resample_script}', 'w')
    chimera_scripts.write('open ' + simulated_map + '\n'
                            'open ' + experimental_map + '\n'
                            'volume #1 step 1' + '\n'
                            'volume #2 step 1' + '\n'
                            'volume resample #1 onGrid #2\n'
                            'fitmap #3 inMap #2\n'
                            'volume #3 step 1' + '\n'
                            'save '+ path_save_name + ' model #3 \n'
                            'exit')
    
    chimera_scripts.close()
    script_finished = False
    while not script_finished:
        try:
            subprocess.run([chimera_path, '--nogui', chimera_scripts.name])
            script_finished = True
        except FileNotFoundError as error:
            raise error
        
    os.remove(chimera_resample_script)


def delete_directory(path):
    """Deletes the specified directory if it exists."""
    if os.path.exists(path) and os.path.isdir(path):
        shutil.rmtree(path)
        print(f"Directory '{path}' has been deleted.")
    else:
        print(f"Directory '{path}' does not exist.")


def delete_file(file_path):
    """Deletes a file if it exists and is a file."""
    if os.path.exists(file_path) and os.path.isfile(file_path):
        os.remove(file_path)

def create_directory(name):
    """Creates a directory with the given name if it doesn't exist."""
    if not os.path.exists(name):
        os.makedirs(name)
        print(f"Directory '{name}' has been created.")
    else:
        print(f"Directory '{name}' already exists.")


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: pass in atleast one density map")
        sys.exit(1)
    else:
        m = sys.argv[1]
  

    label_save_path = "/denoise_cryoem"

    
    experimental_map = f"/experimental_pdb/{m}/{m}.mrc"
    simulated_map = f"/situs_simulated_pdb/{m}/{m}_situs_simulated.mrc"
    pdb_structure = f"/experimental_pdb/{m}/{m}.pdb1"
    # path_save_name = f"/cluster/pixstor/chengji-lab/nabin/denoise_cryoem/simulated_map/{m}_sim_fixed_3.mrc"

    regression_label_path = f"{label_save_path}/regression_labels_situs/{m}_regression_situs.mrc"
    classification_label_path = f"{label_save_path}/classification_labels_situs/{m}_classification_situs.mrc"
    classification_types_label_path = f"{label_save_path}/classification_types_labels_situs/{m}_classification_types_situs.mrc"

    delete_file(regression_label_path)
    delete_file(classification_label_path)
    delete_file(classification_types_label_path)
    

    
    # run_chimera(experimental_map, simulated_map, path_save_name, label_save_path, m)

    # simulated_map = path_save_name
    # get corresponding indices from pdb_structure 
    experimental_map_cord_idx = get_pdb_indices(experimental_map, pdb_structure)
    # print(experimental_map_cord_idx)
    print("Total Experimental map atoms :",   len(experimental_map_cord_idx))
    simulated_map_cord_idx = get_simulated_map_values(simulated_map, pdb_structure)
    
    print("Total Simulated map atoms :",  len(simulated_map_cord_idx))
    print("Total type1 atoms :",  len(simulated_map_cord_idx))



    generate_label_type1(regression_label_path, classification_label_path, experimental_map, experimental_map_cord_idx, simulated_map_cord_idx )
    generate_label_type2(regression_label_path, classification_label_path, simulated_map, experimental_map_cord_idx, simulated_map_cord_idx)
    generate_label_classification_types(experimental_map, pdb_structure, classification_types_label_path)

    # os.remove(path_save_name)

    print("Done", m)
    print("COMPLETED!")



