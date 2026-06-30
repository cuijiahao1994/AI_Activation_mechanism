import hashlib
import os
from ase import io
from ase.io import read
from ase import Atoms
from collections import defaultdict
import itertools

co_replacement = list(sorted(set(['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Al', 'Ga', 'Sc', 'In', 'Tl', 'Si', 'Ge', 'Sn', 'Pb', 'As', 'Sb', 'Bi', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac'])))
fe_replacement = list(sorted(set(['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Al', 'Ga', 'Sc', 'In', 'Tl', 'Si', 'Ge', 'Sn', 'Pb', 'As', 'Sb', 'Bi', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac'])))
c_replacement = list(sorted(set(['N', 'C', 'O', 'S', 'B', 'P'])))
s_replacement = list(sorted(set(['B'])))
p_replacement = list(sorted(set(['N', 'C', 'O', 'S', 'B', 'P'])))
o_replacement = list(sorted(set(['N'])))

pseudo_dir = "/data/home/cuijiahao/potpaw_PBE.64"
potcar_set = "POT_PBE"


transition_metals = {
    
    'Sc': 4.0, 'Ti': 4.5, 'V': 3.5, 'Cr': 3.7, 'Mn': 4.5, 
    'Fe': 4.0, 'Co': 3.5, 'Ni': 6.0, 'Cu': 5.0, 'Zn': 0.0,
    
   
    'Y': 3.0, 'Zr': 4.0, 'Nb': 0, 'Mo': 4.0, 'Tc': 0,
    'Ru': 3.0, 'Rh': 3.0, 'Pd': 3.0, 'Ag': 3.0, 'Cd': 0.0,
    
   
    'Hf': 4.0, 'Ta': 0, 'W': 0.0, 'Re': 0, 'Os': 3.0,
    'Ir': 3.0, 'Pt': 3.0, 'Au': 3.0,
    
  
    'La': 6.0, 'Ce': 5.0, 'Pr': 5.5, 'Nd': 6.0, 'Pm': 0.0,
    'Sm': 6.5, 'Eu': 0.0, 'Gd': 0, 'Tb': 0.0, 'Dy': 0,
    'Ho': 6.0, 'Er': 0, 'Tm': 0.0, 'Yb': 0, 'Lu': 4.0
}


rare_earth_metals = ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
                    'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac']


atoms = io.read('POSCAR')


element_indices = {
    'Co': [i for i, atom in enumerate(atoms) if atom.symbol == 'Co'],
    'Fe': [i for i, atom in enumerate(atoms) if atom.symbol == 'Fe'],
    'C':  [i for i, atom in enumerate(atoms) if atom.symbol == 'C'],
    'S':  [i for i, atom in enumerate(atoms) if atom.symbol == 'S'],
    'P':  [i for i, atom in enumerate(atoms) if atom.symbol == 'P'],
    'O':  [i for i, atom in enumerate(atoms) if atom.symbol == 'O']
}


active_elements = []
for symbol in ['Co', 'Fe', 'C', 'S', 'P', 'O']:
    if element_indices[symbol]:
        replacement = globals()[f"{symbol.lower()}_replacement"]
        active_elements.append((symbol, replacement, element_indices[symbol]))


generated_hashes = set()

def get_structure_hash(atoms):
    atoms_data = [(atom.symbol, tuple(atom.position)) for atom in atoms]
    atoms_str = str(sorted(atoms_data))
    return hashlib.sha256(atoms_str.encode()).hexdigest()

def save_poscar_grouped(atoms, filename):
    elements = sorted(set(atoms.get_chemical_symbols()))
    element_counts = {e: atoms.get_chemical_symbols().count(e) for e in elements}
    
    with open(filename, 'w') as f:
        f.write("ASE Generated\n1.0\n")
        for vec in atoms.get_cell():
            f.write(f"{vec[0]:20.16f} {vec[1]:20.16f} {vec[2]:20.16f}\n")
        f.write(" ".join(elements) + "\n")
        f.write(" ".join(map(str, [element_counts[e] for e in elements])) + "\n")
        f.write("Cartesian\n")
        for elem in elements:
            for atom in atoms:
                if atom.symbol == elem:
                    f.write(f"{atom.position[0]:20.16f} {atom.position[1]:20.16f} {atom.position[2]:20.16f}\n")
    return elements


replace_lists = [elem[1] for elem in active_elements]
for combination in itertools.product(*replace_lists):
    atoms_copy = atoms.copy()
    folder_parts = []
    

    for (sym, _, indices), new_elem in zip(active_elements, combination):
        for idx in indices:
            atoms_copy[idx].symbol = new_elem
        folder_parts.append(new_elem)
    

    folder_name = "_".join(folder_parts)
    

    struct_hash = get_structure_hash(atoms_copy)
    if struct_hash in generated_hashes:
        print(f"Skip repetitive structures: {folder_name}")
        continue
    generated_hashes.add(struct_hash)
    

    os.makedirs(folder_name, exist_ok=True)
    

    poscar_path = os.path.join(folder_name, "POSCAR")
    elements_order = save_poscar_grouped(atoms_copy, poscar_path)
    

    potcar_path = os.path.join(folder_name, "POTCAR")
    with open(potcar_path, 'w') as potcar:
        missing = False
        for elem in elements_order:
            pseudo_path = os.path.join(pseudo_dir, potcar_set, elem, "POTCAR")
            if os.path.exists(pseudo_path):
                with open(pseudo_path) as f:
                    potcar.write(f.read())
            else:
                print(f"Warning: {elem}的POTCAR not exist in{pseudo_path}")
                missing = True
        if missing:
            os.remove(potcar_path)  
    

    incar_path = os.path.join(folder_name, "INCAR")
    with open(incar_path, 'w') as incar:
        incar_content = """SYSTEM = Generated
ENCUT = 400
ISTART = 1
ICHARG = 1
LCHARG = F
LWAVE =  F
LVTOT = F
LVHAR = F
LELF = F
ISMEAR = 0
SIGMA = 0.05
NSW = 500
IBRION = 2
LREAL = Auto
ISIF = 2
EDIFF = 1E-7
EDIFFG = -5E-2
ISPIN = 2
IVDW = 11
LDAU = T
LDAUTYPE = 2
AMIX = 0.2
BMIX = 0.0001
AMIX_MAG = 0.8
BMIX_MAG = 0.0001
NELMIN = 5 
NELM = 300
GGA = PE 
LREAL = Auto 
ALGO = Fast 
POTIM = 0.2
LDAU = .TRUE.
LDAUTYPE = 2
"""
  
        ldaul = []
        ldauu = []
        ldauj = []
        for elem in elements_order:
            if elem in rare_earth_metals:
            
                ldaul.append(3)
                ldauu.append(transition_metals.get(elem, 6.0))  
                ldauj.append(0)
            elif elem in transition_metals and elem not in rare_earth_metals:
                
                ldaul.append(2)
                ldauu.append(transition_metals[elem])
                ldauj.append(0)
            else:
         
                ldaul.append(-1)
                ldauu.append(0)
                ldauj.append(0)

        incar_content += f"LDAUL = {' '.join(map(str, ldaul))}\n"
        incar_content += f"LDAUU = {' '.join(map(str, ldauu))}\n"
        incar_content += f"LDAUJ = {' '.join(map(str, ldauj))}\n"
        incar.write(incar_content)
        

    kpoints_path = os.path.join(folder_name, "KPOINTS")
    with open(kpoints_path, 'w') as kpoints:
        kpoints.write("Automatic mesh\n0\nMonkhorst-Pack\n1 1 1\n0 0 0\n")
    
    
    vasp_script = os.path.join(folder_name, "vasp.sh")
    with open(vasp_script, 'w') as script:
        script.write("""#!/bin/bash
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --ntasks-per-node=56
#SBATCH --partition=partition
#SBATCH --time=24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

source /srv/nfs/share/opt/intel/oneapi-2025.0/setvars.sh
export PATH=/srv/nfs/share/app/vasp.6.4.3/bin:$PATH
ulimit -s unlimited
mpirun -np $SLURM_NPROCS vasp_gam
""")
    
    print(f"Successfully generated: {folder_name}")

print("All the structures have been generated successfully.")