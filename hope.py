import argparse
from Bio.PDB import PDBParser, DSSP
import numpy as np
import math
import pymol
import copy
import sys

# Constantes pour la liste des acides aminés hydrophobes
HYDROPHOBICS_AMINO_ACIDS = ['PHE', 'GLY', 'ILE', 'LEU', 'MET', 'VAL', 'TRP', 'TYR']

def create_amino_acid(code, id_aa, x, y, z):
    """Créer un acide aminé et vérifier s'il est hydrophobe."""
    is_hydrophobic = code in HYDROPHOBICS_AMINO_ACIDS
    point = np.array([x, y, z])
    return {
        'id': id_aa,
        'code': code,
        'is_hydrophobic': is_hydrophobic,
        'point': point
    }

def create_protein(name, mass_center, amino_acid_sequence, full_sequence, best_positions):
    return {
        "name": name, 
        "mass_center": None, 
        "amino_acid_sequence": amino_acid_sequence, 
        "full_sequence": full_sequence, 
        "best_positions": best_positions
    }

def find_best_axis(protein):
    best_axis_val = 0
    best_axis_found = None
    for axis in protein["best_positions"]:
        if axis["best_hydrophobicity"] > best_axis_val:
            best_axis_val = axis["best_hydrophobicity"]
            best_axis_found = copy.deepcopy(axis)
    return best_axis_found

def compute_mass_center(chain):
    result = chain.center_of_mass()
    x, y, z = result[0], result[1], result[2]
    return np.array([x, y, z])

def check_input_file(input_file):
    """Vérifier si le fichier d'entrée est un fichier PDB valide."""
    if not input_file.lower().endswith(".pdb"):
        print("Wrong extension. A PDB file is required.")
        sys.exit(0)

    try:
        parser = PDBParser()
        structure = parser.get_structure("temp", input_file)
        for model in structure:
            for chain in model:
                if chain.id == 'A':
                    return True
    except Exception:
        print("Issue when attempting to read the PDB file. Please check the content of your input file.")
        sys.exit(0)


def calculate_solvant_accessibility(structure, input_file):
    """Calculer l'accessibilité au solvant des résidus avec DSSP."""
    print("Computing solvent accessibility...")
    model = structure[0]
    dssp = DSSP(model, input_file, dssp='dssp')
    return dssp


def parse_pdb(input_file, chain):
    """Analyser le fichier PDB pour construire les informations de la protéine."""
    print("Parsing the PDB file ...")
    p = PDBParser()
    structure = p.get_structure("structure", input_file)

    id_amino_acid = 1  # Counting the accessible residues
    id_full_amino_acid = 1  # Counting all residues
    # Parsing the name
    name = input_file[-8:]
    name = name[:-4]
    protein = create_protein(name=name, mass_center=None,  # Ce sera calculé plus tard
        amino_acid_sequence=[],  # Initialiser la séquence vide
        full_sequence=[],  # Initialiser la séquence complète vide
        best_positions=[])

    # Calculer l'accessibilité au solvant
    dssp_res = calculate_solvant_accessibility(structure, input_file= input_file)

    print("Trimming to get only exposed residues...")
    exposed_residues = []
    for res_id, _, _, asa, _, _, _, _, _, _, _, _, _, _ in dssp_res:
        if asa > 0.3:
            exposed_residues.append(res_id)
    print(f"Found {len(exposed_residues)} exposed residues.")

    model = structure[0]
    chain_selected = model["A"]

    # Calcul du centre de masse
    protein['mass_center'] = compute_mass_center(chain_selected)
    print(f"Mass center is {protein['mass_center']}")

    # Création des acides aminés
    for residue in chain_selected:
        if residue.has_id("CA"):
            atom = residue['CA']
            x, y, z = atom.get_coord()
            new_amino_acid = create_amino_acid(code=residue.get_resname(), id_aa=residue.get_id()[1], x=x, y=y, z=z)
            protein['full_sequence'].append(new_amino_acid)
            id_full_amino_acid += 1
            if residue.get_id()[1] in exposed_residues:
                protein['amino_acid_sequence'].append(new_amino_acid)
                id_amino_acid += 1
    return protein


def get_x(point):
    """Renvoie la coordonnée x d'un point"""
    return point[0]

def get_y(point):
    """Renvoie la coordonnée y d'un point"""
    return point[1]

def get_z(point):
    """Renvoie la coordonnée z d'un point"""
    return point[2]

def find_points(n_points, mass_center):
    """Générer des points répartis uniformément sur une demi-sphère."""
    points = []
    phi = 0
    for k in range(1, n_points + 1):
        h = -1 + (2 * (k - 1) / (n_points - 1))
        if h != 1:
            theta = math.acos(h)
            if k != 1 and k != n_points:
                phi = (phi + (3.6 / math.sqrt(n_points) * (1 / math.sqrt(1 - h * h)))) % (2 * math.pi)

            x = mass_center[0] + math.sin(phi) * math.sin(theta)
            y = mass_center[1] + math.cos(theta)
            z = mass_center[2] + math.cos(phi) * math.sin(theta)
            points.append(np.array([x, y, z]))
    # Get only the points from half circle :
    above_x_axis_points = [point for point in points if point[2] > mass_center[2]]
    return above_x_axis_points


def find_director_vector(point, center_coordinate):
    """Trouver le vecteur directeur reliant deux points."""
    return np.array([center_coordinate[0] - point[0],
                     center_coordinate[1] - point[1],
                     center_coordinate[2] - point[2]])


def define_plane(point, normal_vector):
    """Définir un plan à partir d'un point et d'un vecteur normal."""
    a = normal_vector[0]
    b = normal_vector[1]
    c = normal_vector[2]
    d = -(a * point[0] + b * point[1] + c * point[2])  # https://mathworld.wolfram.com/Plane.html
    return np.array([a, b, c, d])

def get_plane_equation(plane):
    """Renvoie l'équation d'un plan"""
    return f"{plane[0]:.3f}x + {plane[1]:.3f}y + {plane[2]:.3f}z + {plane[3]:.3f} = 0"

def complementary_plane(plane, gap):
    """Crée un plan parallèle décalé d'un certain écart"""
    new_plane = list(plane)
    new_plane[3] += gap
    return tuple(new_plane)

def slide_plane(plane, sliding_window):
    """Décale un plan d'une certaine valeur"""
    plane = list(plane)
    plane[3] += sliding_window
    return tuple(plane)

def is_point_above_plane(point, plane):
    """Vérifier si un point est au-dessus d'un plan."""
    return np.dot(plane[:3], point) + plane[3] > 0


def is_point_below_plane(point, plane):
    """Vérifier si un point est en dessous d'un plan."""
    return np.dot(plane[:3], point) + plane[3] < 0

def create_axis(plane1, plane2):
    """Crée un axe"""
    axis = {
        'plane1': plane1,
        'plane2': plane2,
        'best_hydrophobicity': -1000,  # Hydrophobicité initialement faible
        'best_ratio_of_atoms': 0
    }
    return axis

def explore_axis(amino_acid_sequence, plane1, plane2, ref_axe):
    """Explorer un axe et trouver la meilleure hydrophobicité relative."""
    in_between_planes = []
    n_total_hydrophobic = 0
    n_total_hydrophile = 0
    n_hydrophobe_in_plan = 0
    nb_hydrophile_out_of_plan = 0

    for aa in amino_acid_sequence:
        if is_point_below_plane(aa['point'], plane1) and is_point_above_plane(aa['point'], plane2):
            in_between_planes.append(aa)

        if aa['is_hydrophobic']:
            n_total_hydrophobic += 1
        else:
            n_total_hydrophile += 1

    for aa in amino_acid_sequence:
        if aa not in in_between_planes and not aa['is_hydrophobic']:
            nb_hydrophile_out_of_plan += 1

    for aa in in_between_planes:
        if aa['is_hydrophobic']:
            n_hydrophobe_in_plan += 1

    if not in_between_planes or len(in_between_planes) >= len(amino_acid_sequence):
        return False

    hydrophobicity = (nb_hydrophile_out_of_plan / n_total_hydrophile) + (n_hydrophobe_in_plan / n_total_hydrophobic)
    if hydrophobicity > ref_axe['best_hydrophobicity']:
        # Updating the "best" match
        ref_axe['best_ratio_of_atoms'] = len(in_between_planes)
        ref_axe['best_hydrophobicity'] = hydrophobicity
        ref_axe['plane1'] = copy.deepcopy(axis['plane1'])
        ref_axe['plane2'] = copy.deepcopy(axis['plane2'])
        return True
    else :      
        return False


def show_in_pymol(plane1, plane2, pdb_file, mass_center):
    """Afficher la protéine et les plans dans PyMol."""
    pymol.finish_launching(['pymol', '-q'])
    pymol.cmd.load(pdb_file, "protein")
    pymol.cmd.remove("solvent")

    x_min, x_max = float("inf"), float("-inf")
    y_min, y_max = float("inf"), float("-inf")

    for atom in pymol.cmd.get_model("protein").atom:
        x = atom.coord[0]
        if x < x_min:
            x_min = x
        if x > x_max:
            x_max = x

    for atom in pymol.cmd.get_model("protein").atom:
        y = atom.coord[1]
        if y < y_min:
            y_min = y
        if y > y_max:
            y_max = y

    step = 3

    points_on_plane1 = []
    for x in np.arange(x_min, x_max + step, step):
        for y in np.arange(y_min, y_max + step, step):
            # Calculate z coordinate using the plane equation for plane 1
            z1 = (-plane1[0] * x - plane1[1] * y - plane1[3]) / plane1[2]
            points_on_plane1.append((x, y, z1))

    # Generate points on plane 2
    points_on_plane2 = []
    for x in np.arange(x_min, x_max + step, step):
        for y in np.arange(y_min, y_max + step, step):
            # Calculate z coordinate using the plane equation for plane 2
            z2 = (-plane2[0] * x - plane2[1] * y - plane2[3]) / plane2[2]
            points_on_plane2.append((x, y, z2))

    # Create pseudoatoms for points on plane 1
    for idx, point in enumerate(points_on_plane1):
        x, y, z = point
        atom_name = f"plane1_{idx}"
        pymol.cmd.pseudoatom(atom_name, pos=point, color="white")
        pymol.cmd.show("spheres", f"plane1_{idx}")

    # Create pseudoatoms for points on plane 2
    for idx, point in enumerate(points_on_plane2):
        x, y, z = point
        atom_name = f"plane2_{idx}"
        pymol.cmd.pseudoatom(atom_name, pos=point, color="white")
        pymol.cmd.show("spheres", f"plane2_{idx}")

    # Mass center
    if mass_center is not None:
        atom_name = "mass_center"
        pymol.cmd.pseudoatom(atom_name, pos=[mass_center[0], mass_center[1], mass_center[2]],
                             color="magenta")
        
    # Show the protein structure
    pymol.cmd.show("spheres", atom_name)
    pymol.cmd.show("cartoon", "protein")

def optimizing_width_membrane(gap_membrane, axis_init, amino_acid_sequence, plane_to_consider):
    """Optimiser la largeur de la membrane en explorant au-dessus et en dessous d'un plan."""
    best_axis = copy.deepcopy(axis_init)
    if plane_to_consider == 1:
        axis_init['plane1'] = slide_plane(axis_init['plane1'], gap_membrane)
    else:
        axis_init['plane2'] = slide_plane(axis_init['plane2'], gap_membrane)

    while True:
        hydrophobicity = explore_axis(amino_acid_sequence, axis_init['plane1'], axis_init['plane2'], ref_axe=axis_init)
        if hydrophobicity > best_axis['best_hydrophobicity']:
            best_axis['best_hydrophobicity'] = hydrophobicity
            best_axis['plane1'] = copy.deepcopy(axis_init['plane1'])
            best_axis['plane2'] = copy.deepcopy(axis_init['plane2'])

        # Slide the chosen plane again
        if plane_to_consider == 1:
            axis_init['plane1'] = slide_plane(axis_init['plane1'], gap_membrane)
        else:
            axis_init['plane2'] = slide_plane(axis_init['plane2'], gap_membrane)

        # Exit condition if no improvement
        if not hydrophobicity:
            break

    return best_axis

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='hope.py',
        description='This program locates the membrane in a transmembrane protein and detects transmembrane segments')
    parser.add_argument('filename', help="A PDB file")
    
    # Optionnal arguments :
    parser.add_argument('-n', help="Number of points to place on the sphere. (default is 15)", type=int, dest="n", default=15) 
    parser.add_argument('-w', help="Initial width of the membrane. (default is 14 A)", type=float, dest="width", default=14) 
    parser.add_argument('-g', help="Gap of sliding membrane along an axis. (default is 1 A)", type=float, dest="gap", default=1) 
    parser.add_argument('-m', help="Gap of optimising membrane's width. (default is 1 A)", type=float, dest="gap_membrane", default=1) 
    
    # Parsing : 
    args = parser.parse_args()
    filename = args.filename
    n = args.n
    chain = "A"
    width = args.width
    gap = args.gap
    gap_membrane = args.gap_membrane
      
    print(f"Command : python hope.py {filename} -n {n} -w {width} -g {gap} -m {gap_membrane}")

    # Vérification du fichier d'entrée
    check_input_file(filename)

    # Parsing du fichier PDB
    protein = parse_pdb(filename, chain)
    
    # Génération des points sur une demi-sphère
    directions = find_points(2 * n, protein['mass_center'])
    
    print("Calculating the planes... ")
    for d in directions:
        point = copy.deepcopy(d)
        normal = find_director_vector(point=point, center_coordinate=protein['mass_center'])
        plane1 = define_plane(point=point, normal_vector=normal)
        plane2 = complementary_plane(plane1, gap=width)

        # Create an axis for exploration
        axis = create_axis(plane1, plane2)
        best_axis_tmp = copy.deepcopy(axis)
        
        # Exploration de l'axe
        while explore_axis(protein['amino_acid_sequence'], axis["plane1"], axis["plane2"], ref_axe=best_axis_tmp) :
            axis = create_axis(slide_plane(axis["plane1"], gap), slide_plane(axis["plane2"], gap))

        # Save the best axis position
        protein['best_positions'].append(best_axis_tmp)
    
    # Trouver le meilleur axe
    best_axis = find_best_axis(protein)
    if best_axis is None:
        print("Something went wrong... Try to modify the values of parameters.")
        sys.exit(1)
    
    print(f"Best axis found: {best_axis}")

    print("Optimising membrane width...")
    best_axis_tmp = optimizing_width_membrane(gap_membrane=gap_membrane, axis_init=best_axis, amino_acid_sequence=protein['amino_acid_sequence'], plane_to_consider=2)
    best_axis_tmp2 = optimizing_width_membrane(gap_membrane=-gap_membrane, axis_init=best_axis_tmp, amino_acid_sequence=protein['amino_acid_sequence'], plane_to_consider=2)
    best_axis_tmp3 = optimizing_width_membrane(gap_membrane=gap_membrane, axis_init=best_axis_tmp2, amino_acid_sequence=protein['amino_acid_sequence'], plane_to_consider=1)
    best_axis_tmp4 = optimizing_width_membrane(gap_membrane=-gap_membrane, axis_init=best_axis_tmp3, amino_acid_sequence=protein['amino_acid_sequence'], plane_to_consider=1)

    print("The membrane width is ", abs(best_axis_tmp4['plane1'][3] - best_axis_tmp4['plane2'][3]))
    print("Best axis found overall is", best_axis_tmp4)
    
    # Visualiser dans PyMol
    show_in_pymol(best_axis_tmp4['plane1'], best_axis_tmp4['plane2'], filename, mass_center=protein['mass_center'])