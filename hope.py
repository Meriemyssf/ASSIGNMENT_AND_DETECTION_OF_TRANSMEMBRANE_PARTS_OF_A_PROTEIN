import argparse
from Bio.PDB import PDBParser, DSSP
import numpy as np
import pymol
import copy
import sys

# Constantes pour la liste des acides aminés hydrophobes
HYDROPHOBICS_AMINO_ACIDS = ['PHE', 'GLY', 'ILE', 'LEU', 'MET', 'VAL', 'TRP', 'TYR']


# ======================== Fonctions de Protein.py =========================

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


def compute_mass_center(chain):
    """Calculer le centre de masse d'une chaîne de protéines à partir des atomes CA."""
    total_mass = 0
    weighted_sum = np.zeros(3)

    atomic_masses = {'C': 12.01}

    for residue in chain:
        if residue.has_id("CA"):
            atom = residue['CA']
            mass = atomic_masses.get(atom.element, 1.0)
            total_mass += mass
            weighted_sum += mass * atom.get_coord()

    if total_mass == 0:
        raise ValueError("Aucun atome Cα trouvé dans la chaîne, impossible de calculer le centre de masse.")

    center_of_mass = weighted_sum / total_mass
    return center_of_mass


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
    dssp = DSSP(model, input_file, dssp='mkdssp')
    return dssp


def parse_pdb(input_file, chain='A'):
    """Analyser le fichier PDB pour construire les informations de la protéine."""
    print("Parsing the PDB file ...")
    p = PDBParser()
    structure = p.get_structure("structure", input_file)

    protein = {
        'name': input_file.split('/')[-1].split('.')[0],
        'mass_center': None,
        'amino_acid_sequence': [],
        'full_sequence': [],
        'best_positions': []
    }

    # Calculer l'accessibilité au solvant
    dssp_res = calculate_solvant_accessibility(structure, input_file)

    print("Trimming to get only exposed residues...")
    exposed_residues = [res_id for res_id, _, _, asa, *_ in dssp_res if asa > 0.3]
    print(f"Found {len(exposed_residues)} exposed residues.")

    model = structure[0]
    chain_selected = model[chain]

    # Calcul du centre de masse
    protein['mass_center'] = compute_mass_center(chain_selected)
    print(f"Mass center is {protein['mass_center']}")

    # Création des acides aminés
    for residue in chain_selected:
        if residue.has_id("CA"):
            atom = residue['CA']
            x, y, z = atom.get_coord()
            new_amino_acid = create_amino_acid(residue.get_resname(), residue.get_id()[1], x, y, z)
            protein['full_sequence'].append(new_amino_acid)
            
            if residue.get_id()[1] in exposed_residues:
                protein['amino_acid_sequence'].append(new_amino_acid)

    return protein


# ===================== Fonctions de Geometry.py ===========================

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

    return [point for point in points if point[2] > mass_center[2]]


def find_director_vector(point, center_coordinate):
    """Trouver le vecteur directeur reliant deux points."""
    return np.array([center_coordinate[0] - point[0],
                     center_coordinate[1] - point[1],
                     center_coordinate[2] - point[2]])


def define_plane(point, normal_vector):
    """Définir un plan à partir d'un point et d'un vecteur normal."""
    a, b, c = normal_vector
    d = -(a * point[0] + b * point[1] + c * point[2])
    return np.array([a, b, c, d])


def is_point_above_plane(point, plane):
    """Vérifier si un point est au-dessus d'un plan."""
    return np.dot(plane[:3], point) + plane[3] > 0


def is_point_below_plane(point, plane):
    """Vérifier si un point est en dessous d'un plan."""
    return np.dot(plane[:3], point) + plane[3] < 0


def explore_axis(amino_acid_sequence, plane1, plane2, best_hydrophobicity):
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
        return best_hydrophobicity, False

    hydrophobicity = (nb_hydrophile_out_of_plan / n_total_hydrophile) + (n_hydrophobe_in_plan / n_total_hydrophobic)
    if hydrophobicity > best_hydrophobicity:
        return hydrophobicity, True

    return best_hydrophobicity, False


def show_in_pymol(plane1, plane2, pdb_file, mass_center):
    """Afficher la protéine et les plans dans PyMol."""
    pymol.finish_launching(['pymol', '-q'])
    pymol.cmd.load(pdb_file, "protein")
    pymol.cmd.remove("solvent")

    x_min, x_max = float("inf"), float("-inf")
    y_min, y_max = float("inf"), float("-inf")

    for atom in pymol.cmd.get_model("protein").atom:
        x, y = atom.coord[:2]
        x_min, x_max = min(x_min, x), max(x_max, x)
        y_min, y_max = min(y_min, y), max(y_max, y)

    step = 3
    points_on_plane1 = [(x, y, (-plane1[0] * x - plane1[1] * y - plane1[3]) / plane1[2])
                        for x in np.arange(x_min, x_max + step, step)
                        for y in np.arange(y_min, y_max + step, step)]
    
    points_on_plane2 = [(x, y, (-plane2[0] * x - plane2[1] * y - plane2[3]) / plane2[2])
                        for x in np.arange(x_min, x_max + step, step)
                        for y in np.arange(y_min, y_max + step, step)]

    # Create pseudoatoms for points on plane 1
    for idx, point in enumerate(points_on_plane1):
        atom_name = f"plane1_{idx}"
        pymol.cmd.pseudoatom(atom_name, pos=point, color="yellow")
        pymol.cmd.show("spheres", f"plane1_{idx}")

    # Create pseudoatoms for points on plane 2
    for idx, point in enumerate(points_on_plane2):
        atom_name = f"plane2_{idx}"
        pymol.cmd.pseudoatom(atom_name, pos=point, color="yellow")
        pymol.cmd.show("spheres", f"plane2_{idx}")

    # Mass center
    atom_name = "mass_center"
    pymol.cmd.pseudoatom(atom_name, pos=mass_center, color="magenta")
    pymol.cmd.show("spheres", atom_name)

    pymol.cmd.show("cartoon", "protein")
    pymol.cmd.save(f"output_{pdb_file.split('/')[-1].split('.')[0]}.pse")


def optimizing_width_membrane(gap_membrane, axis_init, amino_acid_sequence, plane_to_consider):
    """Optimiser la largeur de la membrane en explorant au-dessus et en dessous d'un plan."""
    best_axis = copy.deepcopy(axis_init)
    if plane_to_consider == 1:
        axis_init.plane1[3] += gap_membrane  # Modifier le terme constant d du plan
    else:
        axis_init.plane2[3] += gap_membrane

    while explore_axis(amino_acid_sequence, axis_init.plane1, axis_init.plane2, best_axis.best_hydrophobicity)[1]:
        if plane_to_consider == 1:
            axis_init.plane1[3] += gap_membrane
        else:
            axis_init.plane2[3] += gap_membrane

    return best_axis


# ============================= Main Program ==============================

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
      
    print(f"Command : python TM_detect.py {filename} -n {n} -w {width} -g {gap} -m {gap_membrane}")

    # Vérification du fichier d'entrée
    check_input_file(filename)

    # Parsing du fichier PDB
    protein = parse_pdb(filename, chain)
    
    # Génération des points sur une demi-sphère
    directions = find_points(2 * n, protein['mass_center'])
    
    print("Calculating the planes... ")
    for d in directions:
        point = copy.deepcopy(d)
        normal = find_director_vector(point, protein['mass_center'])
        plane1 = define_plane(point, normal)
        plane2 = define_plane(point, normal)
        plane2[3] += width  # Ajuster pour la largeur initiale de la membrane

        best_axis_tmp = {'plane1': plane1, 'plane2': plane2, 'best_hydrophobicity': -1000}
        
        # Exploration de l'axe
        while explore_axis(protein['amino_acid_sequence'], plane1, plane2, best_axis_tmp['best_hydrophobicity'])[1]:
            plane1[3] += gap
            plane2[3] += gap
        
        plane1 = define_plane(point, normal)
        plane2 = define_plane(point, normal)
        plane2[3] += width
        
        while explore_axis(protein['amino_acid_sequence'], plane1, plane2, best_axis_tmp['best_hydrophobicity'])[1]:
            plane1[3] -= gap
            plane2[3] -= gap
        
        protein['best_positions'].append(best_axis_tmp)
    
    # Trouver le meilleur axe
    best_axis = max(protein['best_positions'], key=lambda x: x['best_hydrophobicity'], default=None)
    if best_axis is None:
        print("Something went wrong... Try to modify the values of parameters.")
        sys.exit(1)

    print("Optimising membrane width...")
    best_axis_tmp = optimizing_width_membrane(gap_membrane, best_axis, protein['amino_acid_sequence'], 2)
    best_axis_tmp2 = optimizing_width_membrane(-gap_membrane, best_axis_tmp, protein['amino_acid_sequence'], 2)
    best_axis_tmp3 = optimizing_width_membrane(gap_membrane, best_axis_tmp2, protein['amino_acid_sequence'], 1)
    best_axis_tmp4 = optimizing_width_membrane(-gap_membrane, best_axis_tmp3, protein['amino_acid_sequence'], 1)

    print("The membrane width is ", abs(best_axis_tmp4['plane1'][3] - best_axis_tmp4['plane2'][3]))
    print("Best axis found overall is", best_axis_tmp4)
    
    # Visualiser dans PyMol
    show_in_pymol(best_axis_tmp4['plane1'], best_axis_tmp4['plane2'], filename, mass_center=protein['mass_center'])
