import argparse
from Bio.PDB import PDBParser, DSSP
import numpy as np
import math
import pymol
import copy
import sys
import os


def verifier_fichier_entree(fichier_entree):
    # Vérifier si le fichier d'entrée est bien un fichier PDB valide.
    if not fichier_entree.lower().endswith(".pdb"):
        print("Extension incorrecte. Un fichier PDB est requis.")
        sys.exit(0)

    try:
        parser = PDBParser()
        structure = parser.get_structure("structure_identifiant", fichier_entree)
        for modele in structure:
            for chaine in modele:
                if chaine.id == 'A':
                    return True
    except Exception:
        print("Problème lors de la lecture du fichier PDB. Veuillez vérifier son contenu.")
        sys.exit(0)

# Constantes pour la liste des acides aminés hydrophobes
ACIDES_AMINES_HYDROPHOBES = ['PHE', 'GLY', 'ILE', 'LEU', 'MET', 'VAL', 'TRP', 'TYR']

def creer_acide_amine(code, id_aa, x, y, z):
    #Créer un acide aminé avec ses coordonnées, et vérifier s'il est hydrophobe.
    est_hydrophobe = code in ACIDES_AMINES_HYDROPHOBES
    point = np.array([x, y, z])
    return {
        'id': id_aa,
        'code': code,
        'est_hydrophobe': est_hydrophobe,
        'point': point
    }

def creer_proteine(nom, centre_masse, sequence_aa, sequence_complete, meilleures_positions):
    # Créer un objet protéine avec son nom, son centre de masse, sa séquence d'acides aminés, et ses meilleures positions.
    return {
        "nom": nom, 
        "centre_masse": None, 
        "sequence_aa": sequence_aa, 
        "sequence_complete": sequence_complete, 
        "meilleures_positions": meilleures_positions
    }

def calculer_accessibilite_solvant(structure, fichier_entree):
    # Calculer l'accessibilité au solvant des résidus dans une protéine via DSSP.
    print("Calcul de l'accessibilité au solvant...")
    modele = structure[0]
    dssp = DSSP(modele, fichier_entree, dssp='dssp')
    return dssp

def calculer_centre_masse(chaine):
    # Calculer le centre de masse d'une chaîne de la protéine.
    resultat = chaine.center_of_mass()
    x, y, z = resultat[0], resultat[1], resultat[2]
    return np.array([x, y, z])

def analyser_pdb(fichier_entree, chaine):
    # Analyser le fichier PDB pour extraire la séquence d'acides aminés et leur accessibilité.
    print("Analyse du fichier PDB en cours...")
    p = PDBParser()
    structure = p.get_structure("structure", fichier_entree)

    id_acide_amine = 1  
    id_acide_amine_complet = 1  
    nom = fichier_entree[-8:-4]  # Extraction du nom à partir du nom de fichier
    proteine = creer_proteine(nom=nom, centre_masse=None, sequence_aa=[], sequence_complete=[], meilleures_positions=[])

    # Calcul de l'accessibilité au solvant via DSSP
    resultat_dssp = calculer_accessibilite_solvant(structure, fichier_entree=fichier_entree)

    print("Réduction pour obtenir uniquement les résidus exposés...")
    residus_exposes = []
    for res_id, chaine, residu, asa, sec_struc, phi, psi, acc, torsion, kappa, alpha, zeta, chirality, beta in resultat_dssp:
        if asa > 0.3: # Seuil pour déterminer les résidus exposés
            residus_exposes.append(res_id)

    print(f"{len(residus_exposes)} résidus exposés trouvés.")

    modele = structure[0]
    chaine_selectionnee = modele["A"]

    # Calcul du centre de masse
    proteine['centre_masse'] = calculer_centre_masse(chaine_selectionnee)
    print(f"Le centre de masse est {proteine['centre_masse']}")

    # Extraction des coordonées de carbonne alpha :
    for residu in chaine_selectionnee:
        if residu.has_id("CA"):
            atome = residu['CA']
            x, y, z = atome.get_coord()
            nouvel_acide_amine = creer_acide_amine(code=residu.get_resname(), id_aa=residu.get_id()[1], x=x, y=y, z=z)
            proteine['sequence_complete'].append(nouvel_acide_amine)
            id_acide_amine_complet += 1
            if residu.get_id()[1] in residus_exposes:
                atome = residu['CA']
                x, y, z = atome.get_coord()
                nouvel_acide_amine = creer_acide_amine(code=residu.get_resname(), id_aa=residu.get_id()[1], x=x, y=y, z=z)
                proteine['sequence_aa'].append(nouvel_acide_amine)
                id_acide_amine += 1
    return proteine

def generer_points(n_points, centre_masse):
    # Générer des points répartis uniformément sur une demi-sphère autour du centre de masse.
    points = []
    phi = 0
    for k in range(1, n_points + 1):           #Algorithme de Saff et Kuijlaars
        h = -1 + (2 * (k - 1) / (n_points - 1))
        if h != 1:
            theta = math.acos(h)
            if k != 1 and k != n_points:
                phi = (phi + (3.6 / math.sqrt(n_points) * (1 / math.sqrt(1 - h * h)))) % (2 * math.pi)

            x = centre_masse[0] + math.sin(phi) * math.sin(theta)
            y = centre_masse[1] + math.cos(theta)
            z = centre_masse[2] + math.cos(phi) * math.sin(theta)
            points.append(np.array([x, y, z]))

    # Obtenir uniquement les points au-dessus de l'axe X (demi-sphère)
    points_au_dessus = []
    # Boucle pour filtrer les points 
    for point in points:
        if point[2] > centre_masse[2]:
            points_au_dessus.append(point)
    return points_au_dessus

def trouver_vecteur_directeur(point, coord_centre):
    # Trouver le vecteur directeur reliant un point et le centre de masse
    return np.array([coord_centre[0] - point[0],
                     coord_centre[1] - point[1],
                     coord_centre[2] - point[2]])

def definir_plan(point, vecteur_normal):
    # Définir un plan à partir d'un point et d'un vecteur normal.
    a = vecteur_normal[0]
    b = vecteur_normal[1]
    c = vecteur_normal[2]
    d = -(a * point[0] + b * point[1] + c * point[2])  # Calcul du plan avec l'équation ax + by + cz + d = 0
    return np.array([a, b, c, d])

def plan_parallele(plan, ecart):
    # Créer un plan parallèle décalé d'une certaine distance.
    nouveau_plan = list(plan)
    nouveau_plan[3] += ecart
    return tuple(nouveau_plan)

def glisser_plan(plan, fenetre_glissante):
    # Glisser un plan d'une certaine valeur le long de son axe normal.
    plan = list(plan)
    plan[3] += fenetre_glissante
    return tuple(plan)

def est_point_au_dessus_plan(point, plan):
    # Vérifier si un point est au-dessus d'un plan donné.
    return np.dot(plan[:3], point) + plan[3] > 0

def est_point_en_dessous_plan(point, plan):
    # Vérifier si un point est en dessous d'un plan donné.
    return np.dot(plan[:3], point) + plan[3] < 0

def creer_axe(plan1, plan2):
    # Créer un axe entre deux plans avec une hydrophobicité initiale faible.
    axe = {
        'plan1': plan1,
        'plan2': plan2,
        'meilleure_hydrophobicite': -1000,  # Hydrophobicité initialement faible
    }
    return axe

def trouver_meilleur_axe(proteine):
    # Trouver le meilleur axe de la protéine en fonction de l'hydrophobicité.
    meilleure_valeur_axe = 0
    meilleur_axe_trouve = None
    for axe in proteine["meilleures_positions"]:
        if axe["meilleure_hydrophobicite"] > meilleure_valeur_axe:
            meilleure_valeur_axe = axe["meilleure_hydrophobicite"]
            meilleur_axe_trouve = copy.deepcopy(axe)
    return meilleur_axe_trouve

def hydrophobicite_relative(sequence_aa, plan1, plan2, ref_axe):
    # Explorer un axe et trouver la meilleure hydrophobicité relative entre deux plans.
    entre_plans = []
    n_total_hydrophobes = 0
    n_total_hydrophiles = 0
    n_hydrophobes_entre_plan = 0
    nb_hydrophiles_hors_plan = 0

    # dentification des acides aminés entre les plans
    for aa in sequence_aa:
        if (est_point_en_dessous_plan(aa['point'], plan1) and est_point_au_dessus_plan(aa['point'], plan2)) or (est_point_en_dessous_plan(aa['point'], plan2) and est_point_au_dessus_plan(aa['point'], plan1)) :
            entre_plans.append(aa)

    # Comptage des acides aminés hydrophobes et hydrophiles 
        if aa['est_hydrophobe']:
            n_total_hydrophobes += 1
        else:
            n_total_hydrophiles += 1

    # Comptage des hydrophiles hors des plans 
    for aa in sequence_aa:
        if aa not in entre_plans and not aa['est_hydrophobe']:
            nb_hydrophiles_hors_plan += 1

    # Comptage des hydrophobes entre les plans 
    for aa in entre_plans:
        if aa['est_hydrophobe']:
            n_hydrophobes_entre_plan += 1

    if not entre_plans or len(entre_plans) >= len(sequence_aa):
        return False

    # Calcul de l'hydrophobicité
    hydrophobicite = (nb_hydrophiles_hors_plan / n_total_hydrophiles) + (n_hydrophobes_entre_plan / n_total_hydrophobes)
    
    # Mise à jour de l'axe si l'hydrophobicité est meilleure
    if hydrophobicite > ref_axe['meilleure_hydrophobicite']:
        ref_axe['meilleure_hydrophobicite'] = hydrophobicite
        ref_axe['plan1'] = copy.deepcopy(plan1)
        ref_axe['plan2'] = copy.deepcopy(plan2)
    return len(entre_plans) >0

def hydrophobicite_relative_verif(sequence_aa, plan1, plan2, ref_axe):
    # Vérifier la meilleure hydrophobicité relative entre deux plans.
    entre_plans = []
    n_total_hydrophobes = 0
    n_total_hydrophiles = 0
    n_hydrophobes_entre_plan = 0
    nb_hydrophiles_hors_plan = 0

    # dentification des acides aminés entre les plans
    for aa in sequence_aa:
        if (est_point_en_dessous_plan(aa['point'], plan1) and est_point_au_dessus_plan(aa['point'], plan2)) or (est_point_en_dessous_plan(aa['point'], plan2) and est_point_au_dessus_plan(aa['point'], plan1)) :
            entre_plans.append(aa)

    # Comptage des acides aminés hydrophobes et hydrophiles 
        if aa['est_hydrophobe']:
            n_total_hydrophobes += 1
        else:
            n_total_hydrophiles += 1

    # Comptage des hydrophiles hors des plans 
    for aa in sequence_aa:
        if aa not in entre_plans and not aa['est_hydrophobe']:
            nb_hydrophiles_hors_plan += 1

    # Comptage des hydrophobes entre les plans 
    for aa in entre_plans:
        if aa['est_hydrophobe']:
            n_hydrophobes_entre_plan += 1

    if not entre_plans or len(entre_plans) >= len(sequence_aa):
        return False

    # Calcul de l'hydrophobicité
    hydrophobicite = (nb_hydrophiles_hors_plan / n_total_hydrophiles) + (n_hydrophobes_entre_plan / n_total_hydrophobes)
    
    # Mise à jour de l'axe si l'hydrophobicité est meilleure
    if hydrophobicite > ref_axe['meilleure_hydrophobicite']:
        ref_axe['meilleure_hydrophobicite'] = hydrophobicite
        ref_axe['plan1'] = copy.deepcopy(plan1)
        ref_axe['plan2'] = copy.deepcopy(plan2)
        return True
    return False

def optimiser_largeur_membrane(ecart_membrane, axe_initial, sequence_aa, plan_a_considerer):
    # Optimiser la largeur de la membrane en glissant les plans pour maximiser l'hydrophobicité.
    meilleur_axe = copy.deepcopy(axe_initial)
    if plan_a_considerer == 1:
        axe_initial['plan1'] = glisser_plan(axe_initial['plan1'], ecart_membrane)
    else:
        axe_initial['plan2'] = glisser_plan(axe_initial['plan2'], ecart_membrane)

    while hydrophobicite_relative_verif(sequence_aa, axe_initial['plan1'], axe_initial['plan2'], ref_axe=axe_initial) :
        if plan_a_considerer == 1:
            axe_initial['plan1'] = glisser_plan(axe_initial['plan1'], ecart_membrane)
        else:
            axe_initial['plan2'] = glisser_plan(axe_initial['plan2'], ecart_membrane)
    return meilleur_axe

def afficher_dans_pymol(plan1, plan2, fichier_pdb, centre_masse, points_sphere):
    # Afficher la protéine, les plans et les points sur la demi-sphère dans PyMol. 
    pymol.finish_launching(['pymol', '-q'])
    pymol.cmd.load(fichier_pdb, "protein")
    pymol.cmd.remove("solvent")

    x_min = float("inf")
    x_max = float("-inf")
    y_min = float("inf")
    y_max = float("-inf")

    # Trouver les limites x et y de la protéine pour définir la taille des plans (pour s'assurer que les plans couvrent toute la région d'intérêt)
    # Déterminer les dimensions de la molécule sur l'axe des X
    for atom in pymol.cmd.get_model("protein").atom:
        x = atom.coord[0]
        if x < x_min:
            x_min = x
        if x > x_max:
            x_max = x
    # Déterminer les dimensions de la molécule sur l'axe des Y
    for atom in pymol.cmd.get_model("protein").atom:
        y = atom.coord[1]
        if y < y_min:
            y_min = y
        if y > y_max:
            y_max = y

    ecart = 1 # Espacement entre les points du plan

    points_on_plane1 = []
    for x in np.arange(x_min, x_max + ecart, ecart):
        for y in np.arange(y_min, y_max + ecart, ecart):
            z1 = (-plan1[0] * x - plan1[1] * y - plan1[3]) / plan1[2]  # Calculer les coordonnées de Z en utilisant l'équation du plan 2
            points_on_plane1.append((x, y, z1))

    # Generate points on plane 2
    points_on_plane2 = []
    for x in np.arange(x_min, x_max + ecart, ecart):
        for y in np.arange(y_min, y_max + ecart, ecart):
            z2 = (-plan2[0] * x - plan2[1] * y - plan2[3]) / plan2[2]  # Calculer les coordonnées de Z en utilisant l'équation du plan 2
            points_on_plane2.append((x, y, z2))

    # Créer les pseudo-atoms pour les points sur le plane 1
    for idx, point in enumerate(points_on_plane1):
        x, y, z = point
        atom_name = f"plane1_{idx}"
        pymol.cmd.pseudoatom(atom_name, pos=point, color="blue")
        pymol.cmd.show("spheres", f"plane1_{idx}")

    # Créer les pseudo-atoms pour les points sur le plane 2
    for idx, point in enumerate(points_on_plane2):
        x, y, z = point
        atom_name = f"plane2_{idx}"
        pymol.cmd.pseudoatom(atom_name, pos=point, color="red")
        pymol.cmd.show("licorice", f"plane2_{idx}")

    # Créer des Pseudo-Atoms pour les Points de la Sphère 
    for idx, point in enumerate(points_sphere):
        pymol.cmd.pseudoatom(f"point_{idx}", pos=list(point), color= "white")
        pymol.cmd.show("spheres", f"point_{idx}")

    # Afficher le centre de la sphère
    pymol.cmd.pseudoatom("centre_masse", pos=list(centre_masse), color="red")
    pymol.cmd.show("dots", "centre_masse")

    # Afficher le centre de masse
    atom_name = "mass_center"
    pymol.cmd.pseudoatom(atom_name, pos=[centre_masse[0], centre_masse[1], centre_masse[2]], color="yellow")
    pymol.cmd.set("sphere_transparency", 0.5)
    pymol.cmd.show("spheres", "mass_center")

    # Afficher la structure de la proteine
    pymol.cmd.show("cartoon", "protein")

    # Utilisation de os.path.basename pour extraire le nom du fichier sans son chemin
    nom_fichier_sans_chemin = os.path.basename(fichier_pdb)  
    nom_fichier_sans_extension = nom_fichier_sans_chemin[:-4]  

    # Enregistrer dans un fichier pdb
    pymol.cmd.save(f"output_{nom_fichier_sans_extension}.pse")

# Programme principal 
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='code_TM_detect.py')
    parser.add_argument('fichier', help="Un fichier PDB")
    
    # Arguments optionnels :
    parser.add_argument('-n', help="Nombre de points à placer sur la sphère. (par défaut 20)", type=int, dest="n", default=20) 
    parser.add_argument('-w', help="Largeur initiale de la membrane. (par défaut 15 A)", type=float, dest="largeur", default=15) 
    parser.add_argument('-g', help="Écart de glissement de la membrane le long d'un axe. (par défaut 1 A)", type=float, dest="ecart", default=1) 
    parser.add_argument('-m', help="Écart d'optimisation de la largeur de la membrane. (par défaut 1 A)", type=float, dest="ecart_membrane", default=1) 
    
    # Parsing :
    args = parser.parse_args()
    fichier = args.fichier
    n = args.n
    chaine = "A"
    largeur = args.largeur
    ecart = args.ecart
    ecart_membrane = args.ecart_membrane
      
    print(f"Commande : python code_TM_detect.py {fichier} -n {n} -w {largeur} -g {ecart} -m {ecart_membrane}")

    # Vérification du fichier d'entrée
    verifier_fichier_entree(fichier)

    # Analyse du fichier PDB
    proteine = analyser_pdb(fichier, chaine)
    
    # Génération des points sur une demi-sphère
    directions = generer_points(n, proteine['centre_masse'])
    
    # Pour chaque direction générée 
    for d in directions:
        point = copy.deepcopy(d)
        normal = trouver_vecteur_directeur(point=point, coord_centre=proteine['centre_masse'])
        plan1 = definir_plan(point=point, vecteur_normal=normal)
        plan2 = plan_parallele(plan1, ecart=largeur)

        # Créer un axe pour exploration avec plan 1 et plan 2
        axe = creer_axe(plan1, plan2)
        meilleur_axe_tmp = copy.deepcopy(axe)
        
        # Exploration de l'axe pour optimiser l'hydrophobicité
        while hydrophobicite_relative(proteine['sequence_aa'], axe["plan1"], axe["plan2"], ref_axe=meilleur_axe_tmp) is True :
            # Déplace les plans le long de l'axe
            axe = creer_axe(glisser_plan(axe["plan1"], ecart), glisser_plan(axe["plan2"], ecart))

        # Recréer plan1 et plan2 
        plan1 = definir_plan(point=point, vecteur_normal=normal)
        plan2 = plan_parallele(plan1, ecart=largeur)

        # Exploration dans la direction opposée pour optimiser l'hydrophobicité
        while hydrophobicite_relative(proteine['sequence_aa'], axe["plan1"], axe["plan2"], ref_axe=meilleur_axe_tmp) is True :
            axe = creer_axe(glisser_plan(axe["plan1"], -ecart), glisser_plan(axe["plan2"], -ecart))

        # Sauvegarder la meilleure position d'axe
        proteine['meilleures_positions'].append(meilleur_axe_tmp)
    
    # Trouver le meilleur axe parmi toutes les positions testées
    meilleur_axe = trouver_meilleur_axe(proteine)
    if meilleur_axe is None:
        print("Un problème est survenu...")
        sys.exit(1)
    
    # Affiche le meilleur axe trouvé
    print(f"Meilleur axe trouvé : {meilleur_axe}")

    # Optimisation de la largeur de la membrane dans différentes directions (plan inf/sup, au-dessus/au dessous)
    meilleur_axe_tmp = optimiser_largeur_membrane(ecart_membrane=ecart_membrane, axe_initial=meilleur_axe, sequence_aa=proteine['sequence_aa'], plan_a_considerer=2)
    meilleur_axe_tmp2 = optimiser_largeur_membrane(ecart_membrane=-ecart_membrane, axe_initial=meilleur_axe_tmp, sequence_aa=proteine['sequence_aa'], plan_a_considerer=2)
    meilleur_axe_tmp3 = optimiser_largeur_membrane(ecart_membrane=ecart_membrane, axe_initial=meilleur_axe_tmp2, sequence_aa=proteine['sequence_aa'], plan_a_considerer=1)
    meilleur_axe_tmp4 = optimiser_largeur_membrane(ecart_membrane=-ecart_membrane, axe_initial=meilleur_axe_tmp3, sequence_aa=proteine['sequence_aa'], plan_a_considerer=1)

    # Afficher la largeur de la membrane optimisée et le meilleur axe trouvé globalement
    print("La largeur de la membrane est ", abs(meilleur_axe_tmp4['plan1'][3] - meilleur_axe_tmp4['plan2'][3]))
    print("Meilleur axe trouvé globalement : ", meilleur_axe_tmp4)
    
    # Visualisation du résultat dans PyMol
    afficher_dans_pymol(meilleur_axe_tmp4['plan1'], meilleur_axe_tmp4['plan2'], fichier, centre_masse=proteine['centre_masse'], points_sphere=directions)
