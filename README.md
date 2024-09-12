# Projet court - ASSIGNATION ET DETECTION DES PARTIES TRANSMEMBRANAIRES D'UNE PROTEINE

### Description

Ce projet vise à identifier et détecter les régions transmembranaires d'une protéine à partir de son fichier PDB. L'outil utilise des calculs basés sur la géométrie et l'hydrophobicité des résidus pour localiser les membranes transmembranaires.

Le programme a été réalisé dans le cadre d'un projet universitaire, en se basant sur l'article suivant :

[Référence](https://pubmed.ncbi.nlm.nih.gov/15180935/) : 
    Tusnády GE, Dosztányi Z, Simon I. Transmembrane proteins in the Protein Data Bank: identification and classification. Bioinformatics. 2004 Nov 22;20(17):2964-72. Epub 2004 Jun 4. PubMed PMID: 15180935.
### Installation 
#### Méthode 1 : Cloner le dépôt avec Git
Assurez-vous d'être dans le répertoire où vous souhaitez cloner le dépôt.
```bash
git clone https://github.com/Meriemyssf/projet_court
```
#### Méthode 2 : Télécharger et décompresser le dépôt
Téléchargez le dépôt au format ZIP. Dans un dossier approprié, décompressez-le en sélectionnant **Extraire** ou en utilisant la commande suivante :
```bash
unzip chemin/vers/le/fichier.zip -d chemin/vers/le/dossier
```
Changez ensuite votre répertoire actuel vers la racine du projet :
```bash
cd ./Projet_court
```
Vous devriez voir les fichiers suivants avec la commande ls :
```
code_TM_detect.py  input_code   README.md
environment.yml    output_code  YOUSSEF_Meriem_rapport_projet_court.pdf
```
### Création de l'environnement Conda
Exécutez la commande suivante pour créer un nouvel environnement Conda :
```bash
conda env create -f environment.yml
```
Une fois l'environnement créé, activez-le avec :
```bash
conda activate transmembrane_detection_PC
```
### Exécution du programme 
Après avoir activé l'environnement virtuel, exécutez le programme avec la commande suivante :
```
python code_TM_detect.py input_code/fichier_pdb [options]
```
Ce code sera exécuté avec les paramètres optionnels par défaut suivants :

    -n : Nombre de points à placer sur la sphère (par défaut 20).
    -w : Largeur initiale de la membrane (par défaut 15 Å).
    -g : Écart de déplacement de la membrane le long d'un axe (par défaut 1 Å).
    -m : Écart d'optimisation de la largeur de la membrane (par défaut 1 Å).

Vous pouvez modifier ces paramètres selon vos besoins.

### Résultats
Le programme ouvre PyMol avec la protéine, les plans prédits, la sphère et son centre, ainsi que le centre de masse de la protéine, tous représentés comme objets hétéroatoms. 
La session PyMol contenant ces éléments est sauvegardée dans un fichier avec l'extension “.pse”.

### Dépannage

En cas de problème :

    - Vérifiez que toutes les dépendances sont installées.
    - Assurez-vous que le fichier PDB est valide et correctement formaté.
    - Consultez les messages d'erreur pour plus d'informations.
### Crédits
  **Auteur :** YOUSSEF Meriem
  
  **Niveau d'Études :** Master 2 Bioinformatique