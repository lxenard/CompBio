Tusnády GE, Dosztányi Z, Simon I. Transmembrane proteins in the Protein Data Bank: identification and classification. Bioinformatics. 2004 Nov 22;20(17):2964-72. Epub 2004 Jun 4. PubMed PMID: 15180935
Lien : https://academic.oup.com/bioinformatics/article/20/17/2964/185954

PTM : protéines transmembranaires
PG : protéines globulaires
AAS : accessible au solvant

1 Å = 10-10m

# Notes article

## Abstract

PDB => peu de protéines membranaires représentées, peu annotées, pas de localisation de la bicouche lipidique
pas de méthode permettant de détecter le plan membranaire à partir des coordonnées atomiques des protéines membranaires

uniquement utilisation de l'information structurale pour différencier protéines transmembranaires et protéines globulaires
identification de la localisation la plus probable de la bicouche lipidique

**Protein Data Bank of Transmembrane Proteins**
PDB_TM database : http://www.enzim.hu/PDB_TM LIEN MORT
=> http://pdbtm.enzim.hu/
Base de donnée des protéines transmembranaires (et fragments)

## Introduction

PTM cible de nombreux médicaments (humain + vétérinaire)
difficile de les cristalliser en environnement aqueux d'où sous-représentation dans la PBD
généralement plus grands que PG donc résolution de la structure compliquée en RMN

parties hydrophobes masquées par des détergents amphiphiles
=> traitement des complexes PTM-détergent comme si c'était des protéines solubles



identification des segments transmembranaires

## Methods

calcul  de la surface accessible au solvant pour chacun des atomes pouvant potentiellement interagir avec la bicouche

Pour chaque tranche, somme de l'hydrophobicité des résidus hydrophiles AAS d'une part et des résidus hydrophobes AAS d'autre part

Fonction objectif : fo()
facteur d'hydrophobicité = surface hydrophobe divisée par surface totale
facteur de structure = produit de 3 facteurs : straightness, turn et end-chain factor

facteur straightness = fréquence relative des résidus droits dans dans une tranche
facteur turn = 1 - fréquence relative des résidus 'turn' dans une tranche
facteur end-chain = 1 - fréquence relative des résidus de fin de chaîne dans une tranche

fo() = mean(facteur d'hydrophobicité * facteur de structure pour chaque slice) 

## Results



## Discussion

# Projet

## Requirements

* python 3.8 ou 3.7

* numpy 1.20.3
* biopython 1.79 conda install -c conda-forge biopython
* dssp 3.0 conda install -c salilab dssp LINUX
* dssp 2.0 conda install -c speleo3 dssp WINDOWS

Visualisation
* mayavi 4.7.1  conda install -c anaconda mayavi 3.7
* pymol 2.5.2  conda install -c schrodinger pymol-bundle 3.8

## Résumé étapes

1) Choisir des protéines membranaires et globulaires de référence (simple à analyser, voir banque OPM)
2) Calcul de la surface accessible au solvant de chaque acide aminé avec NACCESS ou DSSP
3) Lecture du PDB et extraction des coordonnées des C alpha
4) Calcul du centre de masse
5) Détermination des droites passant par le centre de masse et quadrillant avec suffisamment de précision toutes les directions
6) Déplacement d'une tranche de 1 Å normale à une droite et calcul de l'hydrophobicité relative des résidus dans les zones exposées dans la tranche
7) Calcul de la position de la membrane par la moyenne de l'hydrophobicité relative et en comparant ces valeurs selon les différentes droites

Ne garder que les Cα
Classifier les résidus en 2 catégories : hydrophobes et non hydrophobes
On cherche la position de la membrane qui maximise le nombre de résidus hydrophobes présents à l'intérieur de la bicouche
Calcul du centre de gravité, et y centrer le repère
Dessin d'une sphère centrée sur le centre de gravité et échantillonnage de la surface de la sphère. Pour commencer, 50 points dans un hémisphère
Les points échantillonnés donnent les vecteurs normaux aux vecteurs directeurs des tranches initiales
Ces tranches sont toujours centrées et on commence par une tranche d'épaisseur 14Å donc 7Å de part et d'autre
Puis pour chaque vecteur, définition de tranches par fenêtre glissante (déplacement de 2Å pour commencer)
Maximisation de la fonction objectif pour chacune des tranches. Une fois la meilleur tranche identifiée, on l'agrandit progressivement par tranche de 1Å jusqu'à maximiser la fonction objectif.
Possible d'identifier plusieurs membranes, dans ce cas elle doivent toutes être traitées.

Fonction objectif pour commencer : ~~rapport entre les résidus externes et ceux internes à la membrane~~
**nb de résidus hydrophobes dans la membrane / nb de résidus hydrophobes total**

Pour visualiser, utiliser PyMol pour plotter la structure ac un code couleur pour les résidus hydrophobes et non hydrophobes, puis tracer la membrane telle qu'elle a été déterminée


## 1) Choix des protéines

OPM : Orientation of Proteins in Membranes database
https://opm.phar.umich.edu/
Commencer par travailler sur [2n90](https://opm.phar.umich.edu/proteins/3252), ne prendre qu'une seule chaîne (1 seul segment transmembranaire)
Puis prendre des protéines avec plus de segments comme par exemple [6g79](https://opm.phar.umich.edu/proteins/3904), ne pas hésiter à couper la prot lorsque trop longue et en dehors de la membrane.

Lecture du fichier PDB

## 2) Filtre DSSP

https://biopython.org/docs/1.75/api/Bio.PDB.DSSP.html
https://anaconda.org/salilab/dssp

https://swift.cmbi.umcn.nl/gv/dssp/





## Implémentation

### Constantes 

N_DIRECTIONS nb de points pour échantillonner la demi-sphère
SLIDING_WINDOW_SIZE incrément pour la fenêtre glissante
SLICE_INIT_SIZE taille initiale de tranche
SLICE_STEP_SIZE incrément pour le changement de taille de la tranche

### Classes

Point
	x
	y
	z

Vector
	Point start (0, 0, 0) attribut de classe
	Point end

Sphere
	Point center (0, 0, 0) attribut de classe
	int radius (peu importe)
	.
	sample(int nb)
		Échantillonne la surface de la demi-sphère z-positive avec nb nombre de points
	*Équation d'une sphère de centre (0, 0, 0) : x² + y² + z² = r²*
	*Chercher un module de géométrie pour ne pas réinventer la roue*

https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
https://gist.github.com/dinob0t/9597525



Slice
	Vector normal
	Point center
	int thickness
	float score
	list(Residue) residues
	.
	find_residues(list(Residue))
		Calcul les 2 plans délimitant la tranche
		Pour chaque résidu
			Détermine si résidu appartient à la tranche
			Si oui
				Ajouter le résidu à la liste residues
	compute_score()
		Rapport entre les résidus externes et ceux internes à la membrane




Un résidu est assimilé à un Cα

Residue -> Res_hydrophobic
		-> Res_polar
	Point position
	*voir en fonction du besoin si hydrophobic = sous-classe ou bool*

Protein
	list(Residue) residues
	list(Slice) slices
	Sphere sphere
	.
	get_nb_residues()
		Renvoie la taille de la liste residues
	compute_slices(int N_DIRECTIONS)
		Echantillonne la surface de la demi-sphère en N_DIRECTIONS points
		Pour chaque point
			Détermine le vecteur passant par l'origine du repère
			Instancie la Slice
			Ajoute la Slice à slices
			Tant qu'on ne sort pas de la protéine
				Translater le centre de la slice de SLIDING_WINDOW_SIZE dans le sens du vecteur
				Instancier la nouvelle tranche
				Ajouter la Slice à slices
			Tant qu'on ne sort pas de la protéine
				Translater le centre de la slice de SLIDING_WINDOW_SIZE dans le sens opposé du vecteur
				Instancier la nouvelle tranche
				Ajouter la Slice à slices



Critère d'arrêt de sortie de la protéine : lorsque plus aucun résidu hydrophobe détecté

Récupérer les Cα
Calculer le centre de gravité de ces Cα
Centrer le repère sur le centre de gravité
Générer une sphère (ou demi-sphère) centrée de rayon quelconque (mais suffisant pour que l'échantillonnage de la surface se fasse facilement donc mieux vaut prendre trop grand que trop petit)
Échantillonner la surface de la demi-sphère z-positive de manière à obtenir n points
Création de n vecteurs à partir de ces n points
Pour chaque vecteur
	Déterminer la tranche centrée normale
	Déterminer les tranches translatées selon le vecteur normal
		Pour chaque tranche
			Identifier les résidus hydrophobes de la tranche
			Calculer le score de la tranche
Comparer les scores (voir si prendre le meilleur score ou si prendre les scores supérieur à un certain seuil afin de pouvoir détecter des membranes multiples)
Pour chaque tranche à score élévé
	Tant que le score augmente
		Agrandir la tranche dans le sens du vecteur
	Tant que le score augmente
		Agrandir la tranche dans le sens opposé du vecteur
	Retourner la nouvelle tranche obtenue
	Dessiner la tranche (PyMol)

	intégrine
	cadhérine
	CMH classe I
	porine

On ne travaille que sur les résidues considérés comme exposés. La partie *"In order to*
improve the membrane detection algorithm, the ‘water accessible surface’ is considered only for those atoms which could potentially interact with the lipid bilayer. These membrane exposed atoms are selected by the following approximate filtering procedure: the protein is cut into 1 Å wide slices along a predefined axis, and around each slice of atoms, test points are placed on a rectangle which embeds all the atoms within that slice. Those atoms lying closest to any of the test points are defined to be on the outside (i.e. possible membrane exposed) of the surface. For all other atoms the ‘water accessible surface area’ is set to zero." décrite dans l'article n'est donc pas à implementer mais est remplacée par la sélection des résidus dont le ASA est supérieur à un certain seuil (0.1~0.3).

Pour vérifier la position du centre de gravité :

```
pseudoatom tmp, pos=[10.0, 17.0, -3.0]
tmp expand 6
```

trouver plan normal à un vecteur

faire ajustement après avoir trouvé la meilleure direction

ne pas prendre de protéine canal

organisation : conda, POO, github
pseudo code ou schéma en annexe

https://www.quora.com/Given-a-point-and-a-plane-how-would-you-determine-which-side-of-the-plane-the-point-lies