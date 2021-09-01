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

## Étapes

1) Choisir des protéines membranaires et globulaires de référence (simple à analyser, voir banque OPM)
2) Calcul de la surface accessible au solvant de chaque acide aminé avec NACCESS ou DSSP
3) Lecture du PDB et extraction des coordonnées des C alpha
4) Calcul du centre de masse
5) Détermination des droites passant par le centre de masse et quadrillant avec suffisamment de précision toutes les directions
6) Déplacement d'une tranche de 1 Å normale à une droite et calcul de l'hydrophobicité relative des résidus dans les zones exposées dans la tranche
7) Calcul de la position de la membrane par la moyenne de l'hydrophobicité relative et en comparant ces valeurs selon les différentes droites

## 1) Choix des protéines

OPM : Orientation of Proteins in Membranes database
https://opm.phar.umich.edu/


Ne garder que les Cα
Classifier les résidus en 2 catégories : hydrophobes et non hydrophobes
On cherche la position de la membrane qui maximise le nombre de résidus hydrophobes présents à l'intérieur de la bicouche
Calcul du centre de gravité, et y centrer le repère
Dessin d'une sphère centrée sur le centre de gravité et échantillonnage de la surface de la sphère
Les points échantillonnés donnent les vecteurs directeurs des tranches
Les tranches sont toujours centrées et on commence par une tranche d'épaisseur 14Å donc 7Å de part et d'autre
Maximisation de la fonction objectif pour chacune des tranches. Une fois la meilleur tranche identifiée, on l'agrandit progressivement par tranche de 1Å jusqu'à maximiser la fonction objectif.
Possible d'identifier plusieurs membranes, dans ce cas elle doivent toutes être traitées.


Fonction objectif pour commencer : rapport entre les résidus externes et ceux internes à la membrane 

Pour visualiser, utiliser PyMol pour plotter la structure ac un code couleur pour les résidus hydrophobes et non hydrophobes, puis tracer la membrane telle qu'elle a été déterminée
