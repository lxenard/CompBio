`membrane_plane.py` est un programme permettant (parfois) de déterminer la position de la membrane sur une protéine transmembranaire uniquement à partir de l'information de structure de celle-ci telle que fournie par un fichier PDB.

Il s'inspire fortement des algorithmes décrits dans [Transmembrane proteins in the Protein Data
Bank: identification and classification](https://doi.org/10.1093/bioinformatics/bth340) (Tusnády, Dosztányi et Simon 2004) et dans [Membrane positioning for high-and low-resolution protein structures through a binary classification approach](https://doi.org/10.1093/protein/gzv063) (Postic, Ghouzm, Guiraud et Gelly 2016).

# Installation

## Linux

Les dépendances requises sont disponibles à la racine, dans le fichier  `membrane_lin.yml`. Pour installer automatiquement ces dépendances dans un environnement conda, saisir la commande :

```conda env create -f membrane_lin.yml```

Ou pour aller plus vite et si mamba a été installé :

```mamba env create -f membrane_lin.yml```

## Windows

L'installation est identique à celle pour Linux, il suffit d'utiliser le fichier `membrane_win.yml` à la place du fichier  `membrane_lin.yml` .

# Utilisation

## Informations générales

`membrane_plane` supporte le choix du modèle et de la chaîne sur lesquels effectuer la détection de la membrane (paramètres `-m` et `-c`). Il est également possible de ne travailler que sur une certaine sélection de résidus (paramètres `-fr` et `-lr`).

En revanche, ce programme ne supporte pas les chaînes comportant des hétéro-atomes (HETATM) ou des atomes ayant une position alternative (ALTLOC).

Une description complète des options de `membrane_plane` est disponible via la commande :

```python ./scripts/membrane_plane.py -h```

Le seul argument obligatoire est le chemin/nom du fichier PDB à traiter. Le programme produit une sortie console renseignant sur l'épaisseur de la membrane et le nombre de résidus que celle-ci englobe. Un fichier PDB nommé `{mon_pdb}_membrane_{options}.pdb` est enregistré dans le dossier contenant le fichier PDB d'origine. Ce nouveau fichier PDB contient, en plus des résidus de la chaîne initiale, deux résidus représentant chacun une paroi de la membrane. Ouvrir le fichier avec PyMol (par exemple) permet de visualiser la position de la membrane par rapport à la protéine transmembranaire.

Le mode de debug (paramètre `-d`) permet également une visualisation 3D de la localisation de la membrane mais celle-ci est plus basique.

## Exemple d'utilisation

Les fichiers PDB utilisés dans les exemples suivants sont disponibles dans le dossier `data`.

### Minimal

En ne renseignant que le chemin du fichier PDB à traiter, tous les paramètres prendront leur valeur par défaut :

```python ./scripts/membrane_plane.py ./data/2n90.pdb```

Le fichier de sortie sera alors : `2n90_membrane_m0_cA_fr1_lr39.pdb`

### Sélection du modèle et le chaîne à traiter

On peut sélectionner le modèle sur lequel travailler grâce à l'option `-m` (ou `--model`) ainsi que la chaîne avec `-c` (ou `--chain`). Ici on sélectionne le modèle 2 et la chaine B :

```python ./scripts/membrane_plane.py ./data/2n90.pdb -m 2 -c 'B'```

Le fichier de sortie sera alors : `2n90_membrane_m2_cB_fr101_lr139.pdb`
Si le modèle ou la chaîne n'existe pas au sein de la protéine, le premier modèle et la première chaîne de la structure seront traités. Par défaut, le modèle est fixé à 0 et la chaîne A.

### Sous-sélection de résidus

On peut également ne travailler que sur une sous-sélection de résidus grâce aux options `-f` (ou `-first_residue`) et `-l` (ou `-last_residue`). La sortie PDB inclura cependant tous les résidus d'origine afin d'éviter de casser les structures. Ici on sélectionne les résidus 10 à 30 :

```python ./scripts/membrane_plane.py ./data/2n90.pdb -f 10 l 30```

Le fichier de sortie sera alors : `2n90_membrane_m0_cA_fr10_lr30.pdb`
Si ces résidus n'existe pas au sein du modèle et de la chaîne, tous les résidus de la chaîne et du modèle seront traités, comme si les options `-r` et `-l` n'avaient pas été renseignées.

### Modifications des paramètres de calcul

Par défaut, les paramètres de calcul sont fixés de manière à limiter le temps d'exécution. Se référer à la documentation du programme pour le détail de ces options (`-h` ou `--help`).

Afin d'éviter les noms de fichier à rallonge, les paramètres de calcul ne sont pas renseignés dans le nom des fichiers PDB de sortie.

### Finalement

Il est bien entendu possible de combiner toutes ces options.

# Licence

Ce programme est sous GNU General Public License v3.0. Pour plus d'informations, consultez le fichier [LICENSE](LICENSE.txt).

