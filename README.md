# ProjetNum-CY-P2MI5C

**Projet de physique moderne de pré-ing 2**

**MI-5 Groupe C**.

## Effet Ramsauer–Townsend

Lorsque deux particules se rencontrent, elles peuvent interagir et l’une peut être déviée par l’autre : il s’agit de la diffusion. L’effet Ramsauer–Townsend est un phénomène quantique où la probabilité de diffusion d’un électron par un atome de gaz noble devient nulle pour certaines énergies. L’objectif de notre projet est de comprendre cet effet à travers un modèle basé sur un puits de potentiel fini et unidimensionnel.

## Contributeurs

- Romain MICHAUT-JOYEUX
- Nathan CHOUPIN
- Ziyad HADDADI

## Langage utilisé

![PYTHON](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)

## Prérequis

Ce programme nécessite **Python 3**. Pour l'installer, tapez dans votre terminal `sudo apt update` puis `sudo apt install python3 python3-pip`.

Il nécessite également les bibliothèques suivantes :

- **ffmpeg**

- **numpy**

- **matplotlib**

Pour installer **ffmpeg**, entrez `sudo apt install ffmpeg`.

Pour installer **numpy** et **matplotlib**, entrez `pip3 install numpy matplotlib`.

## Installation

Afin de récupérer le repos, entrez `git clone https://github.com/bGlmZWxpbmU/ProjetNum-CY-P2MI5C`.

## Utilisation

Pour éxecuter le programme, lancer la commande suivante dans le répertoire du projet : `./main.py`.

Ensuite, vous aurez la possibilité de rentrer vos propres paramètres dans le terminal. Si vous souhaitez utiliser des valeurs par défaut, il vous suffit d'appuyer sur la touche entrée de votre clavier lors de la saisie.

## Revue des programmes

`Etats_stationnaires.py` calcule les états stationnaires d'une particule dans un potentiel donné.

`Parametres.py` permet de récupérer les entrées utilisateurs (pas de temps, pas spatial, etc.).

`Sch1d_solution.py` est un algorithme de résolution d'équation de Schrödinger en une dimension.

`main.py` est la boucle principale du programme.
