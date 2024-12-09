# README: BE Filtre de Kalman

## Description des fichiers

### Scripts principaux

- **`main.m`** : Script principal qui lance le chargement des données, l'initialisation, et l'exécution du filtre de Kalman.
- **`compute_transition_matrix.m`** : Génère les matrices de transition.
- **`compute_jacobian.m`** : Calcule la matrice Jacobienne utilisée dans la mise à jour du filtre de Kalman.
- **`initialize_filter.m`** : Initialise les paramètres du filtre de Kalman.

### Données

- **`carte.dat`** : Coordonnées des amers dans le repère global (X, Y, Z).
- **`mesure_accelero`** : Mesures d'accélération en fonction du temps.
- **`images/imageXXX`** : Fichiers représentant les images capturées à différents instants.

## Exécution

Exécuter le fichier `main.m` dans Matlab.
