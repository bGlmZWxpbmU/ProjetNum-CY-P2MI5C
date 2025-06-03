# Etats_stationnaires.py : calcul et affichage des états stationnaires

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
import os

def calculer(dx, nx, v0, x_debut, x_fin, E):
    # Fonction pour calculer et imprimer les états stationnaires liés (E < 0, V = -V₀)
    def etats_stationnaires_lies(dx, nx, V, n_states=5):

        x = np.linspace(0, (nx - 1) * dx, nx)
        diag = np.full(nx, -2.0)
        offdiag = np.full(nx - 1, 1.0)
        T = (-1 / dx**2) * (np.diag(diag) + np.diag(offdiag, 1) + np.diag(offdiag, -1))
        H = T + np.diag(V)
        energies, states = eigh(H, subset_by_index=(0, n_states - 1))
        
        plt.figure(figsize=(10,6))
        for i in range(n_states):
            psi = states[:, i]
            psi = psi / np.sqrt(np.sum(psi**2) * dx)
            plt.plot(x, psi**2 + energies[i], label=f"État {i} (E = {energies[i]:.2f})")
        plt.plot(x, V, 'k--', label='Potentiel V(x)')
        plt.title("États stationnaires liés (E < 0, V = -V₀)")
        plt.xlabel("x (m)")
        plt.ylabel("Énergie (eV) / Densité de probabilité")
        plt.legend()
        plt.grid()
        
        # Création du dossier output s'il n'existe pas
        output_dir = "output"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Enregistrer l'image
        file_name = f'etats_stat_lies_e={E/v0:.0f}_v0={v0:.0f}.png'
        file_path = os.path.join(output_dir, file_name)
        plt.savefig(file_path)
        print(f"Graphique des états stationnaires liés exporté dans '{output_dir}'")

    # Fonction pour calculer et imprimer les états stationnaires libres (E > 0, V = 0)
    def etats_stationnaires_libres(dx, nx, n_states=5):

        # Paramètres:
        # dx: pas spatial
        # nx: nombre de points spatiaux
        # n_states: nombre d'états à calculer

        x = np.linspace(0, (nx - 1) * dx, nx)
        V_libre = np.zeros(nx)  # Potentiel nul
        
        diag = np.full(nx, -2.0)
        offdiag = np.full(nx - 1, 1.0)
        T = (-1 / dx**2) * (np.diag(diag) + np.diag(offdiag, 1) + np.diag(offdiag, -1))
        H = T + np.diag(V_libre)
        
        # Pour les états libres, on prend les premiers états avec E > 0
        energies, states = eigh(H)
        
        # Sélectionner les états avec énergies positives
        positive_indices = np.where(energies > 0)[0][:n_states]
        
        plt.figure(figsize=(10,6))
        for i, idx in enumerate(positive_indices):
            psi = states[:, idx]
            psi = psi / np.sqrt(np.sum(psi**2) * dx)
            energy = energies[idx]
            plt.plot(x, psi**2 + energy, label=f"État libre {i} (E = {energy:.2f})")
        
        plt.plot(x, V_libre, 'k--', label='Potentiel V(x) = 0')
        plt.title("États stationnaires libres (E > 0, V = 0)")
        plt.xlabel("x (m)")
        plt.ylabel("Énergie (eV) / Densité de probabilité")
        plt.legend()
        plt.grid()
        
        # Enregistrer l'image
        output_dir = "output"

        file_name = f'etats_stat_libres_e={E/v0:.0f}_v0={v0:.0f}.png'
        file_path = os.path.join(output_dir, file_name)
        plt.savefig(file_path)
        print(f"Graphique des états stationnaires libres exporté dans le dossier '{output_dir}'")

    # Nombre d'états à calculer
    n_states = int(input("\nNombre d'états à calculer [ex: 5]: ") or 5)
    
    # Grille spatiale
    x_grid = np.linspace(0, (nx - 1) * dx, nx)
    
    # Création du potentiel pour les états liés
    V_lie = np.zeros(nx)
    V_lie[(x_grid >= x_debut) & (x_grid <= x_fin)] = -v0
    
    print(f"\nCalcul de {n_states} états stationnaires...")
    
    # Calcul des états stationnaires liés (E < 0, V = -V₀)
    print("\n1. États liés dans le puit de potentiel:")
    etats_stationnaires_lies(dx, nx, V_lie, n_states)
    
    # Calcul des états stationnaires libres (E > 0, V = 0)
    print("\n2. États libres (particule libre):")
    etats_stationnaires_libres(dx, nx, n_states)