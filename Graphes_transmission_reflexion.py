# Graphes_transmission_reflexion.py : tracé des courbes T(E) et R(e)

import numpy as np
import matplotlib.pyplot as plt
import os

# Création du dossier output
output_dir = "output"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def ramsauer_townsend():
    # Paramètres
    v0 = 4  # Pontentiel, réduire pour voir l'effet
    m = 9.11e-31  # Masse de l'électron
    hbar = 1.055e-34
    a = 3e-9  # Largeur du puits, augmenter pour amplifier l'effet
    eV = 1.602e-19
    
    # Plage d'énergies
    E_eV = np.linspace(0.01, 10, 1000)
    
    # Calcul des coefficients de transmission et réflexion
    T = np.zeros_like(E_eV)
    R = np.zeros_like(E_eV)
    
    for i, E_ev in enumerate(E_eV):
        E = E_ev * eV  # Conversion en Joules
        
        # Terme dans le sinus
        k_inside = np.sqrt(2 * m * (E + v0 * eV)) / hbar
        sin_term = np.sin(k_inside * a)**2
        
        # Coefficient de transmission, on utilise la formule trouvée dans l'étude analytique
        denominator = 1 + (v0**2) / (4 * E_ev * (E_ev + v0)) * sin_term
        T[i] = 1 / denominator
        
        # Coefficient de réflexion
        R[i] = 1 - T[i]
    
    # Création des graphiques
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # T(E)
    ax1.plot(E_eV, T, 'b-', linewidth=2, label='T(E)')
    ax1.set_ylabel('Coefficient de transmission T')
    ax1.set_title('Effet Ramsauer-Townsend - Transmission')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 1.1)
    ax1.legend()
    
    # R(E)
    ax2.plot(E_eV, R, 'r-', linewidth=2, label='R(E)')
    ax2.set_xlabel('Énergie E (eV)')
    ax2.set_ylabel('Coefficient de réflexion R')
    ax2.set_title('Effet Ramsauer-Townsend - Réflexion')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, 1.1)
    ax2.legend()
    
    plt.tight_layout()
    filepath_TR = os.path.join(output_dir, f"T(E)_R(E)_v0={v0}.png")
    plt.savefig(filepath_TR, dpi=300, bbox_inches='tight')
    
    # Graphique combiné
    plt.figure(figsize=(10, 6))
    plt.plot(E_eV, T, 'b-', linewidth=2, label='Transmission T(E)')
    plt.plot(E_eV, R, 'r-', linewidth=2, label='Réflexion R(E)')
    plt.xlabel('Énergie E (eV)')
    plt.ylabel('Coefficients T et R')
    plt.title('Effet Ramsauer-Townsend - Transmission et réflexion')
    plt.grid(True, alpha=0.3)
    plt.ylim(0, 1.1)
    plt.legend()
    filepath_combined = os.path.join(output_dir, f"combined_v0={v0}.png")
    plt.savefig(filepath_combined, dpi=300, bbox_inches='tight')
    
    print("Graphiques sauvegardés.")

ramsauer_townsend()