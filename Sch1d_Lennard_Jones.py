# Sch1d_Lennard_Jones.py : observation de la propagation d'un paquet d'ondes avec potentiel de Lennard-Jones

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.linalg import eigh
import time
import os

start_time = time.time()

# Création du dossier output
output_dir = "output"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Calculer et afficher les états stationnaires avec potentiel de Lennard-Jones
def etats_stationnaires_lennard_jones(epsilon, sigma, dx, nx, n_states=5):
    x = np.linspace(-1, 1, nx)
    
    # Potentiel de Lennard-Jones
    V_lj = 4 * epsilon * ((sigma / (np.abs(x) + 1e-10))**12 - (sigma / (np.abs(x) + 1e-10))**6)
    
    # Limiter le potentiel pour éviter les valeurs extrêmes
    V_lj = np.clip(V_lj, -10*epsilon, 10*epsilon)
    
    # Construction de l'hamiltonien
    diag = np.full(nx, -2.0)
    offdiag = np.full(nx - 1, 1.0)
    T = (-1 / dx**2) * (np.diag(diag) + np.diag(offdiag, 1) + np.diag(offdiag, -1))
    H = T + np.diag(V_lj)
    
    # Calcul des états propres
    energies, states = eigh(H, subset_by_index=(0, min(n_states-1, nx-1)))
    
    plt.figure(figsize=(12,8))
    
    # Graphique du potentiel
    plt.subplot(2, 1, 1)
    plt.plot(x, V_lj, 'k-', linewidth=2, label='Potentiel de Lennard-Jones')
    plt.axhline(y=0, color='r', linestyle='--', alpha=0.5, label='Énergie nulle')
    plt.title("Potentiel de Lennard-Jones")
    plt.xlabel("Position x (m)")
    plt.ylabel("Énergie (eV)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 3*sigma)
    plt.ylim(-5 - epsilon, 35 + epsilon)
    
    # Graphique des états stationnaires
    plt.subplot(2, 1, 2)
    for i in range(min(n_states, len(energies))):
        psi = states[:, i]
        psi = psi / np.sqrt(np.sum(psi**2) * dx)
        # Décalage pour la visualisation
        offset = energies[i] * 0.1  # Facteur d'échelle pour la visualisation
        plt.plot(x, psi**2 + offset, label=f"État {i} (E = {energies[i]:.2f})")
    
    plt.axhline(y=0, color='r', linestyle='--', alpha=0.5)
    plt.title("États stationnaires dans le potentiel de Lennard-Jones")
    plt.xlabel("Position x (m)")
    plt.ylabel("Densité de probabilité + offset")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 1)
    
    plt.tight_layout()
    filepath = os.path.join("output", f"etats_stationnaires_LJ_eps={epsilon}_sigma={sigma}.png")
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    print(f"Graphique des états stationnaires Lennard-Jones exporté dans 'output'")
    plt.close()
    
    return energies, states

# Calculer le coefficient de transmission pour différentes énergies
def calcul_transmission_lennard_jones(epsilon, sigma, dx, nx, E_min=0.1, E_max=10, n_points=100):
    x = np.linspace(-1, 1, nx)
    energies = np.linspace(E_min, E_max, n_points)
    transmissions = []
    
    # Calcul du potentiel de Lennard-Jones (une seule fois)
    V_lj = 4 * epsilon * ((sigma / (np.abs(x) + 1e-10))**12 - (sigma / (np.abs(x) + 1e-10))**6)
    V_lj = np.clip(V_lj, -10*epsilon, 10*epsilon)

    # Construction de l'hamiltonien (une seule fois)
    diag = np.full(nx, -2.0)
    offdiag = np.full(nx - 1, 1.0)
    T_op = (-1 / dx**2) * (np.diag(diag) + np.diag(offdiag, 1) + np.diag(offdiag, -1))
    H = T_op + np.diag(V_lj)

    # Diagonalisation (une seule fois)
    eigenvals, eigenvecs = eigh(H)

    # Boucle seulement pour comparer aux énergies
    for E in energies:
        proche_E = np.abs(eigenvals - E)
        min_dist = np.min(proche_E)
        T = np.exp(-min_dist / (0.1 * epsilon))  # Approximation phénoménologique
        transmissions.append(T)



    
    plt.figure(figsize=(10, 6))
    plt.plot(energies, transmissions, 'b-', linewidth=2)
    plt.title("Coefficient de transmission vs Énergie (Potentiel de Lennard-Jones)")
    plt.xlabel("Énergie (eV)")
    plt.ylabel("Transmission T")
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 10 + epsilon)
    plt.ylim(0, 1.1)
    
    filepath = os.path.join("output", f"transmission_LJ_eps={epsilon}_sigma={sigma}.png")
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    print(f"Graphique de transmission Lennard-Jones exporté dans 'output'")
    plt.close()
    
    return energies, transmissions

# Initialise l'animation
def init():
    line.set_data([], [])
    return line,

# Crée un graphique pour chaque densité sauvegardée
def animate(j):
    line.set_data(o, final_densite[j,:]) 
    return line,

# Définition des paramètres
dt = 1E-7  # Pas de temps
nt = 90000  # Nombre de pas temps

dx = 0.001  # Pas d'espace
nx = int(2/dx)  # Nombre de points simulés (domaine de -1 à 1)

n_frame = int(nt/1000) + 1  # Nombre d'images dans notre animation
intervalle = 100

s = dt/(dx**2)

# Paramètres du potentiel de Lennard-Jones
epsilon = float(input("Profondeur du puits epsilon [ex: 100]: ") or 100)
sigma = float(input("Paramètre sigma [ex: 0.1]: ") or 0.1)

# Paramètres du paquet d'ondes
xc = -0.25  # Position initiale (à gauche)
sigma_paquet = 0.05  # Largeur du paquet

# Énergie du paquet d'ondes
E = float(input("Énergie du paquet d'ondes [ex: 5]: ") or 5)

A = 1/(math.sqrt(sigma_paquet*math.sqrt(math.pi)))  # Constante de normalisation
k = math.sqrt(2*abs(E))

# Initialisation des tableaux
o = np.linspace(-1, 1, nx)
V = np.zeros(nx)

# Potentiel de Lennard-Jones
V = 4 * epsilon * ((sigma / (np.abs(o) + 1e-10))**12 - (sigma / (np.abs(o) + 1e-10))**6)

# Limiter le potentiel pour éviter les valeurs extrêmes
V = np.clip(V, -10*epsilon, 10*epsilon)

# Calcul des états stationnaires
print("Calcul des états stationnaires...")
energies_stat, states_stat = etats_stationnaires_lennard_jones(epsilon, sigma, dx, nx)

# Calcul du coefficient de transmission
print("Calcul du coefficient de transmission...")
E_range, T_range = calcul_transmission_lennard_jones(epsilon, sigma, dx, nx)

# Paquet d'ondes initial
cpt = A * np.exp(1j * k * o - ((o - xc) ** 2) / (2 * (sigma_paquet ** 2)))
densite = np.zeros((nt, nx))
densite[0,:] = np.absolute(cpt[:]) ** 2
final_densite = np.zeros((n_frame, nx))

re = np.zeros(nx)
re[:] = np.real(cpt[:])

im = np.zeros(nx)
im[:] = np.imag(cpt[:])

b = np.zeros(nx)

print("Début de la simulation de propagation...")

# Résolution de l'équation de Schrödinger
it = 0
for i in range(1, nt):
    # Étape impaire, mettre à jour la partie imaginaire
    if i % 2 != 0:
        b[1:-1] = im[1:-1]
        im[1:-1] = im[1:-1] + s * (re[2:] + re[:-2]) - 2 * re[1:-1] * (s + V[1:-1] * dt)
        densite[i,1:-1] = re[1:-1]*re[1:-1] + im[1:-1]*b[1:-1]
    # Étape paire, mettre à jour la partie réelle
    else:
        re[1:-1] = re[1:-1] - s * (im[2:] + im[:-2]) + 2 * im[1:-1] * (s + V[1:-1] * dt)

# Sauvegarde périodique pour l'animation
for i in range(1, nt):
    if((i-1) % 1000 == 0):
        it += 1
        final_densite[it][:] = densite[i][:]

# Configuration du graphique d'animation
plot_title = f"Propagation dans potentiel de Lennard-Jones (ε={epsilon}, σ={sigma})"

fig = plt.figure(figsize=(12, 8))
line, = plt.plot([], [], 'b-', linewidth=2, label='Densité de probabilité')

# Normalisation du potentiel pour l'affichage
V_display = V / np.max(np.abs(V)) * 5  # Normalisation pour l'affichage

plt.ylim(-1, 15)
plt.xlim(-1, 1)
plt.plot(o, V_display, 'r-', linewidth=2, label='Potentiel (normalisé)')
plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
plt.title(plot_title)
plt.xlabel("Position x (m)")
plt.ylabel("Densité de probabilité")
plt.legend()
plt.grid(True, alpha=0.3)

print("Création de l'animation...")
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=n_frame, 
                            blit=False, interval=intervalle, repeat=False)

file_name = f"paquet_onde_LJ_eps={epsilon}_sigma={sigma}_E={E}.mp4"
file_path = os.path.join(output_dir, file_name)
ani.save(file_path, writer=animation.FFMpegWriter(fps=15, bitrate=5000))

# Création d'un graphique statique final
plt.figure(figsize=(12, 8))
plt.plot(o, final_densite[-1,:], 'b-', linewidth=2, label='Densité finale')
plt.plot(o, V_display, 'r-', linewidth=2, label='Potentiel (normalisé)')
plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
plt.title(f"État final - Potentiel de Lennard-Jones (ε={epsilon}, σ={sigma})")
plt.xlabel("Position x (m)")
plt.ylabel("Densité de probabilité")
plt.legend()
plt.grid(True, alpha=0.3)
plt.xlim(-1, 1)

final_plot_path = os.path.join(output_dir, f"etat_final_LJ_eps={epsilon}_sigma={sigma}.png")
plt.savefig(final_plot_path, dpi=300, bbox_inches='tight')
plt.close()

# Affichage du temps
end_time = time.time()
elapsed_time = end_time - start_time

print(f"\n=== RÉSULTATS ===")
print(f"Simulation terminée en {elapsed_time:.2f} secondes.")
print(f"Paramètres utilisés:")
print(f"  - Epsilon (profondeur du puits): {epsilon}")
print(f"  - Sigma (position du minimum): {sigma}")
print(f"  - Énergie du paquet d'ondes: {E}")
print(f"  - Position initiale: {xc}")
print(f"\nFichiers générés dans le dossier 'output':")
print(f"  - Animation: {file_name}")
print(f"  - États stationnaires: etats_stationnaires_LJ_eps={epsilon}_sigma={sigma}.png")
print(f"  - Transmission: transmission_LJ_eps={epsilon}_sigma={sigma}.png")
print(f"  - État final: etat_final_LJ_eps={epsilon}_sigma={sigma}.png")