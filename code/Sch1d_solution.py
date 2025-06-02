# Sch1d_solution.py : observation de la propagation d’un paquet d’ondes 

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import os

# Collecte les paramètres de simulation dans le terminal
def get_inputs():
    print("--- Simulation de paquet d'ondes ---\n")
    
    # Paramètres temporels
    print("--- Paramètres temporels ---")
    dt = float(input("Pas de temps (s) [valeur positive, ex: 1e-7]: ") or 1e-7)
    nt = int(input("Nombre d'itérations temporelles [entier positif, ex: 90000]: ") or 90000)

    n_frame = int(nt/250) + 1 # Nombre d'images dans notre animation
    intervalle = 1
    
    # Paramètres spatiaux
    print("\n--- Paramètres spatiaux ---")
    dx = float(input("Pas spatial dx (m) [valeur positive, ex: 0.001]: ") or 0.001)
    L_total = float(input("Longueur totale du domaine L (m) [valeur positive, ex: 2]: ") or 2.0)
    nx = int(L_total / dx)
    
    # Paramètres du paquet d'ondes
    print("\n--- Paquet d'ondes ---")
    xc = float(input("Position initiale du paquet (m) [ex: 0.6]: ") or 0.6)
    sigma = float(input("Largeur du paquet (m) [valeur positive, ex: 0.05]: ") or 0.05)
    
    # Paramètres du potentiel
    v0 = -5000  # Potentiel de référence (eV)
    e = 5  # Rapport E/V0 de référence
    E = e * v0  # Énergie calculée (eV)
    
    print("\n--- Puit de potentiel ---")
    x_debut = float(input("Position de début du potentiel (m) [ex: 0.8]: ") or 0.8)
    x_fin = float(input("Position de fin du potentiel (m) [ex: 0.9]: ") or 0.9)
    
    return dt, nt, dx, nx, xc, sigma, E, v0, x_debut, x_fin, n_frame, intervalle

# Initialise l'animation
def init():
    line.set_data([], [])
    return line,

# Crée un graphique pour chaque densite sauvegarde
def animate(j):
    line.set_data(x_grid, final_densite[j,:])
    return line,

dt, nt, dx, nx, xc, sigma, E, v0, x_debut, x_fin, n_frame, intervalle = get_inputs()

start_time = time.time()

# Calculs des paramètres (simulation de l'équation de Schrödinger dépendante du temps)
s = dt / (dx**2)
k = math.sqrt(2 * abs(E))  # Nombre d'onde (m^-1)
A = 1 / (math.sqrt(sigma * math.sqrt(math.pi)))  # Constante de normalisation (m^-1/2)
e = E / v0 # Rapport E/V0

# Initialisation des grilles
x_grid = np.linspace(0, (nx - 1) * dx, nx)  # Grille spatiale (m)
V = np.zeros(nx)  # Potentiel (eV)

# Définition du potentiel dans la région spécifiée
V[(x_grid >= x_debut) & (x_grid <= x_fin)] = v0

# Fonction d'onde complexe
psi = A * np.exp(1j * k * x_grid - ((x_grid - xc) ** 2) / (2 * (sigma ** 2)))

# Séparation partie réelle et imaginaire
re = np.real(psi)
im = np.imag(psi)

# Stocker les densités de probabilité
densite = np.zeros((nt, nx))
densite[0,:] = np.absolute(psi) ** 2  # Densité initiale |psi|^2

final_densite = np.zeros((n_frame, nx))  # Densités pour l'animation
final_densite[0,:] = densite[0,:]  # Première image

b = np.zeros(nx)  # Tableau temporaire

# Résolution de l'équation
it = 1
for i in range(1, nt):

    # Étape impaire : mettre à jour la partie imaginaire
    if i % 2 != 0:  
        b[1:-1] = im[1:-1]

        # Évolution de Im(psi)
        im[1:-1] = im[1:-1] + s * (re[2:] + re[:-2]) - 2 * re[1:-1] * (s + V[1:-1] * dt)

        # Calcul de la densité de probabilité |psi|^2
        densite[i, 1:-1] = re[1:-1]**2 + im[1:-1] * b[1:-1]

    # Étape paire : mettre à jour la partie réelle
    else:  
        # Évolution de Re(psi)
        re[1:-1] = re[1:-1] - s * (im[2:] + im[:-2]) + 2 * im[1:-1] * (s + V[1:-1] * dt)
    
    # Sauvegarde périodique pour l'animation
    if (i-1) % (nt // n_frame) == 0 and it < n_frame:
        final_densite[it,:] = densite[i,:]
        it += 1

# Titre du graphique
plot_title = f"Marche Ascendante avec E/V₀ = {e:.2f}, E = {E:.0f} eV, V₀ = {v0:.0f} eV"

# Configuration du graphique
fig = plt.figure(figsize=(12, 8))

# Affichage de la densité
line, = plt.plot([], [], 'b-', linewidth=2, label='Densité')

# Limites d'affichage
max_densite = np.max(final_densite)
max_potentiel = max(np.max(V), 0)

plt.ylim(0, max(max_densite * 1.2, max_potentiel * 1.2))
plt.xlim(0, np.max(x_grid))

# Affichage du potentiel
plt.plot(x_grid, V, 'orange', linewidth=3, label=f'Potentiel')

# Labels et titre
plt.title(plot_title)
plt.xlabel('Position x (m)')
plt.ylabel('Densité de probabilité (m⁻¹) / Potentiel V (eV)')
plt.legend()
plt.grid()

# Création de l'animation
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=n_frame, blit=False, interval=intervalle, repeat=True)

# Création du dossier output
output_dir = "../output"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Sauvegarde
file_name = f'simulation_e={e:.0f}_v0={v0:.0f}.mp4'
file_path = os.path.join(output_dir, file_name)
ani.save(file_path, writer = animation.FFMpegWriter(fps=120, bitrate=5000))

# Temps d'exécution
end_time = time.time()
elapsed_time = end_time - start_time
print(f"\nSuccès. Temps d'exécution total: {elapsed_time:.2f} secondes")

plt.show()