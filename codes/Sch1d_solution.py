# Sch1d_solution.py : observation de la propagation d'un paquet d'ondes 

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

def simulation(dt, nt, dx, nx, xc, sigma, E, v0, x_debut, x_fin, n_frame, intervalle):
    # Initialise l'animation
    def init():
        line.set_data([], [])
        return line,

    # Crée un graphique pour chaque densite sauvegarde
    def animate(j):
        line.set_data(x_grid, final_densite[j,:])
        return line,

    # Calculs des paramètres (simulation de l'équation de Schrödinger dépendante du temps)
    s = dt / (dx**2)
    k = math.sqrt(2 * E)  # Nombre d'onde
    A = 1 / (math.sqrt(sigma * math.sqrt(math.pi)))  # Constante de normalisation
    e = E / v0 # Rapport E/V0

    # Initialisation des grilles
    x_grid = np.linspace(0, (nx - 1) * dx, nx)  # Grille spatiale
    V = np.zeros(nx)  # Affichage du potentiel

    # Définition du potentiel dans la région spécifiée
    V[(x_grid >= x_debut) & (x_grid <= x_fin)] = -v0

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

    plt.ylim(-5, max(max_densite * 1.2, max_potentiel * 1.2))
    plt.xlim(0, np.max(x_grid))

    # Affichage du potentiel
    plt.plot(x_grid, V, 'orange', linewidth=3, label=f'Potentiel')

    # Labels et titre
    plt.title(plot_title)
    plt.xlabel('Position x (m)')
    plt.ylabel('Densité de probabilité (m^2) / Potentiel V (eV)')
    plt.legend()
    plt.grid()

    # Création de l'animation
    ani = animation.FuncAnimation(fig, animate, init_func=init, frames=n_frame, blit=False, interval=intervalle, repeat=True)

    # Création du dossier output s'il n'existe pas
    output_dir = "output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Sauvegarde
    file_name = f'simulation_e={e:.0f}_v0={v0:.0f}.mp4'
    file_path = os.path.join(output_dir, file_name)
    ani.save(file_path, writer = animation.FFMpegWriter(fps=120, bitrate=5000))
    print(f"\nSimulation enregistrée avec succès dans le dossier {output_dir}.")