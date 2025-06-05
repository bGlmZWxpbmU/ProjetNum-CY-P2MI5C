# Sch1d_stationnaires.py : observation de la propagation et calcul des états stationnaires

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

# Calculer et afficher les états stationnaires liés (E < 0, V = -V0)
def etats_stationnaires_lies(e, dx, nx, V, n_states=5):
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
    plt.title("États stationnaires liés (E < 0, V = -V0)")
    plt.xlabel("x (m)")
    plt.ylabel("Énergie (eV) / Densité de probabilité")
    plt.legend()
    plt.grid()
    filepath = os.path.join("output", f"etats_stationnaires_lies_e={e}.png")
    plt.savefig(filepath)
    print(f"Graphique des états stationnaires liés exporté dans 'output'")
    plt.close()

# Calculer et afficher les états stationnaires libres (E > 0, V = 0)
def etats_stationnaires_libres(e, dx, nx, n_states=5):
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
    plt.xlabel("Position x (m)")
    plt.ylabel("Énergie (eV) / Densité de probabilité")
    plt.legend()
    plt.grid()
    filepath = os.path.join("output", f"etats_stationnaires_libres_e={e}.png")
    plt.savefig(filepath)
    print(f"Graphique des états stationnaires libres exporté dans 'output'")
    plt.close()

# Initialise l'animation
def init():
    line.set_data([], [])
    return line,

#Crée un graphique pour chaque densite sauvegarde
def animate(j):
    line.set_data(o, final_densite[j,:]) 
    return line,

# Définition des paramètres

dt=1E-7 # Pas de temps
nt=90000 # Nombre de pas temps

dx=0.001 # Pas d'espace
nx=int(1/dx)*2 # Nombre de points simulés

n_frame=int(nt/1000)+1 # Nombre d'images dans notre animation
intervalle=100

s=dt/(dx**2)

xc=0.6 # Position initiale
sigma=0.05

v0=-4000 # Potentiel
e=float(input("Rapport E/V0 [ex: 5]: ") or 5) # Saisie utilisateur du rapport e = E/V0
E=e*v0 # Énergie

A=1/(math.sqrt(sigma*math.sqrt(math.pi))) # Constante de normalisation
k=math.sqrt(2*abs(E))

o=np.zeros(nx)
V=np.zeros(nx)

# Initialisation des tableaux
o=np.linspace(0, (nx - 1) * dx, nx)
V=np.zeros(nx)

V[(o >= 0.8) & (o<=0.9)] = v0  # Potentiel dans la région spécifiée

# Calcul des états stationnaires liés (E < 0, V = -V0)
etats_stationnaires_lies(e, dx, nx, V)

# Calcul des états stationnaires libres (E > 0, V = 0)
etats_stationnaires_libres(e, dx, nx)

cpt = A * np.exp(1j * k * o - ((o - xc) ** 2) / (2 * (sigma ** 2))) # Paquet d'ondes cpt = A * e^{i k x}
densite = np.zeros((nt,nx))
densite[0,:] = np.absolute(cpt[:]) ** 2
final_densite = np.zeros((n_frame,nx))

re = np.zeros(nx)
re[:]=np.real(cpt[:])

im=np.zeros(nx)
im[:]=np.imag(cpt[:])

b=np.zeros(nx)

# Résolution de l'équation
it = 0
for i in range(1, nt):
    # Étape impaire, mettre à jour la partie imaginaire
    if i % 2 != 0:
        b[1:-1]=im[1:-1]
        im[1:-1] = im[1:-1] + s * (re[2:] + re[:-2]) - 2 * re[1:-1] * (s + V[1:-1] * dt)
        densite[i,1:-1] = re[1:-1]*re[1:-1] + im[1:-1]*b[1:-1]

    # Étape paire, mettre à jour la partie imaginaire
    else:
        re[1:-1] = re[1:-1] - s * (im[2:] + im[:-2]) + 2 * im[1:-1] * (s + V[1:-1] * dt)

# Sauvegarde périodique pour l'animation
for i in range(1,nt):
    if((i-1)%1000==0):
        it+=1
        final_densite[it][:]=densite[i][:]

# Configuration du graphique

plot_title = "Marche Ascendante avec E/Vo = "+str(e)

fig = plt.figure() # Initialise la figure principale
line, = plt.plot([], [])
plt.ylim(-3,13)
plt.xlim(0,2)
plt.plot(o,V,label="Potentiel")
plt.title(plot_title)
plt.xlabel("x (m)")
plt.ylabel("Densité de probabilité de présence (m^-1)")
plt.legend()

ani = animation.FuncAnimation(fig,animate,init_func=init, frames=n_frame, blit=False, interval=intervalle, repeat=False)
file_name = f"paquet_onde_e={e}.mp4"
file_path = os.path.join(output_dir, file_name)
ani.save(file_path, writer = animation.FFMpegWriter(fps=15, bitrate=5000))

# Affichage du temps
end_time = time.time()
elapsed_time = end_time - start_time

print(f"\nSimulation enregistrée avec succès en {elapsed_time:.2f} secondes.")

plt.show()
