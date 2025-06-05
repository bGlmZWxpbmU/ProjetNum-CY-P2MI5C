# Sch1d.py : observation de la propagation d'un paquet d'ondes 

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import os

start_time = time.time()

# Création du dossier output
output_dir = "output"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

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
o = np.linspace(0, (nx - 1) * dx, nx)
V = np.zeros(nx)

# V[o >= 1] = v0
V[(o >= 0.8) & (o<=0.9)] = v0  # Potentiel dans la région spécifiée

cpt = A * np.exp(1j * k * o - ((o - xc) ** 2) / (2 * (sigma ** 2))) # Paquet d'ondes cpt = A * e^{i k x}
densite=np.zeros((nt,nx))
densite[0,:] = np.absolute(cpt[:]) ** 2
final_densite=np.zeros((n_frame,nx))

re=np.zeros(nx)
re[:]=np.real(cpt[:])

im=np.zeros(nx)
im[:]=np.imag(cpt[:])

b=np.zeros(nx)

# Résolution de l'équation
it=0
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