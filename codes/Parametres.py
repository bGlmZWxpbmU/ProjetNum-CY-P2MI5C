# Parametres.py : collecte les paramètres de simulation dans le terminal

def get_parametres():
    # Paramètres temporels
    print("--- Paramètres temporels ---")
    dt = float(input("Pas de temps (s) [ex: 2e-7]: ") or 2e-7)
    nt = int(input("Nombre d'itérations temporelles [ex: 90000]: ") or 90000)


    # Paramètres spatiaux
    print("\n--- Paramètres spatiaux ---")
    dx = float(input("Pas spatial (m) [ex: 0.001]: ") or 0.001)
    L_total = float(input("Longueur totale du domaine (m) [ex: 2]: ") or 2.0)
    nx = int(L_total / dx)
    
    # Paramètres du paquet d'ondes
    print("\n--- Paquet d'ondes ---")
    xc = float(input("Position initiale du paquet (m) [ex: 0.6]: ") or 0.6)
    sigma = float(input("Largeur du paquet (m) [ex: 0.05]: ") or 0.05)
    
    # Paramètres du potentiel
    print("\n--- Potentiel ---")
    v0 = float(input("Profondeur du potentiel V0 (eV) [ex: 4000]: ") or 4000.0)
    e = float(input("Rapport E / V0 [ex: 0.9]: ") or 0.9)
    E = e * v0  # Énergie calculée (eV)
    
    print("\n--- Puit de potentiel ---")
    x_debut = float(input("Position de début du potentiel (m) [ex: 0.8]: ") or 0.8)
    x_fin = float(input("Position de fin du potentiel (m) [ex: 0.9]: ") or 0.9)

    n_frame = int(nt/100)+1 # Nombre d'images dans notre animation
    intervalle = 1000
    
    return dt, nt, dx, nx, xc, sigma, E, v0, x_debut, x_fin, n_frame, intervalle