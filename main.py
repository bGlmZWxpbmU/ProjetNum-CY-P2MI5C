#!/usr/bin/env python3

# main.py : boucle principale de la simulation

import time

def main():
    start_time = time.time()

    from codes.Parametres import get_parametres

    # Récupère les entrées utilisateur
    print("--- Configuration du paquet d'ondes ---\n")
    dt, nt, dx, nx, xc, sigma, E, v0, x_debut, x_fin, n_frame, intervalle = get_parametres()

    choix = input("\nSouhaitez-vous modéliser la propagation de l'onde ? (o/n) : ")
    if choix == 'o':
        from codes import Sch1d_solution
        Sch1d_solution.simulation(dt, nt, dx, nx, xc, sigma, E, v0, x_debut, x_fin, n_frame, intervalle)

    choix = input("\nSouhaitez-vous calculer les états stationnaires ? (o/n) : ")
    if choix == 'o':
        from codes import Etats_stationnaires
        Etats_stationnaires.calculer(dx, nx, v0, x_debut, x_fin, E)

    # Temps d'exécution
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"\nProgramme exécuté avec succès en {elapsed_time:.2f} secondes")

if __name__ == "__main__":
    main()