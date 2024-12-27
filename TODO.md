# TODO
- [x] Créer un Iplotclaw plus facile à utiliser: [arrowkeyplot.py](https://github.com/giboul/visclaw/blob/Iplot_mpl/src/python/visclaw/arrowkeyplot.py)
- [x] Écrire le fichier `qinit.xyz` par un algorithme de remplissage: fonction [`flood_mask`](https://github.com/giboul/TriftGeoClaw/blob/main/Tsunami/makeqinit.py#L26) dans `makeqinit.py`.
- [x] Enregistrer les flux de neige entrants.
- [x] Interpoler les flux entrants et les amortir par $\rho_\text{neige}/\rho_\text{eau}$ pour créer la bonne vague d'impulsion.
- [x] Ajouter des graphes de l’intumescence au barrage (~~ligne de jauges~~ script python).
- [x] Simuler toutes les avalanches proposées.
- [x] Il faudrait que les avalanches soient introduites plus près du la rétention, la grande distance parcourue entre le bord et le lac est parfois un problème. => Interpolation de $q$ au-dessus du lac.
- [x] Faciliter l'utilisation avec un fichier `config.yaml` (chemin vers le ficier topo, paramètres redondants, préférences...)
- [x] Faire un diagramme expliquant le fonctionnement de GeoClaw. => [AMR flowcharts](https://www.clawpack.org/dev/amrclaw_flowcharts.html)
- [x] Rajouter une région de maille fine autour du barrage (éviter la rupture).
- [ ] Comparer les résultats à ceux des expériences VAW+LHE (test des conditions de bord).

# À explorer 
- code Kyle Mandli : two layer shallow water (water and avalanche?) `rpn2_layered` (clawpack/dev)

