# multiplex-graphlet-analysis

Library for performing graphlet analysis on multiplex networks.
Details and theoretical basis in the accompanying [arXiv paper](https://arxiv.org/abs/2106.13011).

To use, run (in order):
- pipeline.make_networks(...)
- pipeline.make_orbits(...)
- pipeline.make_gcds(...)
- pipeline.make_figures(...)

The set of parameters required differs by what you want to use as the network test set. The functions make_orbits, make_gcds, make_figures use the same parameters as make_networks PLUS the parameter allowed_aspects_orbits.

To run the multiplex analysis presented in the arXiv paper, use:
- different models, constant degree: `test_set_type='random',n_nets=30,n_n=1000,n_l=3,m=2,use_simple_conf=False,use_simple_conf_plex=True` and `allowed_aspects_orbits='all'`
- different models, degree progression: `test_set_type='random_deg_progression',n_nets=30,n_n=1000,n_l=3,m=[1,2,3,4,5,6],use_simple_conf=False,use_simple_conf_plex=True` and `allowed_aspects_orbits='all'`
- graphlet insertion 4 nodes 2 layers: `test_set_type='graphlet_insertion',n_nets=30,n_n=1000,n_l=3,m=2,allowed_aspects_graphlets='all',n_classes=5,n_different_graphlets=20,graphlet_frequency=3,graphlet_size=(4,2)` and `allowed_aspects_orbits='all'`
- graphlet insertion 3 nodes 3 layers: `test_set_type='graphlet_insertion',n_nets=30,n_n=1000,n_l=3,m=2,allowed_aspects_graphlets='all',n_classes=5,n_different_graphlets=10,graphlet_frequency=3,graphlet_size=(3,3)` and `allowed_aspects_orbits='all'`

