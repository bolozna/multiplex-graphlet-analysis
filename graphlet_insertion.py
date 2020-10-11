# insert random graphlets into multiplex networks
import random
import math
import pymnet

def insert_graphlet(M,graphlet,amount):
    ntot = len(M.slices[0])
    ltot = len(M.slices[1])
    n = len(graphlet.slices[0])
    l = len(graphlet.slices[1])
    if isinstance(amount,int):
        insertions = amount
    elif isinstance(amount,float):
        f = math.factorial
        total_node_layer_combinations = (f(ntot)/f(n)/f(ntot-n))*(f(ltot)/f(l)/f(ltot-l))
        insertions = int(amount*total_node_layer_combinations)
    rng = random.Random()
    for insertion_loc in sample_insertion_locs(n,l,insertions,M,rng):
        insert_graphlet_into_loc(insertion_loc,M,graphlet,rng)
    
def insert_graphlet_into_loc((sn,sl),M,graphlet,rng):
    # make sure this works with any couplings
    edges_init = set(pymnet.subnet(M,sn,sl).edges)
    mapping_nodes = dict()
    mapping_layers = dict()
    for n in graphlet.iter_nodes():
        dest_node = rng.choice(tuple(sn))
        mapping_nodes[n] = dest_node
        sn.remove(dest_node)
    for l in graphlet.iter_layers():
        dest_layer = rng.choice(tuple(sl))
        mapping_layers[l] = dest_layer
        sl.remove(dest_layer)
    relabeled_graphlet = pymnet.transforms.relabel(graphlet,mapping_nodes,mapping_layers)
    edges_target = set(relabeled_graphlet.edges)
    edges_to_remove = edges_init-edges_target
    edges_to_add = edges_target-edges_init
    for e in edges_to_remove:
        M[e[0],e[1],e[2],e[3]] = 0
    for e in edges_to_add:
        M[e[0],e[1],e[2],e[3]] = 1

def sample_insertion_locs(nnodes,nlayers,nsamples,M,rng):
    nodes = list(M.iter_nodes())
    layers = list(M.iter_layers())
    sample = set()
    while len(sample) < nsamples:
        sn = set()
        sl = set()
        while len(sn) < nnodes:
            sn.add(rng.choice(nodes))
        while len(sl) < nlayers:
            sl.add(rng.choice(layers))
        sample.add((sn,sl))
    return sample

  