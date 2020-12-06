# insert random graphlets into multiplex networks
# currently only works with 'none' couplings
import random
import math
import pymnet
import itertools

def insert_random_graphlets(M_list,nnodes,nlayers,number_of_graphlets,amounts,allowed_aspects='all',balance=True,print_progress=True):
    # NB! currently only works for integer amounts because of location precomputing
    all_graphlets,invs = pymnet.graphlets.graphlets(nnodes,list(range(nlayers)),couplings=None,allowed_aspects=allowed_aspects)
    random_graphlets = set()
    while len(random_graphlets) < number_of_graphlets:
        random_graphlets.add(random.choice(all_graphlets[nnodes]))
    if print_progress:
        print('Inserting graphlets...')
    count = 0
    total_count = len(M_list)
    for M in M_list:
        if balance:
            forbidden_edge_locs = dict()
            init_edge_numbers = dict()
            for l in M.iter_layers():
                forbidden_edge_locs[l] = set()
                init_edge_numbers[l] = len(M.A[l].edges)
        else:
            forbidden_edge_locs = None
        # precompute non-overlapping locations for each graphlet
        rng = random.Random()
        all_locs = list(sample_insertion_locs(nnodes,nlayers,sum(amounts),M,rng))
        rng.shuffle(all_locs)
        # insert graphlets
        for ii,graphlet in enumerate(random_graphlets):
            precomputed_insertion_locs = all_locs[0:amounts[ii]]
            insert_graphlet(M,graphlet,amounts[ii],forbidden_edge_locs=forbidden_edge_locs,precomputed_insertion_locs=precomputed_insertion_locs)
            del all_locs[0:amounts[ii]]
        if balance:
            for l in M.iter_layers():
                curr_edge_number = len(M.A[l].edges)
                if curr_edge_number > init_edge_numbers[l]:
                    curr_edges = list(M.A[l].edges)
                    deleted_edges = set()
                    while curr_edge_number > init_edge_numbers[l]:
                        candidate_delet = random.choice(curr_edges)
                        if (candidate_delet[0],candidate_delet[1]) not in forbidden_edge_locs[l] and (candidate_delet[1],candidate_delet[0]) not in forbidden_edge_locs[l] and candidate_delet not in deleted_edges:
                            M[candidate_delet[0],candidate_delet[1],l,l] = 0
                            deleted_edges.add(candidate_delet)
                            curr_edge_number -= 1
                elif curr_edge_number < init_edge_numbers[l]:
                    M_nodes = list(M.iter_nodes())
                    curr_edges = set([(e[0],e[1]) for e in M.A[l].edges])
                    added_edges = set()
                    while curr_edge_number < init_edge_numbers[l]:
                        candidate_add = set()
                        while len(candidate_add) < 2:
                            candidate_add.add(random.choice(M_nodes))
                        tca = tuple(candidate_add)
                        if (tca[0],tca[1]) not in curr_edges and (tca[1],tca[0]) not in curr_edges and (tca[0],tca[1]) not in forbidden_edge_locs[l] and (tca[1],tca[0]) not in forbidden_edge_locs[l] and tca not in added_edges:
                            M[tca[0],tca[1],l,l] = 1
                            added_edges.add(tca)
                            curr_edge_number += 1
        count += 1
        if print_progress:
                print('+'*count+'-'*(total_count-count))

def insert_graphlet(M,graphlet,amount,forbidden_edge_locs=None,precomputed_insertion_locs=None):
    ntot = len(M.slices[0])
    ltot = len(M.slices[1])
    n = len(graphlet.slices[0])
    l = len(graphlet.slices[1])
    if isinstance(amount,int):
        insertions = amount
    elif isinstance(amount,float):
        f = math.factorial
        # this amount not the exact total amount because of forbidden overlaps
        total_node_layer_combinations = (f(ntot)/f(n)/f(ntot-n))*(f(ltot)/f(l)/f(ltot-l))
        insertions = int(amount*total_node_layer_combinations)
    rng = random.Random()
    if precomputed_insertion_locs:
        insertion_locs = precomputed_insertion_locs
    else:
        insertion_locs = sample_insertion_locs(n,l,insertions,M,rng)
    for insertion_loc in insertion_locs:
        insert_graphlet_into_loc(insertion_loc,M,graphlet,rng)
        if forbidden_edge_locs:
            for l_insert in insertion_loc[1]:
                for e in itertools.combinations(insertion_loc[0],2):
                    forbidden_edge_locs[l_insert].add(e)
    
def insert_graphlet_into_loc((sn,sl),M,graphlet,rng):
    # unfreeze sn and sl
    sn = set(sn)
    sl = set(sl)
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
        # for...else structure
        for existing_loc in sample:
            if sn.intersection(existing_loc[0]) > 1 and sl.intersection(existing_loc[1]) > 0:
                break
        else:
            sample.add((frozenset(sn),frozenset(sl)))
    return sample

