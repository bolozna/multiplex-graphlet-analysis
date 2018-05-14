import pymnet
import os
import itertools
#import pipeline
import pipeline
import pandas as pd
import re

def graphlet_degree_distributions(network,nnodes,nlayers,allowed_aspects='all',save_name=None):
    
    if save_name != None:
        directory = 'Results/'+save_name+'_'+str(nlayers)
        if not os.path.exists(directory):
            os.makedirs(directory)
        
    #net_nlayers = sum(1 for _ in network.iter_layers())
    if nlayers == 1:
        net_layers = [0]
    else:
        net_layers = list(network.iter_layers())
    net_nlayers = len(net_layers)
    net_nodes = list(network.iter_nodes())
    
    graphlets,invs = pymnet.graphlets.graphlets(nnodes,net_layers,nlayers,couplings=network.couplings[0][0],allowed_aspects=allowed_aspects)
    auts = pymnet.graphlets.automorphism_orbits(graphlets,allowed_aspects=allowed_aspects)
    orbit_is = pymnet.graphlets.orbit_numbers(nnodes,graphlets,auts)
    orbit_list = pymnet.graphlets.ordered_orbit_list(orbit_is)

    if nlayers == 1:
        agg_net = pymnet.MultiplexNetwork()
        for node in net_nodes:
            agg_net.add_node(node)
        for edge in network.edges:
            agg_net[edge[0],edge[1],0] = 1
        network = agg_net
    
    for layer_comb in itertools.combinations(net_layers,nlayers):
        sub_net = pymnet.subnet(network,net_nodes,layer_comb)
        orbits = pymnet.graphlets.orbit_counts_all(sub_net,nnodes,graphlets,invs,auts,orbit_list,allowed_aspects=allowed_aspects)
        
        if save_name != None:
            f_name = directory+'/'+save_name
            for layer in layer_comb:
                f_name += '_'+str(layer)
            f_name += '.txt'
            pipeline.write_orbit_counts(orbits,f_name,net_nodes,orbit_list)
        
    return orbit_list

def read_graphlet_degree_distribution_folder(folder):
    '''
    Reads GDDs from a folder containing files for individual layer combinations.
    Every file must have the same first line (same orbits).
    '''
    orbits = pd.DataFrame()
    header = None
    for orbit_file in os.listdir(folder):
        f_name = folder+'/'+orbit_file
        with open(f_name) as f:
            file_header = f.readline()
        if header == None:
            header = file_header
        else:
            assert file_header == header, 'Files need to have the same first line!'
            
        orbit_list = re.findall('\(.*?\)',header)
        
        if orbits.empty:
            orbits = pd.read_csv(f_name,index_col=0,skiprows=[0],header=None,names=orbit_list)
            orbits.index.name = None
        else:
            orbits2 = pd.read_csv(f_name,index_col=0,skiprows=[0],header=None,names=orbit_list)
            orbits2.index.name = None
            print(orbits2)
            orbits = orbits.add(orbits2, fill_value=0)
    return orbits

