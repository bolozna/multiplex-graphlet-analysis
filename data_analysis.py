import os
import numpy as np
import pandas as pd
from .orbits import redundants
from .graphlet_measures import GCM, GCD_matrix
from .independent_equations import independent_equations, redundant_orbits

def precision_recall(distances, groups):
    '''
    Computes precision-recall curves
    
    Parameters
    ----------
    distances : 2-d array, list of lists
        Distances between networks
    groups : dict
        Group labels of networks, keys are the network indices
        
    Returns
    -------
    pres : np.array
        Precision
    recs : np.array
        Recall
    '''
    
    n = len(distances)
    dists = np.unique(distances)
    u = len(dists)
    
    TP = np.zeros(u)
    FP = np.zeros(u)
    FN = np.zeros(u)
    
    inds = {}
    for i, d in enumerate(dists):
        inds[d] = i
        
    for i in range(n):
        for j in range(i+1, n):
            d = distances[i][j]
            k = inds[d]
            if groups[i] == groups[j]:
                TP[k:] += 1.0
                FN[:k] += 1.0
                
            else:
                FP[k:] += 1.0
                
    if TP[0] + FP[0] == 0:
        TP = TP[1:]
        FP = FP[1:]
        FN = FN[1:]
    
    pres = TP / (TP + FP)
    recs = TP / (TP + FN)
    
    return pres, recs
    
    
def area_under_precision_recall(pres, recs):
    '''
    Compute AUPRs
    
    Parameters
    ----------
    pres : np.array, list
        Precision
    recs : np.array, list
        Recall
        
    Returns
    -------
    aupr : float
    '''
    
    aupr = 0
    pre0 = 1.0
    rec0 = 0.0
    for pre1, rec1 in zip(pres, recs):
        aupr += (rec1 - rec0) * (pre0 + pre1) / 2.0
        pre0 = pre1
        rec0 = rec1
    
    return aupr
    
    
def GCDs(nets, n, n_l, layers, res_dir, orbit_is, orbit_list, no_reds=False, allowed_aspects='all'):
    '''
    Graphlet correlation distances between networks computed using graphlets 
    with n_l layers and up to n nodes.
    
    Parameters
    ----------
    nets : list of strs
        Network names used in the file names
    n : int
        Max number of nodes
    n_l : int
        Number of layers
    layers : list
        Layers used to generate the graphlets
    res_dir : str
        Directory where orbit counts are stored
    orbit_list : list of orbits
        All orbits with n_l layers and max n nodes in the order they are saved
        in the files.
    no_reds : boolean
        If True redundant orbits are removed
        
    Returns
    -------
    gcds : list of lists
        Graphlet correlation distances
    '''
    
    gcms = []
    
    if no_reds:
        inds, eqs = independent_equations(n, n_l, layers, allowed_aspects=allowed_aspects)
        reds = redundant_orbits(inds, eqs, orbit_is, orbit_list) #redundants[n_l]
    
    for net in nets:
        o_dir = res_dir + net + "_" + str(n_l)
        orbits = sum_orbit_counts(o_dir, orbit_list)
        if n_l <= 1:
            cols = get_column_names(orbit_list)
            orbits = orbits.rename(columns=cols)
        
        if no_reds:
            orbits = orbits.drop(reds, axis=1)
            
        col_names = list(orbits)
        for col in col_names:
            if col[1] > str(n):
                orbits = orbits.drop([col], axis=1)
        
        gcm = GCM(orbits)
        gcms.append(gcm)
        
    gcds = GCD_matrix(gcms)
    
    return gcds
    
    
def sum_orbit_counts(o_dir, orbit_list):
    '''
    Sums the orbit counts from all the files in o_dir
    
    Parameters
    ----------
    o_dir : str
        Directory of the orbit counts
    orbit_list :
    
    Returns
    -------
    orbits : pandas dataframe
        Summed orbit counts
    '''
    
    dtypes = {}
    col_names = []
    for orbit in orbit_list:
        dtypes[str(orbit)] = int
        col_names.append(str(orbit))
        
    dtypes['n'] = str
    
    orbits = pd.DataFrame()
    for file in os.listdir(o_dir):
        f_name = o_dir + '/' + file
        if orbits.empty:
            orbits = pd.read_csv(f_name, dtype=dtypes, header=None, skiprows=[0], 
                                 names=col_names, usecols=range(1, len(col_names) + 1)) #index_col=0
        else:
            orbits2 = pd.read_csv(f_name, dtype=dtypes, header=None, skiprows=[0], 
                                  names=col_names, usecols=range(1, len(col_names) + 1)) #index_col=0
            orbits = orbits.add(orbits2, fill_value=0)
            
    return orbits
    
    
def get_column_names(orbit_list):
    '''
    for matching the order of the monoplex orbits with literature
    '''
    
    cols = {}
    for orbit in orbit_list:
        i = orbit[0]
        j = orbit[1]
        k = orbit[2]
        if i == 4 and j == 0:
            cols[str((i,j,k))] = str((i,1,k))
            
        elif i == 4 and j == 1:
            cols[str((i,j,k))] = str((i,0,k))
            
        elif i == 4 and j == 2:
            cols[str((i,j,k))] = str((i,3,k))
            
        elif i == 4 and j == 3:
            cols[str((i,j,k))] = str((i,2,k))
            
    return cols