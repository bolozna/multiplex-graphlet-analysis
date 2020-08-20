import data_analysis, visualization
import pymnet
from pymnet import graphlets
import os
import time
import itertools
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import manifold
from matplotlib.ticker import NullFormatter
from collections import defaultdict as dd
import interface

def main():
    
    n_nets = 5
    n_n = 10
    n_l = 3
    m = 2
    allowed_aspects = 'all'
    
    networks, net_names, boundaries, labels = example_networks(n_nets, n_n, n_l, m)
    
    directory = 'Results'
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    layers = list(range(n_l))
    orbit_lists = {}
    orbit_is_d = {}
    start = time.time()
    nn_nl = [(1,4), (2,4), (3,3)]
    for n_l, n in nn_nl:
        if n_l == 1:
            net_layers = [0]
        else:
            net_layers = layers
            
        nets, invs = graphlets.graphlets(n, net_layers, n_l, allowed_aspects=allowed_aspects)
        auts = graphlets.automorphism_orbits(nets, allowed_aspects=allowed_aspects)
        orbit_is = graphlets.orbit_numbers(n, nets, auts)
        orbit_is_d[n_l] = orbit_is
        orbit_list = graphlets.ordered_orbit_list(orbit_is)
        orbit_lists[n_l] = orbit_list
        
        for net, name in zip(networks, net_names):
            
            interface.graphlet_degree_distributions(net,n,n_l,save_name='interface_'+name)            
            
            o_dir = directory + '/' + name + '_' + str(n_l)
            if not os.path.exists(o_dir):
                os.makedirs(o_dir)
            nodes = net.slices[0]
            if n_l == 1:
                agg_net = pymnet.MultiplexNetwork()
                for node in nodes:
                    agg_net.add_node(node)
                    
                for e in net.edges:
                    agg_net[e[0], e[1], 0] = 1
                    
                net = agg_net
                
            for layer_comb in itertools.combinations(net_layers, n_l):
                sub_net = pymnet.subnet(net, nodes, layer_comb)
                orbits = graphlets.orbit_counts_all(sub_net, n, nets, invs, auts, orbit_list, allowed_aspects=allowed_aspects)
                f_name = o_dir + '/' + name
                for layer in layer_comb:
                    f_name += "_" + str(layer)
                f_name += '.txt'
                write_orbit_counts(orbits, f_name, nodes, orbit_list)
        
    end = time.time()
    print(end - start)
    
    start = time.time()
    res_dir = directory + '/'
    all_gcds = {}
    nn_nls = [(1,4), (1,3), (2,4), (2,3), (3,3)]
    rs = ['', 'R']
    for nn_nl, r in itertools.product(nn_nls, rs):
        n_l, n = nn_nl
        if r == 'R':
            no_reds = True
        else:
            no_reds = False
        
        orbit_list = orbit_lists[n_l]
        orbit_is = orbit_is_d[n_l]
        
        if n_l == 1:
            net_layers = [0]
        elif allowed_aspects == [0]:
            net_layers = layers
        else:
            net_layers = list(range(n_l))
        gcds = data_analysis.GCDs(net_names, n, n_l, net_layers, res_dir, orbit_is, orbit_list, no_reds=no_reds, allowed_aspects=allowed_aspects)
        all_gcds[(n_l, n, r)] = gcds
        
    end = time.time()
    print(end - start)
    
    dist_name = 'GCD'
    fig = precision_recall_plot(all_gcds, boundaries, dist_name)
    for n_l, n, r in all_gcds:
        gcds = all_gcds[(n_l, n, r)]
        title = dist_name + '-' + str(n_l) + '-' + str(n) + r
        fig = MDS_plot(gcds, boundaries, labels, title)
        auprs = pairwise_auprs(gcds, boundaries, labels, title)
        fig = plot_AUPRs(auprs, labels=labels, title=title)
    
    
def example_networks(n_nets, n_n, n_l, m):
    '''
    Generates a test set of networks.
    
    Parameters
    ----------
    n_nets : int
        Number of networks generated from each model.
    n_n : int
        Number of nodes in each network
    n_l : int
        Number of layers in each network
    m : int
        Number of edges added to each new node in each layer in the BA model
        
    Returns
    -------
    networks : list of MultiplexNetworks
        Generated networks
    net_names : list of strs
        Names of the networks
    boundaries : list of ints
        The indices where a new network group begins. If there are three groups
        and five networks in each group, the boundaries are [5,10,15].
    labels : list of strs
        Group labels
    '''
    
    ms = [m] * n_l
    
    ba = []
    ba_plex = []
    conf = []
    conf_plex = []
    er_0 = []
    er_20 = []
    ws = []
    
    ba_names = []
    ba_plex_names = []
    conf_names = []
    conf_plex_names = []
    er_0_names = []
    er_20_names = []
    ws_names = []
    
    # ba
    for i in range(n_nets):
        net = ba_independent_multiplex(n_n, ms, couplings=None)
        ba.append(net)
        ba_names.append('ba_'+str(i))
    
    # ba plex
    for i in range(n_nets):
        net = pymnet.models.ba_total_degree(n_n, ms)
        ba_plex.append(net)
        ba_plex_names.append('ba_plex_' + str(i))
    
    # conf
    for i in range(n_nets):
        net = conf_independent_multiplex(ba[i])
        conf.append(net)
        conf_names.append('conf_'+str(i))
    
    # conf plex
    for i in range(n_nets):
        net = simple_conf_overlaps(ba_plex[i])
        conf_plex.append(net)
        conf_plex_names.append('conf_plex_'+str(i))
    
    # er 0 and er 20
    net0 = ba_plex[0]
    agg_net = pymnet.aggregate(net0, 1)
    n_e = int(round(len(agg_net.edges) / float(n_l)))
    ps0 = {}
    ps20 = {}
    layers = net0.slices[1]
    for nl in range(1, n_l):
        for layer_comb in itertools.combinations(layers, nl):
            ps0[layer_comb] = 0.0
            ps20[layer_comb] = 0.2
            
    for i in range(n_nets):
        net = pymnet.models.er_overlaps_match_aggregated(n_n, n_e, ps0)
        er_0.append(net)
        er_0_names.append('er_0_' + str(i))
        
        net = pymnet.models.er_overlaps_match_aggregated(n_n, n_e, ps20)
        er_20.append(net)
        er_20_names.append('er_20_' + str(i))
    
    # geo
    for i in range(n_nets):
        geo_edge_number = n_n*m # approximate number of edges
        # TODO: finish
    
    # ws: check what the number of edges should be! Now its made so that m = number of connected nearest neighbors parameter of ws
    for i in range(n_nets):
        ws.append(pymnet.models.ws(n_n,[(m/2.0)*n_n]*n_l))
        ws_names.append('ws_'+str(i))
    
    # still missing: conf plex, geo
    # conf degs should be based on ba and conf plex degs should be based on ba plex degs (?) -> normal ba needs to be implemented
    # what is the approximate edge density in geo?
    # what is the starting edge number in ws? The same m as in ba?
    # !!! what should the couplings be? !!!
    # are the er edges per layer (n_e) calculated correctly? n_e seems to be somehow average number of edges per layer in first ba plex net (except some edges get squished by aggregation so its less actually)
    # check edge numbers per layer: should these be the same?
    #for n in nets:
    #    print len([x for x in list(n.edges) if x[2] == x[3]])

    
    networks = ba + ba_plex + conf + conf_plex + er_0 + er_20 + ws
    net_names = ba_names + ba_plex_names + conf_names + conf_plex_names + er_0_names + er_20_names + ws_names
    boundaries = [n_nets, n_nets*2, n_nets*3, n_nets*4, n_nets*5, n_nets*6, n_nets*7]
    labels = ['BA', 'BA-plex', 'conf', 'conf-plex', 'ER$_{0,0}$', 'ER$_{20,20}$', 'WS']
    
    return networks, net_names, boundaries, labels


def ba_independent_multiplex(n, ms, couplings=None):
    net = pymnet.MultiplexNetwork(couplings=couplings)
    for ii in range(len(ms)):
        net.add_layer(ii)
        net.A[ii] = pymnet.nx.barabasi_albert_graph(n, ms[ii]) # apparently the seed network has no edges so average degs are not accurate
    return net

def conf_independent_multiplex(M):
    # apply configuration model to every layer of multiplex network M independently
    # uses 'bag of stubs' approach and multiedges and self-loops are simply not added at all
    M_conf = pymnet.MultiplexNetwork(couplings=M.couplings)
    for l in M.iter_layers():
        M_conf.add_layer(l)
        stubs = reduce(lambda bag,edge:bag+[edge[0]]+[edge[1]], M.A[l].edges, [])
        for sampled_edge in np.random.choice(stubs, (len(stubs)/2,2), replace=False):
            if sampled_edge[0] != sampled_edge[1]:
                M_conf[sampled_edge[0],l][sampled_edge[1],l] = 1
    return M_conf

def get_overlap_degs(net):
    # fixed version of the one in pymnet.diagnostics, function references working
    ol_degs = {}
    nodes = net.slices[0]
    layers = net.slices[1]
    
    net0 = pymnet.subnet(net, nodes, layers)
    
    for n_l in range(len(layers), 0, -1):
        for layer_comb in itertools.combinations(layers, n_l):
            sub_net = pymnet.subnet(net0, nodes, layer_comb)
            agg_net = pymnet.aggregate(sub_net, 1)
            thr_net = pymnet.transforms.threshold(agg_net, n_l)
            ol_degs[layer_comb] = pymnet.degs(thr_net, degstype='nodes')
            
            if n_l > 1:
                for e in thr_net.edges:
                    for layer in layer_comb:
                        net0[e[0], e[1], layer] = 0
                        
    return ol_degs

def simple_conf_overlaps(M):
    # simplified version of conf_overlaps
    # runs a configuration model for each type of edge overlap, but doesn't remove accidental extra overlaps, self-loops or multiedges
    ol_degs = get_overlap_degs(M)
    M_conf_overlap = pymnet.MultiplexNetwork(couplings=M.couplings)
    for layer_comb in ol_degs:
        for l in layer_comb:
            M_conf_overlap.add_layer(l)
        stubs = reduce(lambda bag,node:bag+[node]*ol_degs[layer_comb][node], ol_degs[layer_comb], [])
        for sampled_edge in np.random.choice(stubs, (len(stubs)/2,2), replace=False):
            if sampled_edge[0] != sampled_edge[1]:
                for ll in layer_comb:
                    M_conf_overlap[sampled_edge[0],ll][sampled_edge[1],ll] = 1
    return M_conf_overlap


def precision_recall_plot(all_dists, boundaries, dist_name=''):
    '''
    Plots the Precision-Recall curves and computes the AUPR values
    
    Parameters
    ----------
    all_dists : dict
        All distance matrices, keys specify the number of layers and nodes used,
        and whether redundant orbits are included (n_layers, n_nodes, r).
    boundaries : list of ints
        The indices where a new network group begins. If there are three groups
        and five networks in each group, the boundaries are [5,10,15].
    dist_name : str
        Name of the distance
        
    Returns
    -------
    pr_fig : matplotlib figure
        Figure with Precision-Recall curves
    '''
    
    c = {(1,4) : 'blue', (1,3) : 'cyan', (2,4) : 'purple', (2,3) : 'magenta', (3,3) : 'orange'}
    ls = {'' : '-', 'R' : '--'}
    types = [(1,4,''), (1,4,'R'), (1,3,''), (1,3,'R'), (2,4,''), (2,4,'R'), (2,3,''), (2,3,'R'), (3,3,''), (3,3,'R')]
    
    groups = group_labels(boundaries)
    
    pr_fig = plt.figure()
    ax = plt.subplot(111)
    for n_l, n, r in types:
        if (n_l, n, r) in all_dists:
            dists = all_dists[(n_l, n, r)]
            pres, recs = data_analysis.precision_recall(dists, groups)
            color = c[(n_l, n)]
            linestyle = ls[r]
            label = dist_name + '-' + str(n_l) + '-' + str(n) + r
            ax.plot(recs, pres, color=color, linestyle=linestyle, label=label, linewidth=1.5)
            aupr = data_analysis.area_under_precision_recall(pres, recs)
            print('AUPR ' + label + ': ' + str(aupr))
        
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.legend(loc=7, bbox_to_anchor=(1.42, 0.5))
    
    return pr_fig
    
    
def MDS_plot(dists, boundaries, labels, title=''):
    '''
    Embeds the networks in 3-dimensional space using MDS.
    
    Parameters
    ----------
    dists : list of lists, 2d-array
        Distance matrix
    boundaries : list of ints
        The indices where a new network group begins. If there are three groups
        and five networks in each group, the boundaries are [5,10,15].
    labels : list of strs
        Group labels
    title : str
        Title for the figure
        
    Returns
    -------
    mds_fig : matplotlib figure
    '''
    
    mds = manifold.MDS(n_components=3, dissimilarity='precomputed', max_iter=5000)
    colors = ['yellow', 'magenta', 'cyan', 'orange', 'purple', 'lightgreen', 'white', 'gray']
    
    Y = mds.fit_transform(dists)
    mds_fig = plt.figure()
    ax_mds = plt.subplot(111, projection='3d')
    ms = 75
    b0 = 0
    for b1, color in zip(boundaries, colors):
        ax_mds.scatter(Y[b0:b1, 0], Y[b0:b1, 1], Y[b0:b1, 2], c=color, s=ms)
        b0 = b1
    
    ax_mds.xaxis.set_major_formatter(NullFormatter())
    ax_mds.yaxis.set_major_formatter(NullFormatter())
    ax_mds.zaxis.set_major_formatter(NullFormatter())
    ax_mds.legend(labels, loc=7, bbox_to_anchor=(1.3, 0.5))
    plt.title(title)
    plt.axis('tight')
    
    return mds_fig
    
    
def pairwise_auprs(dists, boundaries, labels, title=''):
    '''
    Computes pairwise AUPRs between different groups.
    
    Parameters
    ----------
    dists : list of lists, 2d-array
        Distance matrix
    boundaries : list of ints
        The indices where a new network group begins. If there are three groups
        and five networks in each group, the boundaries are [5,10,15].
        
    Returns
    -------
    auprs : list of lists
        AUPR values between groups
    '''
    
    n_models = len(boundaries)
    boundaries = [0] + boundaries
    
    auprs = [None] * n_models
    for k in range(n_models):
        auprs[k] = [None] * n_models
        
    for i in range(1, n_models+1):
        b_i0 = boundaries[i-1]
        b_i1 = boundaries[i]
        n_i = b_i1 - b_i0
        for j in range(i, n_models+1):
            if j == i:
                auprs[i-1][j-1] = 1.0
                continue
            
            b_j0 = boundaries[j-1]
            b_j1 = boundaries[j]
            n_j = b_j1 - b_j0
            dists_sub = []
            for k in range(b_i0, b_i1) + range(b_j0, b_j1):
                sub = dists[k][b_i0:b_i1] + dists[k][b_j0:b_j1]
                dists_sub.append(sub)
                
            bounds = [n_i, n_i + n_j]
            groups = group_labels(bounds)
            pres, recs = data_analysis.precision_recall(dists_sub, groups)
            aupr = data_analysis.area_under_precision_recall(pres, recs)
            auprs[i-1][j-1] = aupr
            auprs[j-1][i-1] = aupr
    
    return auprs
    
    
def plot_AUPRs(auprs, labels, title=''):
    '''
    Visualizes pairwise AUPRs (or other measures in symmetric matrix form)
    
    Parameters
    ----------
    auprs : list of lists
        AUPR values between groups
    labels : list of strs
        Group labels
    title : str
        Title for the figure
        
    Returns
    -------
    fig : matplotlib figure
    '''
    
    auprs = auprs[1:]
    for i in range(len(auprs)):
        auprs[i] = auprs[i][:-1]
    mask = np.zeros_like(auprs)
    mask[np.triu_indices_from(mask)] = True
    for i in range(len(auprs)):
        mask[i][i] = False
    
    with sns.axes_style("white"):
        fig, ax = plt.subplots()
        sns.heatmap(auprs, ax=ax, xticklabels=labels[:-1], yticklabels=labels[1:], vmin=0.5, mask=mask)
        plt.title(title)
        
    return fig
    
    
def group_labels(boundaries):
    '''
    Gives group labels for network indices.
    
    Parameters
    ----------
    boundaries : list of ints
        The indices where a new network group begins. If there are three groups
        and five networks in each group, the boundaries are [5,10,15].
        
    Returns
    -------
    groups : dict
        Group labels, keys are network indices and value is the group label. Can
        be given as a parameter to precision_recall.
    '''
    
    groups = {}
    
    b0 = 0
    for k, b in enumerate(boundaries):
        for i in range(b0, b):
            groups[i] = k
            
        b0 = b
        
    return groups


    
    
def write_orbit_counts(orbits, file_name, nodes, orbit_list):
    '''
    Saves the orbit counts into a file.
    
    Parameters
    ----------
    orbits : dd (key: (node, orbit), value: count)
        Orbit counts for all the nodes
    file_name : str
        Name of the file where orbit counts are saved
    nodes : iterable
        Nodes for which the counts are computed / saved
    orbit_list : list of orbits
    '''
    
    file = open(file_name, 'w')
    line = "n"
    for orbit in orbit_list:
        line += "," + str(orbit)
    
    line += "\n"    
    file.write(line)
    
    for node in nodes:
        line = str(node)
        for orbit in orbit_list:
            line += "," + str(orbits[node, orbit])
                
        line += "\n"
        file.write(line)
        
    file.close()
    
    
def write_equations(n, n_l, allowed_aspects='all'):
    '''
    Writes the equations in LaTeX
    
    Parameters
    ----------
    n : int
        maximum number of nodes
    n_l : int
        Number of layers in the generated graphlets
    allowed_aspects : list, string
        the aspects that can be permutated when computing isomorphisms
    '''
    
    layers = list(range(n_l))
    nets, invs = graphlets.graphlets(n, layers, allowed_aspects=allowed_aspects)
    if n_l == 1:
        nets = visualization.order_nets(nets)
        invs = visualization.order_invs(invs)
    auts = graphlets.automorphism_orbits(nets, allowed_aspects)
    eqs = graphlets.orbit_equations(n, nets, auts, invs, allowed_aspects)
    
    subs = dd()
    for eq in eqs:
        if len(eq[0]) != 3:
            orbit1 = eq[0][0]
            orbit2 = eq[1][0]
            sub = graphlets.subtrahend(orbit1, orbit2, nets, auts, invs, allowed_aspects)
            subs[eq] = sub
                
    orbit_is = graphlets.orbit_numbers(n, nets, auts)
    if n_l == 1:
        orbit_is = visualization.order_orbit_is(orbit_is)
    
    for eq in eqs:
        eq_tex = ""
        if eq in subs:
            orbit1 = eq[0][0]
            o1 = orbit_is[orbit1]
            k1 = eq[0][1]
            orbit2 = eq[1][0]
            o2 = orbit_is[orbit2]
            k2 = eq[1][1]
            sub = subs[eq]
            eq_tex += "\\binom{C_{" + str(o1) + "}}{" + str(k1) + "} \\binom{C_{" + str(o2) + "}"
            if sub > 0:
                eq_tex += " - " + str(sub)
                
            eq_tex += "}{" + str(k2) + "}"
            
        else:
            orbit = eq[0]
            o = orbit_is[orbit]
            k = eq[1]
            eq_tex += "\\binom{C_{" + str(o) + "}}{" + str(k) + "}"
            
        eq_tex += " & = & "
        orbits = set()
        coefs = {}
        for orbit in eqs[eq]:
            o = orbit_is[orbit]
            orbits.add(o)
            coef = eqs[eq][orbit]
            coefs[o] = coef
                 
        while len(orbits) > 1:
            o = min(orbits)
            orbits.remove(o)
            coef = coefs[o]
            if coef > 1:
                eq_tex += str(coef) + " "
            
            eq_tex += "C_{" + str(o) + "} + "
            
        o = min(orbits)
        coef = coefs[o]
        if coef > 1:
            eq_tex += str(coef) + " "
        eq_tex += "C_{" + str(o) + "} "
            
        eq_tex += "\\\\ \n"
        print(eq_tex)
    
    
if __name__ == "__main__":
    main()
