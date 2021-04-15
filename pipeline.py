import data_analysis, visualization
import pymnet
from pymnet import graphlets
import os
import time
import itertools
import numpy as np
import random
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import manifold
from matplotlib.ticker import NullFormatter
from collections import defaultdict as dd
from scipy.stats import spearmanr
import interface
import graphlet_insertion
import cPickle

# try importing Dimitrova, Petrovski, Kocarev method
try:
    import Multiplex_Graph
    DPK_available = True
except ImportError:
    DPK_available = False


def main():
    
    n_nets = 10
    n_n = 1000
    n_l = 3
    m = 2
    use_simple_conf=False
    use_simple_conf_plex=True
    print_progress=True
    allowed_aspects = 'all'
    
    networks, net_names, boundaries, labels = example_networks(n_nets, n_n, n_l, m, use_simple_conf, use_simple_conf_plex, print_progress)
    
    directory = 'Results_3'
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
        if print_progress:
            print('Orbit list '+str(n_l)+' layers, '+str(n)+' nodes done')
        
        count = 0
        for net, name in zip(networks, net_names):
            
            # below: interface function for doing all of this, should work but maybe needs a test
            # TODO: test and use the interface instead to simplify this code
            #interface.graphlet_degree_distributions(net,n,n_l,save_name='interface_'+name)
            
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
            
            count += 1
            if print_progress:
                print('+'*count+'-'*(len(networks)-count))
        
        if print_progress:
            print('Orbit counts '+str(n_l)+' layers, '+str(n)+' nodes done')
    
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
        
        if print_progress:
            print('GCDs '+str(n_l)+' layers, '+str(n)+' nodes done')
    
    if DPK_available:
        all_gcds[('DPK','','')] = gcds_for_Dimitrova_Petrovski_Kocarev_method(networks)
        if print_progress:
            print('GCDs DPK method done')
            
    cPickle.dump(all_gcds,open('all_gcds_3.pickle','wb'))
    
    end = time.time()
    print(end - start)
    
    dist_name = 'GCD'
    fig,lgd = precision_recall_plot(all_gcds, boundaries, dist_name)
    fig.savefig('precision_recall_3.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')
    for n_l, n, r in all_gcds:
        pass
        #gcds = all_gcds[(n_l, n, r)]
        #title = dist_name + '-' + str(n_l) + '-' + str(n) + r
        #fig,lgd = MDS_plot(gcds, boundaries, labels, title)
        #fig.savefig('mds_1_'+title+'.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')
        #auprs = pairwise_auprs(gcds, boundaries, labels, title)
        #fig = plot_AUPRs(auprs, labels=labels, title=title)
        #fig.savefig('pairwise_auprs_1_'+title+'.pdf',bbox_inches='tight')

def make_networks(test_set_type='random',n_nets=10,n_n=1000,n_l=3,m=2,use_simple_conf=False,use_simple_conf_plex=True,allowed_aspects_graphlets='all',n_classes=5,n_different_graphlets=5,graphlet_frequency=0.05,graphlet_size=(4,2),print_progress=True):
    '''
    test_set_type == 'random':
        parameters used: n_nets, n_n, n_l, m, use_simple_conf, use_simple_conf_plex, print_progress
    test set type == 'random_deg_progression':
        parameters used: n_nets, n_n, n_l, m, use_simple_conf, use_simple_conf_plex, print_progress
    test_set_type == 'graphlet_insertion':
        parameters used: n_nets, n_n, n_l, m, n_classes, n_different_graphlets, graphlet_frequency, graphlet_size, print_progress
    See function example_networks
    '''
    networks, net_names, boundaries, labels = example_networks(test_set_type,n_nets,n_n,n_l,m,use_simple_conf,use_simple_conf_plex,allowed_aspects_graphlets,n_classes,n_different_graphlets,graphlet_frequency,graphlet_size,print_progress)
    netdir = 'Nets/'
    if not os.path.exists(netdir):
        os.makedirs(netdir)
    if test_set_type == 'random':
        file_prefix = '_'.join([str(n_nets),str(n_n),str(n_l),str(m),str(use_simple_conf),str(use_simple_conf_plex)])
    elif test_set_type == 'random_deg_progression':
        file_prefix = '_'.join(['deg_progression',str(n_nets),str(n_n),str(n_l),str(m).replace(' ',''),str(use_simple_conf),str(use_simple_conf_plex)])
    elif test_set_type == 'graphlet_insertion':
        file_prefix = '_'.join(['graphlet_insertion',str(n_nets),str(n_n),str(n_l),str(m),str(allowed_aspects_graphlets),str(n_classes),str(n_different_graphlets),str(graphlet_frequency),'-'.join([str(a) for a in graphlet_size])])
    with open(netdir+'networks_'+file_prefix+'.pickle','wb') as f:
        cPickle.dump(networks,f)
    with open(netdir+'netnames_'+file_prefix+'.pickle','wb') as g:
        cPickle.dump(net_names,g)
    with open(netdir+'boundaries_'+file_prefix+'.pickle','wb') as h:
        cPickle.dump(boundaries,h)
    with open(netdir+'labels_'+file_prefix+'.pickle','wb') as j:
        cPickle.dump(labels,j)
    return networks,net_names,boundaries,labels

def make_orbits(test_set_type='random',n_nets=10,n_n=1000,n_l=3,m=2,use_simple_conf=False,use_simple_conf_plex=True,allowed_aspects_graphlets='all',n_classes=5,n_different_graphlets=5,graphlet_frequency=0.05,graphlet_size=(4,2),allowed_aspects_orbits='all',print_progress=True):
    netdir = 'Nets/'
    # prefix depending on test_set_type
    if test_set_type == 'random':
        file_prefix = '_'.join([str(n_nets),str(n_n),str(n_l),str(m),str(use_simple_conf),str(use_simple_conf_plex)])
    elif test_set_type == 'random_deg_progression':
        file_prefix = '_'.join(['deg_progression',str(n_nets),str(n_n),str(n_l),str(m).replace(' ',''),str(use_simple_conf),str(use_simple_conf_plex)])
    elif test_set_type == 'graphlet_insertion':
        file_prefix = '_'.join(['graphlet_insertion',str(n_nets),str(n_n),str(n_l),str(m),str(allowed_aspects_graphlets),str(n_classes),str(n_different_graphlets),str(graphlet_frequency),'-'.join([str(a) for a in graphlet_size])])
    with open(netdir+'networks_'+file_prefix+'.pickle','rb') as f:
        networks = cPickle.load(f)
    with open(netdir+'netnames_'+file_prefix+'.pickle','rb') as g:
        net_names = cPickle.load(g)
    if print_progress:
        print('Nets loaded')
    orbit_dir = 'Orbits/'+file_prefix+'_'+str(allowed_aspects_orbits)+'/'
    if not os.path.exists(orbit_dir):
        os.makedirs(orbit_dir)
    orbit_aux_dir = 'Orbits_aux/'+file_prefix+'_'+str(allowed_aspects_orbits)+'/'
    if not os.path.exists(orbit_aux_dir):
        os.makedirs(orbit_aux_dir)
    layers = list(range(n_l))
    orbit_lists = {}
    orbit_is_d = {}
    start = time.time()
    nn_nl = [(1,4), (2,4), (3,3)]
    for n_l_orbit, n in nn_nl:
        if n_l_orbit == 1:
            net_layers = [0]
        else:
            net_layers = layers
        nets, invs = graphlets.graphlets(n, net_layers, n_l_orbit, allowed_aspects=allowed_aspects_orbits)
        auts = graphlets.automorphism_orbits(nets, allowed_aspects=allowed_aspects_orbits)
        orbit_is = graphlets.orbit_numbers(n, nets, auts)
        orbit_is_d[n_l_orbit] = orbit_is
        orbit_list = graphlets.ordered_orbit_list(orbit_is)
        orbit_lists[n_l_orbit] = orbit_list
        if print_progress:
            print('Orbit list '+str(n_l_orbit)+' layers, '+str(n)+' nodes done')
        count = 0
        for net, name in zip(networks, net_names):
            # below: interface function for doing all of this, should work but maybe needs a test
            # TODO: test and use the interface instead to simplify this code, fix save_name
            #interface.graphlet_degree_distributions(net,n,n_l_orbit,save_name='interface_'+name)
            o_dir = orbit_dir + '/' + name + '_' + str(n_l_orbit)
            if not os.path.exists(o_dir):
                os.makedirs(o_dir)
            nodes = net.slices[0]
            if n_l_orbit == 1:
                agg_net = pymnet.MultiplexNetwork()
                for node in nodes:
                    agg_net.add_node(node)
                for e in net.edges:
                    agg_net[e[0], e[1], 0] = 1
                net = agg_net
            for layer_comb in itertools.combinations(net_layers, n_l_orbit):
                f_name = o_dir + '/' + name
                for layer in layer_comb:
                    f_name += "_" + str(layer)
                f_name += '.txt'
                if not os.path.exists(f_name):
                    sub_net = pymnet.subnet(net, nodes, layer_comb)
                    orbits = graphlets.orbit_counts_all(sub_net, n, nets, invs, auts, orbit_list, allowed_aspects=allowed_aspects_orbits)
                    write_orbit_counts(orbits, f_name, nodes, orbit_list)
            count += 1
            if print_progress:
                print('+'*count+'-'*(len(networks)-count))
        if print_progress:
            print('Orbit counts '+str(n_l_orbit)+' layers, '+str(n)+' nodes done')
    with open(orbit_aux_dir+'orbit_is_d','wb') as f:
        cPickle.dump(orbit_is_d,f)
    with open(orbit_aux_dir+'orbit_lists','wb') as g:
        cPickle.dump(orbit_lists,g)
    end = time.time()
    print(end - start)

def make_gcds(test_set_type='random',n_nets=10,n_n=1000,n_l=3,m=2,use_simple_conf=False,use_simple_conf_plex=True,allowed_aspects_graphlets='all',n_classes=5,n_different_graphlets=5,graphlet_frequency=0.05,graphlet_size=(4,2),allowed_aspects_orbits='all',print_progress=True):
    start = time.time()
    if test_set_type == 'random':
        file_prefix = '_'.join([str(n_nets),str(n_n),str(n_l),str(m),str(use_simple_conf),str(use_simple_conf_plex)])
    elif test_set_type == 'random_deg_progression':
        file_prefix = '_'.join(['deg_progression',str(n_nets),str(n_n),str(n_l),str(m).replace(' ',''),str(use_simple_conf),str(use_simple_conf_plex)])
    elif test_set_type == 'graphlet_insertion':
        file_prefix = '_'.join(['graphlet_insertion',str(n_nets),str(n_n),str(n_l),str(m),str(allowed_aspects_graphlets),str(n_classes),str(n_different_graphlets),str(graphlet_frequency),'-'.join([str(a) for a in graphlet_size])])
    orbit_dir = 'Orbits/'+file_prefix+'_'+str(allowed_aspects_orbits)+'/'
    orbit_aux_dir = 'Orbits_aux/'+file_prefix+'_'+str(allowed_aspects_orbits)+'/'
    netdir = 'Nets/'
    with open(netdir+'netnames_'+file_prefix+'.pickle','rb') as f:
        net_names = cPickle.load(f)
    with open(orbit_aux_dir+'orbit_is_d','rb') as g:
        orbit_is_d = cPickle.load(g)
    with open(orbit_aux_dir+'orbit_lists','rb') as h:
        orbit_lists = cPickle.load(h)
    gcd_dir = 'GCDs/'
    gcd_aux_dir = 'GCDs_aux/'+file_prefix+'_'+str(allowed_aspects_orbits)+'/'
    if not os.path.exists(gcd_dir):
        os.makedirs(gcd_dir)
    if not os.path.exists(gcd_aux_dir):
        os.makedirs(gcd_aux_dir)
    layers = list(range(n_l))
    all_gcds = {}
    nn_nls = [(1,4), (1,3), (2,4), (2,3), (3,3)]
    rs = ['', 'R']
    for nn_nl, r in itertools.product(nn_nls, rs):
        n_l_orbit, n = nn_nl
        if r == 'R':
            no_reds = True
        else:
            no_reds = False
        orbit_list = orbit_lists[n_l_orbit]
        orbit_is = orbit_is_d[n_l_orbit]
        if n_l_orbit == 1:
            net_layers = [0]
        elif allowed_aspects_orbits == [0]:
            net_layers = layers
        else:
            net_layers = list(range(n_l_orbit))
        gcd_aux_savename = gcd_aux_dir+'_'.join(filter(None,[str(n_l_orbit),str(n),str(r)]))+'.pickle'
        if os.path.exists(gcd_aux_savename):
            with open(gcd_aux_savename,'rb') as aux_f:
                gcds = cPickle.load(aux_f)
        else:
            gcds = data_analysis.GCDs(net_names, n, n_l_orbit, net_layers, orbit_dir, orbit_is, orbit_list, no_reds=no_reds, allowed_aspects=allowed_aspects_orbits)
            with open(gcd_aux_savename,'wb') as aux_g:
                cPickle.dump(gcds,aux_g)
        all_gcds[(n_l_orbit, n, r)] = gcds
        if print_progress:
            print('GCDs '+str(n_l_orbit)+' layers, '+str(n)+' nodes done')
    if DPK_available:
        gcd_aux_savename = gcd_aux_dir+'DPK'+'.pickle'
        dpk_aux_dir = gcd_aux_dir+'DPK_GCMs/'
        if os.path.exists(gcd_aux_savename):
            with open(gcd_aux_savename,'rb') as aux_h:
                gcds = cPickle.load(aux_h)
        else:
            with open(netdir+'networks_'+file_prefix+'.pickle','rb') as j:
                networks = cPickle.load(j)
            gcds = gcds_for_Dimitrova_Petrovski_Kocarev_method(networks,dpk_aux_dir)
            with open(gcd_aux_savename,'wb') as aux_j:
                cPickle.dump(gcds,aux_j)
        all_gcds[('DPK','','')] = gcds
        if print_progress:
            print('GCDs DPK method done')
    with open(gcd_dir+'all_gcds_'+file_prefix+'_'+str(allowed_aspects_orbits)+'.pickle','wb') as k:
        cPickle.dump(all_gcds,k)
    end = time.time()
    print(end - start)

def make_figures(test_set_type='random',n_nets=10,n_n=1000,n_l=3,m=2,use_simple_conf=False,use_simple_conf_plex=True,allowed_aspects_graphlets='all',n_classes=5,n_different_graphlets=5,graphlet_frequency=0.05,graphlet_size=(4,2),allowed_aspects_orbits='all',print_progress=True):
    netdir = 'Nets/'
    if test_set_type == 'random':
        file_prefix = '_'.join([str(n_nets),str(n_n),str(n_l),str(m),str(use_simple_conf),str(use_simple_conf_plex)])
    elif test_set_type == 'random_deg_progression':
        file_prefix = '_'.join(['deg_progression',str(n_nets),str(n_n),str(n_l),str(m).replace(' ',''),str(use_simple_conf),str(use_simple_conf_plex)])
    elif test_set_type == 'graphlet_insertion':
        file_prefix = '_'.join(['graphlet_insertion',str(n_nets),str(n_n),str(n_l),str(m),str(allowed_aspects_graphlets),str(n_classes),str(n_different_graphlets),str(graphlet_frequency),'-'.join([str(a) for a in graphlet_size])])
    with open(netdir+'boundaries_'+file_prefix+'.pickle','rb') as f:
        boundaries = cPickle.load(f)
    with open(netdir+'labels_'+file_prefix+'.pickle','rb') as g:
        labels = cPickle.load(g)
    gcd_dir = 'GCDs/'
    with open(gcd_dir+'all_gcds_'+file_prefix+'_'+str(allowed_aspects_orbits)+'.pickle','rb') as h:
        all_gcds = cPickle.load(h)
    fig_dir = 'Figures/'+file_prefix+'_'+str(allowed_aspects_orbits)+'/'
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    dist_name = 'GCD'
    fig,lgd = precision_recall_plot(all_gcds, boundaries, dist_name, fig_dir+dist_name+'_AUPRs.txt')
    fig.savefig(fig_dir+'precision_recall.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')
    plt.close(fig)
    fig_mds_combined = plt.figure(figsize=(15,9))
    fig_auprs_combined,axs_auprs_combined = plt.subplots(nrows=2,ncols=3,sharey='row',sharex='col',figsize=(15,9))
    fig_auprs_combined.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.02, hspace=0.02)
    auprs_combined_cbar_ax = fig_auprs_combined.add_axes([0.83, 0.1, 0.02, 0.8])
    mds_combined_index = {(1,3):1,(2,3):2,(3,3):3,(1,4):4,(2,4):5,('DPK',''):6}
    auprs_combined_index = {(1,3):(0,0),(2,3):(0,1),(3,3):(0,2),(1,4):(1,0),(2,4):(1,1),('DPK',''):(1,2)}
    for n_l, n, r in all_gcds:
        gcds = all_gcds[(n_l, n, r)]
        title = dist_name + '-' + str(n_l) + '-' + str(n) + r
        if r == '':
            mds_subax = fig_mds_combined.add_subplot(2,3,mds_combined_index[(n_l,n)],projection='3d')
            if mds_combined_index[(n_l,n)] == 5:
                mds_legend_subax = mds_subax
            fig,lgd = MDS_plot(gcds, boundaries, labels, title, additional_ax=mds_subax)
        else:
            fig,lgd = MDS_plot(gcds, boundaries, labels, title)
        fig.savefig(fig_dir+'mds_'+title+'.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')
        plt.close(fig)
        auprs = pairwise_auprs(gcds, boundaries, labels, title)
        if r == '':
            add_ax = axs_auprs_combined[auprs_combined_index[(n_l,n)][0],auprs_combined_index[(n_l,n)][1]]
            fig = plot_AUPRs(auprs, labels=replace_labels(labels), title=title, additional_ax=add_ax, additional_cbar_ax=auprs_combined_cbar_ax)
        else:
            fig = plot_AUPRs(auprs, labels=replace_labels(labels), title=title)
        fig.savefig(fig_dir+'pairwise_auprs_'+title+'.pdf',bbox_inches='tight')
        plt.close(fig)
    lgd_mds = mds_legend_subax.legend(replace_labels(labels), loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=8, fontsize=15, handletextpad=0.05)
    fig_mds_combined.subplots_adjust(hspace=-0.05,wspace=-0.1)
    fig_mds_combined.savefig(fig_dir+'mds_combined.pdf',bbox_extra_artists=(lgd_mds,),bbox_inches='tight')
    plt.close(fig_mds_combined)
    auprs_combined_cbar_ax.tick_params(axis='both', which='both', length=0, labelsize=15)
    fig_auprs_combined.savefig(fig_dir+'pairwise_auprs_combined.pdf',bbox_inches='tight')
    plt.close(fig_auprs_combined)

def replace_labels(labels):
    new_label_dict = {'BA':'BA-ind','BA-plex':'BA-dep','conf':'Conf-ind','conf-plex':'Conf-dep','ER$_{0,0}$':'ER-0','ER$_{20,20}$':'ER-20','geo':'GEO','WS':'WS','graphlet_set1':'1','graphlet_set2':'2','graphlet_set3':'3','graphlet_set4':'4','graphlet_set5':'5'}
    return [new_label_dict.get(s, s) for s in labels]



def gcds_for_Dimitrova_Petrovski_Kocarev_method(networks,aux_folder):
    gcm_file_names = []
    if not os.path.exists(aux_folder):
        os.makedirs(aux_folder)
    for ii,M in enumerate(networks):
        gcm_file_name = aux_folder+str(ii)+'.pickle'
        if os.path.exists(gcm_file_name):
            gcm_file_names.append(gcm_file_name)
        else:
            assert isinstance(M,pymnet.MultiplexNetwork)
            graph_edges = [(e[0],e[1],e[2]) for e in M.edges if e[2]==e[3]]
            graphF = Multiplex_Graph.GraphFunc(doDirected=False)
            graphF.make_graph_table(graph_edges)
            graphF.make_direct_neighbours(subOne=False)
            graphF.make_zero_orbit()
            graphF.count_tri_graphs()
            orb_mat = graphF.return_orbits_Mat().values
            orb_mat_with_dummy = np.row_stack((orb_mat,[1]*orb_mat.shape[1]))
            gcm = spearmanr(orb_mat_with_dummy)[0]
            with open(gcm_file_name,'wb') as f:
                cPickle.dump(gcm,f)
            gcm_file_names.append(gcm_file_name)
            del gcm
    return GCD_matrix_save_memory(gcm_file_names)

def GCD_matrix_save_memory(gcm_file_names):
    gcds = []
    for gcm_file_1 in gcm_file_names:
        gcds_t = []
        with open(gcm_file_1,'rb') as f:
            gcm1 = cPickle.load(f)
        for gcm_file_2 in gcm_file_names:
            with open(gcm_file_2,'rb') as g:
                gcm2 = cPickle.load(g)
            gcd = graphlets.GCD(gcm1,gcm2)
            gcds_t.append(gcd)
            del gcm2
        gcds.append(gcds_t)
        del gcm1
    return gcds

def example_networks(test_set_type, n_nets, n_n, n_l, m,
                     use_simple_conf=False, use_simple_conf_plex=False,
                     allowed_aspects='all', n_classes=5, n_different_graphlets=5, graphlet_frequency=0.05, graphlet_size=(4,2),
                     print_progress=False):
    '''
    Generates a test set of networks.
    
    Parameters
    ----------
    test_set_type : 'random' or 'random_deg_progression' or 'graphlet_insertion'
        Determines whether random network models (static or progressing avg degs) or the graphlet insertion method is used in creating the test set
    n_nets : int
        Number of networks generated from each model / number of networks in each graphlet insertion set
        Used in all test set types
    n_n : int
        Number of nodes in each network
        Used in all test set types
    n_l : int
        Number of layers in each network
        Used in all test set types
    m : int or list
        Number of edges added to each new node in each layer in the BA model, int (test set type random)
        Same as above except a list defining the progression of the parameter, list (test set type random_deg_progression)
        n_n*m is the number of edges per layer in the baseline ER model, int (test set type graphlet_insertion)
        Used in all test set types
    use_simple_conf : bool
        Is the simplified conf net model used
        Used only in test set type random
    use_simple_conf_plex : bool
        Is the simplified conf plex net model used
        Used only in test set type random
    allowed_aspects : 'all' or list of aspects
        Allowed permutable aspects for determining the set of all possible graphlets
        Used only in test set type graphlet insertion
    n_classes : int
        Number of classes of nets with different inserted graphlets
        Used only in test set type graphlet insertion
    n_different_graphlets : int
        Number of inserted graphlets
        Used only in test set type graphlet insertion
    graphlet_frequency : float
        Proportion of inserted graphlets out of all possible locations
        Used only in test set type graphlet insertion
    graphlet_size : tuple
        Size of inserted graphlets (nnodes,nlayers)
        Used only in test set type graphlet insertion
        
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
    
    
    if test_set_type == 'random':
        
        ms = [m] * n_l
        
        ba = []
        ba_plex = []
        conf = []
        conf_plex = []
        er_0 = []
        er_20 = []
        geo = []
        ws = []
        
        ba_names = []
        ba_plex_names = []
        conf_names = []
        conf_plex_names = []
        er_0_names = []
        er_20_names = []
        geo_names = []
        ws_names = []
        
        # ba
        for i in range(n_nets):
            net = ba_independent_multiplex(n_n, ms, couplings=None)
            ba.append(net)
            ba_names.append('ba_'+str(i))
        if print_progress:
            print('ba done')
        
        # ba plex
        for i in range(n_nets):
            net = pymnet.models.ba_total_degree(n_n, ms)
            ba_plex.append(net)
            ba_plex_names.append('ba_plex_' + str(i))
        if print_progress:
            print('ba plex done')
        
        # conf
        for i in range(n_nets):
            if use_simple_conf:
                net = conf_independent_multiplex(ba[i])
            else:
                net = pymnet.models.conf(ba[i],couplings=None)
                # a small hack to properly set couplings
                net.couplings = [('none',)]
            conf.append(net)
            conf_names.append('conf_'+str(i))
        if print_progress:
            print('conf done')
        
        # conf plex
        for i in range(n_nets):
            if use_simple_conf_plex:
                net = simple_conf_overlaps(ba_plex[i])
            else:
                net = pymnet.models.conf_overlaps(get_overlap_degs(ba_plex[i]))
            conf_plex.append(net)
            conf_plex_names.append('conf_plex_'+str(i))
        if print_progress:
            print('conf plex done')
        
        # er 0 and er 20, aggregated number of edges should be n_n*m*n_l,
        # so the edge parameter should be n_n*m
        n_e_0 = n_n*m
        n_e_20 = n_n*m
        ps0 = {}
        ps20 = {}
        layers = list(range(n_l))
        for nl in range(1, n_l):
            for layer_comb in itertools.combinations(layers, nl):
                ps0[layer_comb] = 0.0
                ps20[layer_comb] = 0.2
        for i in range(n_nets):
            net = pymnet.models.er_overlaps_match_aggregated(n_n, n_e_0, ps0)
            er_0.append(net)
            er_0_names.append('er_0_' + str(i))
            
            net = pymnet.models.er_overlaps_match_aggregated(n_n, n_e_20, ps20)
            er_20.append(net)
            er_20_names.append('er_20_' + str(i))
        if print_progress:
            print('er 0 and er 20 done')
        
        # geo
        for i in range(n_nets):
            geo_edge_number = n_n*m # approximate number of edges, works better for large n_n
            net = pymnet.models.geo(n_n, [geo_edge_number]*n_l)
            geo.append(net)
            geo_names.append('geo_'+str(i))
        if print_progress:
            print('geo done')
        
        # ws
        for i in range(n_nets):
            ws_edge_number = n_n*m # number of edges
            net = pymnet.models.ws(n_n,[ws_edge_number]*n_l,p=0.3)
            ws.append(net)
            ws_names.append('ws_'+str(i))
        if print_progress:
            print('ws done')
        
        # conf degs should be based on ba and conf plex degs should be based on ba plex degs (?) -> normal ba needs to be implemented
        # what is the approximate edge density in geo?
        # what is the starting edge number in ws? The same m as in ba?
        # !!! what should the couplings be? !!!
        # are the er edges per layer (n_e) calculated correctly? n_e seems to be somehow average number of edges per layer in first ba plex net (except some edges get squished by aggregation so its less actually)
        # check edge numbers per layer: should these be the same?
        #for n in nets:
        #    print len([x for x in list(n.edges) if x[2] == x[3]])
    
        
        networks = ba + ba_plex + conf + conf_plex + er_0 + er_20 + geo + ws
        net_names = ba_names + ba_plex_names + conf_names + conf_plex_names + er_0_names + er_20_names + geo_names + ws_names
        boundaries = [n_nets, n_nets*2, n_nets*3, n_nets*4, n_nets*5, n_nets*6, n_nets*7, n_nets*8]
        labels = ['BA', 'BA-plex', 'conf', 'conf-plex', 'ER$_{0,0}$', 'ER$_{20,20}$', 'geo', 'WS']
        
    elif test_set_type == 'graphlet_insertion':
        
        networks = []
        net_names = []
        boundaries = []
        labels = []
        for ii in range(n_classes):
            curr_nets = []
            curr_net_names = []
            for jj in range(n_nets):
                base_net = pymnet.models.er(n_n,edges=[n_n*m]*n_l)
                base_net.couplings = [('none',)]
                curr_nets.append(base_net)
                curr_net_names.append('set'+str(ii)+'_'+str(jj))
            graphlet_insertion.insert_random_graphlets(curr_nets,graphlet_size[0],graphlet_size[1],n_different_graphlets,[graphlet_frequency]*n_different_graphlets,allowed_aspects,balance=True,print_progress=print_progress)
            networks.extend(curr_nets)
            net_names.extend(curr_net_names)
            boundaries.append(n_nets*(ii+1))
            labels.append('graphlet_set'+str(ii+1))
            if print_progress:
                print(labels[-1]+' done')
    
    elif test_set_type == 'random_deg_progression':
        n_diff_degs = len(m)
        assert n_nets % n_diff_degs == 0
        n_per_deg = n_nets // n_diff_degs
        all_res = list()
        for m_instance in m:
            all_res.append(example_networks(test_set_type='random',n_nets=n_per_deg,n_n=n_n,n_l=n_l,m=m_instance,use_simple_conf=use_simple_conf,use_simple_conf_plex=use_simple_conf_plex,print_progress=print_progress))
        boundaries_0 = all_res[0][2]
        labels_0 = all_res[0][3]
        networks = []
        net_names = []
        for ii in range(len(labels_0)):
            for jj in range(n_diff_degs):
                if ii == 0:
                    previous = 0
                else:
                    previous = boundaries_0[ii-1]
                for kk in range(previous,boundaries_0[ii]):
                    networks.append(all_res[jj][0][kk])
                    curr_name = all_res[jj][1][kk].split('_')
                    number = str(int(curr_name[-1])+jj*n_per_deg)
                    curr_name[-1] = number
                    net_names.append('_'.join(curr_name))
        boundaries = [b*n_diff_degs for b in boundaries_0]
        labels = labels_0
    
    return networks, net_names, boundaries, labels


def ba_independent_multiplex(n, ms, couplings=None):
    net = pymnet.MultiplexNetwork(couplings=couplings)
    for ii in range(len(ms)):
        net.add_layer(ii)
        #net.A[ii] = pymnet.nx.barabasi_albert_graph(n, ms[ii]) # apparently the seed network has no edges so average degs are not accurate
        net.A[ii] = ba_single_layer(n, ms[ii])
    # independent random shuffle nodes on each layer to remove iteration order effect
    net = pymnet.transforms.randomize_nodes_by_layer(net)
    return net

def ba_single_layer(n,m):
    # seed network is full net with m+1 nodes
    # the number of edges is thus n*m-0.5*m-0.5*m**2
    # the nodes are iterated in order, so smaller node numbers get larger degs
    # this must be taken into account if made into multiplex
    net = pymnet.full(nodes=m+1,layers=None)
    rng = random.Random()
    stubs = list(range(m+1))*m
    for node in range(m+1,n):
        sample = sample_stubs(stubs,m,rng)
        for target in sample:
            net[node,target] = 1
        stubs.extend([node]*m)
        stubs.extend(sample)
    return net

def sample_stubs(stubs,m,rng):
    sample = set()
    while len(sample) < m:
        sample.add(rng.choice(stubs))
    return sample

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
        if stubs:
            for sampled_edge in np.random.choice(stubs, (len(stubs)/2,2), replace=False):
                if sampled_edge[0] != sampled_edge[1]:
                    for ll in layer_comb:
                        M_conf_overlap[sampled_edge[0],ll][sampled_edge[1],ll] = 1
    return M_conf_overlap


def precision_recall_plot(all_dists, boundaries, dist_name='', AUPR_writefilename=None):
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
    
    if DPK_available:
        c = {(1,4) : 'blue', (1,3) : 'cyan', (2,4) : 'purple', (2,3) : 'magenta', (3,3) : 'orange', ('DPK','') : 'green'}
        ls = {'' : '-', 'R' : '--'}
        types = [(1,4,''), (1,4,'R'), (1,3,''), (1,3,'R'), (2,4,''), (2,4,'R'), (2,3,''), (2,3,'R'), (3,3,''), (3,3,'R'), ('DPK','','')]
    else:
        c = {(1,4) : 'blue', (1,3) : 'cyan', (2,4) : 'purple', (2,3) : 'magenta', (3,3) : 'orange'}
        ls = {'' : '-', 'R' : '--'}
        types = [(1,4,''), (1,4,'R'), (1,3,''), (1,3,'R'), (2,4,''), (2,4,'R'), (2,3,''), (2,3,'R'), (3,3,''), (3,3,'R')]
    
    groups = group_labels(boundaries)
    
    if AUPR_writefilename != None:
        f = open(AUPR_writefilename,'w')
    
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
            if AUPR_writefilename != None:
                f.write('AUPR ' + label + ': ' + str(aupr)+'\n')
    if AUPR_writefilename != None:
        f.close()
        
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    lgd = ax.legend(loc=7, bbox_to_anchor=(1.42, 0.5))
    
    return pr_fig,lgd
    
    
def MDS_plot(dists, boundaries, labels, title='', additional_ax=None):
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
    #colors = ['yellow', 'magenta', 'cyan', 'orange', 'purple', 'lightgreen', 'white', 'gray']
    #colors = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494','#b3b3b3']
    #colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf']
    colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#cccc29','#a65628','#f781bf']
    markers = ['^','v','s','D','P','X','o','*']
    
    Y = mds.fit_transform(dists)
    mds_fig = plt.figure()
    ax_mds = plt.subplot(111, projection='3d')
    ms = 75
    b0 = 0
    for b1, color, marker in zip(boundaries, colors, markers):
        ax_mds.scatter(Y[b0:b1, 0], Y[b0:b1, 1], Y[b0:b1, 2], c=color, s=ms, marker=marker)
        if additional_ax:
            additional_ax.scatter(Y[b0:b1, 0], Y[b0:b1, 1], Y[b0:b1, 2], c=color, s=ms, marker=marker)
        b0 = b1
    if additional_ax:
        additional_ax.xaxis.set_major_formatter(NullFormatter())
        additional_ax.yaxis.set_major_formatter(NullFormatter())
        additional_ax.zaxis.set_major_formatter(NullFormatter())
    ax_mds.xaxis.set_major_formatter(NullFormatter())
    ax_mds.yaxis.set_major_formatter(NullFormatter())
    ax_mds.zaxis.set_major_formatter(NullFormatter())
    lgd = ax_mds.legend(labels, loc=7, bbox_to_anchor=(1.3, 0.5))
    plt.title(title)
    plt.axis('tight')
    
    return mds_fig,lgd
    
    
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
    
    
def plot_AUPRs(auprs, labels, title='', additional_ax=None, additional_cbar_ax=None):
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
        if additional_ax is not None and additional_cbar_ax is not None:
            sns.heatmap(auprs, ax=additional_ax, xticklabels=labels[:-1], yticklabels=labels[1:], vmin=0.5, mask=mask, cbar_ax=additional_cbar_ax)
            additional_ax.tick_params(axis='both', which='both', length=0, labelsize=15, rotation=0)
            additional_ax.tick_params(axis='x', rotation=45)
        sns.heatmap(auprs, ax=ax, xticklabels=labels[:-1], yticklabels=labels[1:], vmin=0.5, mask=mask)
        plt.tick_params(axis='both', which='both', length=0, labelsize=15, rotation=0)
        plt.tick_params(axis='x', rotation=45)
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
