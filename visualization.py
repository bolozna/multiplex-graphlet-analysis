import pymnet
import itertools
import matplotlib.pyplot as plt
import matplotlib.patches
from collections import defaultdict as dd
from pymnet import graphlets
#import main as functions
#import graphlets

def main():
    
    name = 'graphlets/graphlet_2_5_'
    end = '.txt'
    #for i in range(5286):
    #i = 4781
    i = 3389
    net = pymnet.MultilayerNetwork(aspects = 1)
    file_name = name + str(i) + end
    file = open(file_name, 'r')
    for row in file:
        row = row.rstrip()
        row = row.split()
        net[int(row[0]), int(row[1]), int(row[2]), int(row[3])] = 1
    
    #mono_iso_example()
    #multi_layers2_nodes2_part()
    #layers1_nodes4()
    #layers3_nodes3()
    #simple_example()
    layers2_nodes4_a()
    #multi_layers2_nodes2()
    #equation_intermediate()
    #orbit_x3_example()
    '''
    figure = plt.figure()
    k = 0
    #figure, axes = plt.subplots(3, subplot_kw=dict(projection='3d'))
    for i in nets:
        for j in range(len(nets[i])):
            #figure.add_subplot(2,2,k)
            #ax = axes[k]#plt.gca(projection='3d')
            ax = figure.add_axes([(k//2) * 0.5,0.5-(k % 2) * 0.5,0.5,0.5], projection='3d')
            net = nets[i][j]
            nodeColors, nodeLabels = node_colors_and_labels(net, i, j, auts, orbit_is)
            pymnet.draw(net, ax=ax, azim=0, elev=50, camera_dist=7.5, autoscale=False, layout='spectral', defaultLayerAlpha = 0.7, layergap = 0.15, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels)
            k += 1
            
    figure.savefig("fig_test.png")
    '''
    
    nodeCoords = {1 : (0.6, 0.01), 4 : (0.1, 0.9), 3 : (0.9, 0.3)} #{'4' : (0.3, 0.9), '0' : (0.4, 0.5), '3' : (0.8, 0.1)}
    #nodeColors = {('4', '0') : '#fa5418', '0' : 'black', '1' : 'black', '2' : 'black', '3' : 'black'}
    #nodeColors = {(0,0) : '#99ff23', (1,0) : '#fdff23', (2,0) : '#ff9933', (3,0) : '#ff2387', (4,0) : '#5723ff'}
    nodeColorRule = {'colormap' : 'plasma', 'rule' : 'order', 'sequence' : list(range(10)), 'scaleby' : 30}
    
    '''
    figure = plt.figure()
    
    for j in range(len(nets[5])):
        net5 = nets[5][j]
        nodeColors = node_colors(net5, auts, 5, j)
        #figure.add_subplot(2,2,j+1)
        #plt.plot([1,2,3], [4,5,6])
        #ax = plt.gca()
    
    figure = pymnet.draw(net, layout = 'spectral', defaultLayerAlpha = 0.5, layergap = 1.5, nodeCoords = nodeCoords, nodeColorDict = nodeColors)#nodeColorRule = nodeColorRule)
        #fig.axes[0].set_figure(figure)
        #figure.add_subplot(fig.axes[0])
        
    plt.show()
    figure.savefig("fig_test.png")
    '''
    net1 = pymnet.MultilayerNetwork(aspects=1, directed=True)
    net1[1,2,'a']=1
    
    iso = pymnet.get_isomorphism(net, net, backend = "bliss")
    
    #print(aut)


def multi_automorphism_orbits(nets, allowed_aspects='all'):
    '''
    note: not general function, only for 2 nodes 2 layers
    '''
    
    multi_auts = dd()
    for i in nets:
        for j in range(len(nets[i])):
            net = nets[i][j]
            aut_gen = pymnet.get_automorphism_generators(net, allowed_aspects=allowed_aspects)
            n_gen = len(aut_gen)
            if n_gen == 0: # kaikki eri orbiteilla
                for (node, layer) in net.iter_node_layers():
                    multi_auts[i, j, node, layer] = (node, layer)
                        
            elif n_gen == 2: # kaikki samalla orbitilla
                aut = None
                for (node, layer) in net.iter_node_layers():
                    if not aut:
                        aut = (node, layer)
                    
                    multi_auts[i, j, node, layer] = aut
                        
            else:
                n_gen_n = len(aut_gen[0][0])
                n_gen_l = len(aut_gen[0][1])
                if (n_gen_n > 0) and (n_gen_l > 0): # swap both nodes and layers
                    auts = set()
                    for (node, layer) in net.iter_node_layers():
                        aut = None
                        for a in auts:
                            if (a[0] != node) and (a[1] != layer):
                                aut = a
                                break
                            
                        if not aut:
                            aut = (node, layer)
                            auts.add(aut)
                            
                        multi_auts[i, j, node, layer] = aut
                            
                elif n_gen_n > 0: # swap only nodes
                    auts = set()
                    for (node, layer) in net.iter_node_layers():
                        aut = None
                        for a in auts:
                            if a[0] != node:
                                aut = a
                                break
                            
                        if not aut:
                            aut = (node, layer)
                            auts.add(aut)
                            
                        multi_auts[i, j, node, layer] = aut
                            
                else: # swap only layers
                    auts = set()
                    for (node, layer) in net.iter_node_layers():
                        aut = None
                        for a in auts:
                            if a[0] == node:
                                aut = a
                                break
                            
                        if not aut:
                            aut = (node, layer)
                            auts.add(aut)
                            
                        multi_auts[i, j, node, layer] = aut
                            
    return multi_auts
    
    
def multi_orbit_numbers(n, nets, multi_auts):
    
    multi_orbit_is = dd()
    for i in range(2, n+1):
        for j in range(len(nets[i])):
            net = nets[i][j]
            for (node, layer) in net.iter_node_layers():
                aut = multi_auts[i, j, node, layer]
                if not (i, j, aut) in multi_orbit_is:
                    multi_orbit_is[(i, j, aut)] = len(multi_orbit_is)
                        
    return multi_orbit_is
    
    
def multi_node_colors_and_labels(net, i, j, multi_auts, multi_orbit_is, color_ids):
    
    nodeColors = {}
    nodeLabels = {}
    colors = ['#fdff23', '#ff2387', '#5723ff', '#00a324', '#68a300', '#23ff9b', '#ffc264', '#00cc0e', '#37ec00', '#ff9933', '#99ff23']
    used_orbits = set()
    for (node, layer) in net.iter_node_layers():
        aut = multi_auts[i, j, node, layer]
        color = colors[color_ids[i, j, aut]]
        nodeColors[(node, layer)] = color
        orbit_i = multi_orbit_is[i, j, aut]
        if orbit_i in used_orbits:
            nodeLabels[(node, layer)] = ''
        else:
            nodeLabels[(node, layer)] = str(orbit_i)
            used_orbits.add(orbit_i)
                
    return nodeColors, nodeLabels
    
    
def node_colors_and_labels(net, net_i, net_j, auts, orbit_is, color_ids):
    
    nodeColors = {}
    nodeLabels = {}
    colors = ['#fdff23', '#ff2387', '#5723ff', '#00a324', '#37ec00', '#ff9933', '#99ff23'] # dippa l2_n4
    #colors = ['#fdff23', '#ff2387', '#5723ff', '#00d200', '#37ec00', '#ff9933', '#99ff23']
    #colors = ['#c166c1', '#5b005b', '#eaccea', '#990099'] #'#2d002d'
    #colors = ['#ff2387', '#5723ff', '#00a324', '#ffa500']
    #colors = ['#984ea3', '#ff7f00', '#4daf4a', '#377eb8']
    #colors = ['#990099', '#6cbf09', '#ff2387', '#5723ff'] #'#00a324'
    #colors = ['#ff5a36', '#000a96', '#e7c0f2', '#6cbf09'] #'#ff541b' '#f96959'
    #colors = ['#6cbf09', '#cfc5d6', '#f96959', '#000a96']
    #colors = ['#6cbf09', '#990099', '#fedf00', '#cfc5d6'] # colors for mono5
    #colors = ['#990099', '#cfc5d6', '#ffa500', '#6cbf09', '#fedf00'] # colors for 3l 3n
    nodes = net.slices[0]
    layers = net.slices[1]
    used_orbits = set()
    for node in nodes:
        aut = auts[net_i, net_j, node]
        color = colors[color_ids[net_i, net_j, aut]]
        for layer in layers:
            nodeColors[(node, layer)] = color #colors[aut]
            orbit_i = orbit_is[(net_i, net_j, aut)]
            if orbit_i in used_orbits:
                nodeLabels[(node, layer)] = ''
            else:
                nodeLabels[(node, layer)] = str(orbit_i)
                used_orbits.add(orbit_i)
            
    return nodeColors, nodeLabels
    
    
def node_colors_and_labels_nl(net, i, j, auts, orbit_is, color_ids):
    
    nodeColors = {}
    nodeLabels = {}
    colors = ['#fdff23', '#ff2387', '#5723ff', '#00a324', '#ffa500', '#990099', '#ffffff', '#732e1d']
    
    used_orbits = set()
    for nl in net.iter_node_layers():
        aut = auts[i, j, nl]
        color = colors[color_ids[i, j, aut]]
        nodeColors[nl] = color
        orbit_i = orbit_is[(i, j, aut)]
        if orbit_i in used_orbits:
            nodeLabels[nl] = ''
        else:
            nodeLabels[nl] = str(orbit_i)
            used_orbits.add(orbit_i)
            
    return nodeColors, nodeLabels
    
    
def orbit_numbers(n, nets, auts):
    
    orbit_is = dd()
    for k in range(2, n+1):
        for j in range(len(nets[k])):
            net = nets[k][j]
            for node in net.slices[0]:
                aut = auts[(k,j,node)]
                if not (k,j,aut) in orbit_is:
                    orbit_is[(k,j,aut)] = len(orbit_is)
            
    return orbit_is
    
    
def orbit_numbers_nl(n, nets, auts):
    
    orbit_is = dd()
    for k in range(2, n+1):
        for j in range(len(nets[k])):
            net = nets[k][j]
            for nl in net.iter_node_layers():
                aut = auts[(k,j,nl)]
                if not (k,j,aut) in orbit_is:
                    orbit_is[(k,j,aut)] = len(orbit_is)
                    
    return orbit_is
    
    
def node_color_ids(n, nets, auts):
    
    color_ids = dd()
    for i in range(2, n+1):
        for j in range(len(nets[i])):
            ids = set()
            net = nets[i][j]
            for (node, layer) in net.iter_node_layers():
                if (i,j,node) in auts:
                    aut = auts[(i,j,node)]
                else:
                    aut = auts[i,j,node,layer]
                    
                if not (i,j,aut) in color_ids:
                    color_ids[(i,j,aut)] = len(ids)
                    ids.add((i,j,aut))
                    
    return color_ids
    
    
def node_color_ids_nl(n, nets, auts):
    
    color_ids = dd()
    for i in range(2, n+1):
        for j in range(len(nets[i])):
            ids = set()
            net = nets[i][j]
            for nl in net.iter_node_layers():
                aut = auts[i,j,nl]
                if not (i,j,aut) in color_ids:
                    color_ids[(i,j,aut)] = len(ids)
                    ids.add((i,j,aut))
                    
    return color_ids
                
    
    
def co_node_colors_and_labels(node1, both_orbit_nodes, net, i, j, nets, auts, orbit_is):
    
    nodeColors = {}
    nodeLabels = {}
    colors = ['#ff2387', '#5723ff', '#990099']
    nodes = net.slices[0]
    layers = net.slices[1]
    used_orbits = set()
    for node in nodes:
        for layer in layers:
            if node == node1:
                iso = pymnet.get_isomorphism(net, nets[i][j])
                if node in iso[0]:
                    iso_node = iso[0][node]
                else:
                    iso_node = node
                aut = auts[(i,j,iso_node)]
                orbit_i = orbit_is[(i,j,aut)]
                if orbit_i in used_orbits:
                    nodeLabels[(node, layer)] = ''
                else:
                    nodeLabels[(node, layer)] = str(orbit_i)
                    used_orbits.add(orbit_i)
            elif node in both_orbit_nodes:
                nodeColors[(node, layer)] = colors[2]
                nodeLabels[(node, layer)] = ''
            elif int(node) >= 0:
                nodeColors[(node, layer)] = colors[1]
                nodeLabels[(node, layer)] = ''
            else:
                nodeColors[(node, layer)] = colors[0]
                nodeLabels[(node, layer)] = ''
                    
    return nodeColors, nodeLabels, orbit_i
    
    
def edge_colors_and_widths(net, both_orbit_nodes):
    
    edgeColors = {}
    edgeWidths = {}
    edges = net.edges
    for edge in edges:
        n1 = edge[0]
        n2 = edge[1]
        l1 = edge[2]
        l2 = edge[3]
        if not (n1 in both_orbit_nodes or n2 in both_orbit_nodes) and ((n1 < 0 and n2 > 0) or (n1 > 0 and n2 < 0)):
            edgeColors[((n1, l1), (n2, l2))] = 'black'
            edgeWidths[((n1, l1), (n2, l2))] = 4
            
    return edgeColors, edgeWidths
    
    
def edge_colors_and_widths_sub(net, both_orbit_nodes):
    
    colors = ['#ff2387', '#5723ff', '#990099']
    
    edgeColors = {}
    edgeWidths = {}
    edges = net.edges
    for edge in edges:
        n1 = edge[0]
        n2 = edge[1]
        l1 = edge[2]
        l2 = edge[3]
        if n1 in both_orbit_nodes and n2 in both_orbit_nodes:
            edgeColors[((n1, l1), (n2, l2))] = colors[2]
            edgeWidths[((n1, l1), (n2, l2))] = 2
        elif n1 in both_orbit_nodes:
            if n2 > 0:
                edgeColors[((n1, l1), (n2, l2))] = colors[1]
                edgeWidths[((n1, l1), (n2, l2))] = 2
            else:
                edgeColors[((n1, l1), (n2, l2))] = colors[0]
                edgeWidths[((n1, l1), (n2, l2))] = 2
        elif n1 < 0:
            if n2 < 0 or n2 in both_orbit_nodes:
                edgeColors[((n1, l1), (n2, l2))] = colors[0]
                edgeWidths[((n1, l1), (n2, l2))] = 2
        else:
            if n2 > 0 or n2 in both_orbit_nodes:
                edgeColors[((n1, l1), (n2, l2))] = colors[1]
                edgeWidths[((n1, l1), (n2, l2))] = 2
                
    return edgeColors, edgeWidths
    
    
def edge_colors_plex(net):
    
    colors = ['#c166c1', '#5b005b', '#eaccea', '#990099']
    colors = ['#ff2387', '#5723ff', '#40e0d0', '#00a324']
    
    edgeColors = {}
    edges = net.edges
    for edge in edges:
        n1 = edge[0]
        n2 = edge[1]
        l1 = edge[2]
        l2 = edge[3]
        if l1 == l2:
            edgeColors[((n1, l1), (n2, l2))] = colors[l1]
            
    return edgeColors
    
    
def layer_labels(layers):
    
    layerLabels = {}
    for layer in layers:
        layerLabels[layer] = ''
        
    return layerLabels
    

def order_nets(nets):
    
    nets[4][0], nets[4][1] = nets[4][1], nets[4][0]
    nets[4][2], nets[4][3] = nets[4][3], nets[4][2]
    if 5 in nets:
        nets[5][0], nets[5][2] = nets[5][2], nets[5][0]
        nets[5][0], nets[5][5] = nets[5][5], nets[5][0]
        nets[5][0], nets[5][11] = nets[5][11], nets[5][0]
        nets[5][0], nets[5][10] = nets[5][10], nets[5][0]
        nets[5][0], nets[5][6] = nets[5][6], nets[5][0]
        nets[5][0], nets[5][13] = nets[5][13], nets[5][0]
        nets[5][0], nets[5][15] = nets[5][15], nets[5][0]
        nets[5][0], nets[5][14] = nets[5][14], nets[5][0]
        nets[5][0], nets[5][9] = nets[5][9], nets[5][0]
        nets[5][0], nets[5][4] = nets[5][4], nets[5][0]
        nets[5][0], nets[5][8] = nets[5][8], nets[5][0]
        nets[5][0], nets[5][3] = nets[5][3], nets[5][0]
        nets[5][0], nets[5][7] = nets[5][7], nets[5][0]
    
    return nets
    
    
def order_invs(invs):
    
    for i in invs:
        if invs[i] == (4, 0):
            invs[i] = (4, 1)
        elif invs[i] == (4, 1):
            invs[i] = (4, 0)
        elif invs[i] == (4, 2):
            invs[i] = (4, 3)
        elif invs[i] == (4, 3):
            invs[i] = (4, 2)
        elif invs[i] == (5, 0):
            invs[i] = (5, 2)
        elif invs[i] == (5, 2):
            invs[i] = (5, 5)
        elif invs[i] == (5, 5):
            invs[i] = (5, 11)
        elif invs[i] == (5, 11):
            invs[i] = (5, 10)
        elif invs[i] == (5, 10):
            invs[i] = (5, 6)
        elif invs[i] == (5, 6):
            invs[i] = (5, 13)
        elif invs[i] == (5, 13):
            invs[i] = (5, 15)
        elif invs[i] == (5, 15):
            invs[i] = (5, 14)
        elif invs[i] == (5, 14):
            invs[i] = (5, 9)
        elif invs[i] == (5, 9):
            invs[i] = (5, 4)
        elif invs[i] == (5, 4):
            invs[i] = (5, 8)
        elif invs[i] == (5, 8):
            invs[i] = (5, 3)
        elif invs[i] == (5, 3):
            invs[i] = (5, 7)
        elif invs[i] == (5, 7):
            invs[i] = (5, 0)
    
    return invs
    
    
def order_orbit_is(orbit_is):
    
    for i in orbit_is:
        if orbit_is[i] == 1:
            orbit_is[i] = 2
        elif orbit_is[i] == 2:
            orbit_is[i] = 1
        elif orbit_is[i] == 4:
            orbit_is[i] = 5
        elif orbit_is[i] == 5:
            orbit_is[i] = 4
        elif orbit_is[i] == 6:
            orbit_is[i] = 7
        elif orbit_is[i] == 7:
            orbit_is[i] = 6
        elif orbit_is[i] == 9:
            orbit_is[i] = 11
        elif orbit_is[i] == 11:
            orbit_is[i] = 9
        elif orbit_is[i] == 12:
            orbit_is[i] = 13
        elif orbit_is[i] == 13:
            orbit_is[i] = 12
        elif orbit_is[i] == 15:
            orbit_is[i] = 17
        elif orbit_is[i] == 17:
            orbit_is[i] = 15
        elif orbit_is[i] == 18:
            orbit_is[i] = 21
        elif orbit_is[i] == 21:
            orbit_is[i] = 18
        elif orbit_is[i] == 19:
            orbit_is[i] = 20
        elif orbit_is[i] == 20:
            orbit_is[i] = 19
        elif orbit_is[i] == 22:
            orbit_is[i] = 23
        elif orbit_is[i] == 23:
            orbit_is[i] = 22
        elif orbit_is[i] == 24:
            orbit_is[i] = 26
        elif orbit_is[i] == 25:
            orbit_is[i] = 24
        elif orbit_is[i] == 26:
            orbit_is[i] = 25
        elif orbit_is[i] == 27:
            orbit_is[i] = 30
        elif orbit_is[i] == 30:
            orbit_is[i] = 27
        elif orbit_is[i] == 31:
            orbit_is[i] = 33
        elif orbit_is[i] == 33:
            orbit_is[i] = 31
        elif orbit_is[i] == 35:
            orbit_is[i] = 38
        elif orbit_is[i] == 36:
            orbit_is[i] = 37
        elif orbit_is[i] == 37:
            orbit_is[i] = 35
        elif orbit_is[i] == 38:
            orbit_is[i] = 36
        elif orbit_is[i] == 39:
            orbit_is[i] = 42
        elif orbit_is[i] == 41:
            orbit_is[i] = 39
        elif orbit_is[i] == 42:
            orbit_is[i] = 41
        elif orbit_is[i] == 43:
            orbit_is[i] = 44
        elif orbit_is[i] == 44:
            orbit_is[i] = 43
        elif orbit_is[i] == 45:
            orbit_is[i] = 48
        elif orbit_is[i] == 48:
            orbit_is[i] = 45
        elif orbit_is[i] == 46:
            orbit_is[i] = 47
        elif orbit_is[i] == 47:
            orbit_is[i] = 46
        elif orbit_is[i] == 49:
            orbit_is[i] = 50
        elif orbit_is[i] == 50:
            orbit_is[i] = 49
        elif orbit_is[i] == 51:
            orbit_is[i] = 53
        elif orbit_is[i] == 52:
            orbit_is[i] = 51
        elif orbit_is[i] == 53:
            orbit_is[i] = 52
        elif orbit_is[i] == 54:
            orbit_is[i] = 55
        elif orbit_is[i] == 55:
            orbit_is[i] = 54
        elif orbit_is[i] == 56:
            orbit_is[i] = 58
        elif orbit_is[i] == 58:
            orbit_is[i] = 56
        elif orbit_is[i] == 59:
            orbit_is[i] = 60
        elif orbit_is[i] == 60:
            orbit_is[i] = 59
        elif orbit_is[i] == 62:
            orbit_is[i] = 63
        elif orbit_is[i] == 63:
            orbit_is[i] = 64
        elif orbit_is[i] == 64:
            orbit_is[i] = 62
        elif orbit_is[i] == 65:
            orbit_is[i] = 67
        elif orbit_is[i] == 67:
            orbit_is[i] = 65
        elif orbit_is[i] == 70:
            orbit_is[i] = 71
        elif orbit_is[i] == 71:
            orbit_is[i] = 70
    
    return orbit_is
    
    
def node_coords(k):
    
    if k == 1:
        nodeCoords = {0 : (0.9, 0.1), 1 : (0.1, 0.9)}
        
    elif k == 2:
        nodeCoords = {0 : (0.5, 0.5), 1 : (0.9, 0.1), 2 : (0.1, 0.9)}
        
    elif k == 3:
        nodeCoords = {0 : (0.8, 0.5), 1 : (0.2, 0.1), 2 : (0.2, 0.9)}
        
    elif k == 4:
        nodeCoords = {0 : (0.65, 0.35), 1 : (0.35, 0.65), 2 : (0.9, 0.1), 3 : (0.1, 0.9)}
        
    elif k == 5:
        nodeCoords = {0 : (0.55, 0.5), 1 : (0.1, 0.5), 2 : (0.9, 0.1), 3 : (0.9, 0.9)}
        
    elif k == 6:
        nodeCoords = {1 : (0.8, 0.2), 0 : (0.8, 0.8), 3 : (0.2, 0.2), 2 : (0.2, 0.8)}
        
    elif k == 7:
        nodeCoords = {0 : (0.5, 0.5), 3 : (0.1, 0.3), 1 : (0.1, 0.7), 2 : (0.9, 0.5)}
        
    elif k == 8:
        nodeCoords = {0 : (0.9, 0.5), 3 : (0.1, 0.5), 1 : (0.5, 0.1), 2 : (0.5, 0.9)}
        
    elif k == 9:
        nodeCoords = {1 : (0.9, 0.1), 0 : (0.9, 0.9), 2 : (0.1, 0.5), 3 : (0.6, 0.5)}
        
    elif k == 10:
        nodeCoords = {0 : (0.5, 0.5), 1 : (0.7, 0.3), 2 : (0.3, 0.7), 3 : (0.9, 0.1), 4 : (0.1, 0.9)}
        
    elif k == 11:
        nodeCoords = {1 : (0.65, 0.5), 0 : (0.4, 0.5), 3 : (0.1, 0.4), 2 : (0.1, 0.6), 4 : (0.9, 0.5)}
        
    elif k == 12:
        nodeCoords = {0 : (0.5, 0.5), 1 : (0.8, 0.2), 2 : (0.2, 0.8), 3 : (0.8, 0.8), 4 : (0.2, 0.2)}
        
    elif k == 13:
        nodeCoords = {4 : (0.1, 0.5), 1 : (0.5, 0.3), 0 : (0.5, 0.7), 3 : (0.9, 0.1), 2 : (0.9, 0.9)}
        
    elif k == 14:
        nodeCoords = {1 : (0.65, 0.5), 0 : (0.4, 0.5), 4 : (0.1, 0.4), 2 : (0.1, 0.6), 3 : (0.9, 0.5)}
        
    elif k == 15:
        nodeCoords = {0 : (0.5, 0.5), 1 : (0.3, 0.9), 2 : (0.3, 0.1), 3 : (0.7, 0.1), 4 : (0.7, 0.9)}
        
    elif k == 16:
        nodeCoords = {0 : (0.1, 0.5), 1 : (0.45, 0.9), 2 : (0.45, 0.1), 4 : (0.9, 0.25), 3 : (0.9, 0.75)}
        
    elif k == 17:
        nodeCoords = {0 : (0.6, 0.5), 3 : (0.9, 0.5), 2 : (0.35, 0.3), 1 : (0.35, 0.7), 4 : (0.1, 0.5)}
        
    elif k == 18:
        nodeCoords = {0 : (0.6, 0.5), 3 : (0.9, 0.5), 2 : (0.35, 0.3), 1 : (0.35, 0.7), 4 : (0.1, 0.5)}
        
    elif k == 19:
        nodeCoords = {0 : (0.5, 0.5), 4 : (0.3, 0.1), 2 : (0.7, 0.1), 3 : (0.3, 0.9), 1 : (0.7, 0.9)}
        
    elif k == 20:
        nodeCoords = {1 : (0.6, 0.5), 3 : (0.9, 0.5), 4 : (0.35, 0.3), 0 : (0.35, 0.7), 2 : (0.1, 0.5)}
        
    elif k == 21:
        nodeCoords = {0 : (0.9, 0.5), 1 : (0.5, 0.1), 2 : (0.5, 0.5), 3 : (0.5, 0.9), 4 : (0.1, 0.5)}
        
    elif k == 22:
        nodeCoords = {0 : (0.45, 0.65), 1 : (0.9, 0.65), 2 : (0.1, 0.5), 4 : (0.45, 0.35), 3 : (0.9, 0.35)}
        
    elif k == 23:
        nodeCoords = {2 : (0.6, 0.5), 3 : (0.9, 0.5), 4 : (0.35, 0.3), 0 : (0.35, 0.7), 1 : (0.1, 0.5)}
        
    elif k == 24:
        nodeCoords = {0 : (0.65, 0.7), 2 : (0.9, 0.7), 4 : (0.4, 0.5), 3 : (0.65, 0.3), 1 : (0.1, 0.5)}
        
    elif k == 25:
        nodeCoords = {4 : (0.5, 0.4), 3 : (0.9, 0.4), 1 : (0.7, 0.6), 2 : (0.1, 0.4), 0 : (0.3, 0.6)}
        
    elif k == 26:
        nodeCoords = {0 : (0.9, 0.5), 3 : (0.5, 0.1), 1 : (0.5, 0.5), 2 : (0.5, 0.9), 4 : (0.1, 0.5)}
        
    elif k == 27:
        nodeCoords = {3 : (0.6, 0.5), 1 : (0.9, 0.5), 4 : (0.35, 0.3), 0 : (0.35, 0.7), 2 : (0.1, 0.5)}
        
    elif k == 28:
        nodeCoords = {4 : (0.5, 0.5), 3 : (0.8, 0.2), 0 : (0.2, 0.8), 1 : (0.8, 0.8), 2 : (0.2, 0.2)}
        
    elif k == 29:
        nodeCoords = {3 : (0.6, 0.5), 1 : (0.9, 0.5), 4 : (0.1, 0.2), 0 : (0.1, 0.8), 2 : (0.35, 0.5)}
        
    elif k == 30:
        nodeCoords = {0 : (0.1, 0.5), 1 : (0.45, 0.9), 2 : (0.45, 0.1), 4 : (0.9, 0.25), 3 : (0.9, 0.75)}
        
    return nodeCoords
    
    
def layers1_nodes5():
    
    n = 5
    n_layers = 1
    layers = list(range(n_layers))
    nets, invs = graphlets.graphlets(n, layers)
    nets = order_nets(nets)
    
    auts = graphlets.automorphism_orbits(nets)
    orbit_is = orbit_numbers(n, nets, auts)
    orbit_is = order_orbit_is(orbit_is)
    
    color_ids = node_color_ids(n, nets, auts)
    layerLabels = layer_labels(layers)
    figure = plt.figure(figsize=(15,18))
    k = 1
    for i in nets:
        for j in range(len(nets[i])):
            nodeCoords = node_coords(k)
            ax = figure.add_subplot(6,5,k, projection='3d')
            net = nets[i][j]
            nodeColors, nodeLabels = node_colors_and_labels(net, i, j, auts, orbit_is, color_ids)
            pymnet.draw(net, ax=ax, azim=0, elev=90, camera_dist=6, autoscale=False, defaultNodeLabelSize=16, layout='spring', defaultLayerAlpha = 0.0, layergap = 0.15, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels, nodeSizeRule={'scalecoeff': 0.3, 'rule': 'scaled'}, nodeCoords=nodeCoords)
            plt.title(r'$\mathbf{G_{' + str(k-1) + '}}$', fontsize=24, ha='left', position=(0,1), fontweight='heavy')
            k += 1
            
    figure.savefig("figs/mono5.pdf", bbox_inches='tight')


def layers1_nodes4():
    
    n = 4
    n_layers = 1
    layers = list(range(n_layers))
    nets, invs = graphlets.graphlets(n, layers)
    nets = order_nets(nets)
    
    auts = graphlets.automorphism_orbits(nets)
    orbit_is = orbit_numbers(n, nets, auts)
    orbit_is = order_orbit_is(orbit_is)
    
    color_ids = node_color_ids(n, nets, auts)
    layerLabels = layer_labels(layers)
    figure = plt.figure(figsize=(9,9))
    k = 1
    for i in nets:
        for j in range(len(nets[i])):
            nodeCoords = node_coords(k)
            ax = figure.add_subplot(3,3,k, projection='3d')
            net = nets[i][j]
            nodeColors, nodeLabels = node_colors_and_labels(net, i, j, auts, orbit_is, color_ids)
            pymnet.draw(net, ax=ax, azim=0, elev=90, camera_dist=6, autoscale=False, defaultNodeLabelSize=16, layout='spring', defaultLayerAlpha = 0.0, layergap = 0.15, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels, nodeSizeRule={'scalecoeff': 0.3, 'rule': 'scaled'}, nodeCoords=nodeCoords)
            plt.title(r'$\mathbf{G_{' + str(k-1) + '}}$', fontsize=24, ha='left', position=(0,1), fontweight='heavy')
            k += 1
            
    figure.savefig("figs/mono4.pdf", bbox_inches='tight')
    
    
def layers2_nodes3():
    
    n = 3
    n_layers = 2
    layers = list(range(n_layers))
    nets, invs = graphlets.graphlets(n, layers)
    auts = graphlets.automorphism_orbits(nets)
    orbit_is = orbit_numbers(n, nets, auts)
    color_ids = node_color_ids(n, nets, auts)
    layerLabels = layer_labels(layers)
    figure = plt.figure(figsize=(10,8))
    k = 1
    for i in nets:
        for j in range(len(nets[i])):
            ax = figure.add_subplot(4,3,k, projection='3d')
            k += 1
            net = nets[i][j]
            nodeColors, nodeLabels = node_colors_and_labels(net, i, j, auts, orbit_is, color_ids)
            pymnet.draw(net, ax=ax, elev=60, camera_dist=9.5, layout='circular', layergap = 2, defaultNodeLabelSize=16, defaultLayerAlpha = 0.3, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels)
            
    figure.savefig("figs/l2_n3.pdf", bbox_inches='tight')


def layers2_nodes3_nl():
    
    n = 3
    n_layers = 2
    layers = list(range(n_layers))
    nets, invs = graphlets.graphlets(n, layers)
    auts = graphlets.automorphism_orbits_nl(nets)
    orbit_is = orbit_numbers_nl(n, nets, auts)
    color_ids = node_color_ids_nl(n, nets, auts)
    layerLabels = layer_labels(layers)
    figure = plt.figure(figsize=(10,8))
    k = 1
    for i in nets:
        for j in range(len(nets[i])):
            ax = figure.add_subplot(4,3,k, projection='3d')
            k += 1
            net = nets[i][j]
            nodeColors, nodeLabels = node_colors_and_labels_nl(net, i, j, auts, orbit_is, color_ids)
            pymnet.draw(net, ax=ax, elev=60, camera_dist=9.5, layout='circular', layergap = 2, defaultNodeLabelSize=16, defaultLayerAlpha = 0.3, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels)
            
    figure.savefig("figs/l2_n3_nl.pdf", bbox_inches='tight')
    
    
def layers2_nodes3_a():
    
    n = 3
    n_layers = 2
    layers = list(range(n_layers))
    nets, invs = graphlets.graphlets(n, layers, couplings='categorical')
    auts = graphlets.automorphism_orbits(nets)
    orbit_is = orbit_numbers(n, nets, auts)
    color_ids = node_color_ids(n, nets, auts)
    layerLabels = layer_labels(layers)
    figure = plt.figure(figsize=(12,4))
    k = 1
    for i in nets:
        for j in range(len(nets[i])):
            ax = figure.add_subplot(2,6,k, projection='3d')
            k += 1
            net = nets[i][j]
            nodeColors, nodeLabels = node_colors_and_labels(net, i, j, auts, orbit_is, color_ids)
            pymnet.draw(net, ax=ax, elev=40, camera_dist=9.0, layout='circular', layergap = 1.0, defaultNodeLabelSize=16, defaultLayerAlpha = 0.3, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels, layerPadding=0.1)
            
    figure.savefig("/u/26/sallmes1/unix/Documents/Article/TheArticle/figs/l2_n3.pdf", bbox_inches='tight')


def layers2_nodes4_a():
    
    n = 4
    n_layers = 2
    layers = list(range(n_layers))
    nets, invs = graphlets.graphlets(n, layers, couplings='categorical')
    auts = graphlets.automorphism_orbits(nets)
    orbit_is = orbit_numbers(n, nets, auts)
    color_ids = node_color_ids(n, nets, auts)
    layerLabels = layer_labels(layers)
    figure = plt.figure(figsize=(10,12))
    k = 1
    i = 4
    m = 1
    #for i in nets:
    for j in range(len(nets[i])):
        if k > 48:
            k = 1
            f_name = "/u/26/sallmes1/unix/Documents/Article/TheArticle/figs/l2_n4_" + str(m) + ".pdf"
            figure.savefig(f_name, bbox_inches='tight')
            figure = plt.figure(figsize=(10,12))
            m +=1
        ax = figure.add_subplot(8,6,k, projection='3d')
        k += 1
        net = nets[i][j]
        nodeColors, nodeLabels = node_colors_and_labels(net, i, j, auts, orbit_is, color_ids)
        pymnet.draw(net, ax=ax, elev=50, azim=-20, camera_dist=9.0, layout='circular', layerPadding=0.2, layergap = 1.5, defaultNodeLabelAlpha=0.0, defaultNodeLabelSize=16, defaultLayerAlpha = 0.3, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels, nodeSizeRule={'scalecoeff': 0.25, 'rule': 'scaled'})
    
    f_name = "/u/26/sallmes1/unix/Documents/Article/TheArticle/figs/l2_n4_" + str(m) + ".pdf"    
    figure.savefig(f_name, bbox_inches='tight')
    print(k)
    
    
def layers2_nodes4():
    
    n = 4
    n_layers = 2
    layers = list(range(n_layers))
    nets, invs = graphlets.graphlets(n, layers)
    auts = graphlets.automorphism_orbits(nets)
    orbit_is = orbit_numbers(n, nets, auts)
    color_ids = node_color_ids(n, nets, auts)
    layerLabels = layer_labels(layers)
    figure = plt.figure(figsize=(10,14))
    k = 1
    i = 4
    m = 1
    #for i in nets:
    for j in range(len(nets[i])):
        if k > 35:
            k = 1
            f_name = "figs/l2_n4_" + str(m) + ".pdf"
            figure.savefig(f_name, bbox_inches='tight')
            figure = plt.figure(figsize=(10,14))
            m +=1
        ax = figure.add_subplot(7,5,k, projection='3d')
        k += 1
        net = nets[i][j]
        nodeColors, nodeLabels = node_colors_and_labels(net, i, j, auts, orbit_is, color_ids)
        pymnet.draw(net, ax=ax, elev=50, azim=-20, camera_dist=9.0, layout='circular', layerPadding=0.2, layergap = 1.5, defaultNodeLabelAlpha=0.0, defaultNodeLabelSize=16, defaultLayerAlpha = 0.3, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels)
    
    f_name = "figs/l2_n4_" + str(m) + ".pdf"    
    figure.savefig(f_name, bbox_inches='tight')
    print(k)
    

def layers2_nodes4_nl():
    
    n = 4
    n_layers = 2
    layers = list(range(n_layers))
    nets, invs = graphlets.graphlets(n, layers)
    auts = graphlets.automorphism_orbits_nl(nets)
    orbit_is = orbit_numbers_nl(n, nets, auts)
    color_ids = node_color_ids_nl(n, nets, auts)
    layerLabels = layer_labels(layers)
    figure = plt.figure(figsize=(10,14))
    k = 1
    i = 4
    m = 1
    #for i in nets:
    for j in range(len(nets[i])):
        if k > 35:
            k = 1
            f_name = "figs/l2_n4_nl" + str(m) + ".pdf"
            figure.savefig(f_name, bbox_inches='tight')
            figure = plt.figure(figsize=(10,14))
            m +=1
        ax = figure.add_subplot(7,5,k, projection='3d')
        k += 1
        net = nets[i][j]
        nodeColors, nodeLabels = node_colors_and_labels_nl(net, i, j, auts, orbit_is, color_ids)
        pymnet.draw(net, ax=ax, elev=50, azim=-20, camera_dist=9.0, layout='circular', layerPadding=0.2, layergap = 1.5, defaultNodeLabelAlpha=0.0, defaultNodeLabelSize=16, defaultLayerAlpha = 0.3, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels)
    
    f_name = "figs/l2_n4_nl" + str(m) + ".pdf"    
    figure.savefig(f_name, bbox_inches='tight')
    
    
def layers3_nodes3():
    
    n = 3
    n_layers = 3
    layers = list(range(n_layers))
    nets, invs = graphlets.graphlets(n, layers)
    auts = graphlets.automorphism_orbits(nets)
    orbit_is = orbit_numbers(n, nets, auts)
    color_ids = node_color_ids(n, nets, auts)
    layerLabels = layer_labels(layers)
    nodeSizes = {(0, 0) : 0.3, (1, 0) : 0.3, (2, 0) : 0.3, (0, 1) : 0.0, (1, 1) : 0.0, (2, 1) : 0.0, (0, 2) : 0.0, (1, 2) : 0.0, (2, 2) : 0.0}
    figure = plt.figure(figsize=(10,14))
    k = 1
    for i in nets:
        for j in range(len(nets[i])):
            ax = figure.add_subplot(7,5,k, projection='3d')
            k += 1
            net = nets[i][j]
            nodeColors, nodeLabels = node_colors_and_labels(net, i, j, auts, orbit_is, color_ids)
            edgeColors = edge_colors_plex(net)
            pymnet.draw(net, ax=ax, azim=-60, elev=80, camera_dist=10, layout='circular', layergap = 0.8, defaultNodeLabelSize=13, defaultLayerAlpha = 0.0, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels, edgeColorDict=edgeColors, nodeSizeDict=nodeSizes)
            plt.title(r'$\mathbf{P^3_{' + str(k-1) + '}}$', fontsize=20, ha='left', position=(0,1), fontweight='heavy')
            
    figure.savefig("figs/l3_n3.pdf", bbox_inches='tight')
    
   
def layers2_nodes2():
    
    n = 2
    n_layers = 2
    layers = list(range(n_layers))
    layerLabels = layer_labels(layers)
    layers = [set(layers)]
    nets, invs = graphlets_general(n, layers, 1)
    nets[2][0].add_layer(1)
    auts = graphlets.automorphism_orbits(nets, [1])
    orbit_is = orbit_numbers(n, nets, auts)
    color_ids = node_color_ids(n, nets, auts)
    figure = plt.figure(figsize=(12,16))
    k = 1
    for i in nets:
        for j in range(len(nets[i])):
            ax = figure.add_subplot(7,5,k, projection='3d')
            k += 1
            net = nets[i][j]
            nodeColors, nodeLabels = node_colors_and_labels(net, i, j, auts, orbit_is, color_ids)
            pymnet.draw(net, ax=ax, layout='circular', camera_dist=8, layergap = 0.7, defaultLayerAlpha = 0.3, edgeStyleRule={'inter': '-', 'intra': '-', 'rule': 'edgetype'}, defaultNodeLabelSize=16, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels)
            
    figure.savefig("figs/l2_n2_layer.pdf", bbox_inches='tight')
    
    
def multi_layers2_nodes2():
    
    n = 2
    n_layers = 2
    layers = list(range(n_layers))
    layerLabels = layer_labels(layers)
    layers = [set(layers)]
    nets, invs = graphlets_general(n, layers, 1)
    nets[2][0].add_layer(1)
    multi_auts = multi_automorphism_orbits(nets, allowed_aspects=[1])
    multi_orbit_is = multi_orbit_numbers(n, nets, multi_auts)
    color_ids = node_color_ids(n, nets, multi_auts)
    figure = plt.figure(figsize=(12,16))
    k = 1
    for i in nets:
        for j in range(len(nets[i])):
            ax = figure.add_subplot(7,5,k, projection='3d')
            k += 1
            net = nets[i][j]
            nodeColors, nodeLabels = multi_node_colors_and_labels(net, i, j, multi_auts, multi_orbit_is, color_ids)
            pymnet.draw(net, ax=ax, layout='circular', camera_dist=8, layergap = 0.7, defaultLayerAlpha = 0.3, layerPadding=0.3, edgeStyleRule={'inter': '-', 'intra': '-', 'rule': 'edgetype'}, defaultNodeLabelSize=16, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels)
            
    figure.savefig("figs/multi_l2_n2_layer.pdf", bbox_inches='tight')


def multi_layers2_nodes2_part():
    
    n = 2
    n_layers = 2
    layers = list(range(n_layers))
    layerLabels = layer_labels(layers)
    layers = [set(layers)]
    nets, invs = graphlets_general(n, layers, 1)
    nets[2][0].add_layer(1)
    multi_auts = multi_automorphism_orbits(nets, allowed_aspects=[0,1])
    multi_orbit_is = multi_orbit_numbers(n, nets, multi_auts)
    color_ids = node_color_ids(n, nets, multi_auts)
    figure = plt.figure(figsize=(11,8))
    k = 1
    for i in nets:
        for j in range(len(nets[i])):
            ax = figure.add_subplot(3,4,k, projection='3d')
            k += 1
            net = nets[i][j]
            nodeColors, nodeLabels = multi_node_colors_and_labels(net, i, j, multi_auts, multi_orbit_is, color_ids)
            pymnet.draw(net, ax=ax, layout='circular', camera_dist=8, layergap = 0.7, defaultLayerAlpha = 0.3, layerPadding=0.3, edgeStyleRule={'inter': '-', 'intra': '-', 'rule': 'edgetype'}, defaultNodeLabelSize=16, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels)
            
            if k > 12:
                break
            
    figure.savefig("figs/multi_l2_n2_12.pdf", bbox_inches='tight')
    

def graphlets_general(n, layers, aspects, allowed_aspects='all'):
    '''
    Function for creating all the general multilayer graphlets up to n nodes
    
    Parameters
    ----------
    n: int
        maximum number of nodes
    layers: list of iterables
        each iterable contains the layers of one aspect
    aspects: int
        number of aspects
    allowed_aspects: list, string
        the aspects that can be permutated when computing isomorphisms
    
    Returns
    -------
    nets: dict (key: n_nodes, value: list of networks)
        graphlets
    invariants: dict (key: str(complete invariant), value: tuple(n_nodes, net_index in nets))
        complete invariants of the graphlets
        
    Notes
    -----
    the aggregated networks of the graphlets are connected,
    only undirected graphlets atm
    '''
    
    nets = {}
    invariants = {}
    nets[1] = []
    invariants[1] = {}
    
    layers_a = list(itertools.product(*layers))
    layer_combs = graphlets.layer_combinations(layers_a)
    for layer_comb in layer_combs:
        self_links = list(itertools.combinations(layer_comb, 2))
        n_links = len(self_links)
        for n_l in range(n_links + 1):
            for self_comb in itertools.combinations(self_links, n_l):
                net0 = pymnet.MultilayerNetwork(aspects=aspects, fullyInterconnected=False)
                for layer in layers_a:
                    for k in range(aspects):
                        net0.add_layer(layer[k], k)
                for layer in layer_comb:
                    if aspects < 2:
                        layer = layer[0]
                    net0.add_node(0, layer)
                    
                for e in self_comb:
                    edge = [0, 0]
                    for j in range(aspects):
                        edge += [e[0][j], e[1][j]]
                        
                    net0[tuple(edge)] = 1
                
                ci = pymnet.get_complete_invariant(net0, allowed_aspects)
                ci_s = str(ci)
                if not ci_s in invariants[1]:
                    invariants[1][ci_s] = (0, len(nets[1]))
                    nets[1].append(net0)
            
    for i in range(1, n):
        nets0 = nets[i]
        nets1 = []
        nodes = list(range(i))
        for net in nets0:
            for layer_comb in layer_combs:
                node_layers = net.iter_node_layers()
                edges = list(itertools.product(layer_comb, node_layers))
                for n_e in range(1, len(edges)+1):
                    for edge_comb in itertools.combinations(edges, n_e):
                        self_links = list(itertools.combinations(layer_comb, 2))
                        n_links = len(self_links)
                        for n_l in range(n_links + 1):
                            for self_comb in itertools.combinations(self_links, n_l):
                                new_net = pymnet.subnet(net, nodes, *layers)
                                for layer in layer_comb:
                                    if aspects < 2:
                                        layer = layer[0]
                                    new_net.add_node(i, layer)
                                for e in edge_comb:
                                    edge = [i, e[1][0]]
                                    for j in range(aspects):
                                        edge += [e[0][j], e[1][j+1]]
                                        
                                    new_net[tuple(edge)] = 1
                                    
                                for e in self_comb:
                                    edge = [i, i]
                                    for j in range(aspects):
                                        edge += [e[0][j], e[1][j]]
                                        
                                    new_net[tuple(edge)] = 1
                                    
                                ci = pymnet.get_complete_invariant(new_net, allowed_aspects)
                                ci_s = str(ci)
                                if not ci_s in invariants:
                                    invariants[ci_s] = (i+1, len(nets1))
                                    nets1.append(new_net)
                            
        nets[i+1] = nets1
        
    del nets[1]
    del invariants[1]
    
    return nets, invariants
    
    
def equation_intermediate():
    
    n = 5
    n_layers = 1
    layers = list(range(n_layers))
    nets, invs = graphlets.graphlets(n, layers)
    nets = order_nets(nets)
    invs = order_invs(invs)
    auts = graphlets.automorphism_orbits(nets)
    orbit_is = orbit_numbers(n, nets, auts)
    orbit_is = order_orbit_is(orbit_is)
    
    layerLabels = layer_labels(layers)
    #nodeCoords = {1 : (0.1, 0.1), 0 : (0.5, 0.5), -3 : (0.9, 0.9), -2 : (0.1,0.9), 2 : (0.9,0.1)}
    used_orbits = set()
    orbit1 = (3,0,0)
    orbit2 = (3,1,0)
    nets_co = graphlets.combine_orbits(orbit1, orbit2, nets)
    for net in nets_co:
        ci = str(pymnet.get_complete_invariant(net))
        i = invs[ci][0]
        j = invs[ci][1]
        k = 1
        figure = plt.figure(figsize=(12,2.5))
        i1 = orbit1[0]
        j1 = orbit1[1]
        net1 = nets[i1][j1]
        nodeColors, nodeLabels, orbit_i = co_node_colors_and_labels(orbit1[2], set(), net1, i1, j1, nets, auts, orbit_is)
        edgeColors, edgeWidths = edge_colors_and_widths_sub(net1, set([orbit1[2]]))
        nodeCoords = {0 : (0.7, 0.5), 1 : (0.3, 0.7), 2 : (0.3, 0.3)}
        ax = figure.add_subplot(1,4,k,projection='3d')
        pymnet.draw(net1, ax=ax, azim=0, elev=90, camera_dist=5.1, layout='circular', defaultNodeLabelSize=18, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels, nodeCoords=nodeCoords, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths)
        k += 1
        i2 = orbit2[0]
        j2 = orbit2[1]
        net2 = pymnet.subnet(net, [orbit1[2], -2, -3], layers)
        nodeColors, nodeLabels, orbit_i = co_node_colors_and_labels(orbit2[2], set(), net2, i2, j2, nets, auts, orbit_is)
        edgeColors, edgeWidths = edge_colors_and_widths_sub(net2, set([orbit2[2]]))
        nodeCoords = {0 : (0.7, 0.5), -3 : (0.3, 0.7), -2 : (0.3, 0.3)}
        ax = figure.add_subplot(1,4,k,projection='3d')
        pymnet.draw(net2, ax=ax, azim=0, elev=90, camera_dist=5.1, layout='circular', defaultNodeLabelSize=18, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels, nodeCoords=nodeCoords, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths)
        k += 1
        ax = figure.add_subplot(1,4,k,projection='3d')
        ax1tr = ax.transData
        pymnet.draw(net2, ax=ax, defaultEdgeAlpha=0, defaultNodeColor='white', defaultNodeLabelAlpha=0, azim=0, elev=90, camera_dist=5.1, layout='circular', defaultNodeLabelSize=18, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeLabelDict=nodeLabels, nodeCoords=nodeCoords, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths)
        k += 1
        nodeColors, nodeLabels, orbit_i = co_node_colors_and_labels(orbit2[2], set(), net, i, j, nets, auts, orbit_is)
        edgeColors, edgeWidths = edge_colors_and_widths_sub(net, set([orbit2[2]]))
        nodeCoords1 = {0 : (0.9, 0.5), -3 : (0.1, 0.25), 1 : (0.6, 0.1), -2 : (0.1, 0.75), 2 : (0.6, 0.9)}
        ax = figure.add_subplot(1,4,k,projection='3d')
        pymnet.draw(net, ax=ax, azim=0, elev=90, camera_dist=5.1, layout='circular', defaultNodeLabelSize=18, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels, nodeCoords=nodeCoords1, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths)
        figtr = figure.transFigure.inverted()
        ptB = figtr.transform(ax1tr.transform((-0.055, 0.0)))
        ptE = figtr.transform(ax1tr.transform((0.055, 0.0)))
        arrow = matplotlib.patches.FancyArrowPatch(ptB, ptE, fc = "k", transform=figure.transFigure, arrowstyle='simple', mutation_scale = 40)
        figure.patches.append(arrow)
        figure.savefig("figs/co_mono2_2.pdf", bbox_inches='tight')
        k = 1
        figure = plt.figure(figsize=(9,2.5))
        ax = figure.add_subplot(1,3,k,projection='3d')
        pymnet.draw(net, ax=ax, azim=0, elev=90, camera_dist=5.1, layout='circular', defaultNodeLabelSize=18, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels, nodeCoords=nodeCoords1, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths)
        k += 1
        ax = figure.add_subplot(1,3,k,projection='3d')
        pymnet.draw(net, ax=ax, defaultEdgeAlpha=0, defaultNodeColor='white', defaultNodeLabelAlpha=0, azim=0, elev=90, camera_dist=5.1, layout='circular', defaultNodeLabelSize=18, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeLabelDict=nodeLabels, nodeCoords=nodeCoords1, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths)
        ax1tr = ax.transData
        figtr = figure.transFigure.inverted()
        ptB = figtr.transform(ax1tr.transform((-0.055, 0.0)))
        ptE = figtr.transform(ax1tr.transform((0.055, 0.0)))
        arrow = matplotlib.patches.FancyArrowPatch(ptB, ptE, fc = "k", transform=figure.transFigure, arrowstyle='simple', mutation_scale = 40)
        figure.patches.append(arrow)
        k += 1
        nets_n_nodes = graphlets.merge_nodes([orbit1[2]], net)
        net_m = nets_n_nodes[0][0]
        both_m = nets_n_nodes[0][1]
        ci_m = str(pymnet.get_complete_invariant(net_m))
        im = invs[ci_m][0]
        jm = invs[ci_m][1]
        nodeColors, nodeLabels, orbit_i = co_node_colors_and_labels(orbit2[2], both_m, net_m, im, jm, nets, auts, orbit_is)
        edgeColors, edgeWidths = edge_colors_and_widths_sub(net_m, both_m)
        nodeCoords2 = {0 : (0.9, 0.5), -2 : (0.1, 0.5), 2 : (0.5, 0.1), 1 : (0.5, 0.9)}
        ax = figure.add_subplot(1,3,k,projection='3d')
        pymnet.draw(net_m, ax=ax, azim=0, elev=90, camera_dist=5.1, layout='circular', defaultNodeLabelSize=18, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels, nodeCoords=nodeCoords2, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths)
        nets_ad = graphlets.add_possible_edges([orbit1[2]], net)
        nets_ad += graphlets.add_possible_edges(both_m, net_m)
        figure.savefig("figs/mn_mono2.pdf", bbox_inches='tight')
        
        figure = plt.figure(figsize=(9,9))
        k = 1
        for net_a in nets_ad:
            ci = str(pymnet.get_complete_invariant(net_a))
            i = invs[ci][0]
            j = invs[ci][1]
            if k == 7:
                both_orbits = both_m
                nodeCoords = nodeCoords2
            else:
                both_orbits = set([orbit1[2]])
                nodeCoords = nodeCoords1
            nodeColors, nodeLabels, orbit_i = co_node_colors_and_labels(orbit1[2], both_orbits, net_a, i, j, nets, auts, orbit_is)
            if not orbit_i in used_orbits:
                ax = figure.add_subplot(3,3,k,projection='3d')
                k += 1
                edgeColors, edgeWidths = edge_colors_and_widths(net_a, both_orbits)
                pymnet.draw(net_a, ax=ax, azim=0, elev=90, camera_dist=5.1, layout='circular', defaultNodeLabelSize=18, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeColorDict = nodeColors, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths, nodeLabelDict=nodeLabels, nodeCoords=nodeCoords)
                used_orbits.add(orbit_i)
                
        figure.savefig("figs/ad_mono2.pdf", bbox_inches='tight')
        
        figure = plt.figure(figsize=(10,3))
        net = nets_ad[4]
        ci = str(pymnet.get_complete_invariant(net))
        i = invs[ci][0]
        j = invs[ci][1]
        nodeColors, nodeLabels, orbit_i = co_node_colors_and_labels(orbit2[2], set(), net, i, j, nets, auts, orbit_is)
        edgeColors, edgeWidths = edge_colors_and_widths_sub(net, set([orbit1[2]]))
        ax = figure.add_subplot(1,3,1,projection='3d')
        pymnet.draw(net, ax=ax, azim=0, elev=90, camera_dist=5.1, layout='circular', defaultNodeLabelSize=18, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels, nodeCoords=nodeCoords1, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths)
        nodeCoords1 = {0 : (0.9, 0.5), 1 : (0.1, 0.25), -3 : (0.6, 0.1), -2 : (0.1, 0.75), 2 : (0.6, 0.9)}
        ax = figure.add_subplot(1,3,2,projection='3d')
        pymnet.draw(net, ax=ax, azim=0, elev=90, camera_dist=5.1, layout='circular', defaultNodeLabelSize=18, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels, nodeCoords=nodeCoords1, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths)
        nodeCoords1 = {0 : (0.9, 0.5), -3 : (0.1, 0.25), -2 : (0.6, 0.1), 1 : (0.1, 0.75), 2 : (0.6, 0.9)}
        ax = figure.add_subplot(1,3,3,projection='3d')
        pymnet.draw(net, ax=ax, azim=0, elev=90, camera_dist=5.1, layout='circular', defaultNodeLabelSize=18, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels, nodeCoords=nodeCoords1, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths)
        figure.savefig("figs/sub_mono2.pdf", bbox_inches='tight')
        
                
def equation_intermediate2():
    
    n = 4
    n_layers = 3
    layers = list(range(n_layers))
    nets, invs = graphlets.graphlets(n, layers)   
    auts = graphlets.automorphism_orbits(nets)
    orbit_is = orbit_numbers(n, nets, auts)
    
    layerLabels = layer_labels(layers)
    nodeCoords = {-3 : (0.1, 0.1), 1 : (0.1, 0.9), -2 : (0.9, 0.1), 0 : (0.9, 0.9)}
    i1 = 2
    j1 = 0
    i2 = 3
    j2 = 15
    orbit1 = (i1,j1,0)
    orbit2 = (i2,j2,0)
    nets_co = graphlets.combine_orbits(orbit1, orbit2, nets)
    figure = plt.figure(figsize=(6.3,4))
    net1 = nets[i1][j1]
    ax = figure.add_subplot(1,2,1,projection='3d')
    nodeColors, nodeLabels, orbit_i = co_node_colors_and_labels(orbit1[2], set(), net1, i1, j1, nets, auts, orbit_is)
    pymnet.draw(net1, ax=ax, layout='circular', defaultLayerAlpha=0.3, defaultNodeLabelSize=18, layergap=0.7, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels)
    net2 = pymnet.subnet(nets_co[0], [orbit1[2], -2, -3], layers)
    ax = figure.add_subplot(1,2,2,projection='3d')
    nodeColors, nodeLabels, orbit_i = co_node_colors_and_labels(orbit1[2], set(), net2, i2, j2, nets, auts, orbit_is)
    pymnet.draw(net2, ax=ax, layout='circular', defaultLayerAlpha=0.3, defaultNodeLabelSize=18, layergap=0.7, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels)
    figure.savefig("figs/co_plex0.pdf", bbox_inches='tight')
    figure = plt.figure(figsize=(10,4))
    k = 1
    for net in nets_co:
        ci = str(pymnet.get_complete_invariant(net))
        i = invs[ci][0]
        j = invs[ci][1]
        nodeColors, nodeLabels, orbit_i = co_node_colors_and_labels(orbit2[2], set(), net, i, j, nets, auts, orbit_is)
        ax = figure.add_subplot(1,3,k,projection='3d')
        pymnet.draw(net, ax=ax, layout='spectral', defaultLayerAlpha=0.3, defaultNodeLabelSize=18, layergap=0.7, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels, nodeCoords=nodeCoords)
        k += 1
        
    figure.savefig("figs/co_plex.pdf", bbox_inches='tight')
    
    
def equation():
    
    n = 5
    n_layers = 1
    layers = list(range(n_layers))
    nets, invs = graphlets.graphlets(n, layers)   
    auts = graphlets.automorphism_orbits(nets)
    orbit_is = orbit_numbers(n, nets, auts)
    orbit_eqs = graphlets.orbit_equations(n, layers, nets, auts, invs)
    
    layerLabels = layer_labels(layers)
    orbit0 = (3,1,0)
    eq = orbit_eqs[(orbit0, 2)]
    figure = plt.figure(figsize=(11,5))
    k = 1
    for orbit in eq:
        i = orbit[0]
        j = orbit[1]
        node = orbit[2]
        net = nets[i][j]
        nodeColors, nodeLabels, _ = co_node_colors_and_labels(node, set(), net, i, j, nets, auts, orbit_is)
        
        if k == 1:
            m = 6
        elif k == 2:
            m = 5
        elif k == 3:
            m = 4
        elif k == 4:
            m = 1
        elif k == 5:
            m = 3
        elif k == 6:
            m = 8
        elif k == 7:
            m = 2
        elif k == 8:
            m = 7
        else:
            m = k
        
        ax = figure.add_subplot(2,4,m,projection='3d')
        k += 1
        pymnet.draw(net, ax=ax, azim=0, elev=90, camera_dist=5.1, layout='spring', defaultNodeLabelSize=18, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels)


def labels_colors_simple(net, node0, node_comb, nets, auts, invs, orbit_is):

    nodeLabels = {}
    nodeColors = {}
    edgeColors = {}
    edgeWidths = {}
    edgeStyles = {}
    
    n = len(node_comb)
    layers = net.slices[1]
    sub = pymnet.subnet(net, node_comb, layers)
    ci_sub = str(pymnet.get_complete_invariant(sub))
    iso = invs[ci_sub]
    i = iso[0]
    j = iso[1]
    iso_net = nets[i][j]
    iso = pymnet.get_isomorphism(sub, iso_net)
    k = iso[0][node0]
    
    nodes = net.slices[0]
    for node in nodes:
        for l in layers:
            if node == node0:
                aut = auts[(i,j,k)]
                nodeLabels[(node, l)] = orbit_is[(i,j,aut)]
            else:
                nodeLabels[(node, l)] = ''
                
            if node in node_comb:
                nodeColors[(node, l)] = '#ff2387'
            else:
                nodeColors[(node, l)] = 'grey'
            
    for e in net.edges:
        n1 = e[0]
        n2 = e[1]
        l = e[2]
        if n1 in node_comb and n2 in node_comb:
            edgeColors[((n1,l), (n2, l))] = 'black'
            edgeWidths[((n1,l), (n2, l))] = 2
        else:
            edgeStyles[((n1,l), (n2, l))] = ':'
            
    return nodeLabels, nodeColors, edgeColors, edgeWidths, edgeStyles
    
    
def simple_example():
    
    net = pymnet.MultilayerNetwork(aspects=1)
    net[1,2,'a'] = 1
    net[1,3,'a'] = 1
    net[1,4,'a'] = 1
    net[1,5,'a'] = 1
    net[2,3,'a'] = 1
    net[4,5,'a'] = 1
    
    layerLabels = {'a' : ""}
    nodeCoords = {1 : (0.5, 0.5), 2 : (0.1, 0.6), 3 : (0.1, 0.4), 4 : (0.9, 0.6), 5 : (0.9, 0.4)} 
    
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(1,1,1,projection='3d')
    pymnet.draw(net, ax=ax, azim=-90, elev=90, camera_dist=5.1, defaultLayerAlpha=0, defaultNodeLabelSize=18, layerLabelDict = layerLabels, nodeCoords=nodeCoords)
    fig.savefig("figs/simple_example.pdf", bbox_inches='tight')
    
    n = 5
    nets, invs = graphlets.graphlets(n, net.slices[1])
    nets = order_nets(nets)
    invs = order_invs(invs)
    auts = graphlets.automorphism_orbits(nets)
    orbit_is = orbit_numbers(n, nets, auts)
    orbit_is = order_orbit_is(orbit_is)
    
    node0 = 1
    
    figure = plt.figure(figsize=(15,9))
    k = 1
    
    for i in range(1,5):
        for comb in itertools.combinations(range(2,6), i):
            node_comb = set(comb) | set([node0])
            nodeLabels, nodeColors, edgeColors, edgeWidths, edgeStyles = labels_colors_simple(net, node0, node_comb, nets, auts, invs, orbit_is)
            
            ax = figure.add_subplot(3,5,k,projection='3d')
            k += 1
            pymnet.draw(net, ax=ax, azim=-90, elev=90, camera_dist=5.1, defaultLayerAlpha=0, defaultNodeLabelSize=18, layerLabelDict = layerLabels, nodeCoords=nodeCoords, nodeColorDict=nodeColors, nodeLabelDict=nodeLabels, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths, edgeStyleDict=edgeStyles)
            
    figure.savefig("figs/graphlet_degree.pdf", bbox_inches='tight')


def mono_iso_example():
    
    n = 5
    n_layers = 1
    layers = list(range(n_layers))
    nets, invs = graphlets.graphlets(n, layers)
    nets = order_nets(nets)
    
    auts = graphlets.automorphism_orbits(nets)
    orbit_is = orbit_numbers(n, nets, auts)
    orbit_is = order_orbit_is(orbit_is)
    
    color_ids = node_color_ids(n, nets, auts)
    layerLabels = layer_labels(layers)
    figure = plt.figure(figsize=(15,18))
    k = 1
    i = 5
    j = 9
    coord_k = 19
    nodeLabels = {0 : '0', 1 : '1', 2 : '2', 3 : '3', 4 : '4'}
    for l in range(3):
        nodeCoords = node_coords(coord_k)
        ax = figure.add_subplot(6,5,k, projection='3d')
        net = nets[i][j]
        #nodeColors, nodeLabels = node_colors_and_labels(net, i, j, auts, orbit_is, color_ids)
        pymnet.draw(net, ax=ax, azim=0, elev=90, camera_dist=6, autoscale=False, defaultNodeLabelSize=16, layout='spring', defaultLayerAlpha = 0.0, layergap = 0.15, layerLabelDict=layerLabels, nodeLabelDict=nodeLabels, nodeSizeRule={'scalecoeff': 0.3, 'rule': 'scaled'}, nodeCoords=nodeCoords)
        plt.title(r'$\mathbf{G_{' + str(l) + '}}$', fontsize=24, ha='left', position=(0,1), fontweight='heavy')
        k += 1
        if l == 0:
            nodeLabels = {(0,0) : '3', (1,0) : '4', (2,0) : '1', (3,0) : '0', (4,0) : '2'}
            
        else:
            nodeLabels = {(0,0) : '0', (1,0) : '4', (2,0) : '1', (3,0) : '2', (4,0) : '3'}
            
    figure.savefig("figs/iso_mono.pdf", bbox_inches='tight')
    
    
def isomorphism_example():
    
    net0 = pymnet.MultilayerNetwork(aspects=1, fullyInterconnected=False)
    net0.add_node(1, 'a')
    net0[1,2,'b','a'] = 1
    
    net1 = pymnet.MultilayerNetwork(aspects=1, fullyInterconnected=False)
    net1.add_node(2, 'a')
    net1[1,2,'a','b'] = 1
    
    net2 = pymnet.MultilayerNetwork(aspects=1, fullyInterconnected=False)
    net2.add_node(1, 'b')
    net2[1,2,'a','b'] = 1
    
    net3 = pymnet.MultilayerNetwork(aspects=1, fullyInterconnected=False)
    net3.add_node(2, 'b')
    net3[1,2,'b','a'] = 1
    
    nets = [net0, net1, net2, net3]
    
    nodeCoords = [{2 : (0.9, 0.9), 1 : (0.1, 0.1)}, {1 : (0.9, 0.9), 2 : (0.1, 0.1)}, {2 : (0.9, 0.9), 1 : (0.1, 0.1)}, {1 : (0.9, 0.9), 2 : (0.1, 0.1)}]
    layerOrders = [{'a' : 1, 'b' : 2}, {'a' : 1, 'b' : 2}, {'a' : 2, 'b' : 1}, {'a' : 2, 'b' : 1}]
    
    figure = plt.figure(figsize=(10,8))
    for i in range(len(nets)):
        k = i + 1
        ax = figure.add_subplot(2,2,k,projection='3d')
        pymnet.draw(nets[i], ax=ax, layout="spectral", camera_dist=8.2, defaultLayerAlpha=0.3, defaultNodeLabelSize=18, defaultLayerLabelSize=18, defaultLayerLabelLoc=(-0.2,0), defaultEdgeColor='black', defaultEdgeWidth=3, layerPadding=0.2, nodeCoords=nodeCoords[i], layerOrderDict=layerOrders[i])
    
    figure.savefig("figs/iso_ex.pdf", bbox_inches='tight')
    
    
def isomorphism_example2():
    
    net0 = pymnet.MultilayerNetwork(aspects=1, fullyInterconnected=False)
    net0[1,2,'a'] = 1
    net0[2,3,'a'] = 1
    net0[1,1,'a','b'] = 1
    net0[1,3,'b'] = 1
    
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(1,1,1,projection='3d')
    pymnet.draw(net0, ax=ax, layout="circular", elev=30, camera_dist=8.6, defaultLayerAlpha=0.3, defaultNodeLabelSize=18, defaultLayerLabelSize=18, layerPadding=0.2, defaultLayerLabelLoc=(-0.2,0))
    plt.title(r'$\mathbf{M_{' + str(0) + '}}$', fontsize=24, ha='left', position=(0,1), fontweight='heavy')
    fig.savefig("figs/iso_ex2_1.pdf", bbox_inches='tight')
    
    nodeLabels1 = {(1,'a') : 2, (1,'b') : 2, (2,'a') : 3, (3,'a') : 1, (3, 'b') : 1}
    nodeLabels2 = {(1,'a') : 1, (1,'b') : 1, (2,'a') : 2, (3,'a') : 3, (3, 'b') : 3}
    nodeLabels3 = {(1,'a') : 2, (1,'b') : 2, (2,'a') : 3, (3,'a') : 1, (3, 'b') : 1}
    nodeLabels4 = {(1,'a') : 2, (1,'b') : 3, (2,'a') : 3, (3,'a') : 1, (3, 'b') : 1}
    nodeLabels = [nodeLabels1, nodeLabels2, nodeLabels3, nodeLabels4]
    
    layerLabels = [{'a' : 'a', 'b' : 'b'}, {'a' : 'b', 'b' : 'a'}, {'a' : 'b', 'b' : 'a'}, {'a' : 'a', 'b' : 'b'}]
    
    figure = plt.figure(figsize=(16,4))
    for i in range(4):
        k = i + 1
        ax = figure.add_subplot(1,4,k,projection='3d')
        pymnet.draw(net0, ax=ax, layout="circular", elev=30, camera_dist=8.6, defaultLayerAlpha=0.3, defaultNodeLabelSize=18, defaultLayerLabelSize=18, layerPadding=0.2, defaultLayerLabelLoc=(-0.2,0), nodeLabelDict=nodeLabels[i], layerLabelDict=layerLabels[i])
        plt.title(r'$\mathbf{M_{' + str(k) + '}}$', fontsize=24, ha='left', position=(0,1), fontweight='heavy')
        
    figure.savefig("figs/iso_ex2_2.pdf", bbox_inches='tight')
    
    
def orbit_x3_example():
    
    n = 5
    n_layers = 1
    layers = list(range(n_layers))
    nets, invs = graphlets.graphlets(n, layers)
    nets = order_nets(nets)
    invs = order_invs(invs)
    auts = graphlets.automorphism_orbits(nets)
    orbit_is = orbit_numbers(n, nets, auts)
    orbit_is = order_orbit_is(orbit_is)
    
    layerLabels = layer_labels(layers)
    orbit1 = (3,0,0)
    orbit2 = (3,0,0)
    nets_co = graphlets.combine_orbits(orbit1, orbit2, nets)
    
    net = nets_co[0]
    ci = str(pymnet.get_complete_invariant(net))
    i = invs[ci][0]
    j = invs[ci][1]
    k = 1
    figure = plt.figure(figsize=(12,5))
    i1 = orbit1[0]
    j1 = orbit1[1]
    net1 = nets[i1][j1]
    #nodeColors, nodeLabels, orbit_i = co_node_colors_and_labels(orbit1[2], set(), net1, i1, j1, nets, auts, orbit_is)
    edgeColors, edgeWidths = edge_colors_and_widths_sub(net1, set([orbit1[2]]))
    nodeCoords = {0 : (0.7, 0.5), 1 : (0.3, 0.7), 2 : (0.3, 0.3)}
    nodeLabels = {(0, 0) : '2', (1, 0) : '', (2, 0) : '', (-2, 0) : '', (-3, 0) : ''}
    ax = figure.add_subplot(2,4,k,projection='3d')
    pymnet.draw(net1, ax=ax, azim=0, elev=90, camera_dist=5.1, layout='circular', defaultNodeLabelSize=18, defaultNodeLabelAlpha=1.0, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeCoords=nodeCoords, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths, nodeLabelDict=nodeLabels) #nodeColorDict = nodeColors
    k += 1
    i2 = orbit2[0]
    j2 = orbit2[1]
    net2 = pymnet.subnet(net, [orbit1[2], -2, -3], layers)
    #nodeColors, nodeLabels, orbit_i = co_node_colors_and_labels(orbit2[2], set(), net2, i2, j2, nets, auts, orbit_is)
    edgeColors, edgeWidths = edge_colors_and_widths_sub(net2, set([orbit2[2]]))
    nodeCoords = {0 : (0.7, 0.5), -3 : (0.3, 0.7), -2 : (0.3, 0.3)}
    ax = figure.add_subplot(2,4,k,projection='3d')
    pymnet.draw(net2, ax=ax, azim=0, elev=90, camera_dist=5.1, layout='circular', defaultNodeLabelSize=18, defaultNodeLabelAlpha=1.0, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeCoords=nodeCoords, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths, nodeLabelDict=nodeLabels) # nodeColorDict = nodeColors, nodeLabelDict=nodeLabels
    k += 1
    ax = figure.add_subplot(2,4,k,projection='3d')
    ax1tr = ax.transData
    pymnet.draw(net2, ax=ax, defaultEdgeAlpha=0, defaultNodeColor='white', defaultNodeLabelAlpha=0, azim=0, elev=90, camera_dist=5.1, layout='circular', defaultNodeLabelSize=18, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeCoords=nodeCoords, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths) #nodeLabelDict=nodeLabels
    k += 1
    #nodeColors, nodeLabels, orbit_i = co_node_colors_and_labels(orbit2[2], set(), net, i, j, nets, auts, orbit_is)
    edgeColors, edgeWidths = edge_colors_and_widths_sub(net, set([orbit2[2]]))
    #nodeCoords1 = {0 : (0.9, 0.5), -3 : (0.1, 0.25), 1 : (0.6, 0.1), -2 : (0.1, 0.75), 2 : (0.6, 0.9)}
    ax = figure.add_subplot(2,4,k,projection='3d')
    nodeCoords2 = {0 : (0.5, 0.5), 1 : (0.8, 0.8), 2 : (0.2, 0.8), -2 : (0.2, 0.2), -3 : (0.8, 0.2)}
    nodeLabels2 = {(0, 0) : '716', (1, 0) : '', (2, 0) : '', (-2, 0) : '', (-3, 0) : ''}
    pymnet.draw(net, ax=ax, azim=0, elev=90, camera_dist=5.1, layout='spring', defaultNodeLabelSize=18, defaultNodeLabelAlpha=1.0, defaultLayerAlpha=.0, layerLabelDict=layerLabels, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths, nodeCoords=nodeCoords2, nodeLabelDict=nodeLabels2)#, nodeColorDict = nodeColors
    figtr = figure.transFigure.inverted()
    ptB = figtr.transform(ax1tr.transform((-0.055, 0.0)))
    ptE = figtr.transform(ax1tr.transform((0.055, 0.0)))
    arrow = matplotlib.patches.FancyArrowPatch(ptB, ptE, fc = "k", transform=figure.transFigure, arrowstyle='simple', mutation_scale = 40)
    figure.patches.append(arrow)
    
    k += 1
    i1 = orbit1[0]
    j1 = orbit1[1]
    net1 = nets[i1][j1]
    #nodeColors, nodeLabels, orbit_i = co_node_colors_and_labels(orbit1[2], set(), net1, i1, j1, nets, auts, orbit_is)
    edgeColors, edgeWidths = edge_colors_and_widths_sub(net2, set([orbit1[2]]))
    #nodeCoords = {0 : (0.7, 0.5), 1 : (0.3, 0.7), 2 : (0.3, 0.3)}
    ax = figure.add_subplot(2,4,k,projection='3d')
    pymnet.draw(net2, ax=ax, azim=0, elev=90, camera_dist=5.1, layout='circular', defaultNodeLabelSize=18, defaultNodeLabelAlpha=1.0, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeCoords=nodeCoords, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths, nodeLabelDict=nodeLabels) #nodeColorDict = nodeColors, nodeLabelDict=nodeLabels
    k += 1
    i1 = orbit1[0]
    j1 = orbit1[1]
    net1 = nets[i1][j1]
    #nodeColors, nodeLabels, orbit_i = co_node_colors_and_labels(orbit1[2], set(), net1, i1, j1, nets, auts, orbit_is)
    edgeColors, edgeWidths = edge_colors_and_widths_sub(net2, set([orbit1[2]]))
    #nodeCoords = {0 : (0.7, 0.5), 1 : (0.3, 0.7), 2 : (0.3, 0.3)}
    ax = figure.add_subplot(2,4,k,projection='3d')
    pymnet.draw(net2, ax=ax, azim=0, elev=90, camera_dist=5.1, layout='circular', defaultNodeLabelSize=18, defaultNodeLabelAlpha=1.0, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeCoords=nodeCoords, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths, nodeLabelDict=nodeLabels) #nodeColorDict = nodeColors, nodeLabelDict=nodeLabels
    k += 1
    ax = figure.add_subplot(2,4,k,projection='3d')
    ax1tr = ax.transData
    pymnet.draw(net2, ax=ax, defaultEdgeAlpha=0, defaultNodeColor='white', defaultNodeLabelAlpha=0, azim=0, elev=90, camera_dist=5.1, layout='circular', defaultNodeLabelSize=18, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeCoords=nodeCoords, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths) #nodeLabelDict=nodeLabels
    figtr = figure.transFigure.inverted()
    ptB = figtr.transform(ax1tr.transform((-0.055, 0.0)))
    ptE = figtr.transform(ax1tr.transform((0.055, 0.0)))
    arrow = matplotlib.patches.FancyArrowPatch(ptB, ptE, fc = "k", transform=figure.transFigure, arrowstyle='simple', mutation_scale = 40)
    figure.patches.append(arrow)
    k += 1
    ax = figure.add_subplot(2,4,k,projection='3d')
    nodeLabels2 = {(0, 0) : '412', (1, 0) : '', (2, 0) : '', (-2, 0) : '', (-3, 0) : ''}
    pymnet.draw(net, ax=ax, azim=0, elev=90, camera_dist=5.1, layout='spring', defaultNodeLabelSize=18, defaultNodeLabelAlpha=1.0, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, defaultEdgeColor='#ff2387', edgeWidthDict=edgeWidths, nodeCoords=nodeCoords2, nodeLabelDict=nodeLabels2) #nodeCoords=nodeCoords1, nodeColorDict = nodeColors, nodeLabelDict=nodeLabels
    
    figure.savefig("figs/co_plex3.pdf", bbox_inches='tight')
    
    figure = plt.figure(figsize=(4,4))
    ax = figure.add_subplot(1,1,1,projection='3d')
    edgeColors = {((0,0),(-2,0)) : '#ff2387', ((0,0),(-3,0)) : '#5723ff'}
    nodeLabels = {(0, 0) : '4', (1, 0) : '', (2, 0) : '', (-2, 0) : '', (-3, 0) : ''}
    pymnet.draw(net2, ax=ax, azim=0, elev=90, camera_dist=5.1, layout='circular', defaultNodeLabelSize=26, defaultNodeLabelAlpha=1.0, defaultLayerAlpha=0.0, layerLabelDict=layerLabels, nodeCoords=nodeCoords, edgeColorDict=edgeColors, edgeWidthDict=edgeWidths, nodeLabelDict=nodeLabels) #nodeColorDict = nodeColors, nodeLabelDict=nodeLabels
    
    figure.savefig("figs/co_plex3_0.pdf", bbox_inches='tight')
    
    
    
if __name__ == "__main__":
    main()