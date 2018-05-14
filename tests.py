import filecmp
import pymnet
import interface
import data_analysis
import pandas as pd

def test_example_network_interface_equivalence():
    '''
    Tests if example network results are the same for pipeline functions
    and interface functions. Results need to be in results folder (interface
    results appended with "interface_")
    '''
    directory = 'Results/'
    total = 0
    successes = 0
    for name in ['ba_plex','er_0','er_20']:
        for network_num in range(0,5):
            for nlayers in [1,2,3]:
                fullname = name+'_'+str(network_num)
                folder = directory+fullname+'_'+str(nlayers)+'/'
                interface_folder = directory+'interface_'+fullname+'_'+str(nlayers)+'/'
                if nlayers == 1:
                    total += 1
                    if filecmp.cmp(folder+fullname+'_0.txt',interface_folder+'interface_'+fullname+'_0.txt',shallow=False):
                        successes += 1
                elif nlayers == 2:
                    total += 1
                    if (filecmp.cmp(folder+fullname+'_0_1.txt',interface_folder+'interface_'+fullname+'_0_1.txt',shallow=False)
                        and filecmp.cmp(folder+fullname+'_0_2.txt',interface_folder+'interface_'+fullname+'_0_2.txt',shallow=False)
                        and filecmp.cmp(folder+fullname+'_1_2.txt',interface_folder+'interface_'+fullname+'_1_2.txt',shallow=False)):
                            successes += 1
                elif nlayers == 3:
                    total += 1
                    if filecmp.cmp(folder+fullname+'_0_1_2.txt',interface_folder+'interface_'+fullname+'_0_1_2.txt',shallow=False):
                            successes += 1
    print('The number of networks is '+str(total))
    print('The number of successes is '+str(successes))
    return total,successes
    
def test_orbit_count_sum_functions():
    '''
    Loads all example results and checks that they sum to equivalent data frames.
    Checks all three loading functions.
    '''
    # what's causing the difference in PR curves?
    directory = 'Results/'
    allowed_aspects = 'all'
    
    total = 0
    successes = 0
    for name in ['ba_plex','er_0','er_20']:
        for network_num in range(0,5):
            for nlayers in [1,2,3]:
                fullname = name+'_'+str(network_num)
                folder = directory+fullname+'_'+str(nlayers)
                
                if nlayers == 1 or nlayers == 2:
                    n = 4
                else:
                    n = 3
                    
                if nlayers == 1:
                    net_layers = [0]
                else:
                    net_layers = [0,1,2]
                
                graphlets,invs = pymnet.graphlets.graphlets(n,net_layers,nlayers,allowed_aspects=allowed_aspects)
                auts = pymnet.graphlets.automorphism_orbits(graphlets,allowed_aspects=allowed_aspects)
                orbit_is = pymnet.graphlets.orbit_numbers(n,graphlets,auts)
                orbit_list = pymnet.graphlets.ordered_orbit_list(orbit_is)
                
                df1 = data_analysis.sum_orbit_counts(folder,orbit_list)
                df2 = data_analysis.sum_orbit_counts_v2(folder,orbit_list)
                df3 = interface.read_graphlet_degree_distribution_folder(folder)
                
                total += 1
                if df1.equals(df2) and df2.equals(df3) and df3.equals(df1):
                    successes += 1
                
    print('The number of networks is '+str(total))
    print('The number of successes is '+str(successes))
    return total,successes
    
def test_nonnumeric_labels():
    net = pymnet.MultiplexNetwork(couplings='categorical')
    net['node1','layer1']['node2','layer1'] = 1
    net['node1','layer2']['node2','layer2'] = 1
    net['node2','layer2']['node3','layer2'] = 1
    orbit_list = interface.graphlet_degree_distributions(net,3,2,save_name='nonnumeric')
    df = interface.read_graphlet_degree_distribution_folder('Results/nonnumeric_2')
    df.sort_index(inplace=True)
    colnames = [str(s) for s in orbit_list]
    comparison_df = pd.DataFrame(0,index=['node1','node2','node3'],columns=colnames)
    comparison_df.at['node1','(2, 1, 0)'] = 1
    comparison_df.at['node1','(3, 2, 2)'] = 1
    comparison_df.at['node3','(2, 0, 0)'] = 1
    comparison_df.at['node3','(3, 2, 1)'] = 1
    comparison_df.at['node2','(2, 0, 0)'] = 1
    comparison_df.at['node2','(2, 1, 0)'] = 1
    comparison_df.at['node2','(3, 2, 0)'] = 1
    comparison_df.sort_index(inplace=True)
    assert df.equals(comparison_df), "Nonnumeric labels not read correctly"
    print('Nonnumeric read successful')
    
if __name__ == '__main__':
    print('Interface:')
    test_example_network_interface_equivalence()
    print('Orbit sum/load functions:')
    test_orbit_count_sum_functions()
    print('Nonnumeric labels:')
    test_nonnumeric_labels()
