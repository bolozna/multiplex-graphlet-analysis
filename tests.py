import filecmp

def test_example_network_interface_equivalence():
    '''
    Tests if example network results are the same for pipeline functions
    and interface functions. Results need to be in results folder.
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
    
if __name__ == '__main__':
    test_example_network_interface_equivalence()
