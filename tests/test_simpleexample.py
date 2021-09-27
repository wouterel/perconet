
import sys
sys.path.append("..")
import perconet as pn
import numpy as np

  
number_of_nodes = 4
max_coordination = 6
testnet = pn.PeriodicNetwork(number_of_nodes, max_coordination,verbose=False)

testnet.add_edge(0,1,np.array([0,0,0])) #add a bond between nodes 0 and 1 that crosses no boundary
testnet.add_edge(2,1,np.array([0,0,0])) #add a bond between nodes 2 and 1 that crosses no boundary
testnet.add_edge(2,3,np.array([0,0,0])) #add a bond between nodes 2 and 3 that crosses no boundary
testnet.add_edge(0,3,np.array([0,0,0])) #add a bond between nodes 0 and 3 that crosses no boundary

#we now have a square that doesn't do anything

if not testnet.crosses_boundaries:
    print ("Network does not cross periodic boundaries in any direction, therefore it does not percolate")


testnet.add_edge(1,0,np.array([1,0,0])) #add a bond between nodes 1 and 0 that crosses the x-boundary
testnet.add_edge(2,3,np.array([1,0,0])) #add a bond between nodes 1 and 0 that crosses the x-boundary
testnet.add_edge(3,0,np.array([0,-1,0])) #add a bond between nodes 3 and 0 that crosses the y-boundary in the negative direction


if testnet.needs_reducing:
    print ("Reducing network to simpler form...")
    my_reduced_network = testnet.get_reduced_network() #think of a way to return it already in the right format you would get with add_edges
    #print(my_reduced_network)
    myloops=pn.LoopFinder(my_reduced_network,verbose=False)
else:          
    myloops=pn.LoopFinder(testnet,verbose=False)
loops=myloops.get_independent_loops()


print ("number of loops \n" , len(loops[0]), " \n independent loops: \n", loops[0])


