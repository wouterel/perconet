# This file is part of the PercoNet package
# (c) 2022 Eindhoven University of Technology
# Released under the BSD 3-clause license
# See LICENSE file for details
# Authors: Chiara Raffaelli, Wouter G. Ellenbroek
# Contributing: https://github.com/wouterel/perconet

import numpy as np


class PeriodicNetwork:
    """Store and analyze the topology of a periodic net.

    Periodic nets are graphs embedded in a periodic topology. This class stores
    the topology of such a graph for the case of a threedimensional periodic
    box (a 3-torus) and allows to determine its percolation properties.
    The key question to determine percolation is whether the graph contains any
    loops that are topologically nontrivial in that they allow to travel from a
    node to itself with a net nonzero number of boundary crossings.

    Args:
        n (int):
            The number of nodes of the graph.
        max_degree (int):
            The largest number of edges coming out of any node.
    """
    def __init__(self, n: int, max_degree=6, verbose=True):
        if n < 1:
            raise ValueError("Number of nodes must be a positive integer.")
        self.number_of_nodes = n
        self.max_degree = max_degree
        # allocate arrays for per-node info
        self.boundary_crossing = np.zeros((n, max_degree, 3), dtype=int)
        self.neighbors = -1 * np.ones((n, max_degree), dtype=int)
        self.edges_list = -1 * np.ones((n, max_degree), dtype=int)
        self.neighbors_counter = np.zeros(n, dtype=int)
        self.edges_counter = 0  # keeps track of the total number of edges while building edge list
        self.simple_edges_list = []  # is this duplicate info?
        self.simple_boundary_crossing = []
        self.bond_is_across_boundary = []
        self.n_internal_edges = 0
        self.crosses_boundaries = 0  # should I set anything? is there a way to keep ot empty?
        self.verbose = verbose

    def get_number_of_nodes(self):
        return self.number_of_nodes

    def add_edge(self, node1: int, node2: int, boundary_vector):
        # by default(for now, boundary crossing information is taken from boundary vector)
        """
        Add an edge to the periodic network

        Args:
            node1 (int):
                The index of the first node of the pair that defines this edge.
                Must satisfy 0 <= node1 < number_of_nodes (index is 0-based)
            node2 (int):
                The index of the second node of the pair that defines this edge.
                Must satisfy 0 <= node2 < number_of_nodes (index is 0-based)
            boundary_vector ((3) int):
                Vector of three integers (bvx,bvy,bvz) denoting the number of times the edge
                wraps around the x, y, and z boundaries, respectively.
                The sign indicates the wrapping direction (e.g. (-1,0,0) indicates that the edge
                goes around the x-boundary in the negative x-direction
                when going from node1 to node2.

        Returns:
            (bool): True if succesful. False if an error occurred.
        """
        if(node1 >= self.number_of_nodes):
            print(f"Error in add_edge(): node {node1} does not exist (N = {self.number_of_nodes})")
            return False
        if(node2 >= self.number_of_nodes):
            print(f"Error in add_edge(): node {node2} does not exist (N = {self.number_of_nodes})")
            return False
        if(self.neighbors_counter[node1] == self.max_degree):
            print(f"Cannot add edge: node {node1} already has {self.max_degree} edges")
            return False
        if(self.neighbors_counter[node2] == self.max_degree):
            print(f"Cannot add edge: node {node2} already has {self.max_degree} edges")
            return False
        if self.verbose:
            print("boundary vector is: ", boundary_vector, self.neighbors_counter)

        # Now we will add the bond data by appending to the various lists
        # Convert to dataclass in the future to be more robust
        if any(boundary_vector):
            # This bond crosses the boundary
            self.bond_is_across_boundary.append(True)
            self.crosses_boundaries += 1
        else:
            # This bond does not cross the boundary
            self.bond_is_across_boundary.append(False)
            self.n_internal_edges += 1
        self.simple_edges_list.append([node1, node2])
        # this list constructor makes the next line work regardless of whether boundary_vector
        # is passed as a Python list or numpy array.
        self.simple_boundary_crossing.append([x for x in boundary_vector])

        # Next, add bond data to some node-based lists
        self.neighbors[node1, self.neighbors_counter[node1]] = node2
        self.edges_list[node1, self.neighbors_counter[node1]] = self.edges_counter
        self.boundary_crossing[node1, self.neighbors_counter[node1], :] = boundary_vector
        self.neighbors_counter[node1] += 1

        if node2 != node1:
            # do the same inverting nodes
            self.neighbors[node2, self.neighbors_counter[node2]] = node1
            self.edges_list[node2, self.neighbors_counter[node2]] = self.edges_counter
            self.boundary_crossing[node2, self.neighbors_counter[node2], :] = \
                -np.asarray(boundary_vector)  # fix this
            self.neighbors_counter[node2] += 1

        # update edges_counter
        self.edges_counter += 1
        return True

    def get_number_of_neighbors(self, i):
        """
        Get the number of bonds of node i

        Args:
            i (int): node number

        Returns:
            int: The number of edges (bonds) involving node i
        """

        return np.count_nonzero(self.neighbors[i, :] != -1)

    def get_neighbors(self, i, padded=True):
        """
        Get array of neighbor indices of node i.

        Args:
            i (int): node number
            padded (bool, optional): If true (the default), the list will be padded
                with values -1 to the value of maximum_neighbors_per_node passed to
                the constructor. Otherwise the length will be the number of neighbors of i.

        Returns:
            :obj:`list` of int: list of neighbors (along bonds) of node i

        """
        neighbors_of_i = self.neighbors[i, :]
        if padded:
            return neighbors_of_i
        if not padded:
            stripped_neighbors_of_i = neighbors_of_i[neighbors_of_i != -1]
            return stripped_neighbors_of_i

    def get_neighbor(self, i, n_index):
        """
        Get n_index'th neighbor of node i

        Args:
            i (int): node number
            n_index (int): the position of the neighbor in the list returned by get_neighbors()

        Returns:
            int: The index of that neighbor (the value of get_neighbors(i)[n_index])

        """
        return self.neighbors[i, n_index]

    def get_edges(self, i):
        """ Returns edge corresponding to """
        return self.edges_list[i, :]

    def get_number_of_edges(self):
        """ Returns number of edges """
        return self.edges_counter

    def get_edge(self, node1, k):
        """
        Get the edge number of the k'th edge of node1.

        Args:
            node1 (int): node number
            k (int): the position of the edge in the list of edges of node1.

        Returns:
            int: The edge number of that edge (to be used as an index in arrays of edge properties)
        """
        return self.edges_list[node1, k]

    def get_boundary_crossing(self, i, n_index):
        """
        Get the boundary crossing vector of the n_index'th neighbor of node i.

        Args:
            i (int): node number
            n_index (int): index of neighbor in neighbor list of i

        Returns:
            int[3]: The vector of integers [bx, by, bz] denoting the number of
            times each boundary is crossed by this edge.
        """
        return self.boundary_crossing[i, n_index, :]

    def needs_reducing(self):
        """
        Determine if the network could be reduced using internal connected component decomposition.

        LoopFinder will perform this reduction automatically so there will not usually be a need
        for the user to call this function themselves.

        Returns:
            bool: True if the network has any edges that do not cross any boundary.
        """
        return (self.n_internal_edges > 0)

    def __label_component(self, start,  current_label, labels):
        """
        Label the entire connected component to which node start belongs with label current_label.

        Args:
            start (int): node number
            current_label (int): the label for this connected component
            labels (:obj:`numpy.ndarray`): numpy array (dtype=int) containing the label
                of each node/vertex. This array is updated by this recursive routine.
        """
        labels[start] = current_label
        # assert len(self.neighbors[1]) == self.max_degree  # using this to simplify next line
        for index in range(self.max_degree):
            neigh = self.neighbors[start, index]
            edge_is_outside = self.bond_is_across_boundary[self.edges_list[start, index]]
            if not edge_is_outside and neigh != -1 and labels[neigh] == -1:  # empty: not connected
                self.__label_component(neigh, current_label, labels)

    def decompose(self):
        """
        Obtain the cluster decomposition of the network (using only internal bonds).

        Returns:
            Tuple[:obj:`List` of int, int]: A list with the cluster ID of each node
            and the number of clusters
        """
        # initiate with first cluster label
        current_label = 0
        # initialise list of labels (-1 == unlabeled)
        labels = -1 * np.ones(self.number_of_nodes, dtype=int)
        for node in range(self.number_of_nodes):
            if labels[node] == -1:
                # this node is still unlabeled. Start recursion to label its connected component.
                self.__label_component(node, current_label, labels)
                current_label += 1
        n_labels = np.amax(labels) + 1
        # At this point, all elements of labels are >=0 and < n_labels
        # assert np.amin(labels) == 0
        return labels, n_labels

    def nodeid_to_clusterid(self, clusterlabels):
        if self.verbose:
            print(self.simple_edges_list, self.simple_boundary_crossing)
        reduced_network = []
        for edge, newnodes in enumerate(self.simple_edges_list):
            if self.bond_is_across_boundary[edge]:
                connected_labels = clusterlabels[newnodes].tolist()
                # parameters are now in two python lists which we can concatenate
                # using the addition operator
                edgedata = connected_labels + self.simple_boundary_crossing[edge]
                reduced_network.append(edgedata)
        reduced_network = np.asarray(reduced_network)
        # The main reason to collect the edge data in a 5-column numpy array is
        # that we can now use np.unique to get rid of duplicate edges
        reduced_network = np.unique(reduced_network, axis=0)
        return reduced_network

    # methods that reduce network
    def get_reduced_network(self):
        """
        Generate the reduced network with identical boundary crossing properties
        but no internal edges.

        Returns:
            :obj:`PeriodicNetwork`: The reduced network
        """
        if self.n_internal_edges == 0:
            # The network has no edges that do not cross the boundary.
            # Therefore the result of reducing would be identical to the current network.
            # Returning original network.
            return self

        # First find a cluster decomposition of the network excluding the boundary-crossing edges
        clusterlabels, n_labels = self.decompose()
        if self.verbose:
            print("labels:", clusterlabels, n_labels)

        # each cluster in this decomposition becomes a node in the reduced network
        # the boundary crossing edges will be put back in, now with cluster IDs
        # rather than node IDs to indicate what they connect
        reduced_network_list = self.nodeid_to_clusterid(clusterlabels)
        if self.verbose:
            print("reduced network list:", reduced_network_list)
        if len(reduced_network_list) == 0:
            # return network without any bonds, but with the proper number of nodes
            return PeriodicNetwork(n_labels, 1)
        # next line finds the number of occurrences of the most-occuring number in the
        # first two columns of reduced_network_list
        largest_functionality = np.max(np.unique(reduced_network_list[:, 0:2],
                                                 return_counts=True)[1])
        # construct new PeriodicNetwork object with n_labels nodes
        reduced_network = PeriodicNetwork(n_labels, largest_functionality, verbose=self.verbose)
        # add the edges of the reduced network
        for edge_info in reduced_network_list:
            if self.verbose:
                print(edge_info)
            reduced_network.add_edge(edge_info[0], edge_info[1], edge_info[2:])
        return reduced_network
