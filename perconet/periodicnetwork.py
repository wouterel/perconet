# This file is part of the perconet package
# (c) 2022 Eindhoven University of Technology
# Released under EUPL v1.2
# See LICENSE file for details
# Contributors:
# * Chiara Raffaelli
# * Wouter G. Ellenbroek
# Contributing: https://github.com/wouterel/perconet

import numpy as np

# TO DO: Make the documentation reflect that the package works for n-tori, not just 3-tori


class PeriodicNetwork:
    """Store and analyze the topology of a periodic net.

    Periodic nets are graphs embedded in a periodic topology. This class stores
    the topology of such a graph for the case of a `d`-dimensional periodic
    box (a `d`-torus). The dimensionaly defaults to 3 for use in contexts where
    the box represents physical space, but the :obj:`PeriodicNetwork` and
    :obj:`LoopFinder` classes work for arbitrary dimension.

    The class stores, for every edge in the graph, a d-dimensional vector
    indicating the boundary-wrapping properties of that edge. See
    :obj:`PeriodicNetwork.add_edge` for details. This
    information is then used by :obj:`LoopFinder` to determine the percolation
    properties.

    Args:
        n (int):
            The number of nodes of the graph.
        max_degree (int):
            The largest number of edges coming out of any node.
        verbose (bool, optional):
            Print debugging information to stdout. Defaults to False.
        dim (int, optional):
            Spatial dimension. Defaults to 3.
    """
    def __init__(self, n: int, max_degree=6, verbose=False, dim=3):
        if n < 1:
            raise ValueError("Number of nodes must be a positive integer.")
        self.number_of_nodes = n
        self.max_degree = max_degree
        self.verbose = verbose
        self.dimension = dim
        # allocate arrays for per-node info
        self.boundary_crossing = np.zeros((n, max_degree, dim), dtype=int)
        self.neighbors = -1 * np.ones((n, max_degree), dtype=int)
        self.edges_list = -1 * np.ones((n, max_degree), dtype=int)
        self.neighbors_counter = np.zeros(n, dtype=int)
        self.n_total_edges = 0  # keeps track of the total number of edges while building edge list
        self.simple_edges_list = []  # is this duplicate info?
        self.simple_boundary_crossing = []
        self.bond_is_across_boundary = []
        self.n_internal_edges = 0
        self.n_boundary_edges = 0

    def get_number_of_nodes(self):
        return self.number_of_nodes

    def get_dimension(self):
        return self.dimension

    def add_edge(self, node1: int, node2: int, boundary_vector):
        """
        Add an edge to the periodic network

        Args:
            node1 (int):
                The index of the first node of the pair that defines this edge.
                Valid values range from 0 up to (but not including) the number of nodes of
                the network (node indices are 0-based).
            node2 (int):
                The index of the second node of the pair that defines this edge.
                Valid values range from 0 up to (but not including) the number of nodes of
                the network (node indices are 0-based).
            boundary_vector (:obj:`List` of int):
                List (or numpy array) of integers denoting the number of times the edge
                wraps around each boundary, respectively. The length of this
                list must be equal to the dimensionality of the network (which defaults to 3 but
                can be overridden during initialization).
                The sign indicates the wrapping direction (e.g. (-1,0,0) indicates that the edge
                goes around the `x`-boundary in the negative `x`-direction
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
        if len(boundary_vector) != self.dimension:
            print(f"Incorrect boundary vector: {len(boundary_vector)} elements given " +
                  f"but {self.dimension} expected")
            return False

        # Now we will add the bond data by appending to the various lists
        # Convert to dataclass in the future to be more robust
        if any(boundary_vector):
            # This bond crosses the boundary
            self.bond_is_across_boundary.append(True)
            self.n_boundary_edges += 1
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
        self.edges_list[node1, self.neighbors_counter[node1]] = self.n_total_edges
        self.boundary_crossing[node1, self.neighbors_counter[node1], :] = boundary_vector
        self.neighbors_counter[node1] += 1

        if node2 != node1:
            # do the same inverting nodes
            self.neighbors[node2, self.neighbors_counter[node2]] = node1
            self.edges_list[node2, self.neighbors_counter[node2]] = self.n_total_edges
            self.boundary_crossing[node2, self.neighbors_counter[node2], :] = \
                -np.asarray(boundary_vector)  # fix this
            self.neighbors_counter[node2] += 1

        # update n_total_edges
        self.n_total_edges += 1
        return True

    def get_number_of_neighbors(self, node):
        """
        Get the number of bonds of node.

        Args:
            node (int): node number

        Returns:
            int: The number of edges (bonds) involving this node
        """

        return np.count_nonzero(self.neighbors[node, :] != -1)

    def get_neighbors(self, node, padded=True):
        """
        Get array of neighbor indices of node.

        Args:
            node (int): node number
            padded (bool, optional): If true (the default), the list will be padded
                with values -1 to the value of maximum_neighbors_per_node passed to
                the constructor. Otherwise the length will be the number of neighbors of i.

        Returns:
            :obj:`numpy.ndarray`: Numpy array (dtype=int) containing list of neighbors of node.
        """
        neighbors_of_node = self.neighbors[node, :]
        if padded:
            return neighbors_of_node
        else:
            stripped_neighbors_of_i = neighbors_of_node[neighbors_of_node != -1]
            return stripped_neighbors_of_i

    def get_neighbor(self, node, nb_index):
        """
        Get nb_index'th neighbor of node.

        Args:
            node (int): node number
            nb_index (int): index of neighbor in neighbor list of node

        Returns:
            int:
                The index of that neighbor (the value of `get_neighbors(i)[nb_index])`.
                A return value of -1 indicates that the edge does not exist.

        """
        return self.neighbors[node, nb_index]

    def get_edges(self, node, padded=True):
        """
        Get the list of edges linking to node.

        Args:
            node (int): The index of the node for which to return the edge list
            padded (bool, optional): If true (the default), the list will be padded
                with values -1 to the value of maximum_neighbors_per_node passed to
                the constructor. Otherwise the length will be the number of neighbors of node.

        Returns:
            :obj:`numpy.ndarray`:
            Numpy array (dtype=int) containing the edge numbers of all edges involving node.

        """
        edges_of_node = self.edges_list[node, :]
        if padded:
            return edges_of_node
        return edges_of_node[edges_of_node != -1]

    def get_number_of_edges(self):
        """
        Get total number of edges in network.

        Returns:
            int: Total number of edges (bonds) in the network
        """
        return self.n_total_edges

    def get_edge(self, node, nb_index):
        """
        Get the edge number of the nb_index'th edge of node.

        Args:
            node (int): node number
            nb_index (int): index of neighbor in neighbor list of node

        Returns:
            int:
                The edge number of that edge (to be used as an index in arrays of edge properties).
                A return value of -1 indicates that the edge does not exist.

        """
        return self.edges_list[node, nb_index]

    def get_boundary_crossing(self, node, nb_index):
        """
        Get the boundary crossing vector of the nb_index'th neighbor of node.

        Args:
            node (int): node number
            nb_index (int): index of neighbor in neighbor list of node

        Returns:
            :obj:`numpy.ndarray`: The list of integers denoting the number of
            times each boundary is crossed by this edge. Provided as
            a numpy array with length equal to the dimensionality of the network
            and dtype=int.
        """
        return self.boundary_crossing[node, nb_index, :]

    def crosses_boundaries(self):
        """
        Check if the network contains any edges that cross a boundary.

        Returns:
            bool: True if the network has any edges that cross a boundary.
        """
        return (self.n_boundary_edges > 0)

    def needs_reducing(self):
        """
        Determine if the network could be reduced using internal connected component decomposition.

        LoopFinder will perform this reduction automatically so there will not usually be a need
        for the user to call this function themselves.

        Returns:
            bool: True if the network has any edges that do not cross any boundary.
        """
        return (self.n_internal_edges > 0)

    def __label_component(self, start,  current_label, labels, internal_only=True):
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
            if neigh != -1 and labels[neigh] == -1:
                if not (internal_only and
                        self.bond_is_across_boundary[self.edges_list[start, index]]):
                    self.__label_component(neigh, current_label, labels)

    def decompose(self, internal_only=True):
        """
        Obtain the cluster decomposition of the network. This method is used
        by :obj:`LoopFinder` (using internal bonds only) to reduce the network for faster
        loop finding, but can also be used for generic cluster analysis.

        Args:
            internal_only (bool, optional): Defaults to True.
                If true, use only bonds that do not cross
                any boundary for the cluster decomposition.

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
                self.__label_component(node, current_label, labels, internal_only=internal_only)
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
        reduced_network = PeriodicNetwork(n_labels, largest_functionality,
                                          verbose=self.verbose, dim=self.dimension)
        # add the edges of the reduced network
        for edge_info in reduced_network_list:
            if self.verbose:
                print(edge_info)
            reduced_network.add_edge(edge_info[0], edge_info[1], edge_info[2:])
        return reduced_network
