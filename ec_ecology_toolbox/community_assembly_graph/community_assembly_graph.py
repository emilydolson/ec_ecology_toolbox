from copy import deepcopy
from queue import PriorityQueue
from numpy import dot
from numpy.linalg import matrix_power


def random_walk_matrix_multiplication(graph, start_node, num_iterations=100):
    """
    Function used to calculate probabilities of ending up at a node i.e. node passing probability
    Implements iterative random walk theory using matrix multiplication
    Refer: https://chih-ling-hsu.github.io/2020/05/20/random-walks

    :param graph: manipulated graph with sink nodes on which the random walk algorithm is run
    :type graph: dict {Node: dict {Node: float}}
    :param start_node: the start node of the random walk
    :type start_node: user-defined Node class
    :param num_iterations: number of iterations to run for the random walk. default: 100
    :type num_iterations: int
    :return: tuple including final passing probabilities after convergence and the node directory for index matching
    :rtype: tuple (list [float], dict {Node: int})
    """

    # A directory for tracking the index position in the matrix for each node
    node_dir = {}
    for idx, node in enumerate(graph):
        node_dir[node] = idx

    # initializing the matrix whose length is equal to the number of nodes in the graph
    matrix = [None]*len(graph)

    # list 'v' holds starting probability of all nodes. 0 for all nodes, except, start node will be 1
    v = [0]*len(graph)
    v[node_dir[start_node]] = 1

    for node, edges in graph.items():
        # Creating the row for each node with its respective outgoing edges
        row = [0]*len(graph)
        for adjacent_node, weight in edges.items():
            pos = node_dir[adjacent_node]
            row[pos] = weight

        # Adding the row to matrix into it's correct index position
        # noinspection PyTypeChecker
        matrix[node_dir[node]] = row

    # algorithm to calculate final passing probability
    # dot product of 'v' and the matrix raised to the power of the number of iterations
    rw_probs = dot(v, matrix_power(matrix, num_iterations))

    return rw_probs, node_dir


def update_priorities(graph, leaves, start_node, node_passing_p, num_iterations):
    """
    Updates the probabilities of the discovered nodes in the priority queue
    Called when new paths to visited nodes or cycles are discovered in the explored graph

    :param graph: explored graph with all visited nodes and discovered paths
    :type graph: dict {Node: dict {Node: float}}
    :param leaves: all leaf nodes that have been discovered and are in the priority queue
    :type leaves: set {Node}
    :param start_node: the starting node of exploration
    :type start_node: user-defined Node class
    :param node_passing_p: holds the updated passing probability of all discovered nodes
    :type node_passing_p: dict {Node: float}
    :param num_iterations: number of iterations to run for the random walk. default: 100
    :type num_iterations: int
    :return: new priority queue with updated probabilities of all leaf nodes
    :rtype: queue.PriorityQueue -> tuple (float, Node)
    """
    # New priority queue to hold nodes with updated priorities
    new_pq = PriorityQueue()

    for curr_node in leaves:
        # Creates sink nodes with self edges for all leaf nodes
        graph[curr_node] = {curr_node: 1.0}

    # function call to calculate passing probability of created sink nodes
    pr, idx = random_walk_matrix_multiplication(graph, start_node, num_iterations)

    for curr_node in leaves:
        # Removes the created sink nodes
        del graph[curr_node]

        # retrieves new probabilities and updates the new priority queue
        new_p = pr[idx[curr_node]]

        new_pq.put(tuple([-1*new_p, curr_node]))
        node_passing_p[curr_node] = new_p

    return new_pq


def CAG(start_node, user_func, num_nodes=10, normalize_w_self_edges=False, num_iterations=100):
    """
    Main function that is used to explore the graph based on user function and start node
    Graph is explored based on continuous probability of reaching a node signified by edge weights
    Probability of reaching a node calculated from user given start node
    User defines the bumber of nodes to be discoered and returned

    :param start_node: starting node of exploration
    :type start_node: user-defined Node class
    :param user_func: user defined function the returns all outgoing edges with weights from input node
    :type user_func: function
    :param num_nodes: number of nodes to explore
    :type num_nodes: int
    :param normalize_w_self_edges: normalize nodes by including self edges if True. default: False
    :type normalize_w_self_edges: bool
    :param num_iterations: number of iterations to run for random walk. default: 100
    :type num_iterations: int
    :return: explored graph until required number of nodes reached or full graph is explored
    :rtype: dict {Node: dict {Node: float}}
    """

    # priority queue that holds the discovered nodes based on the probability of reaching it
    p_queue = PriorityQueue()

    # set contains all nodes that have been visited and explored
    visited = set()

    # set contains all nodes discovered but yet to be explored
    discovered = set()

    # set contains all active entries into the priority queue
    # to keep a check on duplicate entries of nodes into the queue
    # only used for optimization
    active_pq_nodes = set()

    # dictionary that contains the explored portion of the graph
    explored_graph = {}

    # same as explored graph
    # contains self edges of nodes while normalizing
    # used only if parameter normalize_w_self_edges is True
    normalized_explored_graph = {}

    # initializing the exploration with the start node
    # queue is min_heap so passing probability negated
    p_queue.put(tuple([-1, start_node]))
    discovered.add(start_node)
    active_pq_nodes.add(start_node)
    node_passing_p = {start_node: 1}

    # traversing the graph until end condition is reached
    while len(visited) < num_nodes:
        if p_queue.empty():
            # graph completely explored. No new nodes to explore
            break

        if discovered == visited:
            # graph completely explored. No new nodes to explore
            break

        # Retrieving node of the highest priority i.e. most likely node to reach next
        curr_p, curr_node = p_queue.get()

        if curr_node not in active_pq_nodes:
            # optimization
            # node has already been handled. most likely duplicate entry in th queue
            continue

        # node will now be explored. removing prom active priority queue nodes
        active_pq_nodes.remove(curr_node)

        # priorities are negated when added since p_queue is based on min_heap
        curr_p *= -1

        if curr_node in visited:
            # new path to an already explored node is found
            # this alters the priorities of nodes already in the queue i.e. leaf nodes
            # therefore, function to update their priorities and the queue is called

            curr_leaves = discovered - visited
            p_queue = update_priorities(explored_graph, curr_leaves, start_node, node_passing_p, num_iterations)

            # active nodes in the new queue are only the leaf nodes
            active_pq_nodes = curr_leaves

            # update handled. proceed to next node in the queue
            continue

        # exploring current node by calling user function
        # which returns all outgoing edges from the nodes with edge weights
        edges = user_func(curr_node)

        # updated the explored graph dict with newly explored edges
        explored_graph[curr_node] = edges

        if normalize_w_self_edges:
            # normalizing node by adding self edge
            edge_weight_sum = round(sum(edges.values()), 4)

            # if sum of edge weights already add up to 1, no changes required
            # if not, add a self edge to ensure sum equals to 1
            # given condition, sum of edge weights will never exceed 1
            if edge_weight_sum < 1:
                normalized_edges = deepcopy(edges)
                normalized_edges[curr_node] = round(1 - edge_weight_sum, 4)
                normalized_explored_graph[curr_node] = normalized_edges
            else:
                normalized_explored_graph[curr_node] = edges

        # handling the discovery of all adjacent nodes
        for adjacent_node, edge_weight in edges.items():
            # calculating the passing probability of the adjacent_node
            adjacent_node_p = edge_weight * curr_p

            # Adding 1 to the passing probability of adjacent node if it is already explored
            # To ensure these nodes will be retrieved from the priority queue first and priorities are updated
            if adjacent_node in visited:
                adjacent_node_p += 1

            # updating the adjacent nodes passing probability in the node_passing_p dict
            if adjacent_node in node_passing_p:
                node_passing_p[adjacent_node] += adjacent_node_p
            else:
                node_passing_p[adjacent_node] = adjacent_node_p

            # Adding the adjacent node into the priority queue
            new_p = node_passing_p[adjacent_node]
            p_queue.put(tuple([-1*new_p, adjacent_node]))

            # updating the discovery of the adjacent node
            discovered.add(adjacent_node)
            active_pq_nodes.add(adjacent_node)

        # adding the current node to the set of visited/explored nodes
        visited.add(curr_node)

    # Returning the explored graph
    if normalize_w_self_edges:
        return normalized_explored_graph
    return explored_graph


if __name__ == '__main__':
    import time
    from ec_ecology_toolbox.community_assembly_graph.custom_user_functions import custom_user_func_3
    from ec_ecology_toolbox.community_assembly_graph.benchmarking.benchmarking import benchmarking_func2
    from ec_ecology_toolbox.community_assembly_graph.example_nodes import Chr_Node

    start_time = time.time()
    N = Chr_Node(1)
    graph = CAG(N, custom_user_func_3, 20, False, 250)
    end_time = time.time()

    print(graph)
    print("Run Time:", end_time - start_time)

    start_time = time.time()
    N = Chr_Node(1)
    graph = CAG(N, benchmarking_func2, 5, False, 250)
    end_time = time.time()

    print(graph)
    print("Run Time:", end_time - start_time)




