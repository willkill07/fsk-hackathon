#!/usr/bin/env python

'''
 Robert Searles
 Department of Computer and Information Sciences
 University of Delaware

 FSK.py
 Python implementation of the Fast Subtree Kernel.
'''

#################################################################
###    Includes
#################################################################

# Needed for system operations
import sys, os, shutil

# Parsing command line arguments
import argparse

# Multiprocessing for threading
import multiprocessing, functools

# Similarity Threshold
DELTA = 0.35

#################################################################
###    Class/Function Declarations
#################################################################

def LoadData(file_list, encoded_graph_list, num_graphs):
    with open(file_list, 'r') as files:
        fl = files.readlines()

        # Open each graph
        for i in xrange(num_graphs):
            curr_file = fl[i]
            with open(os.path.abspath(curr_file[:-1]), 'r') as current_graph:
                encoded_trees = list()

                # Gather encoded feature vectors
                for line in current_graph:
                    feature_vector = map(int, line.split())
                    encoded_trees.append(feature_vector)
                encoded_graph_list.append(encoded_trees)
                current_graph.close()

        files.close()

def SubtreeSimilarity(fv_1, fv_2):
    normalized_distance = 0.0
    fv_length = len(fv_1)

    for k in xrange(fv_length):
        normalized_distance += float(abs(float(fv_1[k]) - float(fv_2[k])) / max(1.0, max(float(fv_1[k]), float(fv_2[k]))))

    normalized_distance = normalized_distance / fv_length
    # print normalized_distance
    return normalized_distance

def PairwiseGraphSimilarity(i, encoded_graph_list, num_graphs, graph_j):
    global_similarity = 0.0

    graph_A = encoded_graph_list[i]
    graph_B = graph_j#encoded_graph_list[j]

    total_subtree_comparisons = len(graph_A) + len(graph_B)

    # Compare graph A to graph B
    for feature_vector_1 in graph_A:
        curr_min = 1.0
        for feature_vector_2 in graph_B:
            curr_min = min(curr_min, SubtreeSimilarity(feature_vector_1, feature_vector_2))
        if curr_min < DELTA:
            global_similarity += 1.0


    # Compare graph B to graph A
    for feature_vector_1 in graph_B:
        curr_min = 1.0
        for feature_vector_2 in graph_A:
            curr_min = min(curr_min, SubtreeSimilarity(feature_vector_1, feature_vector_2))
        if curr_min < DELTA:
            global_similarity += 1.0

    normalized_similarity = 1.0 - (global_similarity / total_subtree_comparisons)
    return normalized_similarity

def ComputeGraphSimilarity(encoded_graph_list, num_graphs):
    kernel_matrix = list()

    for i in xrange(0, num_graphs):
        print "--> Comparing graph " + str(i) + " to the dataset"

        # Pass arguments to function
        func = functools.partial(PairwiseGraphSimilarity, i, encoded_graph_list, num_graphs)

        # Create multiprocessing pool
        pool = multiprocessing.Pool(NUM_THREADS)

        # Map over final argument
        multi_return = pool.map(func, encoded_graph_list)

        # Append result to kernel matrix
        kernel_matrix.append(multi_return)
        
        # for j in xrange(0, num_graphs, NUM_THREADS):
        #     kernel_matrix[i].append(PairwiseGraphSimilarity(i, j, encoded_graph_list, num_graphs))
#            PairwiseGraphSimilarity(i, j, encoded_graph_list, num_graphs, kernel_matrix)
#        kernel_matrix.append(POOL)

    for row in kernel_matrix:
        print " ".join(str("%.2f" % x) for x in row)
    return

#################################################################
###    Script Execution
#################################################################

def main():
    print "--> Loading " + str(NUM_GRAPHS) + " Graphs"

    encoded_graph_list = list()
    LoadData(FILE_LIST, encoded_graph_list, NUM_GRAPHS)

    print "--> Computing Graph Similarities using " + str(NUM_THREADS) + " threads"
    ComputeGraphSimilarity(encoded_graph_list, NUM_GRAPHS)

if __name__ == "__main__":
    # Used to parse options
    parser = argparse.ArgumentParser()

    # Add arguments
    parser.add_argument("file_list", help="File contatining the list of encoded graph files we will process.", type=str)
    parser.add_argument("-g", "--graphs", metavar="num_graphs", help="Number of graphs to process", type=int, default=10)
    parser.add_argument("-t", "--threads", metavar="num_threads", help="Number of CPU threads to use", type=int, default=1)

    # Parse arguments
    args = parser.parse_args()

    # Set globals
    FILE_LIST = args.file_list
    NUM_GRAPHS = args.graphs
    NUM_THREADS = args.threads
    
    # execute
    main()
