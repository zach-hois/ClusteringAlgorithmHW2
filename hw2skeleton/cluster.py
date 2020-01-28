from .utils import Atom, Residue, ActiveSite
import numpy as np
#import matplotlib.pyplot as plt

def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """
    similarity = 0.0
    print("THIS LINE IS TO TEST THAT COMPUTE SIMILARITY IS RUNNING")
    

    # Fill in your code here!
    
    
    return similarity


def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    # Fill in your code here!
    compute_similarity(active_sites[0], active_sites[1])

    k = 5 #arbitrary beginning value
    centroids = {
        i+1: [np.random.randint(0,130), np.random.randint(0,130)]
        for i in range(k)
    }
    
    #def assign(active_sites, centroids):
    print(active_sites[0].name)
    print(active_sites[0].residues)
    print(active_sites[0].residues[0].atoms)


    return []


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []
