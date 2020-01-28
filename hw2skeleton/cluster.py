from .utils import Atom, Residue, ActiveSite
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
import copy

def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """
    similarity = 0.0


    print("THIS LINE IS TO TEST THAT COMPUTE SIMILARITY IS RUNNING")
    x = site_a.residues[0].atoms[0].coords
    y = site_b.residues[0].atoms[0].coords

    #planning on computing all pairs and then summing up

    distance = math.sqrt(sum([(a - b) ** 2 for a, b in zip(x, y)]))
    print("EUCLIDEAN DISTANCE:",distance)
    
    
    
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


    df = pd.DataFrame({ #placeholder, DF will soon be similarity matrix
    'x': [12, 20, 28, 18, 29, 33, 24, 45, 45, 52, 51, 52, 55, 53, 55, 61, 64, 69, 72],
    'y': [39, 36, 30, 52, 54, 46, 55, 59, 63, 70, 66, 63, 58, 23, 14, 8, 19, 7, 24]
    })
    k = 3 #arbitrary beginning value
    centroids = {
        i+1: [np.random.randint(0,80), np.random.randint(0,80)] #select a value from dataframe for centroid
        for i in range(k)
    } #randomly select three of the points to be the arbitrary centroids to start algorithm
    

    ##### BEGIN ASSIGNMENT TO CENTROID #####
    def assignment(df, centroids): #we will assign each value in df to a centroid
        for i in centroids.keys():
        #euclidian distance calculation
            df['distance_from_{}'.format(i)] = ( #distance of each point from the centroid
              np.sqrt((df['x'] - centroids[i][0]) ** 2 + (df['y'] - centroids[i][1]) ** 2)
        )
        centroid_distance_cols = ['distance_from_{}'.format(i) for i in centroids.keys()]
        df['closest'] = df.loc[:, centroid_distance_cols].idxmin(axis=1)
        df['closest'] = df['closest'].map(lambda x: int(x.lstrip('distance_from_')))
       
        return df

    df = assignment(df, centroids) #initialize the new dataframe with each sorted to a centroid

    def update(k): #change the value of the centroids to be more representative of the cluster
        for i in centroids.keys():
            centroids[i][0] = np.mean(df[df['closest'] == i]['x'])
            centroids[i][1] = np.mean(df[df['closest'] == i]['y'])
        return k #k means takes the average x and y and moves the centroid closer 

    centroids = update(centroids) #update the previous arbitrary centroids to be more accurate

    ##### REPEAT ASSIGNMENT TO CENTROID #####

    df = assignment(df, centroids)

    ##### REPEAT UNTIL THERE ARE NO MORE CHANGES TO THE CLUSTERS #####
    while True: 
        closest_centroids = df['closest'].copy(deep=True)
        centroids = update(centroids)
        df = assignment(df, centroids)
        if closest_centroids.equals(df['closest']):
            break


    fig = plt.figure(figsize=(5, 5))
    plt.scatter(df['x'], df['y'])
    for i in centroids.keys():
         plt.scatter(*centroids[i])
    plt.xlim(0, 80)
    plt.ylim(0, 80)
    plt.show()

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
