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
    #print("THIS LINE IS TO TEST THAT COMPUTE SIMILARITY IS RUNNING")
    
    ## my GOAL is to have my similarity matrix compare how many amino acids
    ## residues of each protein have in common and normalize it to total they could have in common
    ##  of the residues to account for the varying length of the residues


    AminoAcids = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                  'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 
                  'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 
                  'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'} 
                  #shorten them to their 1 letter code to make this easier


    def AminoAcidSequence(active_site):
        return [AminoAcids[i.type[0:3]] for i in active_site.residues]
        #convert the input sequences so they are more workable

    seqA = AminoAcidSequence(site_a) #take the inputs to the compute similarity function
    seqB = AminoAcidSequence(site_b) #and convert their sequences using new function
    observedRes = set(seqA + seqB) #these are the residues that are observed in the inputs
    total = 0 #initializations
    common = 0
    for res in observedRes: #for each amino acid residue in observed residues
        inA = seqA.count(res) #count the # of occurrences of the residue in each seq
        inB = seqB.count(res)
        if inA < inB: #res more common in B than A it is added to total and common calc for A
            common += inA
            total += inB
        else: #opposite
            common += inB
            total += inA
    similarity = (common / total)*100 # therefore it returns the # in common vs max possible in common, example:
    # AND vs ANA
    # A - 1 in common out of 2 max
    # N - 1 in common out of 1
    # D - 0 in common out of 1
    # (1 + 1 + 0) / (2 + 1 + 1) = 2/4 = similarity of 0.5 or 50%
    #print("SIMILARITY IS:", similarity)
    return similarity

def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    My k-means clusering algorithm was adapted from: http://benalexkeen.com/k-means-clustering-in-python/
    """
    # Fill in your code here!
    SMatrix = np.asarray([ #compile the similarity matrix by iterating all of the 
                np.asarray([  #protein structures over each other
                    compute_similarity(i,j) for i in active_sites]) 
                         for j in active_sites])


    df = pd.DataFrame({'x': np.asarray(np.mean(SMatrix, axis = 0)),
                       'y': np.asarray(np.mean(SMatrix, axis = 1))})

    k = 3 #arbitrary beginning value
    centroids = {
        i+1: [np.random.randint(0,100), np.random.randint(0,100)] #select a value from dataframe for centroid
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
        centroids = update(centroids) #previously defined update function realigns centroids
        df = assignment(df, centroids) #new dataframe is created with new centroids
        if closest_centroids.equals(df['closest']): 
            break #end the while loop if there are no newly assigned points 

    #example plot for what would be shown in 2d
    fig = plt.figure(figsize=(5, 5))
    plt.scatter(df['x'], df['y'])
    for i in centroids.keys():
         plt.scatter(*centroids[i])
    plt.xlim(0, 80)
    plt.ylim(0, 80)
    #plt.show()



    #####messing around with comparing the active site residues
    """
    for i in range(1):
        #print(len(active_sites[i].residues))
        for j in range(130):
            a = (active_sites[i].residues[i].atoms)
            #print(a)

            b = (active_sites[j].residues[i].atoms)
            #print(b)
            c = [item for item in b if item in a]
            print(c)

    return []
    """ 



def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []
