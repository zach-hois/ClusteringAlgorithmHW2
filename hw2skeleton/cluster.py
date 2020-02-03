from .utils import Atom, Residue, ActiveSite
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
import copy
import umap
import seaborn as sns


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
    df = 100 - pd.DataFrame(
            np.asarray([ #compile the similarity matrix by iterating all of the 
                np.asarray([  #protein structures over each other
                    compute_similarity(i,j) for i in active_sites]) 
                         for j in active_sites]))

    k = 3 #arbitrary beginning value
    centroids = {
        i+1: [i*25] #select a value from dataframe for centroid
        for i in range(k)
    } #select three of the points to be the arbitrary centroids to start algorithm
    #originally this was going to be [np.random.randint(0, 100)] and i was happy to learn that
    #it correctly clustered the points together each time, however the centroid they were 
    #assigned to was not always consistent due to the random draw so the assert statement
    #was difficult to right. as a solution I picked values 25, 50, 75

    ##### BEGIN ASSIGNMENT TO CENTROID #####
    def assignment(df, centroids): #we will assign each value in df to a centroid
        for i in centroids.keys():
        #euclidian distance calculation
            for j in range(len(df.index)):
                df['distance_from_{}'.format(i)] = ( #distance of each point from the centroid
                    np.sqrt((centroids[i][0] - df[j] ) ** 2)
        )
        centroid_distance_cols = ['distance_from_{}'.format(i) for i in centroids.keys()]
        df['closest'] = df.loc[:, centroid_distance_cols].idxmin(axis=1)
        df['closest'] = df['closest'].map(lambda x: int(x.lstrip('distance_from_')))
        return df

    df = assignment(df, centroids) #initialize the new dataframe with each sorted to a centroid

    def update(k): #change the value of the centroids to be more representative of the cluster
        for i in centroids.keys():
            for j in range(len(df.index)):
                centroids[i][0] = np.mean(df[df['closest'] == i][j])
        return k #k means takes the average x and y and moves the centroid closer 

    centroids = update(centroids) #update the previous arbitrary centroids to be more accurate
    ##### REPEAT ASSIGNMENT TO CENTROID #####

    df = assignment(df, centroids)

    ##### REPEAT UNTIL THERE ARE NO MORE CHANGES TO THE CLUSTERS #####
    while True: 
        closest_centroids = df['closest'].copy(deep=True)
        centroids = update(centroids) #previously defined update function realigns centroids
        df = assignment(df, centroids) #new dataframe is created with new centroids
        #print(df.head())
        if closest_centroids.equals(df['closest']): 
            break #end the while loop if there are no newly assigned points 

    #example plot for what would be shown in 2d
    """
    reducer = umap.UMAP() #UMAP to make it easy to visualize
    embedding = reducer.fit_transform(df)


    plt.scatter(embedding[:, 0], embedding[:, 1], c=df['closest']) #make the plot
    plt.title('UMAP projection of df', fontsize=24)
    #for i in centroids.keys():
     #   plt.scatter(*centroids[i], *centroids[i])
    plt.show()
    """
    clusterList = list(range(0, df.shape[0])) #initialize output 

    for i in range(len(df.index)):
        for j in range(k):
            clusterList[j] = list(df.index[df.closest == j+1])
    #save the row # (residue) to the correct column in the output corresponding to its centroid
    #this iteration is similar to the one from the initialization of the centroid

    return clusterList #output is a list with proteins assigned to clusters

def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """
    # Fill in your code here!
    def flatten(flatList): #this is a solution in order to flatten a nested list
        if isinstance(flatList, list): #its used in the clustering loop
            keepGoing = True #go more for list
            while keepGoing:
                lll = []
                keepGoing = False #finite 
                for x in flatList:
                    if isinstance(x, list):
                        lll += x
                        keepGoing = True
                    else:
                        lll.append(x)
                flatList = lll
        return flatList #returns a flat list so i can swap the clusters later on 

    df = 100 - np.asarray([ #compile the similarity matrix by iterating all of the 
                        np.asarray([  #protein structures over each other
                            compute_similarity(i,j) for i in active_sites]) 
                                for j in active_sites]) #same as before

    inCluster = np.array(list(range(0, df.shape[0]))) #which cluster will the point be in

    clusterList = list(range(0, df.shape[0])) #all clusters, will be out output

    for i in range(0, df.shape[0]): #format the df to be useable in future
        df[i][i] = np.inf 
    
    while (len(clusterList)) > 2: #run while there are still points to be assigned
        near = np.argmin(df)
 
        nearI, nearJ = int(near / df.shape[0]), near % df.shape[0]
        if inCluster[nearI] != inCluster[nearJ]:  #this will combine clusters
            old = clusterList[inCluster[nearJ]] #here we join them if this is triggered
            oldNearJ = inCluster[nearJ] #hold this but itll be overwritten 
            clusterList[inCluster[nearI]] = [
                clusterList[inCluster[nearJ]], clusterList[inCluster[nearI]]]

            inCluster[flatten(clusterList[inCluster[nearJ]])] = inCluster[nearI] #set the new location of a member

            clusterList.remove(old) #remove old value 

            for i in range(0, len(inCluster)):
                if inCluster[i] > oldNearJ:
                    inCluster[i] -= 1   #work out way down clusters when one is created
        df[nearI][nearJ] = np.inf #format similar to before this began 
        df[nearJ][nearI] = np.inf

    return clusterList #finish and return

    """
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(df)

    plt.scatter(embedding[:, 0], embedding[:, 1])
    plt.title('UMAP projection of df', fontsize=24)
    #for i in centroids.keys():
     #   plt.scatter(*centroids[i], *centroids[i])
    plt.show()
    """

def comparison(clusterList, active_sites):
    """
    I am going to use a comparison method that averages the intracluster similarities
    and the intercluster dissimilarities to make a score that will be returned and compared
    for both of the algorithms
    """
    index = {active_sites[i]: i for i in range(0,len(active_sites))} #for any amount of active sites used

    df = 100 - np.asarray([ #compile the similarity matrix by iterating all of the 
                        np.asarray([  #protein structures over each other
                            compute_similarity(i,j) for i in active_sites]) 
                                for j in active_sites])

    n = len(clusterList)
    intraCluster = [(np.sum([[df[index[u]][index[v]] for u in x] for v in x
                         ])) / (len(x) ** 2) for x in clusterList] #intracluster similarity computation

    score = 0 #initialization
    for i in range(0, len(clusterList)): #for all submitted clusters
        for j in range(i, len(clusterList)): #nested ticker
            score += np.sqrt(intraCluster[i] * intraCluster[j]) * ( #* intercluster dissimilarity
                1.0 - np.sum([[df[index[u]][index[v]] for u in clusterList[i]] for v in clusterList[j]]
                                ) / (len(clusterList[i]) + len(clusterList[j]))) #add the calculated score and normalize
    
    score = (1 / (n * (n-1))) * score #compute the comparison value normalized to the amount of clusters

    return score 



