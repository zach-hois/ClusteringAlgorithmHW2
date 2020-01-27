# HW2 Skeleton

[![Build
Status](https://travis-ci.org/zach-hois/ClusteringAlgorithmHW2.svg?branch=master)](https://travis-ci.org/zach-hois/ClusteringAlgorithmHW2)

Skeleton for clustering project.

## assignment

1. Implement a similarity metric
2. Implement a clustering method based on a partitioning algorithm
3. Implement a clustering method based on a hierarchical algorithm
4. Answer the questions given in the homework assignment


## structure

The main file that you will need to modify is `cluster.py` and the corresponding `test_cluster.py`. `utils.py` contains helpful classes that you can use to represent Active Sites. `io.py` contains some reading and writing files for interacting with PDB files and writing out cluster info.

```
.
├── README.md
├── data
│   ...
├── hw2skeleton
│   ├── __init__.py
│   ├── __main__.py
│   ├── cluster.py
│   ├── io.py
│   └── utils.py
└── test
    ├── test_cluster.py
    └── test_io.py
```

## usage

To use the package, first run

```
conda install --yes --file requirements.txt
```

to install all the dependencies in `requirements.txt`. Then the package's
main function (located in `hw2skeleton/__main__.py`) can be run as
follows

```
python -m hw2skeleton -P data test.txt
```

## testing

Testing is as simple as running

```
python -m pytest
```

from the root directory of this project.


## contributors

Original design by Scott Pegg. Refactored and updated by Tamas Nagy.


## Notes from 1/27/2020

End goal objectives: cluster partitioning, cluster hierarchy
Class Atom (line in .pdb), Active Site (.pdb file)
pdb file is a protein active site, with each residue being a block 
Utils file is built to read in the .pdb
Goal: design a similarity metric between atoms (Think about the biological meaning!)

idea: bin by percentages of oxygen, carbon, and nitrogen

only compare over the minimum length and add extra length to the comparison 

qualities we need in a distance metric: triangle inequality (a2b + b2c > a2c)
symmetric (a2b = b2a)
distance a2a = 0

after done clustering, figure out a way to visualize it. (reduce the dimensions)

how can we compare the clusterings?
silhouette score (how cohesive is structure) **
	which is better and why? (plot)
Tacard (how similar?)


ActiveSite[0].residues
.type
.number
.atoms
take an active site, unzip all of the structure, align based on coordinates and type, make dataframe for active site, then compare based on stuff

gave clusting algorithm the coordinates and compute similarity
OR
use the similarity metric as the clustering input

