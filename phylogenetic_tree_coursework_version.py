from Bio import SeqIO
import math
import numpy as np
from pathlib import Path

# Specify the file name
file_name = "king_penguins_truncated.fasta"
# Get the absolute path
abs_file_path = Path(__file__).parent / file_name
records = list(SeqIO.parse(abs_file_path, "fasta"))
n = len(records)
m = len(records[0])
distMatrix = np.zeros((n,n))
distMatrixJC = np.zeros((n,n))
for j in range(0,n):
    for k in range(0,j):
        w1 = records[j].seq
        w2 = records[k].seq
        distMatrix[j][k] = sum(c1 != c2 for c1,c2 in zip(w1,w2))/m
        distMatrixJC[j][k] = -3*math.log(1-4*distMatrix[j][k]/3)/4
pairs = []
ind = [list(range(n))]
clusters = [[x] for x in range(n)]
heights = [0]*(n-1)
nodesLeft = n
d = distMatrixJC
while nodesLeft>2:
    minDist = np.amin(d[np.nonzero(d)])
    newClust = np.where(d == minDist)
    a = newClust[0][0]
    b = newClust[1][0]

    pairs.append([a,b])

    newInd = [0]*(n-1)
    newInd[:2] = [x for x in ind[n-nodesLeft] if x != a and x!=b]
    ind.append(newInd)

    heights[n-nodesLeft] = d[a,b]/2
    clustDist = np.zeros(nodesLeft-1)
    clustDist[:nodesLeft-2] = [(len(clusters[a])*d[max(a,x),min(a,x)]+\
    len(clusters[b])*d[max(b,x),min(b,x)])/(len(clusters[a])+len(clusters[b]))\
    for x in range(nodesLeft) if x!=a and x!=b]
    clusters.append(clusters[a]+clusters[b])
    del clusters[min(a,b)]
    del clusters[max(a-1,b-1)]

    d = np.delete(d,[a,b],0)
    d = np.delete(d,[a,b],1)
    d = np.c_[d,np.zeros(nodesLeft-2)]
    d = np.r_[d,[clustDist]]

    nodesLeft = nodesLeft-1

print(clusters)

minDist = np.amin(d[np.nonzero(d)])
newClust = np.where(d == minDist)
a = newClust[0][0]
b = newClust[1][0]

pairs.append([a,b])

newInd = [0]*(n-1)
newInd[:2] = [x for x in ind[2] if x != a and x!=b]
ind.append(newInd)

heights[3] = d[a,b]/2
print(heights)
