import numpy as np
from Bio import SeqIO
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt

class DNASequenceClusterer:
    def __init__(self, file_name):
        self.file_name = file_name
        self.records = list(SeqIO.parse(file_name, "fasta"))
        self.n = len(self.records)
        self.m = len(self.records[0])
        self.dist_matrix = self.compute_distance_matrix()

    def compute_distance_matrix(self):
        dist_matrix = np.zeros((self.n, self.n))
        for j in range(self.n):
            for k in range(j):
                w1 = self.records[j].seq
                w2 = self.records[k].seq
                dist_matrix[j][k] = sum(c1 != c2 for c1, c2 in zip(w1, w2)) / self.m
        return dist_matrix

    def perform_upgma_clustering(self):
        return linkage(self.dist_matrix, method='average')

    def plot_dendrogram(self, Z):
        plt.figure(figsize=(10, 6))
        dendrogram(Z, labels=[str(i) for i in range(self.n)], leaf_rotation=90)
        plt.title('UPGMA Clustering Dendrogram')
        plt.xlabel('Sample Index')
        plt.ylabel('Distance')
        plt.show()

if __name__ == "__main__":
    file_name = "king_penguins_truncated.fasta"
    clusterer = DNASequenceClusterer(file_name)
    Z = clusterer.perform_upgma_clustering()
    clusterer.plot_dendrogram(Z)