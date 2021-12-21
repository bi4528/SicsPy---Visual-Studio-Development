import os
from sicspym import Graph, SicsAlgorithms, read_amalfi, read_ldgraphs_lab

if __name__ == "__main__":

    head, tail = os.path.split(__file__)
    file1 = os.path.join(head, "si2_b03_m200.A00")
    file2 = os.path.join(head, "si2_b03_m200.B00")
    g = read_amalfi("adjacency_listmat", "undirected_tag", repr(file1).strip("'"))
    h = read_amalfi("adjacency_listmat", "undirected_tag", repr(file2).strip("'"))
    s = SicsAlgorithms("adjacency_listmat", "undirected_tag")
    s.backjumping_bitset_degreeprune_ind(g, h, "GCF", mapping=False)
