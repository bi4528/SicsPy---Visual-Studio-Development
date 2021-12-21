
from sicspym import Graph, SicsAlgorithms, read_amalfi, read_ldgraphs_lab

if __name__ == "__main__":

    m = Graph("adjacency_listmat", "undirected_tag", 3)
    m.set_vertex_label(0, "red")
    m.set_vertex_label(1, "blue")
    m.set_vertex_label(2, "green")
    m.add_edge(0, 1)
    m.add_edge(1, 2)
    m.add_edge(2, 0)
    
    n = Graph("adjacency_listmat", "undirected_tag", 4)
    n.set_vertex_label(0, "red")
    n.set_vertex_label(1, "blue")
    n.set_vertex_label(2, "green")
    n.set_vertex_label(3, "green")
    n.add_edge(0, 1)
    n.add_edge(1, 2)
    n.add_edge(2, 0)
    n.add_edge(1, 3)
    n.add_edge(3, 0)

    md = Graph("adjacency_degreesortedlistmat", "undirected_tag", graph=m)
    nd = Graph("adjacency_degreesortedlistmat", "undirected_tag", graph=n)

    c = SicsAlgorithms("adjacency_listmat", "undirected_tag")
    
    c.backjumping_bitset_degreeprune_ind(m, n, "GCF", mapping=True)
    c.forwardchecking_bitset_mrv_degreeprune_ind(m, n, mapping=True)
    c.forwardchecking_bitset_mrv_degreeprune_ind(m, n)
    cd = SicsAlgorithms("adjacency_degreesortedlistmat", "undirected_tag")
    cd.backjumping_bitset_degreeprune_ind(md, nd, "GCF", mapping=True)
