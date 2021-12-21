
from sicspyem import adjacency_list_bidirectional_tag_alg as alg, adjacency_list_bidirectional_tag as graph, ostream_redirect

if __name__ == "__main__":

    m = graph(3)
    m.set_vertex_label(0, "red")
    m.set_vertex_label(1, "blue")
    m.set_vertex_label(2, "green")
    m.add_edge(0, 1)
    m.add_edge(1, 2)
    m.add_edge(2, 0)
    
    n = m = graph(4)
    n.set_vertex_label(0, "red")
    n.set_vertex_label(1, "blue")
    n.set_vertex_label(2, "green")
    n.set_vertex_label(3, "green")
    n.add_edge(0, 1)
    n.add_edge(1, 2)
    n.add_edge(2, 0)
    n.add_edge(1, 3)
    n.add_edge(3, 0)

    a = alg()
    with ostream_redirect(stdout=True, stderr=True):
        a.backjumping_bitset_degreeprune_ind(m, n, "DEG", False)

    