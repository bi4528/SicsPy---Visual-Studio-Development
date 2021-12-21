#include <fstream>
#include <iostream>
#include <string>

#include <adjacency_listmat.h>
#include <adjacency_list.h>
#include <adjacency_degreesortedlistmat.h>
#include <vertex_order.h>
#include <lazyforwardcheckingbackjumping_low_bitset_degreeprune_ind.h>
#include <backjumping_bitset_degreeprune_ind.h>
#include <backjumping_bitset_degreesequenceprune_ind.h>
#include <backjumping_degreeprune_ind.h>
#include <backjumping_degreesequenceprune_ind.h>
#include <backjumping_ind.h>
#include <backtracking_adjacentconsistency_forwardcount_ind.h>
#include <backtracking_adjacentconsistency_ind.h>
#include <backtracking_adjacentconsistency_precount_ind.h>
#include <backtracking_bitset_degreeprune_ind.h>
#include <backtracking_bitset_degreesequenceprune_ind.h>
#include <backtracking_degreeprune_adjacentconsistency_forwardcount_ind.h>
#include <backtracking_degreeprune_adjacentconsistency_ind.h>
#include <backtracking_degreeprune_adjacentconsistency_precount_ind.h>
#include <backtracking_degreeprune_ind.h>
#include <backtracking_degreesequenceprune_ind.h>
#include <backtracking_forwardcount_ind.h>
#include <backtracking_ind.h>
#include <backtracking_parent_adjacentconsistency_forwardcount_ind.h>
#include <backtracking_parent_adjacentconsistency_ind.h>
#include <backtracking_parent_adjacentconsistency_precount_ind.h>
#include <backtracking_parent_degreeprune_adjacentconsistency_forwardcount_ind.h>
#include <backtracking_parent_degreeprune_adjacentconsistency_ind.h>
#include <backtracking_parent_degreeprune_adjacentconsistency_precount_ind.h>
#include <backtracking_parent_degreeprune_ind.h>
#include <backtracking_parent_forwardcount_ind.h>
#include <backtracking_parent_ind.h>
#include <conflictbackjumping_degreeprune_ind.h>
#include <conflictbackjumping_degreesequenceprune_ind.h>
#include <conflictbackjumping_ind.h>
#include <forwardchecking_bitset_degreeprune_ac1_ind.h>
#include <forwardchecking_bitset_degreeprune_countingalldifferent_ind.h>
#include <forwardchecking_bitset_degreeprune_ind.h>
#include <forwardchecking_bitset_degreesequenceprune_ac1_ind.h>
#include <forwardchecking_bitset_degreesequenceprune_countingalldifferent_ind.h>
#include <forwardchecking_bitset_degreesequenceprune_ind.h>
#include <forwardchecking_bitset_mrv_degreeprune_ac1_ind.h>
#include <forwardchecking_bitset_mrv_degreeprune_countingalldifferent_ind.h>
#include <forwardchecking_bitset_mrv_degreeprune_ind.h>
#include <forwardchecking_bitset_mrv_degreesequenceprune_ac1_ind.h>
#include <forwardchecking_bitset_mrv_degreesequenceprune_countingalldifferent_ind.h>
#include <forwardchecking_bitset_mrv_degreesequenceprune_ind.h>
#include <forwardchecking_degreeprune_ind.h>
#include <forwardchecking_degreesequenceprune_ind.h>
#include <forwardchecking_ind.h>
#include <forwardchecking_mrv_degreeprune_ind.h>
#include <lazyforwardchecking_degreeprune_ind.h>
#include <lazyforwardchecking_degreesequenceprune_ind.h>
#include <lazyforwardchecking_ind.h>
#include <lazyforwardchecking_low_bitset_degreeprune_ind.h>
#include <lazyforwardchecking_low_bitset_degreesequenceprune_ind.h>
#include <lazyforwardchecking_low_degreeprune_ind.h>
#include <lazyforwardchecking_low_degreesequenceprune_ind.h>
#include <lazyforwardchecking_low_ind.h>
#include <lazyforwardchecking_low_parent_degreeprune_ind.h>
#include <lazyforwardchecking_low_parent_ind.h>
#include <lazyforwardchecking_parent_degreeprune_ind.h>
#include <lazyforwardchecking_parent_degreesequenceprune_ind.h>
#include <lazyforwardchecking_parent_ind.h>
#include <lazyforwardcheckingbackjumping_low_bitset_degreeprune_ind.h>
#include <lazyforwardcheckingbackjumping_low_bitset_degreesequenceprune_ind.h>
#include <read_amalfi.h>
#include <read_gal.h>
#include <read_galv.h>
#include <read_galve.h>
#include <read_gf.h>
#include <read_ldgraphs.h>

#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>

using adjacency_list_undirected_tag = sics::adjacency_list<uint16_t, sics::undirected_tag, std::string>;
using adjacency_list_bidirectional_tag = sics::adjacency_list<uint16_t, sics::bidirectional_tag, std::string>;

using adjacency_listmat_undirected_tag = sics::adjacency_listmat<uint16_t, sics::undirected_tag, std::string>;
using adjacency_listmat_bidirectional_tag = sics::adjacency_listmat<uint16_t, sics::bidirectional_tag, std::string>;

using adjacency_degreesortedlistmat_undirected_tag = sics::adjacency_degreesortedlistmat<uint16_t, sics::undirected_tag, std::string>;
using adjacency_degreesortedlistmat_bidirectional_tag = sics::adjacency_degreesortedlistmat<uint16_t, sics::bidirectional_tag, std::string>;


template <typename graph_type>
class sics_algorithms
{
public:

    void backjumping_bitset_degreeprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backjumping_bitset_degreeprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backjumping_bitset_degreesequenceprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backjumping_bitset_degreesequenceprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backjumping_degreeprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backjumping_degreeprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backjumping_degreesequenceprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backjumping_degreesequenceprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backjumping_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backjumping_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_adjacentconsistency_forwardcount_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_adjacentconsistency_forwardcount_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_adjacentconsistency_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_adjacentconsistency_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_adjacentconsistency_precount_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_adjacentconsistency_precount_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_bitset_degreeprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_bitset_degreeprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_bitset_degreesequenceprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_bitset_degreesequenceprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_degreeprune_adjacentconsistency_forwardcount_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_degreeprune_adjacentconsistency_forwardcount_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_degreeprune_adjacentconsistency_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_degreeprune_adjacentconsistency_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_degreeprune_adjacentconsistency_precount_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_degreeprune_adjacentconsistency_precount_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_degreeprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_degreeprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_degreesequenceprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_degreesequenceprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_forwardcount_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_forwardcount_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_parent_adjacentconsistency_forwardcount_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_parent_adjacentconsistency_forwardcount_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_parent_adjacentconsistency_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_parent_adjacentconsistency_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_parent_adjacentconsistency_precount_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_parent_adjacentconsistency_precount_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_parent_degreeprune_adjacentconsistency_forwardcount_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_parent_degreeprune_adjacentconsistency_forwardcount_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_parent_degreeprune_adjacentconsistency_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_parent_degreeprune_adjacentconsistency_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_parent_degreeprune_adjacentconsistency_precount_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_parent_degreeprune_adjacentconsistency_precount_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_parent_degreeprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_parent_degreeprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_parent_forwardcount_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_parent_forwardcount_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void backtracking_parent_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::backtracking_parent_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void conflictbackjumping_degreeprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::conflictbackjumping_degreeprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void conflictbackjumping_degreesequenceprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::conflictbackjumping_degreesequenceprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void conflictbackjumping_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::conflictbackjumping_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void forwardchecking_bitset_degreeprune_ac1_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::forwardchecking_bitset_degreeprune_ac1_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void forwardchecking_bitset_degreeprune_countingalldifferent_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::forwardchecking_bitset_degreeprune_countingalldifferent_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void forwardchecking_bitset_degreeprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::forwardchecking_bitset_degreeprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void forwardchecking_bitset_degreesequenceprune_ac1_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::forwardchecking_bitset_degreesequenceprune_ac1_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void forwardchecking_bitset_degreesequenceprune_countingalldifferent_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::forwardchecking_bitset_degreesequenceprune_countingalldifferent_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void forwardchecking_bitset_degreesequenceprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::forwardchecking_bitset_degreesequenceprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void forwardchecking_bitset_mrv_degreeprune_ac1_ind(graph_type g, graph_type h, bool mapping) {

        int count = 0;

        sics::forwardchecking_bitset_mrv_degreeprune_ac1_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            });

        std::cout << count << std::endl;
    };

    void forwardchecking_bitset_mrv_degreeprune_countingalldifferent_ind(graph_type g, graph_type h, bool mapping) {

        int count = 0;

        sics::forwardchecking_bitset_mrv_degreeprune_countingalldifferent_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            });

        std::cout << count << std::endl;
    };

    void forwardchecking_bitset_mrv_degreeprune_ind(graph_type g, graph_type h, bool mapping) {

        int count = 0;

        sics::forwardchecking_bitset_mrv_degreeprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            });

        std::cout << count << std::endl;
    };

    void forwardchecking_bitset_mrv_degreesequenceprune_ac1_ind(graph_type g, graph_type h, bool mapping) {

        int count = 0;

        sics::forwardchecking_bitset_mrv_degreesequenceprune_ac1_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            });

        std::cout << count << std::endl;
    };

    void forwardchecking_bitset_mrv_degreesequenceprune_countingalldifferent_ind(graph_type g, graph_type h, bool mapping) {

        int count = 0;

        sics::forwardchecking_bitset_mrv_degreesequenceprune_countingalldifferent_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            });

        std::cout << count << std::endl;
    };

    void forwardchecking_bitset_mrv_degreesequenceprune_ind(graph_type g, graph_type h, bool mapping) {

        int count = 0;

        sics::forwardchecking_bitset_mrv_degreesequenceprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            });

        std::cout << count << std::endl;
    };

    void forwardchecking_degreeprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::forwardchecking_degreeprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void forwardchecking_degreesequenceprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::forwardchecking_degreesequenceprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void forwardchecking_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::forwardchecking_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void forwardchecking_mrv_degreeprune_ind(graph_type g, graph_type h, bool mapping) {

        int count = 0;

        sics::forwardchecking_mrv_degreeprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            });

        std::cout << count << std::endl;
    };

    void lazyforwardchecking_degreeprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::lazyforwardchecking_degreeprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void lazyforwardchecking_degreesequenceprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::lazyforwardchecking_degreesequenceprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void lazyforwardchecking_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::lazyforwardchecking_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void lazyforwardchecking_low_bitset_degreeprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::lazyforwardchecking_low_bitset_degreeprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void lazyforwardchecking_low_bitset_degreesequenceprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::lazyforwardchecking_low_bitset_degreesequenceprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void lazyforwardchecking_low_degreeprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::lazyforwardchecking_low_degreeprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void lazyforwardchecking_low_degreesequenceprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::lazyforwardchecking_low_degreesequenceprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void lazyforwardchecking_low_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::lazyforwardchecking_low_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void lazyforwardchecking_low_parent_degreeprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::lazyforwardchecking_low_parent_degreeprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void lazyforwardchecking_low_parent_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::lazyforwardchecking_low_parent_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void lazyforwardchecking_parent_degreeprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::lazyforwardchecking_parent_degreeprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void lazyforwardchecking_parent_degreesequenceprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::lazyforwardchecking_parent_degreesequenceprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void lazyforwardchecking_parent_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::lazyforwardchecking_parent_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void lazyforwardcheckingbackjumping_low_bitset_degreeprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::lazyforwardcheckingbackjumping_low_bitset_degreeprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

    void lazyforwardcheckingbackjumping_low_bitset_degreesequenceprune_ind(graph_type g, graph_type h, std::string v_order, bool mapping) {

        int count = 0;
        auto index_order_g = sics::vertex_order_GreatestConstraintFirst(g);

        if (v_order == "DEG") {
            index_order_g = sics::vertex_order_DEG(g);
        }
        else if (v_order == "RDEG") {
            index_order_g = sics::vertex_order_RDEG(g);
        }

        sics::lazyforwardcheckingbackjumping_low_bitset_degreesequenceprune_ind(
            g,
            h,
            [&count, &mapping](auto const& S) {
                ++count;
                if (mapping == true) {
                    std::cout << "{";
                    for (int i = 0; i < S.map.size(); ++i) {
                        std::cout << i << ": " << S.map[i];
                        if (i < S.map.size() - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << "}" << std::endl;
                }
                return true;
            },
            index_order_g);

        std::cout << count << std::endl;
    };

};

template <typename graph_type>
graph_type read_amalfi(std::string myfile_name) {

    std::ifstream myfile(myfile_name, std::ios::in | std::ios::binary);

    std::istream& myfile_add = myfile;
    graph_type g = sics::read_amalfi<graph_type>(myfile_add);

    return g;
};

template <typename graph_type>
graph_type read_gal(std::string myfile_name) {

    std::ifstream myfile(myfile_name, std::ios::in | std::ios::binary);

    std::istream& myfile_add = myfile;
    graph_type g = sics::read_gal<graph_type>(myfile_add);

    return g;
};

template <typename graph_type>
graph_type read_galv(std::string myfile_name) {

    std::ifstream myfile(myfile_name, std::ios::in | std::ios::binary);

    std::istream& myfile_add = myfile;
    graph_type g = sics::read_galv<graph_type>(myfile_add);

    return g;
};

template <typename graph_type>
graph_type read_gf(std::string myfile_name) {

    std::ifstream myfile(myfile_name, std::ios::in | std::ios::binary);

    std::istream& myfile_add = myfile;
    graph_type g = sics::read_gf<graph_type>(myfile_add);

    return g;
};

template <typename graph_type>
graph_type read_ldgraphs_unl(std::string myfile_name) {

    std::ifstream myfile(myfile_name, std::ios::in | std::ios::binary);

    std::istream& myfile_add = myfile;
    graph_type g = sics::read_ldgraphs_unl<graph_type>(myfile_add);

    return g;
};

template <typename graph_type>
graph_type read_ldgraphs_lab(std::string myfile_name) {

    std::ifstream myfile(myfile_name, std::ios::in | std::ios::binary);

    std::istream& myfile_add = myfile;
    graph_type g = sics::read_ldgraphs_lab<graph_type>(myfile_add);

    return g;
};


namespace py = pybind11;

PYBIND11_MODULE(sicspyem, m) {

    py::add_ostream_redirect(m, "ostream_redirect");

    py::class_<adjacency_list_bidirectional_tag>(m, "adjacency_list_bidirectional_tag")
        .def(py::init<std::uint16_t>())
        .def("set_vertex_label", &adjacency_list_bidirectional_tag::set_vertex_label)
        .def("add_edge", &adjacency_list_bidirectional_tag::add_edge);

    py::class_<adjacency_list_undirected_tag>(m, "adjacency_list_undirected_tag")
        .def(py::init<std::uint16_t>())
        .def("set_vertex_label", &adjacency_list_undirected_tag::set_vertex_label)
        .def("add_edge", &adjacency_list_undirected_tag::add_edge);

    py::class_<adjacency_listmat_bidirectional_tag>(m, "adjacency_listmat_bidirectional_tag")
        .def(py::init<std::uint16_t>())
        .def("set_vertex_label", &adjacency_listmat_bidirectional_tag::set_vertex_label)
        .def("add_edge", &adjacency_listmat_bidirectional_tag::add_edge);

    py::class_<adjacency_listmat_undirected_tag>(m, "adjacency_listmat_undirected_tag")
        .def(py::init<std::uint16_t>())
        .def("set_vertex_label", &adjacency_listmat_undirected_tag::set_vertex_label)
        .def("add_edge", &adjacency_listmat_undirected_tag::add_edge);

    py::class_<adjacency_degreesortedlistmat_bidirectional_tag>(m, "adjacency_degreesortedlistmat_bidirectional_tag")
        .def(py::init<adjacency_list_bidirectional_tag>())
        .def(py::init<adjacency_listmat_bidirectional_tag>())
        .def(py::init<adjacency_degreesortedlistmat_bidirectional_tag>());

    py::class_<adjacency_degreesortedlistmat_undirected_tag>(m, "adjacency_degreesortedlistmat_undirected_tag")
        .def(py::init<adjacency_list_undirected_tag>())
        .def(py::init<adjacency_listmat_undirected_tag>())
        .def(py::init<adjacency_degreesortedlistmat_undirected_tag>());

    m.def("read_amalfi_adjacency_list_bidirectional_tag", &read_amalfi<adjacency_list_bidirectional_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_amalfi_adjacency_list_undirected_tag", &read_amalfi<adjacency_list_undirected_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_amalfi_adjacency_listmat_bidirectional_tag", &read_amalfi<adjacency_listmat_bidirectional_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_amalfi_adjacency_listmat_undirected_tag", &read_amalfi<adjacency_listmat_undirected_tag>, py::arg("myfile_name"), py::return_value_policy::move);

    m.def("read_gal_adjacency_list_bidirectional_tag", &read_gal<adjacency_list_bidirectional_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_gal_adjacency_list_undirected_tag", &read_gal<adjacency_list_undirected_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_gal_adjacency_listmat_bidirectional_tag", &read_gal<adjacency_listmat_bidirectional_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_gal_adjacency_listmat_undirected_tag", &read_gal<adjacency_listmat_undirected_tag>, py::arg("myfile_name"), py::return_value_policy::move);

    m.def("read_galv_adjacency_list_bidirectional_tag", &read_galv<adjacency_list_bidirectional_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_galv_adjacency_list_undirected_tag", &read_galv<adjacency_list_undirected_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_galv_adjacency_listmat_bidirectional_tag", &read_galv<adjacency_listmat_bidirectional_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_galv_adjacency_listmat_undirected_tag", &read_galv<adjacency_listmat_undirected_tag>, py::arg("myfile_name"), py::return_value_policy::move);

    /*m.def("read_galve_adjacency_list_bidirectional_tag", &read_galve<adjacency_list_bidirectional_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_galve_adjacency_list_undirected_tag", &read_galve<adjacency_list_undirected_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_galve_adjacency_listmat_bidirectional_tag", &read_galve<adjacency_listmat_bidirectional_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_galve_adjacency_listmat_undirected_tag", &read_galve<adjacency_listmat_undirected_tag>, py::arg("myfile_name"), py::return_value_policy::move);*/

    m.def("read_gf_adjacency_list_bidirectional_tag", &read_gf<adjacency_list_bidirectional_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_gf_adjacency_list_undirected_tag", &read_gf<adjacency_list_undirected_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_gf_adjacency_listmat_bidirectional_tag", &read_gf<adjacency_listmat_bidirectional_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_gf_adjacency_listmat_undirected_tag", &read_gf<adjacency_listmat_undirected_tag>, py::arg("myfile_name"), py::return_value_policy::move);

    m.def("read_ldgraphs_unl_adjacency_list_bidirectional_tag", &read_ldgraphs_unl<adjacency_list_bidirectional_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_ldgraphs_unl_adjacency_list_undirected_tag", &read_ldgraphs_unl<adjacency_list_undirected_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_ldgraphs_unl_adjacency_listmat_bidirectional_tag", &read_ldgraphs_unl<adjacency_listmat_bidirectional_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_ldgraphs_unl_adjacency_listmat_undirected_tag", &read_ldgraphs_unl<adjacency_listmat_undirected_tag>, py::arg("myfile_name"), py::return_value_policy::move);

    m.def("read_ldgraphs_lab_adjacency_list_bidirectional_tag", &read_ldgraphs_lab<adjacency_list_bidirectional_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_ldgraphs_lab_adjacency_list_undirected_tag", &read_ldgraphs_lab<adjacency_list_undirected_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_ldgraphs_lab_adjacency_listmat_bidirectional_tag", &read_ldgraphs_lab<adjacency_listmat_bidirectional_tag>, py::arg("myfile_name"), py::return_value_policy::move);
    m.def("read_ldgraphs_lab_adjacency_listmat_undirected_tag", &read_ldgraphs_lab<adjacency_listmat_undirected_tag>, py::arg("myfile_name"), py::return_value_policy::move);

    py::class_<sics_algorithms<adjacency_list_bidirectional_tag>>(m, "adjacency_list_bidirectional_tag_alg")
        .def(py::init<>())
        .def("backjumping_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backjumping_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backjumping_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_degreeprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backjumping_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backjumping_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backjumping_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_adjacentconsistency_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_degreeprune_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_adjacentconsistency_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_degreeprune_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_degreeprune_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_forwardcount_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_parent_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_adjacentconsistency_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_parent_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_parent_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_parent_degreeprune_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_adjacentconsistency_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_parent_degreeprune_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_parent_degreeprune_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_parent_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_forwardcount_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_parent_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::backtracking_parent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("conflictbackjumping_degreeprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::conflictbackjumping_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("conflictbackjumping_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::conflictbackjumping_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("conflictbackjumping_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::conflictbackjumping_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreeprune_ac1_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::forwardchecking_bitset_degreeprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreeprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::forwardchecking_bitset_degreeprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::forwardchecking_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreesequenceprune_ac1_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::forwardchecking_bitset_degreesequenceprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreesequenceprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::forwardchecking_bitset_degreesequenceprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::forwardchecking_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreeprune_ac1_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::forwardchecking_bitset_mrv_degreeprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreeprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::forwardchecking_bitset_mrv_degreeprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreeprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::forwardchecking_bitset_mrv_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreesequenceprune_ac1_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::forwardchecking_bitset_mrv_degreesequenceprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreesequenceprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::forwardchecking_bitset_mrv_degreesequenceprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::forwardchecking_bitset_mrv_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_degreeprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::forwardchecking_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::forwardchecking_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::forwardchecking_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_mrv_degreeprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::forwardchecking_mrv_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("lazyforwardchecking_degreeprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::lazyforwardchecking_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::lazyforwardchecking_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::lazyforwardchecking_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::lazyforwardchecking_low_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::lazyforwardchecking_low_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_degreeprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::lazyforwardchecking_low_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::lazyforwardchecking_low_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::lazyforwardchecking_low_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_parent_degreeprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::lazyforwardchecking_low_parent_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_parent_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::lazyforwardchecking_low_parent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_parent_degreeprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::lazyforwardchecking_parent_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_parent_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::lazyforwardchecking_parent_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_parent_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::lazyforwardchecking_parent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardcheckingbackjumping_low_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::lazyforwardcheckingbackjumping_low_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardcheckingbackjumping_low_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_bidirectional_tag>::lazyforwardcheckingbackjumping_low_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"));

    py::class_<sics_algorithms<adjacency_list_undirected_tag>>(m, "adjacency_list_undirected_tag_alg")
        .def(py::init<>())
        .def("backjumping_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backjumping_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backjumping_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_degreeprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backjumping_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backjumping_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backjumping_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_adjacentconsistency_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_degreeprune_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_adjacentconsistency_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_degreeprune_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_degreeprune_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_forwardcount_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_parent_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_adjacentconsistency_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_parent_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_parent_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_parent_degreeprune_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_adjacentconsistency_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_parent_degreeprune_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_parent_degreeprune_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_parent_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_forwardcount_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_parent_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::backtracking_parent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("conflictbackjumping_degreeprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::conflictbackjumping_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("conflictbackjumping_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::conflictbackjumping_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("conflictbackjumping_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::conflictbackjumping_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreeprune_ac1_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::forwardchecking_bitset_degreeprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreeprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::forwardchecking_bitset_degreeprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::forwardchecking_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreesequenceprune_ac1_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::forwardchecking_bitset_degreesequenceprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreesequenceprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::forwardchecking_bitset_degreesequenceprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::forwardchecking_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreeprune_ac1_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::forwardchecking_bitset_mrv_degreeprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreeprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::forwardchecking_bitset_mrv_degreeprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreeprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::forwardchecking_bitset_mrv_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreesequenceprune_ac1_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::forwardchecking_bitset_mrv_degreesequenceprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreesequenceprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::forwardchecking_bitset_mrv_degreesequenceprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::forwardchecking_bitset_mrv_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_degreeprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::forwardchecking_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::forwardchecking_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::forwardchecking_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_mrv_degreeprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::forwardchecking_mrv_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("lazyforwardchecking_degreeprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::lazyforwardchecking_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::lazyforwardchecking_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::lazyforwardchecking_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::lazyforwardchecking_low_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::lazyforwardchecking_low_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_degreeprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::lazyforwardchecking_low_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::lazyforwardchecking_low_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::lazyforwardchecking_low_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_parent_degreeprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::lazyforwardchecking_low_parent_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_parent_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::lazyforwardchecking_low_parent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_parent_degreeprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::lazyforwardchecking_parent_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_parent_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::lazyforwardchecking_parent_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_parent_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::lazyforwardchecking_parent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardcheckingbackjumping_low_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::lazyforwardcheckingbackjumping_low_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardcheckingbackjumping_low_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_list_undirected_tag>::lazyforwardcheckingbackjumping_low_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"));

    py::class_<sics_algorithms<adjacency_listmat_bidirectional_tag>>(m, "adjacency_listmat_bidirectional_tag_alg")
        .def(py::init<>())
        .def("backjumping_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backjumping_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backjumping_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backjumping_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backjumping_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backjumping_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_adjacentconsistency_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_degreeprune_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_adjacentconsistency_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_degreeprune_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_degreeprune_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_forwardcount_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_parent_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_adjacentconsistency_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_parent_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_parent_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_parent_degreeprune_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_adjacentconsistency_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_parent_degreeprune_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_parent_degreeprune_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_parent_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_forwardcount_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_parent_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::backtracking_parent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("conflictbackjumping_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::conflictbackjumping_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("conflictbackjumping_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::conflictbackjumping_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("conflictbackjumping_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::conflictbackjumping_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreeprune_ac1_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::forwardchecking_bitset_degreeprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreeprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::forwardchecking_bitset_degreeprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::forwardchecking_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreesequenceprune_ac1_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::forwardchecking_bitset_degreesequenceprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreesequenceprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::forwardchecking_bitset_degreesequenceprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::forwardchecking_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreeprune_ac1_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::forwardchecking_bitset_mrv_degreeprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreeprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::forwardchecking_bitset_mrv_degreeprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::forwardchecking_bitset_mrv_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreesequenceprune_ac1_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::forwardchecking_bitset_mrv_degreesequenceprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreesequenceprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::forwardchecking_bitset_mrv_degreesequenceprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::forwardchecking_bitset_mrv_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::forwardchecking_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::forwardchecking_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::forwardchecking_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_mrv_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::forwardchecking_mrv_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("lazyforwardchecking_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::lazyforwardchecking_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::lazyforwardchecking_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::lazyforwardchecking_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::lazyforwardchecking_low_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::lazyforwardchecking_low_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::lazyforwardchecking_low_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::lazyforwardchecking_low_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::lazyforwardchecking_low_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_parent_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::lazyforwardchecking_low_parent_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_parent_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::lazyforwardchecking_low_parent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_parent_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::lazyforwardchecking_parent_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_parent_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::lazyforwardchecking_parent_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_parent_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::lazyforwardchecking_parent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardcheckingbackjumping_low_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::lazyforwardcheckingbackjumping_low_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardcheckingbackjumping_low_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_bidirectional_tag>::lazyforwardcheckingbackjumping_low_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"));

    py::class_<sics_algorithms<adjacency_listmat_undirected_tag>>(m, "adjacency_listmat_undirected_tag_alg")
        .def(py::init<>())
        .def("backjumping_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backjumping_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backjumping_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backjumping_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backjumping_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backjumping_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_adjacentconsistency_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_degreeprune_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_adjacentconsistency_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_degreeprune_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_degreeprune_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_forwardcount_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_parent_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_adjacentconsistency_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_parent_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_parent_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_parent_degreeprune_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_adjacentconsistency_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_parent_degreeprune_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_parent_degreeprune_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_parent_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_forwardcount_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_parent_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::backtracking_parent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("conflictbackjumping_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::conflictbackjumping_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("conflictbackjumping_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::conflictbackjumping_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("conflictbackjumping_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::conflictbackjumping_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreeprune_ac1_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::forwardchecking_bitset_degreeprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreeprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::forwardchecking_bitset_degreeprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::forwardchecking_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreesequenceprune_ac1_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::forwardchecking_bitset_degreesequenceprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreesequenceprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::forwardchecking_bitset_degreesequenceprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::forwardchecking_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreeprune_ac1_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::forwardchecking_bitset_mrv_degreeprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreeprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::forwardchecking_bitset_mrv_degreeprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::forwardchecking_bitset_mrv_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreesequenceprune_ac1_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::forwardchecking_bitset_mrv_degreesequenceprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreesequenceprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::forwardchecking_bitset_mrv_degreesequenceprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::forwardchecking_bitset_mrv_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::forwardchecking_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::forwardchecking_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::forwardchecking_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_mrv_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::forwardchecking_mrv_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("lazyforwardchecking_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::lazyforwardchecking_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::lazyforwardchecking_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::lazyforwardchecking_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::lazyforwardchecking_low_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::lazyforwardchecking_low_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::lazyforwardchecking_low_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::lazyforwardchecking_low_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::lazyforwardchecking_low_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_parent_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::lazyforwardchecking_low_parent_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_parent_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::lazyforwardchecking_low_parent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_parent_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::lazyforwardchecking_parent_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_parent_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::lazyforwardchecking_parent_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_parent_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::lazyforwardchecking_parent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardcheckingbackjumping_low_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::lazyforwardcheckingbackjumping_low_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardcheckingbackjumping_low_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_listmat_undirected_tag>::lazyforwardcheckingbackjumping_low_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"));

    py::class_<sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>>(m, "adjacency_degreesortedlistmat_bidirectional_tag_alg")
        .def(py::init<>())
        .def("backjumping_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backjumping_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backjumping_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backjumping_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backjumping_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backjumping_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_adjacentconsistency_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_degreeprune_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_adjacentconsistency_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_degreeprune_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_degreeprune_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_forwardcount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_parent_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_adjacentconsistency_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_parent_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_parent_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_parent_degreeprune_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_adjacentconsistency_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_parent_degreeprune_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_parent_degreeprune_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_parent_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_forwardcount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_parent_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::backtracking_parent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("conflictbackjumping_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::conflictbackjumping_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("conflictbackjumping_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::conflictbackjumping_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("conflictbackjumping_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::conflictbackjumping_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreeprune_ac1_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::forwardchecking_bitset_degreeprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreeprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::forwardchecking_bitset_degreeprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::forwardchecking_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreesequenceprune_ac1_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::forwardchecking_bitset_degreesequenceprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreesequenceprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::forwardchecking_bitset_degreesequenceprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::forwardchecking_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreeprune_ac1_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::forwardchecking_bitset_mrv_degreeprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreeprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::forwardchecking_bitset_mrv_degreeprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::forwardchecking_bitset_mrv_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreesequenceprune_ac1_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::forwardchecking_bitset_mrv_degreesequenceprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreesequenceprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::forwardchecking_bitset_mrv_degreesequenceprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::forwardchecking_bitset_mrv_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::forwardchecking_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::forwardchecking_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::forwardchecking_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_mrv_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::forwardchecking_mrv_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("lazyforwardchecking_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::lazyforwardchecking_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::lazyforwardchecking_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::lazyforwardchecking_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::lazyforwardchecking_low_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::lazyforwardchecking_low_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::lazyforwardchecking_low_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::lazyforwardchecking_low_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::lazyforwardchecking_low_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_parent_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::lazyforwardchecking_low_parent_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_parent_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::lazyforwardchecking_low_parent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_parent_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::lazyforwardchecking_parent_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_parent_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::lazyforwardchecking_parent_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_parent_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::lazyforwardchecking_parent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardcheckingbackjumping_low_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::lazyforwardcheckingbackjumping_low_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardcheckingbackjumping_low_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_bidirectional_tag>::lazyforwardcheckingbackjumping_low_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"));

    py::class_<sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>>(m, "adjacency_degreesortedlistmat_undirected_tag_alg")
        .def(py::init<>())
        .def("backjumping_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backjumping_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backjumping_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backjumping_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backjumping_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backjumping_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backjumping_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_adjacentconsistency_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_degreeprune_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_adjacentconsistency_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_degreeprune_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_degreeprune_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_forwardcount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_parent_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_adjacentconsistency_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_parent_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_parent_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_adjacentconsistency_forwardcount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_parent_degreeprune_adjacentconsistency_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_adjacentconsistency_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_parent_degreeprune_adjacentconsistency_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_adjacentconsistency_precount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_parent_degreeprune_adjacentconsistency_precount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_parent_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_forwardcount_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_parent_forwardcount_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("backtracking_parent_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::backtracking_parent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("conflictbackjumping_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::conflictbackjumping_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("conflictbackjumping_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::conflictbackjumping_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("conflictbackjumping_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::conflictbackjumping_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreeprune_ac1_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::forwardchecking_bitset_degreeprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreeprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::forwardchecking_bitset_degreeprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::forwardchecking_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreesequenceprune_ac1_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::forwardchecking_bitset_degreesequenceprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreesequenceprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::forwardchecking_bitset_degreesequenceprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::forwardchecking_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreeprune_ac1_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::forwardchecking_bitset_mrv_degreeprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreeprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::forwardchecking_bitset_mrv_degreeprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::forwardchecking_bitset_mrv_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreesequenceprune_ac1_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::forwardchecking_bitset_mrv_degreesequenceprune_ac1_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreesequenceprune_countingalldifferent_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::forwardchecking_bitset_mrv_degreesequenceprune_countingalldifferent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_bitset_mrv_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::forwardchecking_bitset_mrv_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("forwardchecking_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::forwardchecking_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::forwardchecking_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::forwardchecking_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("forwardchecking_mrv_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::forwardchecking_mrv_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("mapping"))
        .def("lazyforwardchecking_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::lazyforwardchecking_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::lazyforwardchecking_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::lazyforwardchecking_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::lazyforwardchecking_low_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::lazyforwardchecking_low_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::lazyforwardchecking_low_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::lazyforwardchecking_low_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::lazyforwardchecking_low_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_parent_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::lazyforwardchecking_low_parent_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_low_parent_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::lazyforwardchecking_low_parent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_parent_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::lazyforwardchecking_parent_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_parent_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::lazyforwardchecking_parent_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardchecking_parent_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::lazyforwardchecking_parent_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardcheckingbackjumping_low_bitset_degreeprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::lazyforwardcheckingbackjumping_low_bitset_degreeprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"))
        .def("lazyforwardcheckingbackjumping_low_bitset_degreesequenceprune_ind",
            &sics_algorithms<adjacency_degreesortedlistmat_undirected_tag>::lazyforwardcheckingbackjumping_low_bitset_degreesequenceprune_ind,
            py::arg("graph1"),
            py::arg("graph2"),
            py::arg("vertex_order"),
            py::arg("mapping"));

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
