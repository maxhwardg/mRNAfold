#ifndef MWMRNA_LIB_CODON_GRAPH_FOLD_HPP
#define MWMRNA_LIB_CODON_GRAPH_FOLD_HPP

#include <stddef.h>

#include <algorithm>
#include <execution>
#include <stdexcept>
#include <utility>
#include <vector>

#include "codon_graph.hpp"
#include "md_array.hpp"
#include "rna.hpp"

using std::pair;
using std::vector;

namespace mwmrna::fold_codons {

template <typename EnergyModelT>
class GraphTraceIterator;

template <typename EnergyModelT>
struct GraphFoldContext {
    using ViennaParamsT = typename EnergyModelT::params_type;
    using SemiringT = typename ViennaParamsT::semiring_type;
    using EnergyT = typename SemiringT::energy_type;
    EnergyModelT em;
    int max_two_loop_size = 30;
    int min_hairpin_loop_size = 3;
    CodonGraph<SemiringT> graph;
    bool parallel = false;
    bool lyngso = false;
    vector<vector<EnergyT>> pair_bonus{};
};

template <typename EnergyModelT>
class GraphFold {
   private:
    using EnergyT = typename GraphFoldContext<EnergyModelT>::EnergyT;
    using SemiringT = typename GraphFoldContext<EnergyModelT>::SemiringT;

   public:
    GraphFold(GraphFoldContext<EnergyModelT> _ctx) : ctx(_ctx), g(ctx.graph) { FillTables(); }
    GraphFold(GraphFoldContext<EnergyModelT>&& _ctx) : ctx(_ctx), g(ctx.graph) { FillTables(); }

    friend GraphTraceIterator<EnergyModelT>;

    /// @brief Maximum number of nodes any nucleotide index has
    int max_nodes;
    /// @brief Represents an undecided codon
    int tbd_cdn;
    /// @brief Number of bases in the sequence
    int n;

    /// @brief sum_paths[ni,nj] is the sum of all graph paths from node ni to node nj
    /// Assumes ni comes before nj
    MdArray<EnergyT, 2> sum_paths;

    /// @brief DP table. (i,j) -> sum (under semiring) of structues where node i and node j are
    /// paired graph nodes
    MdArray<EnergyT, 2> p_table;
    /// @brief DP table for multiloops. (branch, j, i) -> sum (under semiring) of all
    /// multiloops fragments containing at least branch branches between node i and node j
    /// @note We index by j,i for cache-friendliness when filling the ML table
    MdArray<EnergyT, 3> ml_table;
    /// @brief DP table for external loops (i) -> sum (under semiring) of all external loops
    /// for the suffix starting at node i
    MdArray<EnergyT, 1> el_table;

    // Represents valid internal loops unpaired configurations for Lyngso's optimization
    enum class L_IlState { IL22 = 0, IL23, IL32, ILGood, NUM_STATES };

    /// @brief DP table for a generalised lyngso's internal loop optimization.
    /// (up, ilstate, i,j) -> sum (under semiring) of internal loop parts from node i to node j
    /// assuming i and j are unpaired and with up unpaired nucleotides (not including i and j).
    /// Also, the internal loop is currently in ilstate. We need ilstate to avoid special cases for
    /// 2x2, 2x3, 3x2 loops, which must be deal with separately as they have different energy
    /// functions to general internal loops. Note that this is also true for 1xn (including 1x1)
    /// loops, but we will never use the l_table for these cases, since we always start with a 2x2,
    /// since all valid general internal loops have one side >= 2.
    MdArray<EnergyT, 4> l_table;
    SemiringT sr;

    /// @brief The context under which folding was computed
    GraphFoldContext<EnergyModelT> ctx;

    /// @brief Alias for ctx.graph
    const CodonGraph<SemiringT>& g;

   protected:
    void FillTables() {
        if (ctx.min_hairpin_loop_size < 3)
            throw std::runtime_error("Hairpin loop size < 3 not supported");

        n = ctx.graph.Nts();

        if (ctx.pair_bonus.size()) {
            bool good = static_cast<int>(ctx.pair_bonus.size()) == n;
            for (auto& row : ctx.pair_bonus) {
                good &= static_cast<int>(row.size()) == n;
            }
            if (!good) {
                throw std::runtime_error(
                    "Pair bonus table must either be empty or be size N x N (where is the number "
                    "of nucleotides)");
            }
        }


        int nodes = g.NumNodes();

        ml_table = MdArray<EnergyT, 3>(sr.Zero(), {3, nodes, nodes});
        p_table = MdArray<EnergyT, 2>(sr.Zero(), {nodes, nodes});
        el_table = MdArray<EnergyT, 1>(sr.Zero(), {nodes});
        if (ctx.lyngso)
            l_table = MdArray<EnergyT, 4>(
                sr.Zero(),
                {ctx.max_two_loop_size + 1, static_cast<int>(L_IlState::NUM_STATES), nodes, nodes});

        sum_paths = MdArray<EnergyT, 2>(sr.Zero(), {nodes, nodes});
        for (int i = 0; i < n; ++i) {
            for (int nodei : g.NtNodes(i)) {
                sum_paths[{nodei, nodei}] = sr.One();
            }
        }
        for (int i = n - 1; i >= 0; --i) {
            for (int nodei : g.NtNodes(i)) {
                for (int j = i + 1; j < n; ++j) {
                    for (int nodej : g.NtNodes(j)) {
                        for (auto& e : g.EdgesL2R(nodei)) {
                            sum_paths[{nodei, nodej}] =
                                sr.Add(sum_paths[{nodei, nodej}],
                                       sr.MultMany(sum_paths[{e.to, nodej}], e.weight));
                        }
                    }
                }
            }
        }

        std::vector<int> indexes;
        if (ctx.parallel) {
            indexes = std::vector<int>(n);
            std::iota(indexes.begin(), indexes.end(), 0);
        }

        // Dp loop
        for (int i = n - 1; i >= 0; --i) {
            if (ctx.parallel) {
                std::for_each(std::execution::par_unseq,
                              indexes.begin() + i + ctx.min_hairpin_loop_size + 1, indexes.end(),
                              [&](int j) { FillP(i, j); });
                if (ctx.lyngso) {
                    for (int sz = 0; sz <= ctx.max_two_loop_size; ++sz) {
                        std::for_each(std::execution::par_unseq,
                                      indexes.begin() + i + ctx.min_hairpin_loop_size + 1,
                                      indexes.end(), [&](int j) { FillL(sz, i, j); });
                    }
                }
                std::for_each(std::execution::par_unseq, indexes.begin() + i, indexes.end(),
                              [&](int j) { FillMl(i, j); });
            } else {
                for (int j = i + ctx.min_hairpin_loop_size + 1; j < n; ++j) {
                    FillP(i, j);
                }
                if (ctx.lyngso) {
                    for (int sz = 0; sz <= ctx.max_two_loop_size; ++sz) {
                        for (int j = i + ctx.min_hairpin_loop_size + 1; j < n; ++j) {
                            FillL(sz, i, j);
                        }
                    }
                }
                for (int j = i; j < n; ++j) {
                    FillMl(i, j);
                }
            }
            FillEl(i);
        }
    }

    void FillP(int i, int j) {
        EnergyT e;

        for (int inode : g.NtNodes(i)) {
            for (int jnode : g.NtNodes(j)) {
                auto bi = g.NodeBase(inode);
                auto bj = g.NodeBase(jnode);
                BasePair bp = BasesToPair(bi, bj);
                // Ignore invalid base pairs
                if (bp == BasePair::NN) continue;

                e = sr.Zero();

                // Stacks
                for (auto& [ip1node, bip1, ip1w] : g.EdgesL2R(inode)) {
                    for (auto& [jm1node, bjm1, jm1w] : g.EdgesR2L(jnode)) {
                        BasePair bpinner = BasesToPair(bip1, bjm1);
                        if (bpinner == BasePair::NN) continue;
                        e = sr.Add(e, sr.MultMany(p_table[{ip1node, jm1node}],
                                                  ctx.em.Stack(bp, bpinner), ip1w, jm1w));
                    }
                }

                // Bulges (5')
                for (int k = i + 2; (j - 1) - k - 1 >= ctx.min_hairpin_loop_size &&
                                    k - i - 1 <= ctx.max_two_loop_size;
                     ++k) {
                    for (int knode : g.NtNodes(k)) {
                        for (int lnode : g.NtNodes(j - 1)) {
                            Base bk = g.NodeBase(knode);
                            Base bl = g.NodeBase(lnode);
                            BasePair bpkl = BasesToPair(bk, bl);
                            if (bpkl == BasePair::NN) continue;
                            e = sr.Add(e, sr.MultMany(p_table[{knode, lnode}],
                                                      ctx.em.Bulge(bp, bpkl, k - i - 1),
                                                      sum_paths[{inode, knode}],
                                                      sum_paths[{lnode, jnode}]));
                        }
                    }
                }

                // Bulges (3')
                for (int l = j - 2; j - l - 1 <= ctx.max_two_loop_size &&
                                    l - (i + 1) - 1 >= ctx.min_hairpin_loop_size;
                     --l) {
                    for (int knode : g.NtNodes(i + 1)) {
                        for (int lnode : g.NtNodes(l)) {
                            Base bk = g.NodeBase(knode);
                            Base bl = g.NodeBase(lnode);
                            BasePair bpkl = BasesToPair(bk, bl);
                            if (bpkl == BasePair::NN) continue;
                            e = sr.Add(e, sr.MultMany(p_table[{knode, lnode}],
                                                      ctx.em.Bulge(bp, bpkl, j - l - 1),
                                                      sum_paths[{inode, knode}],
                                                      sum_paths[{lnode, jnode}]));
                        }
                    }
                }

                // Internal loops (1x1)
                if (i + 2 + ctx.min_hairpin_loop_size < j - 2) {
                    for (auto& [ip1node, bip1, ip1w] : g.EdgesL2R(inode)) {
                        for (auto& [ip2node, bip2, ip2w] : g.EdgesL2R(ip1node)) {
                            for (auto& [jm1node, bjm1, jm1w] : g.EdgesR2L(jnode)) {
                                for (auto& [jm2node, bjm2, jm2w] : g.EdgesR2L(jm1node)) {
                                    BasePair bpkl = BasesToPair(bip2, bjm2);
                                    if (bpkl == BasePair::NN) continue;
                                    e = sr.Add(e, sr.MultMany(p_table[{ip2node, jm2node}],
                                                              ctx.em.IL11(bp, bpkl, bip1, bjm1),
                                                              ip1w, ip2w, jm1w, jm2w));
                                }
                            }
                        }
                    }
                }

                // Internal loops (1xn)
                // Ignores 1x1 and handles 1x2 as a special case
                for (int k = j - 3;
                     j - k <= ctx.max_two_loop_size && k - (i + 2) - 1 >= ctx.min_hairpin_loop_size;
                     --k) {
                    // Note j-3 is because 1x1 is handled separately
                    for (auto& [ip1node, bip1, ip1w] : g.EdgesL2R(inode)) {
                        for (auto& [ip2node, bip2, ip2w] : g.EdgesL2R(ip1node)) {
                            for (int knode : g.NtNodes(k)) {
                                Base bk = g.NodeBase(knode);
                                BasePair bpkl = BasesToPair(bip2, bk);
                                if (bpkl == BasePair::NN) continue;
                                for (auto& [jm1node, bjm1, jm1w] : g.EdgesR2L(jnode)) {
                                    for (auto& [kp1node, bkp1, kp1w] : g.EdgesL2R(knode)) {
                                        if (j - k - 1 == 2) {
                                            e = sr.Add(
                                                e,
                                                sr.MultMany(p_table[{ip2node, knode}],
                                                            ctx.em.IL12(bp, bpkl, bip1, bjm1, bkp1),
                                                            ip1w, ip2w, jm1w, kp1w,
                                                            sum_paths[{kp1node, jm1node}]));
                                        } else {
                                            e = sr.Add(e,
                                                       sr.MultMany(p_table[{ip2node, knode}],
                                                                   ctx.em.IL1N(bp, bpkl, bip1, bjm1,
                                                                               bkp1, j - k - 1),
                                                                   ip1w, ip2w, jm1w, kp1w,
                                                                   sum_paths[{kp1node, jm1node}]));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // Internal loops (nx1)
                // Ignores 1x1 and handles 2x1 as a special case
                for (int k = i + 3;
                     (j - 2) - k - 1 >= ctx.min_hairpin_loop_size && k - i <= ctx.max_two_loop_size;
                     ++k) {
                    // Note i+3 is because 1x1 is handled separately
                    for (auto& [jm1node, bjm1, jm1w] : g.EdgesR2L(jnode)) {
                        for (auto& [jm2node, bjm2, jm2w] : g.EdgesR2L(jm1node)) {
                            for (int knode : g.NtNodes(k)) {
                                Base bk = g.NodeBase(knode);
                                BasePair bpkl = BasesToPair(bk, bjm2);
                                if (bpkl == BasePair::NN) continue;
                                for (auto& [ip1node, bip1, ip1w] : g.EdgesL2R(inode)) {
                                    for (auto& [km1node, bkm1, km1w] : g.EdgesR2L(knode)) {
                                        if (k - i - 1 == 2) {
                                            e = sr.Add(
                                                e,
                                                sr.MultMany(p_table[{knode, jm2node}],
                                                            ctx.em.IL21(bp, bpkl, bip1, bjm1, bkm1),
                                                            ip1w, jm1w, jm2w, km1w,
                                                            sum_paths[{ip1node, km1node}]));
                                        } else {
                                            e = sr.Add(e,
                                                       sr.MultMany(p_table[{knode, jm2node}],
                                                                   ctx.em.ILN1(bp, bpkl, bip1, bjm1,
                                                                               bkm1, k - i - 1),
                                                                   ip1w, jm1w, jm2w, km1w,
                                                                   sum_paths[{ip1node, km1node}]));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // Internal loops. Handles 2x2 and 2x3 as special cases. Also handles general
                // internal loops unless lyngso is enabled
                int max_unpaired_on_one_side = ctx.max_two_loop_size;
                if (ctx.lyngso) {
                    // Set this to 3 so that these loops handle the special cases (2x2, 2x3)
                    // But the lyngso DP handles the rest
                    max_unpaired_on_one_side = 3;
                }
                for (int k = i + 3; k - i - 1 <= max_unpaired_on_one_side && k < j - 3; ++k) {
                    for (int l = j - 3; j - l - 1 <= max_unpaired_on_one_side &&
                                        j - l - 1 + k - i - 1 <= ctx.max_two_loop_size &&
                                        l - k - 1 >= ctx.min_hairpin_loop_size;
                         --l) {
                        int lup = k - i - 1, rup = j - l - 1;
                        for (int knode : g.NtNodes(k)) {
                            for (int lnode : g.NtNodes(l)) {
                                Base bk = g.NodeBase(knode);
                                Base bl = g.NodeBase(lnode);
                                BasePair bpkl = BasesToPair(bk, bl);
                                if (bpkl == BasePair::NN) continue;
                                for (auto& [km1node, bkm1, km1w] : g.EdgesR2L(knode)) {
                                    for (auto& [lp1node, blp1, lp1w] : g.EdgesL2R(lnode)) {
                                        for (auto& [ip1node, bip1, ip1w] : g.EdgesL2R(inode)) {
                                            for (auto& [jm1node, bjm1, jm1w] : g.EdgesR2L(jnode)) {
                                                if (lup == 2 && rup == 2) {
                                                    e = sr.Add(
                                                        e,
                                                        sr.MultMany(p_table[{knode, lnode}],
                                                                    ctx.em.IL22(bp, bpkl, bip1,
                                                                                bjm1, bkm1, blp1),
                                                                    ip1w, jm1w, km1w, lp1w,
                                                                    sum_paths[{ip1node, km1node}],
                                                                    sum_paths[{lp1node, jm1node}]));
                                                } else if (lup == 2 && rup == 3) {
                                                    e = sr.Add(
                                                        e,
                                                        sr.MultMany(p_table[{knode, lnode}],
                                                                    ctx.em.IL23(bp, bpkl, bip1,
                                                                                bjm1, bkm1, blp1),
                                                                    ip1w, jm1w, km1w, lp1w,
                                                                    sum_paths[{ip1node, km1node}],
                                                                    sum_paths[{lp1node, jm1node}]));
                                                } else if (lup == 3 && rup == 2) {
                                                    e = sr.Add(
                                                        e,
                                                        sr.MultMany(p_table[{knode, lnode}],
                                                                    ctx.em.IL32(bp, bpkl, bip1,
                                                                                bjm1, bkm1, blp1),
                                                                    ip1w, jm1w, km1w, lp1w,
                                                                    sum_paths[{ip1node, km1node}],
                                                                    sum_paths[{lp1node, jm1node}]));
                                                } else if (!ctx.lyngso) {
                                                    // 3x3 are handled separately by lyngso DP
                                                    e = sr.Add(
                                                        e,
                                                        sr.MultMany(p_table[{knode, lnode}],
                                                                    ctx.em.NormalInternalLoop(
                                                                        bp, bpkl, bip1, bjm1, bkm1,
                                                                        blp1, lup, rup),
                                                                    ip1w, jm1w, km1w, lp1w,
                                                                    sum_paths[{ip1node, km1node}],
                                                                    sum_paths[{lp1node, jm1node}]));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // Lyngso general internal loops. Handles 3x3 and larger internal loops if lyngso is
                // enabled
                // Minimum size is 6 since 1x1, 1x2, 2x2, 2x3 are handled separately
                if (ctx.lyngso) {
                    for (int sz = 6; sz <= ctx.max_two_loop_size; ++sz) {
                        for (auto& [ip1node, bip1, ip1w] : g.EdgesL2R(inode)) {
                            for (auto& [jm1node, bjm1, jm1w] : g.EdgesR2L(jnode)) {
                                for (auto& [ip2node, bip2, ip2w] : g.EdgesL2R(ip1node)) {
                                    for (auto& [jm2node, bjm2, jm2w] : g.EdgesR2L(jm1node)) {
                                        e = sr.Add(
                                            e,
                                            sr.MultMany(
                                                l_table[{sz - 4, static_cast<int>(L_IlState::IL22),
                                                         ip2node, jm2node}],
                                                ctx.em.ILInnerMismatch(bp, bip1, bjm1),
                                                ctx.em.ILInit(sz), ip1w, jm1w, ip2w, jm2w));
                                    }
                                }
                            }
                        }
                    }
                }

                // Hairpin loops
                // ASSUMPTION: Special hairpins are small (e.g., < 6nts)
                auto special_hairpins = ctx.em.SpecialHairpins(j - i - 1);

                if (special_hairpins.empty()) {
                    // No special hairpins, just enumerate dangles
                    for (auto& [ip1node, bip1, ip1w] : g.EdgesL2R(inode)) {
                        for (auto& [jm1node, bjm1, jm1w] : g.EdgesR2L(jnode)) {
                            e = sr.Add(e,
                                       sr.MultMany(ctx.em.HairpinNormal(bp, bip1, bjm1, j - i - 1),
                                                   ip1w, jm1w, sum_paths[{ip1node, jm1node}]));
                        }
                    }
                } else {
                    // Special hairpins
                    // Deal with special hairpins by enumerating all paths
                    // This may be faster as a mutable DFS, but this is unlikely to be a bottleneck
                    vector<std::tuple<Primary, EnergyT, int>> stk;
                    stk.push_back({{bi}, sr.One(), inode});
                    while (!stk.empty()) {
                        Primary pri;
                        EnergyT ecurr;
                        int currnode;
                        tie(pri, ecurr, currnode) = stk.back();
                        stk.pop_back();
                        if (static_cast<int>(pri.size()) == j - i) {
                            ecurr = sr.MultMany(ecurr, sum_paths[{currnode, jnode}]);
                            pri.push_back(bj);
                            auto it =
                                std::find_if(special_hairpins.begin(), special_hairpins.end(),
                                             [&pri](const auto& p) { return p.first == pri; });
                            if (it != special_hairpins.end()) {
                                e = sr.Add(e, sr.MultMany(it->second, ecurr));
                            } else {
                                e = sr.Add(
                                    e, sr.MultMany(ctx.em.HairpinNormal(
                                                       bp, pri[1], pri[pri.size() - 2], j - i - 1),
                                                   ecurr));
                            }
                        } else {
                            for (auto& edge : g.EdgesL2R(currnode)) {
                                auto pri_ = pri;
                                pri_.push_back(edge.base);
                                stk.push_back({pri_, sr.MultMany(ecurr, edge.weight), edge.to});
                            }
                        }
                    }
                }

                // Multi loops
                for (auto& [ip1node, bip1, ip1w] : g.EdgesL2R(inode)) {
                    for (auto& [jm1node, bjm1, jm1w] : g.EdgesR2L(jnode)) {
                        e = sr.Add(e, sr.MultMany(ml_table[{2, jm1node, ip1node}], ip1w, jm1w,
                                                  ctx.em.MultiClosing(bp, Base::N, Base::N)));
                    }
                }

                if (ctx.pair_bonus.size()) {
                    e = sr.MultMany(e, ctx.pair_bonus[i][j]);
                }

                p_table[{inode, jnode}] = e;
            }
        }
    }

    void FillL(int sz, int i, int j) {
        for (int inode : g.NtNodes(i)) {
            for (int jnode : g.NtNodes(j)) {
                for (int ilstate = 0; ilstate < static_cast<int>(L_IlState::NUM_STATES);
                     ++ilstate) {
                    Base bi = g.NodeBase(inode);
                    Base bj = g.NodeBase(jnode);
                    EnergyT e = sr.Zero();

                    // Size 0 is a base case
                    if (sz == 0) {
                        // We may only end in a good state
                        if (static_cast<L_IlState>(ilstate) == L_IlState::ILGood) {
                            for (auto& [ip1node, bip1, ip1w] : g.EdgesL2R(inode)) {
                                for (auto& [jm1node, bjm1, jm1w] : g.EdgesR2L(jnode)) {
                                    BasePair bp = BasesToPair(bip1, bjm1);
                                    if (bp == BasePair::NN) continue;
                                    e = sr.Add(e, sr.MultMany(p_table[{ip1node, jm1node}],
                                                              ctx.em.ILOuterMismatch(bp, bi, bj),
                                                              ip1w, jm1w));
                                }
                            }
                        }
                        l_table[{0, ilstate, inode, jnode}] = e;
                    } else {
                        // Use all sizes > 0 for asymmetry
                        // This is the only place asymmetry is added

                        if (i + sz + 1 + ctx.min_hairpin_loop_size < j - 1) {
                            L_IlState next_state = static_cast<L_IlState>(ilstate);
                            if (sz == 1 && next_state == L_IlState::IL22)
                                next_state = L_IlState::IL32;
                            else
                                next_state = L_IlState::ILGood;
                            for (int knode : g.NtNodes(i + sz)) {
                                e = sr.Add(
                                    e, sr.MultMany(
                                           ctx.em.ILAsymmetry(sz),
                                           l_table[{0, static_cast<int>(next_state), knode, jnode}],
                                           sum_paths[{inode, knode}]));
                            }
                            next_state = static_cast<L_IlState>(ilstate);
                            if (sz == 1 && next_state == L_IlState::IL22)
                                next_state = L_IlState::IL23;
                            else
                                next_state = L_IlState::ILGood;
                            for (int lnode : g.NtNodes(j - sz)) {
                                e = sr.Add(
                                    e, sr.MultMany(
                                           ctx.em.ILAsymmetry(sz),
                                           l_table[{0, static_cast<int>(next_state), inode, lnode}],
                                           sum_paths[{lnode, jnode}]));
                            }
                        }
                        // If sz >= 2 we can also try making i and j unpaired with no asymmetry
                        if (sz >= 2) {
                            // Make i and j unpaired. No asymmetry added
                            for (auto& [ip1node, bip1, ip1w] : g.EdgesL2R(inode)) {
                                // This loop could be avoided by adding a bit to the DP
                                // It's unlikely to be worth it in practice
                                for (auto& [jm1node, bjm1, jm1w] : g.EdgesR2L(jnode)) {
                                    // Zero asymmetry can have a contribution depending on model
                                    // parameters, so we need to deal with it as a special case
                                    EnergyT asymw = sz == 2 ? ctx.em.ILAsymmetry(0) : sr.One();
                                    // We always move to a good state, since adding one unpaired to
                                    // each side at least makes a 3x3
                                    e = sr.Add(
                                        e, sr.MultMany(
                                               l_table[{sz - 2, static_cast<int>(L_IlState::ILGood),
                                                        ip1node, jm1node}],
                                               asymw, ip1w, jm1w));
                                }
                            }
                        }
                        l_table[{sz, ilstate, inode, jnode}] = e;
                    }
                }
            }
        }
    }

    void FillMl(int i, int j) {
        EnergyT e;
        for (int br = 0; br <= 2; ++br) {
            for (int inode : g.NtNodes(i)) {
                for (int jnode : g.NtNodes(j)) {
                    if (i == j && inode != jnode) continue;

                    e = sr.Zero();
                    Base bi = g.NodeBase(inode);
                    Base bj = g.NodeBase(jnode);

                    // Unpaired nucleotide
                    if (i + 1 < j) {
                        // Next nucleotide comes before j
                        for (auto& [ip1node, bip1, ip1w] : g.EdgesL2R(inode)) {
                            e = sr.Add(e, sr.MultMany(ctx.em.MultiUnpaired(), ip1w,
                                                      ml_table[{br, jnode, ip1node}]));
                        }
                    } else if (i + 1 == j) {
                        // Next nucleotide is j
                        e = sr.MultMany(ctx.em.MultiUnpaired(), ml_table[{br, jnode, jnode}],
                                        sum_paths[{inode, jnode}]);
                    } else if (i == j && br == 0) {
                        // i==j. Only allowed to be unpaired if br == 0
                        e = ctx.em.MultiUnpaired();
                    }

                    // Branches
                    // k<j because we deal with the (i,j) branch as a special case
                    for (int k = i + ctx.min_hairpin_loop_size + 1; k < j; ++k) {
                        // Could be optimized by only looping over pairs
                        for (int knode : g.NtNodes(k)) {
                            Base bk = g.NodeBase(knode);
                            BasePair bp = BasesToPair(bi, bk);
                            if (bp == BasePair::NN) continue;
                            // This could be optimized by using more dp
                            for (auto [kp1node, bkp1, kp1w] : g.EdgesL2R(knode)) {
                                e = sr.Add(
                                    e, sr.MultMany(p_table[{inode, knode}],
                                                   ml_table[{std::max(0, br - 1), jnode, kp1node}],
                                                   ctx.em.MultiBranch(Base::N, bp, Base::N), kp1w));
                            }
                        }
                    }

                    BasePair bp = BasesToPair(bi, bj);
                    if (br <= 1 && bp != BasePair::NN) {
                        // Special case for (i,j) branch. I.e., ending on a branch
                        e = sr.Add(e, sr.MultMany(p_table[{inode, jnode}],
                                                  ctx.em.MultiBranch(Base::N, bp, Base::N)));
                    }

                    ml_table[{br, jnode, inode}] = e;
                }
            }
        }
    }

    void FillEl(int i) {
        EnergyT e;
        for (int inode : g.NtNodes(i)) {
            Base bi = g.NodeBase(inode);

            if (i + 1 < n) {
                // Unpaired nucleotide
                e = sr.Zero();
                for (auto& [ip1node, bip1, ip1w] : g.EdgesL2R(inode)) {
                    e = sr.Add(e, sr.MultMany(el_table[{ip1node}], ip1w));
                }
            } else {
                // Base case
                e = sr.One();
            }

            // Branches
            for (int j = i + ctx.min_hairpin_loop_size + 1; j < n; ++j) {
                for (int jnode : g.NtNodes(j)) {
                    Base bj = g.NodeBase(jnode);
                    BasePair bp = BasesToPair(bi, bj);
                    if (bp == BasePair::NN) continue;
                    if (j + 1 < n) {
                        for (auto& [jp1node, bjp1, jp1w] : g.EdgesL2R(jnode)) {
                            e = sr.Add(e, sr.MultMany(el_table[{jp1node}],
                                                      ctx.em.ExtBranch(Base::N, bp, Base::N),
                                                      p_table[{inode, jnode}], jp1w));
                        }
                    } else {
                        e = sr.Add(e, sr.MultMany(ctx.em.ExtBranch(Base::N, bp, Base::N),
                                                  p_table[{inode, jnode}]));
                    }
                }
            }

            el_table[{inode}] = e;
        }
    }
};

}  // namespace mwmrna::fold_codons
#endif  // MWMRNA_LIB_CODON_GRAPH_FOLD_HPP