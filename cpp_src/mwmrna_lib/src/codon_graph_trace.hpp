#ifndef MWMRNA_LIB_CODON_GRAPH_TRACE_HPP
#define MWMRNA_LIB_CODON_GRAPH_TRACE_HPP

#include <random>
#include <set>
#include <cmath>

#include "codon_graph.hpp"
#include "codon_graph_fold.hpp"
#include "energy_semiring.hpp"
#include "rna.hpp"
#include "structure.hpp"

namespace mwmrna::fold_codons {

template <typename EnergyModelT>
struct GraphTraceContext {
    using EnergyT = typename GraphFoldContext<EnergyModelT>::EnergyT;
    using SemiringT = typename GraphFoldContext<EnergyModelT>::SemiringT;
    int max_pq_size = -1;
    double subopt_randomness = 0.0;
    EnergyT threshold = SemiringT().Zero();
    int rng_seed = 0;
    /// @bried Chance to keep a trace
    double keep_chance = 1.0;
};

template <typename EnergyT>
struct CodonGraphTrace {
    Primary pri;
    Matching match;
    EnergyT energy;
};

template <typename EnergyModelT>
class GraphTraceIterator {
   private:
    using EnergyT = typename GraphFoldContext<EnergyModelT>::EnergyT;
    using SemiringT = typename GraphFoldContext<EnergyModelT>::SemiringT;
    static_assert(std::is_base_of<MFESemiring<EnergyT>, SemiringT>::value,
                  "Semiring must be MFE to do tracebacks");

   public:
    GraphTraceIterator(const GraphFold<EnergyModelT>& fld,
                       const GraphTraceContext<EnergyModelT>& ctx = {})
        : fld(fld) {
        Init(ctx);
    }

    GraphTraceIterator(const GraphFold<EnergyModelT>&& fld,
                       const GraphTraceContext<EnergyModelT>& ctx = {})
        : fld(std::move(fld)) {
        Init(ctx);
    }

    bool HasNext() { return !set_pq.empty(); }

    CodonGraphTrace<EnergyT> Next() {
        while (1) {
            PartialTrace pt = *set_pq.begin();
            set_pq.erase(set_pq.begin());
            if (pt.idxs.empty()) {
                CodonGraphTrace<EnergyT> res{pt.pri, pt.match, pt.min_energy};
                return res;
            }
            DpIdx idx = pt.idxs.back();
            pt.idxs.pop_back();
            switch (idx.table) {
                case DpTable::P:
                    TraceP(idx, pt);
                    break;
                case DpTable::L:
                    TraceL(idx, pt);
                    break;
                case DpTable::EL:
                    TraceEl(idx, pt);
                    break;
                case DpTable::ML:
                    TraceMl(idx, pt);
                    break;
                case DpTable::SP:
                    TraceSumPaths(idx, pt);
                    break;
                default:
                    throw new std::runtime_error("Invalid table");
                    break;
            }
        }
    }

   private:
    void Init(const GraphTraceContext<EnergyModelT>& ctx) {
        if (ctx.subopt_randomness < 0 || ctx.subopt_randomness > 1) {
            throw new std::runtime_error("subopt_randomness must be in [0, 1]");
        }
        this->rng.seed(ctx.rng_seed);
        this->subopt_randomness = ctx.subopt_randomness;
        nt_keep_chance = std::pow(ctx.keep_chance, 1.0/fld.g.Nts());

        nt_keep_chance_n.resize(fld.g.Nts() + 1, 1.0);
        for (int i = 1; i <= fld.g.Nts(); ++i) {
            nt_keep_chance_n[i] = nt_keep_chance_n[i - 1] * nt_keep_chance;
        }

        node_to_ntidx.resize(fld.g.NumNodes());
        for (int i = 0; i < fld.g.Nts(); ++i) {
            for (int j : fld.g.NtNodes(i)) {
                node_to_ntidx[j] = i;
            }
        }

        InitTrace(ctx.threshold, ctx.max_pq_size);
    }

    /// @brief Must be called before Next/HasNext
    /// @param threshold only trace sequence-structures with energy < threshold
    void InitTrace(EnergyT threshold, int _max_pq_size = -1) {
        trace_threshold = threshold;
        std::uniform_real_distribution<double> dist(1.0 - this->subopt_randomness, 1.0);
        for (int node : fld.g.NtNodes(0)) {
            if (fld.el_table[{node}] < threshold) {
                EnergyT energy = fld.el_table[{node}];
                double priority = energy * dist(rng);
                set_pq.insert({.pri = Primary(fld.n, Base::N),
                               .match = Matching(fld.n, Matching::kUnmatched),
                               .idxs = {{.table = DpTable::EL, .i = node}},
                               .min_energy = energy,
                               .priority = priority});
            }
        }
        max_pq_size = _max_pq_size;
    }

    // These map onto the fld.ml_table, fld.el_table, fld.l_table, fld.sum_paths
    enum class DpTable { ML, EL, P, L, SP };

    struct DpIdx {
        DpTable table;
        int i{-1}, j{-1}, k{-1}, l{-1};
    };

    struct PartialTrace {
        Primary pri;
        Matching match;
        std::vector<DpIdx> idxs;
        EnergyT min_energy;
        double priority;
        auto operator<=>(const PartialTrace& other) const { return priority <=> other.priority; }
    };

    void RelaxTrace(EnergyT e, std::initializer_list<DpIdx> idxs, EnergyT base_e,
                    const PartialTrace& pt, int nts) {
        if (e >= trace_threshold) return;
        EnergyT min_energy = pt.min_energy + e - base_e;
        std::uniform_real_distribution<double> dist_pri(1.0 - this->subopt_randomness, 1.0);
        double priority = min_energy * dist_pri(rng);
        // Avoid an expensive copy if possible
        if (max_pq_size != -1 && static_cast<int>(set_pq.size()) >= max_pq_size &&
            priority >= set_pq.rbegin()->priority)
            return;

        // Only do keep_chance if it's enabled
        // Don't make the PQ empty after dropping out a trace
        if (nt_keep_chance_n[1] < 1.0 && !set_pq.empty()) {
            double p = std::uniform_real_distribution<double>(0.0, 1.0)(rng);
            if (p > nt_keep_chance_n[nts]) return;
        }

        PartialTrace npt = pt;
        npt.min_energy = min_energy;
        npt.priority = priority;
        for (const auto& nidx : idxs) npt.idxs.push_back(nidx);
        set_pq.insert(std::move(npt));
        while (max_pq_size != -1 && static_cast<int>(set_pq.size()) > max_pq_size) {
            set_pq.erase(std::prev(set_pq.end()));
        }
    }

    void TraceP(const DpIdx& idx, PartialTrace& pt) {
        int inode = idx.i, jnode = idx.j;
        int i = node_to_ntidx[inode], j = node_to_ntidx[jnode];
        Base bi = fld.g.NodeBase(inode), bj = fld.g.NodeBase(jnode);
        BasePair bp = BasesToPair(bi, bj);
        EnergyT base_e = fld.p_table[{inode, jnode}];

        pt.match[i] = j;
        pt.match[j] = i;
        pt.pri[i] = bi;
        pt.pri[j] = bj;

        EnergyT pair_bonus_e = sr.One();
        if (!fld.ctx.pair_bonus.empty()) {
            pair_bonus_e = fld.ctx.pair_bonus[i][j];
        }

        auto relax = [&](EnergyT e, std::initializer_list<DpIdx> idxs, int nts,
                         std::initializer_list<int> unpaired_nodes = {}) {
            // Copy in unpaired nodes that have base changes to track
            for (int node : unpaired_nodes)
                pt.pri[this->node_to_ntidx[node]] = fld.g.NodeBase(node);
            RelaxTrace(sr.MultMany(e, pair_bonus_e), idxs, base_e, pt, nts);
            // Clear out the unpaired nodes
            for (int node : unpaired_nodes) pt.pri[this->node_to_ntidx[node]] = Base::N;
        };

        // Stacks
        for (auto& [ip1node, bip1, ip1w] : fld.g.EdgesL2R(inode)) {
            for (auto& [jm1node, bjm1, jm1w] : fld.g.EdgesR2L(jnode)) {
                BasePair bpinner = BasesToPair(bip1, bjm1);
                if (bpinner == BasePair::NN) continue;
                relax(sr.MultMany(fld.p_table[{ip1node, jm1node}], fld.ctx.em.Stack(bp, bpinner),
                                  ip1w, jm1w),
                      {{DpTable::P, ip1node, jm1node}}, 2);
            }
        }

        // Bulges (5')
        for (int k = i + 2; (j - 1) - k - 1 >= fld.ctx.min_hairpin_loop_size &&
                            k - i - 1 <= fld.ctx.max_two_loop_size;
             ++k) {
            for (int knode : fld.g.NtNodes(k)) {
                for (int lnode : fld.g.NtNodes(j - 1)) {
                    Base bk = fld.g.NodeBase(knode);
                    Base bl = fld.g.NodeBase(lnode);
                    BasePair bpkl = BasesToPair(bk, bl);
                    if (bpkl == BasePair::NN) continue;
                    relax(sr.MultMany(fld.p_table[{knode, lnode}],
                                      fld.ctx.em.Bulge(bp, bpkl, k - i - 1),
                                      fld.sum_paths[{inode, knode}], fld.sum_paths[{lnode, jnode}]),
                          {{DpTable::P, knode, lnode},
                           {DpTable::SP, inode, knode},
                           {DpTable::SP, lnode, jnode}},
                          2 + k - i - 1);
                }
            }
        }

        // Bulges (3')
        for (int l = j - 2; j - l - 1 <= fld.ctx.max_two_loop_size &&
                            l - (i + 1) - 1 >= fld.ctx.min_hairpin_loop_size;
             --l) {
            for (int knode : fld.g.NtNodes(i + 1)) {
                for (int lnode : fld.g.NtNodes(l)) {
                    Base bk = fld.g.NodeBase(knode);
                    Base bl = fld.g.NodeBase(lnode);
                    BasePair bpkl = BasesToPair(bk, bl);
                    if (bpkl == BasePair::NN) continue;
                    relax(sr.MultMany(fld.p_table[{knode, lnode}],
                                      fld.ctx.em.Bulge(bp, bpkl, j - l - 1),
                                      fld.sum_paths[{inode, knode}], fld.sum_paths[{lnode, jnode}]),
                          {{DpTable::P, knode, lnode},
                           {DpTable::SP, inode, knode},
                           {DpTable::SP, lnode, jnode}},
                          2 + j - l - 1);
                }
            }
        }

        // Internal loops (1x1)
        if (i + 2 + fld.ctx.min_hairpin_loop_size < j - 2) {
            for (auto& [ip1node, bip1, ip1w] : fld.g.EdgesL2R(inode)) {
                for (auto& [ip2node, bip2, ip2w] : fld.g.EdgesL2R(ip1node)) {
                    for (auto& [jm1node, bjm1, jm1w] : fld.g.EdgesR2L(jnode)) {
                        for (auto& [jm2node, bjm2, jm2w] : fld.g.EdgesR2L(jm1node)) {
                            BasePair bpkl = BasesToPair(bip2, bjm2);
                            if (bpkl == BasePair::NN) continue;
                            relax(sr.MultMany(fld.p_table[{ip2node, jm2node}],
                                              fld.ctx.em.IL11(bp, bpkl, bip1, bjm1), ip1w, ip2w,
                                              jm1w, jm2w),
                                  {{DpTable::P, ip2node, jm2node}}, 4, {ip1node, jm1node});
                        }
                    }
                }
            }
        }

        // Internal loops (1xn)
        // Ignores 1x1 and handles 1x2 as a special case
        for (int k = j - 3;
             j - k <= fld.ctx.max_two_loop_size && k - (i + 2) - 1 >= fld.ctx.min_hairpin_loop_size;
             --k) {
            // Note j-3 is because 1x1 is handled separately
            for (auto& [ip1node, bip1, ip1w] : fld.g.EdgesL2R(inode)) {
                for (auto& [ip2node, bip2, ip2w] : fld.g.EdgesL2R(ip1node)) {
                    for (int knode : fld.g.NtNodes(k)) {
                        Base bk = fld.g.NodeBase(knode);
                        BasePair bpkl = BasesToPair(bip2, bk);
                        if (bpkl == BasePair::NN) continue;
                        for (auto& [jm1node, bjm1, jm1w] : fld.g.EdgesR2L(jnode)) {
                            for (auto& [kp1node, bkp1, kp1w] : fld.g.EdgesL2R(knode)) {
                                if (j - k - 1 == 2) {
                                    relax(sr.MultMany(fld.p_table[{ip2node, knode}],
                                                      fld.ctx.em.IL12(bp, bpkl, bip1, bjm1, bkp1),
                                                      ip1w, ip2w, jm1w, kp1w,
                                                      fld.sum_paths[{kp1node, jm1node}]),
                                          {{DpTable::P, ip2node, knode},
                                           {DpTable::SP, kp1node, jm1node}}, 5,
                                          {ip1node});
                                } else {
                                    relax(sr.MultMany(fld.p_table[{ip2node, knode}],
                                                      fld.ctx.em.IL1N(bp, bpkl, bip1, bjm1, bkp1,
                                                                      j - k - 1),
                                                      ip1w, ip2w, jm1w, kp1w,
                                                      fld.sum_paths[{kp1node, jm1node}]),
                                          {{DpTable::P, ip2node, knode},
                                           {DpTable::SP, kp1node, jm1node}}, 3+j-k-1,
                                          {ip1node});
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
             (j - 2) - k - 1 >= fld.ctx.min_hairpin_loop_size && k - i <= fld.ctx.max_two_loop_size;
             ++k) {
            // Note i+3 is because 1x1 is handled separately
            for (auto& [jm1node, bjm1, jm1w] : fld.g.EdgesR2L(jnode)) {
                for (auto& [jm2node, bjm2, jm2w] : fld.g.EdgesR2L(jm1node)) {
                    for (int knode : fld.g.NtNodes(k)) {
                        Base bk = fld.g.NodeBase(knode);
                        BasePair bpkl = BasesToPair(bk, bjm2);
                        if (bpkl == BasePair::NN) continue;
                        for (auto& [ip1node, bip1, ip1w] : fld.g.EdgesL2R(inode)) {
                            for (auto& [km1node, bkm1, km1w] : fld.g.EdgesR2L(knode)) {
                                if (k - i - 1 == 2) {
                                    relax(sr.MultMany(fld.p_table[{knode, jm2node}],
                                                      fld.ctx.em.IL21(bp, bpkl, bip1, bjm1, bkm1),
                                                      ip1w, jm1w, jm2w, km1w,
                                                      fld.sum_paths[{ip1node, km1node}]),
                                          {{DpTable::P, knode, jm2node},
                                           {DpTable::SP, ip1node, km1node}}, 5,
                                          {jm1node});
                                } else {
                                    relax(sr.MultMany(fld.p_table[{knode, jm2node}],
                                                      fld.ctx.em.ILN1(bp, bpkl, bip1, bjm1, bkm1,
                                                                      k - i - 1),
                                                      ip1w, jm1w, jm2w, km1w,
                                                      fld.sum_paths[{ip1node, km1node}]),
                                          {{DpTable::P, knode, jm2node},
                                           {DpTable::SP, ip1node, km1node}}, 3+k-i-1,
                                          {jm1node});
                                }
                            }
                        }
                    }
                }
            }
        }

        // Internal loops. Handles 2x2 and 2x3 as special cases. Also handles general
        // internal loops unless lyngso is enabled
        int max_unpaired_on_one_side = fld.ctx.max_two_loop_size;
        if (fld.ctx.lyngso) {
            // Set this to 3 so that these loops handle the special cases (2x2, 2x3)
            // But the lyngso DP handles the rest
            max_unpaired_on_one_side = 3;
        }
        for (int k = i + 3; k - i - 1 <= max_unpaired_on_one_side && k < j - 3; ++k) {
            for (int l = j - 3; j - l - 1 <= max_unpaired_on_one_side &&
                                j - l - 1 + k - i - 1 <= fld.ctx.max_two_loop_size &&
                                l - k - 1 >= fld.ctx.min_hairpin_loop_size;
                 --l) {
                int lup = k - i - 1, rup = j - l - 1;
                for (int knode : fld.g.NtNodes(k)) {
                    for (int lnode : fld.g.NtNodes(l)) {
                        Base bk = fld.g.NodeBase(knode);
                        Base bl = fld.g.NodeBase(lnode);
                        BasePair bpkl = BasesToPair(bk, bl);
                        if (bpkl == BasePair::NN) continue;
                        for (auto& [km1node, bkm1, km1w] : fld.g.EdgesR2L(knode)) {
                            for (auto& [lp1node, blp1, lp1w] : fld.g.EdgesL2R(lnode)) {
                                for (auto& [ip1node, bip1, ip1w] : fld.g.EdgesL2R(inode)) {
                                    for (auto& [jm1node, bjm1, jm1w] : fld.g.EdgesR2L(jnode)) {
                                        if (lup == 2 && rup == 2) {
                                            relax(sr.MultMany(fld.p_table[{knode, lnode}],
                                                              fld.ctx.em.IL22(bp, bpkl, bip1, bjm1,
                                                                              bkm1, blp1),
                                                              ip1w, jm1w, km1w, lp1w,
                                                              fld.sum_paths[{ip1node, km1node}],
                                                              fld.sum_paths[{lp1node, jm1node}]),
                                                  {{DpTable::P, knode, lnode},
                                                   {DpTable::SP, ip1node, km1node},
                                                   {DpTable::SP, lp1node, jm1node}}, 2+lup+rup);
                                        } else if (lup == 2 && rup == 3) {
                                            relax(sr.MultMany(fld.p_table[{knode, lnode}],
                                                              fld.ctx.em.IL23(bp, bpkl, bip1, bjm1,
                                                                              bkm1, blp1),
                                                              ip1w, jm1w, km1w, lp1w,
                                                              fld.sum_paths[{ip1node, km1node}],
                                                              fld.sum_paths[{lp1node, jm1node}]),
                                                  {{DpTable::P, knode, lnode},
                                                   {DpTable::SP, ip1node, km1node},
                                                   {DpTable::SP, lp1node, jm1node}}, 2+lup+rup);
                                        } else if (lup == 3 && rup == 2) {
                                            relax(sr.MultMany(fld.p_table[{knode, lnode}],
                                                              fld.ctx.em.IL32(bp, bpkl, bip1, bjm1,
                                                                              bkm1, blp1),
                                                              ip1w, jm1w, km1w, lp1w,
                                                              fld.sum_paths[{ip1node, km1node}],
                                                              fld.sum_paths[{lp1node, jm1node}]),
                                                  {{DpTable::P, knode, lnode},
                                                   {DpTable::SP, ip1node, km1node},
                                                   {DpTable::SP, lp1node, jm1node}}, 2+lup+rup);
                                        } else if (!fld.ctx.lyngso) {
                                            // 3x3 are handled separately by lyngso DP
                                            relax(sr.MultMany(fld.p_table[{knode, lnode}],
                                                              fld.ctx.em.NormalInternalLoop(
                                                                  bp, bpkl, bip1, bjm1, bkm1, blp1,
                                                                  lup, rup),
                                                              ip1w, jm1w, km1w, lp1w,
                                                              fld.sum_paths[{ip1node, km1node}],
                                                              fld.sum_paths[{lp1node, jm1node}]),
                                                  {{DpTable::P, knode, lnode},
                                                   {DpTable::SP, ip1node, km1node},
                                                   {DpTable::SP, lp1node, jm1node}}, 2+lup+rup);
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
        if (fld.ctx.lyngso) {
            for (int sz = 6; sz <= fld.ctx.max_two_loop_size; ++sz) {
                for (auto& [ip1node, bip1, ip1w] : fld.g.EdgesL2R(inode)) {
                    for (auto& [jm1node, bjm1, jm1w] : fld.g.EdgesR2L(jnode)) {
                        for (auto& [ip2node, bip2, ip2w] : fld.g.EdgesL2R(ip1node)) {
                            for (auto& [jm2node, bjm2, jm2w] : fld.g.EdgesR2L(jm1node)) {
                                relax(
                                    sr.MultMany(
                                        fld.l_table[{sz - 4,
                                                     static_cast<int>(
                                                         GraphFold<EnergyModelT>::L_IlState::IL22),
                                                     ip2node, jm2node}],
                                        fld.ctx.em.ILInnerMismatch(bp, bip1, bjm1),
                                        fld.ctx.em.ILInit(sz), ip1w, jm1w, ip2w, jm2w),
                                    {{DpTable::L, ip2node, jm2node, sz - 4,
                                      static_cast<int>(GraphFold<EnergyModelT>::L_IlState::IL22)}}, 4,
                                    {ip1node, jm1node});
                            }
                        }
                    }
                }
            }
        }

        auto relax_hairpin = [&](EnergyT e, const Primary& pri) {
            for (int k = i + 1; k < j; ++k) pt.pri[k] = pri[k - i];
            RelaxTrace(sr.MultMany(e, pair_bonus_e), {}, base_e, pt, j-i-1);
            for (int k = i + 1; k < j; ++k) pt.pri[k] = Base::N;
        };

        // Hairpin loops
        // ASSUMPTION: Special hairpins are small (e.g., < 6nts)
        auto special_hairpins = fld.ctx.em.SpecialHairpins(j - i - 1);

        if (special_hairpins.empty()) {
            // No special hairpins, just enumerate dangles
            for (auto& [ip1node, bip1, ip1w] : fld.g.EdgesL2R(inode)) {
                for (auto& [jm1node, bjm1, jm1w] : fld.g.EdgesR2L(jnode)) {
                    relax(sr.MultMany(fld.ctx.em.HairpinNormal(bp, bip1, bjm1, j - i - 1), ip1w,
                                      jm1w, fld.sum_paths[{ip1node, jm1node}]),
                          {{DpTable::SP, ip1node, jm1node}}, j - i - 1);
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
                    ecurr = sr.MultMany(ecurr, fld.sum_paths[{currnode, jnode}]);
                    pri.push_back(bj);
                    auto it = std::find_if(special_hairpins.begin(), special_hairpins.end(),
                                           [&pri](const auto& p) { return p.first == pri; });
                    if (it != special_hairpins.end()) {
                        relax_hairpin(sr.MultMany(it->second, ecurr), pri);
                    } else {
                        relax_hairpin(sr.MultMany(fld.ctx.em.HairpinNormal(
                                                      bp, pri[1], pri[pri.size() - 2], j - i - 1),
                                                  ecurr),
                                      pri);
                    }
                } else {
                    for (auto& edge : fld.g.EdgesL2R(currnode)) {
                        auto pri_ = pri;
                        pri_.push_back(edge.base);
                        stk.push_back({pri_, sr.MultMany(ecurr, edge.weight), edge.to});
                    }
                }
            }
        }

        // Multi loops
        for (auto& [ip1node, bip1, ip1w] : fld.g.EdgesL2R(inode)) {
            for (auto& [jm1node, bjm1, jm1w] : fld.g.EdgesR2L(jnode)) {
                relax(sr.MultMany(fld.ml_table[{2, jm1node, ip1node}], ip1w, jm1w,
                                  fld.ctx.em.MultiClosing(bp, Base::N, Base::N)),
                      {{DpTable::ML, ip1node, jm1node, 2}}, 0);
            }
        }
    }

    void TraceL(const DpIdx& idx, PartialTrace& pt) {
        int inode = idx.i, jnode = idx.j, sz = idx.k, ilstate_idx = idx.l;
        int i = node_to_ntidx[inode], j = node_to_ntidx[jnode];
        using IlStateT = typename GraphFold<EnergyModelT>::L_IlState;
        auto ilstate = static_cast<IlStateT>(ilstate_idx);
        Base bi = fld.g.NodeBase(inode);
        Base bj = fld.g.NodeBase(jnode);
        EnergyT base_e = fld.l_table[{sz, ilstate_idx, inode, jnode}];

        pt.pri[i] = bi;
        pt.pri[j] = bj;

        auto relax = [&](EnergyT e, std::initializer_list<DpIdx> idxs, int nts) {
            RelaxTrace(e, idxs, base_e, pt, nts);
        };

        // Size 0 is a base case
        if (sz == 0) {
            // We may only end in a good state
            if (ilstate == IlStateT::ILGood) {
                for (auto& [ip1node, bip1, ip1w] : fld.g.EdgesL2R(inode)) {
                    for (auto& [jm1node, bjm1, jm1w] : fld.g.EdgesR2L(jnode)) {
                        BasePair bp = BasesToPair(bip1, bjm1);
                        if (bp == BasePair::NN) continue;
                        relax(sr.MultMany(fld.p_table[{ip1node, jm1node}],
                                          fld.ctx.em.ILOuterMismatch(bp, bi, bj), ip1w, jm1w),
                              {{DpTable::P, ip1node, jm1node}}, 2);
                    }
                }
            }
        } else {
            // Use all sizes > 0 for asymmetry
            // This is the only place asymmetry is added

            if (i + sz + 1 + fld.ctx.min_hairpin_loop_size < j - 1) {
                IlStateT next_state = static_cast<IlStateT>(ilstate);
                if (sz == 1 && next_state == IlStateT::IL22)
                    next_state = IlStateT::IL32;
                else
                    next_state = IlStateT::ILGood;
                for (int knode : fld.g.NtNodes(i + sz)) {
                    relax(sr.MultMany(fld.ctx.em.ILAsymmetry(sz),
                                      fld.l_table[{0, static_cast<int>(next_state), knode, jnode}],
                                      fld.sum_paths[{inode, knode}]),
                          {{DpTable::L, knode, jnode, 0, static_cast<int>(next_state)},
                           {DpTable::SP, inode, knode}}, sz);
                }
                next_state = static_cast<IlStateT>(ilstate);
                if (sz == 1 && next_state == IlStateT::IL22)
                    next_state = IlStateT::IL23;
                else
                    next_state = IlStateT::ILGood;
                for (int lnode : fld.g.NtNodes(j - sz)) {
                    relax(sr.MultMany(fld.ctx.em.ILAsymmetry(sz),
                                      fld.l_table[{0, static_cast<int>(next_state), inode, lnode}],
                                      fld.sum_paths[{lnode, jnode}]),
                          {{DpTable::L, inode, lnode, 0, static_cast<int>(next_state)},
                           {DpTable::SP, lnode, jnode}}, sz);
                }
            }
            // If sz >= 2 we can also try making i and j unpaired with no asymmetry
            if (sz >= 2) {
                // Make i and j unpaired. No asymmetry added
                for (auto& [ip1node, bip1, ip1w] : fld.g.EdgesL2R(inode)) {
                    // This loop could be avoided by adding a bit to the DP
                    // It's unlikely to be worth it in practice
                    for (auto& [jm1node, bjm1, jm1w] : fld.g.EdgesR2L(jnode)) {
                        // Zero asymmetry can have a contribution depending on model
                        // parameters, so we need to deal with it as a special case
                        EnergyT asymw = sz == 2 ? fld.ctx.em.ILAsymmetry(0) : sr.One();
                        // We always move to a good state, since adding one unpaired to
                        // each side at least makes a 3x3
                        relax(sr.MultMany(fld.l_table[{sz - 2, static_cast<int>(IlStateT::ILGood),
                                                       ip1node, jm1node}],
                                          asymw, ip1w, jm1w),
                              {{DpTable::L, ip1node, jm1node, sz - 2,
                                static_cast<int>(IlStateT::ILGood)}}, 2);
                    }
                }
            }
        }
    }

    void TraceMl(const DpIdx& idx, PartialTrace& pt) {
        int inode = idx.i, jnode = idx.j, br = idx.k;
        int i = node_to_ntidx[inode], j = node_to_ntidx[jnode];
        EnergyT base_e = fld.ml_table[{br, jnode, inode}];

        Base bi = fld.g.NodeBase(inode), bj = fld.g.NodeBase(jnode);

        pt.pri[i] = bi;
        pt.pri[j] = bj;

        auto relax = [&](EnergyT e, std::initializer_list<DpIdx> idxs, int nts) {
            RelaxTrace(e, idxs, base_e, pt, nts);
        };

        if (i == j && inode != jnode) return;

        // Unpaired nucleotide
        if (i + 1 < j) {
            // Next nucleotide comes before j
            for (auto& [ip1node, bip1, ip1w] : fld.g.EdgesL2R(inode)) {
                relax(sr.MultMany(fld.ctx.em.MultiUnpaired(), ip1w,
                                  fld.ml_table[{br, jnode, ip1node}]),
                      {{DpTable::ML, ip1node, jnode, br}}, 1);
            }
        } else if (i + 1 == j) {
            // Next nucleotide is j
            relax(sr.MultMany(fld.ctx.em.MultiUnpaired(), fld.ml_table[{br, jnode, jnode}],
                              fld.sum_paths[{inode, jnode}]),
                  {{DpTable::ML, jnode, jnode, br}, {DpTable::SP, inode, jnode}}, 1);
        } else if (i == j && br == 0) {
            // i==j. Only allowed to be unpaired if br == 0
            relax(fld.ctx.em.MultiUnpaired(), {}, 1);
        }

        // Branches
        // k<j because we deal with the (i,j) branch as a special case
        for (int k = i + fld.ctx.min_hairpin_loop_size + 1; k < j; ++k) {
            for (int knode : fld.g.NtNodes(k)) {
                Base bk = fld.g.NodeBase(knode);
                BasePair bp = BasesToPair(bi, bk);
                if (bp == BasePair::NN) continue;
                for (auto [kp1node, bkp1, kp1w] : fld.g.EdgesL2R(knode)) {
                    relax(sr.MultMany(fld.p_table[{inode, knode}],
                                      fld.ml_table[{std::max(0, br - 1), jnode, kp1node}],
                                      fld.ctx.em.MultiBranch(Base::N, bp, Base::N), kp1w),
                          {{DpTable::P, inode, knode},
                           {DpTable::ML, kp1node, jnode, std::max(0, br - 1)}}, 2);
                }
            }
        }

        BasePair bp = BasesToPair(bi, bj);
        if (br <= 1 && bp != BasePair::NN) {
            // Special case for (i,j) branch. I.e., ending on a branch
            relax(sr.MultMany(fld.p_table[{inode, jnode}],
                              fld.ctx.em.MultiBranch(Base::N, bp, Base::N)),
                  {{DpTable::P, inode, jnode}}, 2);
        }
    }

    void TraceEl(const DpIdx& idx, PartialTrace& pt) {
        int inode = idx.i;
        int i = node_to_ntidx[inode];
        EnergyT base_e = fld.el_table[{inode}];
        Base bi = fld.g.NodeBase(inode);

        pt.pri[i] = bi;

        auto relax = [&](EnergyT e, std::initializer_list<DpIdx> idxs, int nts) {
            RelaxTrace(e, idxs, base_e, pt, nts);
        };

        if (i + 1 < fld.n) {
            // Unpaired nucleotide
            for (auto& [ip1node, bip1, ip1w] : fld.g.EdgesL2R(inode)) {
                relax(sr.MultMany(fld.el_table[{ip1node}], ip1w), {{DpTable::EL, ip1node}}, 1);
            }
        } else {
            // Base case
            relax(sr.One(), {}, 1);
        }

        // Branches
        for (int j = i + fld.ctx.min_hairpin_loop_size + 1; j < fld.n; ++j) {
            for (int jnode : fld.g.NtNodes(j)) {
                Base bj = fld.g.NodeBase(jnode);
                BasePair bp = BasesToPair(bi, bj);
                if (bp == BasePair::NN) continue;
                if (j + 1 < fld.n) {
                    for (auto& [jp1node, bjp1, jp1w] : fld.g.EdgesL2R(jnode)) {
                        relax(sr.MultMany(fld.el_table[{jp1node}],
                                          fld.ctx.em.ExtBranch(Base::N, bp, Base::N),
                                          fld.p_table[{inode, jnode}], jp1w),
                              {{DpTable::P, inode, jnode}, {DpTable::EL, jp1node}}, 2);
                    }
                } else {
                    relax(sr.MultMany(fld.ctx.em.ExtBranch(Base::N, bp, Base::N),
                                      fld.p_table[{inode, jnode}]),
                          {{DpTable::P, inode, jnode}}, 2);
                }
            }
        }
    }

    void TraceSumPaths(const DpIdx& idx, PartialTrace& pt) {
        int inode = idx.i, jnode = idx.j;
        int i = node_to_ntidx[inode], j = node_to_ntidx[jnode];

        EnergyT base_e = fld.sum_paths[{inode, jnode}];
        Base bi = fld.g.NodeBase(inode), bj = fld.g.NodeBase(jnode);

        pt.pri[i] = bi;
        pt.pri[j] = bj;

        auto relax = [&](EnergyT e, std::initializer_list<DpIdx> idxs) {
            RelaxTrace(e, idxs, base_e, pt, 0);
        };

        if (inode == jnode) {
            relax(sr.One(), {});
            return;
        }

        for (auto& e : fld.g.EdgesL2R(inode)) {
            relax(sr.MultMany(fld.sum_paths[{e.to, jnode}], e.weight),
                  {{DpTable::SP, e.to, jnode}});
        }
    }

    vector<int> node_to_ntidx;
    // Precalculated keep_chance^n
    vector<double> nt_keep_chance_n;
    double subopt_randomness, nt_keep_chance;
    int max_pq_size;
    std::multiset<PartialTrace> set_pq;
    EnergyT trace_threshold;
    GraphFold<EnergyModelT> fld;
    SemiringT sr;
    std::default_random_engine rng;
};
}  // namespace mwmrna::fold_codons

#endif  // MWMRNA_LIB_CODON_GRAPH_TRACE_HPP