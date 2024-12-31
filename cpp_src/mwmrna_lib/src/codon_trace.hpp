#ifndef MWMRNA_LIB_TRACE_CODONS_HPP
#define MWMRNA_LIB_TRACE_CODONS_HPP

#include "codon_fold.hpp"
#include "energy_model.hpp"
#include "rna.hpp"
#include "structure.hpp"

namespace mwmrna::fold_codons {

template <typename EnergyModelT>
struct TraceContext {
    using EnergyT = typename FoldContext<EnergyModelT>::EnergyT;
    using SemiringT = typename FoldContext<EnergyModelT>::SemiringT;
    int max_pq_size = -1;
    double subopt_randomness = 0.0;
    /// Use max value/3 as a safe maximum
    /// Some DP recursions rely on adding a valid energy to an "infinity" value to represent an
    /// invalid call. For example, the ML table.
    EnergyT threshold = SemiringT().Zero() / 3;
    int rng_seed = 0;
};

template <typename EnergyT>
struct Trace {
    Primary pri;
    Matching match;
    EnergyT energy;
};

template <typename EnergyModelT>
class TraceIterator {
   private:
    using EnergyT = typename FoldContext<EnergyModelT>::EnergyT;
    using SemiringT = typename FoldContext<EnergyModelT>::SemiringT;
    static_assert(std::is_base_of<MFESemiring<EnergyT>, SemiringT>::value,
                  "Semiring must be MFE to do tracebacks");

   public:
    TraceIterator(const Fold<EnergyModelT>& fld, const TraceContext<EnergyModelT>& ctx = {})
        : fld(fld) {
        Init(ctx);
    }

    TraceIterator(const Fold<EnergyModelT>&& fld, const TraceContext<EnergyModelT>& ctx = {})
        : fld(std::move(fld)) {
        Init(ctx);
    }

    bool HasNext() { return !set_pq.empty(); }

    Trace<EnergyT> Next() {
        while (1) {
            PartialTrace pt = *set_pq.begin();
            set_pq.erase(set_pq.begin());
            if (pt.idxs.empty()) {
                Trace<EnergyT> res{pt.pri, pt.match, pt.min_energy};
                return res;
            }
            DpIdx idx = pt.idxs.back();
            pt.idxs.pop_back();
            switch (idx.table) {
                case DpTable::P:
                    TraceP(idx, pt);
                    break;
                case DpTable::OMM:
                    TraceOmm(idx, pt);
                    break;
                case DpTable::EL:
                    TraceEl(idx, pt);
                    break;
                case DpTable::ML:
                    TraceMl(idx, pt);
                    break;
                case DpTable::IMM:
                    TraceImm(idx, pt);
                    break;
                case DpTable::UpCodons:
                    TraceUnpairedCodons(idx, pt);
                    break;
                case DpTable::SC_SpecialHairpin:
                    TraceScSpecialHairpin(idx, pt);
                    break;
                case DpTable::SC_UpCodon:
                    TraceScUpCodon(idx, pt);
                    break;
                default:
                    throw new std::runtime_error("Invalid table");
                    break;
            }
        }
    }

   private:
    void Init(const TraceContext<EnergyModelT>& ctx) {
        if (ctx.subopt_randomness < 0 || ctx.subopt_randomness > 1) {
            throw new std::runtime_error("subopt_randomness must be in [0, 1]");
        }
        this->rng.seed(ctx.rng_seed);
        this->subopt_randomness = ctx.subopt_randomness;
        InitTrace(ctx.threshold, ctx.max_pq_size);
    }

    /// @brief Must be called before Next/HasNext
    /// @param threshold only trace sequence-structures with energy < threshold
    void InitTrace(EnergyT threshold, int _max_pq_size = -1) {
        trace_threshold = threshold;
        if (fld.el_table[{0, fld.tbd_cdn}] < threshold) {
            EnergyT energy = fld.el_table[{0, fld.tbd_cdn}];
            set_pq.insert({.pri = Primary(fld.n, Base::N),
                           .match = Matching(fld.n, Matching::kUnmatched),
                           .idxs = {{.table = DpTable::EL, .i = 0, .cdni = fld.tbd_cdn}},
                           .min_energy = energy,
                           .priority = energy});
        }
        max_pq_size = _max_pq_size;
    }
    // These map onto the fld.ml_table, fld.el_table, fld.omm_table, fld.p_table,
    // fld.unpaired_codons IMM maps onto inner matches. We don't actually store a table for these,
    // but it is an implicit DP table. SC_ represent special cases for dealing with non-table events
    // These include emitting a specific unpaired codons, or a special hairpin respectively
    enum class DpTable { ML, EL, OMM, P, UpCodons, IMM, SC_UpCodon, SC_SpecialHairpin };

    struct DpIdx {
        DpTable table;
        int br{-1}, i{-1}, j{-1}, cdni{-1}, cdnj{-1};
    };

    struct PartialTrace {
        Primary pri;
        Matching match;
        std::vector<DpIdx> idxs;
        EnergyT min_energy;
        double priority;
        bool operator<(const PartialTrace& other) const { return priority < other.priority; }
        bool operator>(const PartialTrace& other) const { return priority > other.priority; }
    };

    void RelaxTrace(EnergyT e, std::initializer_list<DpIdx> idxs, EnergyT base_e,
                    const PartialTrace& pt) {
        if (e >= trace_threshold) return;
        EnergyT min_energy = pt.min_energy + e - base_e;
        std::uniform_real_distribution<double> dist(1.0 - this->subopt_randomness, 1.0);
        double priority = min_energy * dist(rng);
        // Avoid an expensive copy if possible
        if (max_pq_size != -1 && static_cast<int>(set_pq.size()) >= max_pq_size &&
            priority >= set_pq.rbegin()->priority)
            return;
        PartialTrace npt = pt;
        npt.min_energy = min_energy;
        npt.priority = priority;
        for (const auto& nidx : idxs) npt.idxs.push_back(nidx);
        set_pq.insert(std::move(npt));
        while (max_pq_size != -1 && static_cast<int>(set_pq.size()) > max_pq_size) {
            set_pq.erase(std::prev(set_pq.end()));
        }
    }

    void AddCodon(Primary& pri, int i, int cdni) {
        for (int k = i / 3 * 3; k < i / 3 * 3 + 3; ++k) pri[k] = fld.idx2base[{k, cdni}];
    }

    void TraceP(const DpIdx& idx, PartialTrace& pt) {
        int i = idx.i, j = idx.j, cdni = idx.cdni, cdnj = idx.cdnj;
        BasePair bp = BasesToPair(fld.idx2base[{i, cdni}], fld.idx2base[{j, cdnj}]);
        EnergyT base_e = fld.p_table[{i, j, cdni, cdnj}];

        pt.match[i] = j;
        pt.match[j] = i;
        // These will sometimes lead to redundant adds, but that's fine
        AddCodon(pt.pri, i, cdni);
        AddCodon(pt.pri, j, cdnj);

        auto relax = [&](EnergyT e, std::initializer_list<DpIdx> idxs) {
            RelaxTrace(e, idxs, base_e, pt);
        };

        // Stacks
        // We need to try all codons for i+1 and j-1 if they are in different codons
        // to i and j. Note that we can ignore cases where i+1=j-1, since these
        // are invalid for min_hairpin_loop_size >= 3 and the paired table lookup
        // will be zero.
        if (i / 3 != (i + 1) / 3) {
            for (int cdnip1 = 0; cdnip1 < fld.num_codons[i + 1]; ++cdnip1) {
                if (j / 3 != (j - 1) / 3) {
                    for (int cdnjm1 = 0; cdnjm1 < fld.num_codons[j - 1]; ++cdnjm1) {
                        relax(sr.MultMany(
                                  fld.idx_cdn_score[{i + 1, cdnip1}],
                                  fld.idx_cdn_score[{j - 1, cdnjm1}],
                                  fld.p_table[{i + 1, j - 1, cdnip1, cdnjm1}],
                                  fld.ctx.em.Stack(bp, BasesToPair(fld.idx2base[{i + 1, cdnip1}],
                                                                   fld.idx2base[{j - 1, cdnjm1}]))),
                              {{.table = DpTable::P,
                                .i = i + 1,
                                .j = j - 1,
                                .cdni = cdnip1,
                                .cdnj = cdnjm1}});
                    }
                } else {
                    relax(
                        sr.MultMany(fld.idx_cdn_score[{i + 1, cdnip1}],
                                    fld.p_table[{i + 1, j - 1, cdnip1, cdnj}],
                                    fld.ctx.em.Stack(bp, BasesToPair(fld.idx2base[{i + 1, cdnip1}],
                                                                     fld.idx2base[{j - 1, cdnj}]))),
                        {{.table = DpTable::P,
                          .i = i + 1,
                          .j = j - 1,
                          .cdni = cdnip1,
                          .cdnj = cdnj}});
                }
            }
        } else {
            if (j / 3 != (j - 1) / 3) {
                for (int cdnjm1 = 0; cdnjm1 < fld.num_codons[j - 1]; ++cdnjm1) {
                    relax(sr.MultMany(
                              fld.idx_cdn_score[{j - 1, cdnjm1}],
                              fld.p_table[{i + 1, j - 1, cdni, cdnjm1}],
                              fld.ctx.em.Stack(bp, BasesToPair(fld.idx2base[{i + 1, cdni}],
                                                               fld.idx2base[{j - 1, cdnjm1}]))),
                          {{.table = DpTable::P,
                            .i = i + 1,
                            .j = j - 1,
                            .cdni = cdni,
                            .cdnj = cdnjm1}});
                }
            } else {
                relax(sr.MultMany(fld.p_table[{i + 1, j - 1, cdni, cdnj}],
                                  fld.ctx.em.Stack(bp, BasesToPair(fld.idx2base[{i + 1, cdni}],
                                                                   fld.idx2base[{j - 1, cdnj}]))),
                      {{.table = DpTable::P, .i = i + 1, .j = j - 1, .cdni = cdni, .cdnj = cdnj}});
            }
        }

        // Bulges
        // 5' bulges
        for (int k = i + 2; (j - 1) - k - 1 >= fld.ctx.min_hairpin_loop_size &&
                            k - i - 1 <= fld.ctx.max_two_loop_size;
             ++k) {
            int cdnjm1_lb = 0, cdnjm1_ub = fld.num_codons[j - 1];
            int cdnk_lb = 0, cdnk_ub = fld.num_codons[k];
            if (j / 3 == (j - 1) / 3) {
                cdnjm1_lb = cdnj;
                cdnjm1_ub = cdnj + 1;
            }
            if (k / 3 == i / 3) {
                cdnk_lb = cdni;
                cdnk_ub = cdni + 1;
            }
            EnergyT extra_codons = sr.One();
            if (i / 3 + 1 <= k / 3 - 1) extra_codons = fld.unpaired_codons[{i / 3 + 1, k / 3 - 1}];
            for (int cdnjm1 = cdnjm1_lb; cdnjm1 < cdnjm1_ub; ++cdnjm1) {
                EnergyT cdn_sc_jm1 =
                    j / 3 == (j - 1) / 3 ? sr.One() : fld.idx_cdn_score[{j - 1, cdnjm1}];
                for (int cdnk = cdnk_lb; cdnk < cdnk_ub; ++cdnk) {
                    EnergyT cdn_sc_k = k / 3 == i / 3 ? sr.One() : fld.idx_cdn_score[{k, cdnk}];
                    BasePair bpkl =
                        BasesToPair(fld.idx2base[{k, cdnk}], fld.idx2base[{j - 1, cdnjm1}]);
                    if (bpkl == BasePair::NN) continue;

                    relax(sr.MultMany(extra_codons, cdn_sc_k, cdn_sc_jm1,
                                      fld.p_table[{k, j - 1, cdnk, cdnjm1}],
                                      fld.ctx.em.Bulge(bp, bpkl, k - i - 1)),
                          {{.table = DpTable::P, .i = k, .j = j - 1, .cdni = cdnk, .cdnj = cdnjm1},
                           {.table = DpTable::UpCodons, .i = i / 3 + 1, .j = k / 3 - 1}});
                }
            }
        }
        // 3' bulges
        for (int k = j - 2; j - k - 1 <= fld.ctx.max_two_loop_size &&
                            k - (i + 1) - 1 >= fld.ctx.min_hairpin_loop_size;
             --k) {
            int cdnip1_lb = 0, cdnip1_ub = fld.num_codons[i + 1];
            int cdnk_lb = 0, cdnk_ub = fld.num_codons[k];
            if ((i + 1) / 3 == i / 3) {
                cdnip1_lb = cdni;
                cdnip1_ub = cdni + 1;
            }
            if (k / 3 == j / 3) {
                cdnk_lb = cdnj;
                cdnk_ub = cdnj + 1;
            }
            EnergyT extra_codons = sr.One();
            if (k / 3 + 1 <= j / 3 - 1) extra_codons = fld.unpaired_codons[{k / 3 + 1, j / 3 - 1}];
            for (int cdnip1 = cdnip1_lb; cdnip1 < cdnip1_ub; ++cdnip1) {
                EnergyT cdn_sc_ip1 =
                    (i + 1) / 3 == i / 3 ? sr.One() : fld.idx_cdn_score[{i + 1, cdnip1}];
                for (int cdnk = cdnk_lb; cdnk < cdnk_ub; ++cdnk) {
                    EnergyT cdn_sc_k = k / 3 == j / 3 ? sr.One() : fld.idx_cdn_score[{k, cdnk}];
                    BasePair bpkl =
                        BasesToPair(fld.idx2base[{i + 1, cdnip1}], fld.idx2base[{k, cdnk}]);
                    if (bpkl == BasePair::NN) continue;
                    relax(sr.MultMany(extra_codons, cdn_sc_k, cdn_sc_ip1,
                                      fld.p_table[{i + 1, k, cdnip1, cdnk}],
                                      fld.ctx.em.Bulge(bp, bpkl, j - k - 1)),
                          {{.table = DpTable::P, .i = i + 1, .j = k, .cdni = cdnip1, .cdnj = cdnk},
                           {.table = DpTable::UpCodons, .i = k / 3 + 1, .j = j / 3 - 1}});
                }
            }
        }

        // Internal loops

        // 1x1
        int cdnip2_lb = 0, cdnip2_ub = fld.num_codons[i + 2];
        int cdnjm2_lb = 0, cdnjm2_ub = fld.num_codons[j - 2];
        if ((i + 2) / 3 == i / 3) {
            cdnip2_lb = cdni;
            cdnip2_ub = cdni + 1;
        }
        if ((j - 2) / 3 == j / 3) {
            cdnjm2_lb = cdnj;
            cdnjm2_ub = cdnj + 1;
        }
        for (int cdnip2 = cdnip2_lb; cdnip2 < cdnip2_ub; ++cdnip2) {
            EnergyT cdn_sc_ip2 =
                (i + 2) / 3 == i / 3 ? sr.One() : fld.idx_cdn_score[{i + 2, cdnip2}];
            for (int cdnjm2 = cdnjm2_lb; cdnjm2 < cdnjm2_ub; ++cdnjm2) {
                EnergyT cdn_sc_jm2 =
                    (j - 2) / 3 == j / 3 ? sr.One() : fld.idx_cdn_score[{j - 2, cdnjm2}];
                BasePair bpkl =
                    BasesToPair(fld.idx2base[{i + 2, cdnip2}], fld.idx2base[{j - 2, cdnjm2}]);
                if (bpkl == BasePair::NN) continue;
                Base bip1 = fld.idx2base[{i + 1, (i + 1) / 3 == i / 3 ? cdni : cdnip2}];
                Base bjm1 = fld.idx2base[{j - 1, (j - 1) / 3 == j / 3 ? cdnj : cdnjm2}];
                relax(sr.MultMany(fld.p_table[{i + 2, j - 2, cdnip2, cdnjm2}], cdn_sc_ip2,
                                  cdn_sc_jm2, fld.ctx.em.IL11(bp, bpkl, bip1, bjm1)),
                      {{.table = DpTable::P,
                        .i = i + 2,
                        .j = j - 2,
                        .cdni = cdnip2,
                        .cdnj = cdnjm2}});
            }
        }

        // 1x2
        for (int cdnip2 = cdnip2_lb; cdnip2 < cdnip2_ub; ++cdnip2) {
            EnergyT cdn_sc_ip2 =
                (i + 2) / 3 == i / 3 ? sr.One() : fld.idx_cdn_score[{i + 2, cdnip2}];
            for (int cdnjm3 = 0; cdnjm3 < fld.num_codons[j - 3]; ++cdnjm3) {
                EnergyT cdn_sc_jm3 = fld.idx_cdn_score[{j - 3, cdnjm3}];
                Base bip1 = fld.idx2base[{i + 1, (i + 1) / 3 == i / 3 ? cdni : cdnip2}];
                Base bjm1 = fld.idx2base[{j - 1, (j - 1) / 3 == j / 3 ? cdnj : cdnjm3}];
                Base bjm2 = fld.idx2base[{j - 2, (j - 2) / 3 == j / 3 ? cdnj : cdnjm3}];
                BasePair bpkl =
                    BasesToPair(fld.idx2base[{i + 2, cdnip2}], fld.idx2base[{j - 3, cdnjm3}]);
                if (bpkl == BasePair::NN) continue;
                relax(sr.MultMany(fld.p_table[{i + 2, j - 3, cdnip2, cdnjm3}], cdn_sc_ip2,
                                  cdn_sc_jm3, fld.ctx.em.IL12(bp, bpkl, bip1, bjm1, bjm2)),
                      {{.table = DpTable::P,
                        .i = i + 2,
                        .j = j - 3,
                        .cdni = cdnip2,
                        .cdnj = cdnjm3}});
            }
        }

        // 2x1
        for (int cdnip3 = 0; cdnip3 < fld.num_codons[i + 3]; ++cdnip3) {
            EnergyT cdn_sc_ip3 = fld.idx_cdn_score[{i + 3, cdnip3}];
            Base bip1 = fld.idx2base[{i + 1, (i + 1) / 3 == i / 3 ? cdni : cdnip3}];
            Base bip2 = fld.idx2base[{i + 2, (i + 2) / 3 == i / 3 ? cdni : cdnip3}];
            for (int cdnjm2 = cdnjm2_lb; cdnjm2 < cdnjm2_ub; ++cdnjm2) {
                EnergyT cdn_sc_jm2 =
                    (j - 2) / 3 == j / 3 ? sr.One() : fld.idx_cdn_score[{j - 2, cdnjm2}];
                Base bjm1 = fld.idx2base[{j - 1, (j - 1) / 3 == j / 3 ? cdnj : cdnjm2}];
                BasePair bpkl =
                    BasesToPair(fld.idx2base[{i + 3, cdnip3}], fld.idx2base[{j - 2, cdnjm2}]);
                if (bpkl == BasePair::NN) continue;
                relax(sr.MultMany(fld.p_table[{i + 3, j - 2, cdnip3, cdnjm2}], cdn_sc_ip3,
                                  cdn_sc_jm2, fld.ctx.em.IL21(bp, bpkl, bip1, bjm1, bip2)),
                      {{.table = DpTable::P,
                        .i = i + 3,
                        .j = j - 2,
                        .cdni = cdnip3,
                        .cdnj = cdnjm2}});
            }
        }

        // 2x2
        for (int cdnip3 = 0; cdnip3 < fld.num_codons[i + 3]; ++cdnip3) {
            EnergyT cdn_sc_ip3 = fld.idx_cdn_score[{i + 3, cdnip3}];
            Base bip1 = fld.idx2base[{i + 1, (i + 1) / 3 == i / 3 ? cdni : cdnip3}];
            Base bip2 = fld.idx2base[{i + 2, (i + 2) / 3 == i / 3 ? cdni : cdnip3}];

            for (int cdnjm3 = 0; cdnjm3 < fld.num_codons[j - 3]; ++cdnjm3) {
                EnergyT cdn_sc_jm3 = fld.idx_cdn_score[{j - 3, cdnjm3}];
                Base bjm1 = fld.idx2base[{j - 1, (j - 1) / 3 == j / 3 ? cdnj : cdnjm3}];
                Base bjm2 = fld.idx2base[{j - 2, (j - 2) / 3 == j / 3 ? cdnj : cdnjm3}];
                BasePair bpkl =
                    BasesToPair(fld.idx2base[{i + 3, cdnip3}], fld.idx2base[{j - 3, cdnjm3}]);
                if (bpkl == BasePair::NN) continue;
                relax(sr.MultMany(fld.p_table[{i + 3, j - 3, cdnip3, cdnjm3}], cdn_sc_ip3,
                                  cdn_sc_jm3, fld.ctx.em.IL22(bp, bpkl, bip1, bjm1, bip2, bjm2)),
                      {{.table = DpTable::P,
                        .i = i + 3,
                        .j = j - 3,
                        .cdni = cdnip3,
                        .cdnj = cdnjm3}});
            }
        }

        // 2x3
        for (int cdnip3 = 0; cdnip3 < fld.num_codons[i + 3]; ++cdnip3) {
            EnergyT cdn_sc_ip3 = fld.idx_cdn_score[{i + 3, cdnip3}];
            Base bip1 = fld.idx2base[{i + 1, (i + 1) / 3 == i / 3 ? cdni : cdnip3}];
            Base bip2 = fld.idx2base[{i + 2, (i + 2) / 3 == i / 3 ? cdni : cdnip3}];
            for (int cdnjm4 = 0; cdnjm4 < fld.num_codons[j - 4]; ++cdnjm4) {
                EnergyT cdn_sc_jm4 = fld.idx_cdn_score[{j - 4, cdnjm4}];
                Base bjm1, bjm3;
                BasePair bpkl =
                    BasesToPair(fld.idx2base[{i + 3, cdnip3}], fld.idx2base[{j - 4, cdnjm4}]);
                if (bpkl == BasePair::NN) continue;
                EnergyT p = fld.p_table[{i + 3, j - 4, cdnip3, cdnjm4}];
                // This can only happen when the 3nt unpaired region is a single codon
                if ((j - 1) / 3 == (j - 3) / 3) {
                    for (int cdnjm1 = 0; cdnjm1 < fld.num_codons[j - 1]; ++cdnjm1) {
                        EnergyT cdn_sc_jm1 = fld.idx_cdn_score[{j - 1, cdnjm1}];
                        bjm1 = fld.idx2base[{j - 1, cdnjm1}];
                        bjm3 = fld.idx2base[{j - 3, cdnjm1}];
                        relax(sr.MultMany(p, cdn_sc_ip3, cdn_sc_jm4, cdn_sc_jm1,
                                          fld.ctx.em.IL23(bp, bpkl, bip1, bjm1, bip2, bjm3)),
                              {{.table = DpTable::P,
                                .i = i + 3,
                                .j = j - 4,
                                .cdni = cdnip3,
                                .cdnj = cdnjm4},
                               {.table = DpTable::SC_UpCodon, .i = j - 1, .cdni = cdnjm1}});
                    }
                } else {
                    // TODO: There probably doesn't need to be a branch here
                    bjm1 = fld.idx2base[{j - 1, (j - 1) / 3 == j / 3 ? cdnj : cdnjm4}];
                    bjm3 = fld.idx2base[{j - 3, (j - 3) / 3 == j / 3 ? cdnj : cdnjm4}];

                    relax(sr.MultMany(p, cdn_sc_ip3, cdn_sc_jm4,
                                      fld.ctx.em.IL23(bp, bpkl, bip1, bjm1, bip2, bjm3)),
                          {{.table = DpTable::P,
                            .i = i + 3,
                            .j = j - 4,
                            .cdni = cdnip3,
                            .cdnj = cdnjm4}});
                }
            }
        }

        // 3x2
        for (int cdnjm3 = 0; cdnjm3 < fld.num_codons[j - 3]; ++cdnjm3) {
            EnergyT cdn_sc_jm3 = fld.idx_cdn_score[{j - 3, cdnjm3}];
            Base bjm1 = fld.idx2base[{j - 1, (j - 1) / 3 == j / 3 ? cdnj : cdnjm3}];
            Base bjm2 = fld.idx2base[{j - 2, (j - 2) / 3 == j / 3 ? cdnj : cdnjm3}];
            for (int cdnip4 = 0; cdnip4 < fld.num_codons[i + 4]; ++cdnip4) {
                EnergyT cdn_sc_ip4 = fld.idx_cdn_score[{i + 4, cdnip4}];
                Base bip1, bip3;
                BasePair bpkl =
                    BasesToPair(fld.idx2base[{i + 4, cdnip4}], fld.idx2base[{j - 3, cdnjm3}]);
                if (bpkl == BasePair::NN) continue;
                EnergyT p = fld.p_table[{i + 4, j - 3, cdnip4, cdnjm3}];
                // This can only happen when the 3nt unpaired region is a single codon
                if ((i + 1) / 3 == (i + 3) / 3) {
                    for (int cdnip1 = 0; cdnip1 < fld.num_codons[i + 1]; ++cdnip1) {
                        EnergyT cdn_sc_ip1 = fld.idx_cdn_score[{i + 1, cdnip1}];
                        bip1 = fld.idx2base[{i + 1, cdnip1}];
                        bip3 = fld.idx2base[{i + 3, cdnip1}];

                        relax(sr.MultMany(p, cdn_sc_ip1, cdn_sc_ip4, cdn_sc_jm3,
                                          fld.ctx.em.IL32(bp, bpkl, bip1, bjm1, bip3, bjm2)),
                              {{.table = DpTable::P,
                                .i = i + 4,
                                .j = j - 3,
                                .cdni = cdnip4,
                                .cdnj = cdnjm3},
                               {.table = DpTable::SC_UpCodon, .i = i + 1, .cdni = cdnip1}});
                    }
                } else {
                    // TODO: There probably doesn't need to be a branch here
                    bip1 = fld.idx2base[{i + 1, (i + 1) / 3 == i / 3 ? cdni : cdnip4}];
                    bip3 = fld.idx2base[{i + 3, (i + 3) / 3 == i / 3 ? cdni : cdnip4}];

                    relax(sr.MultMany(p, cdn_sc_ip4, cdn_sc_jm3,
                                      fld.ctx.em.IL32(bp, bpkl, bip1, bjm1, bip3, bjm2)),
                          {{.table = DpTable::P,
                            .i = i + 4,
                            .j = j - 3,
                            .cdni = cdnip4,
                            .cdnj = cdnjm3}});
                }
            }
        }

        // 1xn
        for (int k = j - 4;
             j - k <= fld.ctx.max_two_loop_size && k - (i + 2) - 1 >= fld.ctx.min_hairpin_loop_size;
             --k) {
            // Note for j-4: 1x1 and 1x2 are handled above
            // cdnk cannot be the same as cdnj since we know there are >2 nts
            int cdnk_lb = 0, cdnk_ub = fld.num_codons[k];
            int cdnkp1_lb = 0, cdnkp1_ub = fld.num_codons[k + 1];

            for (int cdnip2 = cdnip2_lb; cdnip2 < cdnip2_ub; ++cdnip2) {
                EnergyT cdn_sc_ip2 =
                    (i + 2) / 3 == i / 3 ? sr.One() : fld.idx_cdn_score[{i + 2, cdnip2}];
                int cdnip1 = i / 3 == (i + 1) / 3 ? cdni : cdnip2;
                for (int cdnk = cdnk_lb; cdnk < cdnk_ub; ++cdnk) {
                    EnergyT cdn_sc_k = fld.idx_cdn_score[{k, cdnk}];
                    BasePair bpkl =
                        BasesToPair(fld.idx2base[{i + 2, cdnip2}], fld.idx2base[{k, cdnk}]);
                    if (bpkl == BasePair::NN) continue;
                    if ((k + 1) / 3 == k / 3) {
                        cdnkp1_lb = cdnk;
                        cdnkp1_ub = cdnk + 1;
                    }
                    for (int cdnkp1 = cdnkp1_lb; cdnkp1 < cdnkp1_ub; ++cdnkp1) {
                        EnergyT cdn_sc_kp1 =
                            (k + 1) / 3 == k / 3 ? sr.One() : fld.idx_cdn_score[{k + 1, cdnkp1}];
                        int cdnjm1_lb = 0, cdnjm1_ub = fld.num_codons[j - 1];
                        if ((j - 1) / 3 == j / 3) {
                            cdnjm1_lb = cdnj;
                            cdnjm1_ub = cdnj + 1;
                        } else if ((j - 1) / 3 == (k + 1) / 3) {
                            cdnjm1_lb = cdnkp1;
                            cdnjm1_ub = cdnkp1 + 1;
                        }
                        for (int cdnjm1 = cdnjm1_lb; cdnjm1 < cdnjm1_ub; ++cdnjm1) {
                            EnergyT cdn_sc_jm1 = (j - 1) / 3 == j / 3 || (j - 1) / 3 == (k + 1) / 3
                                                     ? sr.One()
                                                     : fld.idx_cdn_score[{j - 1, cdnjm1}];
                            Base bip1 = fld.idx2base[{i + 1, cdnip1}];
                            Base bjm1 = fld.idx2base[{j - 1, cdnjm1}];
                            Base bkp1 = fld.idx2base[{k + 1, cdnkp1}];
                            EnergyT extra_codons = sr.One();
                            if ((k + 1) / 3 + 1 <= (j - 1) / 3 - 1)
                                extra_codons =
                                    fld.unpaired_codons[{(k + 1) / 3 + 1, (j - 1) / 3 - 1}];

                            // Note that k+1 and j-1 might be newly added codons, i+1 can never
                            // be So we add the SC_UpCodon for k+1 and j-1 Sometimes this is
                            // redundant, but that's fine
                            relax(
                                sr.MultMany(fld.p_table[{i + 2, k, cdnip2, cdnk}], extra_codons,
                                            cdn_sc_ip2, cdn_sc_k, cdn_sc_kp1, cdn_sc_jm1,
                                            fld.ctx.em.IL1N(bp, bpkl, bip1, bjm1, bkp1, j - k - 1)),
                                {{.table = DpTable::P,
                                  .i = i + 2,
                                  .j = k,
                                  .cdni = cdnip2,
                                  .cdnj = cdnk},
                                 {.table = DpTable::UpCodons,
                                  .i = (k + 1) / 3 + 1,
                                  .j = (j - 1) / 3 - 1},
                                 {.table = DpTable::SC_UpCodon, .i = j - 1, .cdni = cdnjm1},
                                 {.table = DpTable::SC_UpCodon, .i = k + 1, .cdni = cdnkp1}});
                        }
                    }
                }
            }
        }

        // nx1
        for (int k = i + 4;
             (j - 2) - k - 1 >= fld.ctx.min_hairpin_loop_size && k - i <= fld.ctx.max_two_loop_size;
             ++k) {
            int cdnk_lb = 0, cdnk_ub = fld.num_codons[k];
            int cdnkm1_lb = 0, cdnkm1_ub = fld.num_codons[k - 1];
            for (int cdnjm2 = cdnjm2_lb; cdnjm2 < cdnjm2_ub; ++cdnjm2) {
                EnergyT cdn_sc_jm2 =
                    (j - 2) / 3 == j / 3 ? sr.One() : fld.idx_cdn_score[{j - 2, cdnjm2}];
                int cdnjm1 = j / 3 == (j - 1) / 3 ? cdnj : cdnjm2;
                for (int cdnk = cdnk_lb; cdnk < cdnk_ub; ++cdnk) {
                    EnergyT cdn_sc_k = fld.idx_cdn_score[{k, cdnk}];
                    BasePair bpkl =
                        BasesToPair(fld.idx2base[{k, cdnk}], fld.idx2base[{j - 2, cdnjm2}]);
                    if (bpkl == BasePair::NN) continue;
                    if ((k - 1) / 3 == k / 3) {
                        cdnkm1_lb = cdnk;
                        cdnkm1_ub = cdnk + 1;
                    }
                    for (int cdnkm1 = cdnkm1_lb; cdnkm1 < cdnkm1_ub; ++cdnkm1) {
                        EnergyT cdn_sc_km1 =
                            (k - 1) / 3 == k / 3 ? sr.One() : fld.idx_cdn_score[{k - 1, cdnkm1}];
                        int cdnip1_lb = 0, cdnip1_ub = fld.num_codons[i + 1];
                        if ((i + 1) / 3 == i / 3) {
                            cdnip1_lb = cdni;
                            cdnip1_ub = cdni + 1;
                        } else if ((i + 1) / 3 == (k - 1) / 3) {
                            cdnip1_lb = cdnkm1;
                            cdnip1_ub = cdnkm1 + 1;
                        }
                        for (int cdnip1 = cdnip1_lb; cdnip1 < cdnip1_ub; ++cdnip1) {
                            EnergyT cdn_sc_ip1 = (i + 1) / 3 == i / 3 || (i + 1) / 3 == (k - 1) / 3
                                                     ? sr.One()
                                                     : fld.idx_cdn_score[{i + 1, cdnip1}];
                            Base bip1 = fld.idx2base[{i + 1, cdnip1}];
                            Base bjm1 = fld.idx2base[{j - 1, cdnjm1}];
                            Base bkm1 = fld.idx2base[{k - 1, cdnkm1}];
                            EnergyT extra_codons = sr.One();
                            if ((i + 1) / 3 + 1 <= (k - 1) / 3 - 1)
                                extra_codons =
                                    fld.unpaired_codons[{(i + 1) / 3 + 1, (k - 1) / 3 - 1}];

                            // Note that i+1 and k-1 might be newly added codons, j-1 can never
                            // be So we add the SC_UpCodon for k+1 and j-1 Sometimes this is
                            // redundant, but that's fine
                            relax(
                                sr.MultMany(fld.p_table[{k, j - 2, cdnk, cdnjm2}], extra_codons,
                                            cdn_sc_k, cdn_sc_km1, cdn_sc_jm2, cdn_sc_ip1,
                                            fld.ctx.em.ILN1(bp, bpkl, bip1, bjm1, bkm1, k - i - 1)),
                                {{.table = DpTable::P,
                                  .i = k,
                                  .j = j - 2,
                                  .cdni = cdnk,
                                  .cdnj = cdnjm2},
                                 {.table = DpTable::UpCodons,
                                  .i = (i + 1) / 3 + 1,
                                  .j = (k - 1) / 3 - 1},
                                 {.table = DpTable::SC_UpCodon, .i = i + 1, .cdni = cdnip1},
                                 {.table = DpTable::SC_UpCodon, .i = k - 1, .cdni = cdnkm1}});
                        }
                    }
                }
            }
        }

        // arbitrary internal loops

        // Precompute sum for inner mismatches
        // Used when inner mismatch codons are independent of outer mismatch codons
        int cdnip1_lb = 0, cdnip1_ub = fld.num_codons[i + 1];
        int cdnjm1_lb = 0, cdnjm1_ub = fld.num_codons[j - 1];
        if (i / 3 == (i + 1) / 3) {
            cdnip1_lb = cdni;
            cdnip1_ub = cdni + 1;
        }
        if (j / 3 == (j - 1) / 3) {
            cdnjm1_lb = cdnj;
            cdnjm1_ub = cdnj + 1;
        }
        EnergyT e_inner_mismatches = sr.Zero();
        for (int cdnip1 = cdnip1_lb; cdnip1 < cdnip1_ub; ++cdnip1) {
            EnergyT cdn_sc_ip1 =
                (i + 1) / 3 == i / 3 ? sr.One() : fld.idx_cdn_score[{i + 1, cdnip1}];
            for (int cdnjm1 = cdnjm1_lb; cdnjm1 < cdnjm1_ub; ++cdnjm1) {
                EnergyT cdn_sc_jm1 =
                    (j - 1) / 3 == j / 3 ? sr.One() : fld.idx_cdn_score[{j - 1, cdnjm1}];
                Base bip1 = fld.idx2base[{i + 1, cdnip1}];
                Base bjm1 = fld.idx2base[{j - 1, cdnjm1}];
                e_inner_mismatches = sr.Add(
                    e_inner_mismatches, sr.MultMany(cdn_sc_ip1, cdn_sc_jm1,
                                                    fld.ctx.em.ILInnerMismatch(bp, bip1, bjm1)));
            }
        }

        for (int k = i + 3; k - i - 1 <= fld.ctx.max_two_loop_size && k < j - 3; ++k) {
            for (int l = j - 3; j - l - 1 + k - i - 1 <= fld.ctx.max_two_loop_size &&
                                l - k - 1 >= fld.ctx.min_hairpin_loop_size;
                 --l) {
                int lup = k - i - 1, rup = j - l - 1;

                // Ignore 2x2 and 2x3 internal loops as these are handled above
                if ((lup == 3 && rup == 2) || (lup == 2 && rup == 3)) continue;
                if (lup == 2 && rup == 2) continue;

                int cdnip1_lb = 0, cdnip1_ub = fld.num_codons[i + 1];
                int cdnjm1_lb = 0, cdnjm1_ub = fld.num_codons[j - 1];
                if (i / 3 == (i + 1) / 3) {
                    cdnip1_lb = cdni;
                    cdnip1_ub = cdni + 1;
                }
                if (j / 3 == (j - 1) / 3) {
                    cdnjm1_lb = cdnj;
                    cdnjm1_ub = cdnj + 1;
                }

                // If at least one codon is shared, solve by enumerating all codon
                // options
                if ((k - 1) / 3 == (i + 1) / 3 || (l + 1) / 3 == (j - 1) / 3) {
                    int cdnk_lb = 0, cdnk_ub = fld.num_codons[k];
                    int cdnl_lb = 0, cdnl_ub = fld.num_codons[l];
                    for (int cdnk = cdnk_lb; cdnk < cdnk_ub; ++cdnk) {
                        EnergyT cdn_sc_k = fld.idx_cdn_score[{k, cdnk}];
                        for (int cdnl = cdnl_lb; cdnl < cdnl_ub; ++cdnl) {
                            EnergyT cdn_sc_l = fld.idx_cdn_score[{l, cdnl}];
                            BasePair bpkl =
                                BasesToPair(fld.idx2base[{k, cdnk}], fld.idx2base[{l, cdnl}]);
                            if (bpkl == BasePair::NN) continue;
                            if ((i + 1) / 3 == k / 3) {
                                cdnip1_lb = cdnk;
                                cdnip1_ub = cdnk + 1;
                            }
                            if ((j - 1) / 3 == l / 3) {
                                cdnjm1_lb = cdnl;
                                cdnjm1_ub = cdnl + 1;
                            }
                            for (int cdnip1 = cdnip1_lb; cdnip1 < cdnip1_ub; ++cdnip1) {
                                EnergyT cdn_sc_ip1 = (i + 1) / 3 == i / 3 || (i + 1) / 3 == k / 3
                                                         ? sr.One()
                                                         : fld.idx_cdn_score[{i + 1, cdnip1}];
                                for (int cdnjm1 = cdnjm1_lb; cdnjm1 < cdnjm1_ub; ++cdnjm1) {
                                    EnergyT cdn_sc_jm1 =
                                        (j - 1) / 3 == j / 3 || (j - 1) / 3 == l / 3
                                            ? sr.One()
                                            : fld.idx_cdn_score[{j - 1, cdnjm1}];
                                    if ((k - 1) / 3 == (i + 1) / 3 && (l + 1) / 3 == (j - 1) / 3) {
                                        // All mismatches share codons. Only happens for
                                        // 3x3s
                                        Base bip1 = fld.idx2base[{i + 1, cdnip1}];
                                        Base bjm1 = fld.idx2base[{j - 1, cdnjm1}];
                                        Base bkm1 = fld.idx2base[{k - 1, cdnip1}];
                                        Base blp1 = fld.idx2base[{l + 1, cdnjm1}];
                                        // No need for combos since this is only for 3x3
                                        relax(sr.MultMany(
                                                  fld.p_table[{k, l, cdnk, cdnl}], cdn_sc_ip1,
                                                  cdn_sc_jm1, cdn_sc_k, cdn_sc_l,
                                                  fld.ctx.em.NormalInternalLoop(
                                                      bp, bpkl, bip1, bjm1, bkm1, blp1, lup, rup)),
                                              {{.table = DpTable::P,
                                                .i = k,
                                                .j = l,
                                                .cdni = cdnk,
                                                .cdnj = cdnl},
                                               {.table = DpTable::SC_UpCodon,
                                                .i = i + 1,
                                                .cdni = cdnip1},
                                               {.table = DpTable::SC_UpCodon,
                                                .i = j - 1,
                                                .cdni = cdnjm1}});
                                    } else if ((k - 1) / 3 == (i + 1) / 3) {
                                        // Only 5' mismatch shares codons
                                        int cdnlp1_lb = 0, cdnlp1_ub = fld.num_codons[l + 1];
                                        if ((l + 1) / 3 == l / 3) {
                                            cdnlp1_lb = cdnl;
                                            cdnlp1_ub = cdnl + 1;
                                        }
                                        for (int cdnlp1 = cdnlp1_lb; cdnlp1 < cdnlp1_ub; ++cdnlp1) {
                                            EnergyT cdn_sc_lp1 =
                                                (l + 1) / 3 == l / 3
                                                    ? sr.One()
                                                    : fld.idx_cdn_score[{l + 1, cdnlp1}];
                                            Base bip1 = fld.idx2base[{i + 1, cdnip1}];
                                            Base bjm1 = fld.idx2base[{j - 1, cdnjm1}];
                                            Base bkm1 = fld.idx2base[{k - 1, cdnip1}];
                                            Base blp1 = fld.idx2base[{l + 1, cdnlp1}];
                                            EnergyT extra_codons = sr.One();
                                            // Only need to consider extra_codons on the
                                            // 3' side
                                            if ((l + 1) / 3 + 1 <= (j - 1) / 3 - 1)
                                                extra_codons = fld.unpaired_codons[{
                                                    (l + 1) / 3 + 1, (j - 1) / 3 - 1}];
                                            relax(sr.MultMany(fld.p_table[{k, l, cdnk, cdnl}],
                                                              extra_codons, cdn_sc_ip1, cdn_sc_jm1,
                                                              cdn_sc_k, cdn_sc_l, cdn_sc_lp1,
                                                              fld.ctx.em.NormalInternalLoop(
                                                                  bp, bpkl, bip1, bjm1, bkm1, blp1,
                                                                  lup, rup)),
                                                  {{.table = DpTable::P,
                                                    .i = k,
                                                    .j = l,
                                                    .cdni = cdnk,
                                                    .cdnj = cdnl},
                                                   {.table = DpTable::UpCodons,
                                                    .i = (l + 1) / 3 + 1,
                                                    .j = (j - 1) / 3 - 1},
                                                   {.table = DpTable::SC_UpCodon,
                                                    .i = i + 1,
                                                    .cdni = cdnip1},
                                                   {.table = DpTable::SC_UpCodon,
                                                    .i = l + 1,
                                                    .cdni = cdnlp1},
                                                   {.table = DpTable::SC_UpCodon,
                                                    .i = j - 1,
                                                    .cdni = cdnjm1}});
                                        }
                                    } else if ((l + 1) / 3 == (j - 1) / 3) {
                                        // Only 3' mismatch shares codons
                                        int cdnkm1_lb = 0, cdnkm1_ub = fld.num_codons[k - 1];
                                        if ((k - 1) / 3 == k / 3) {
                                            cdnkm1_lb = cdnk;
                                            cdnkm1_ub = cdnk + 1;
                                        }
                                        for (int cdnkm1 = cdnkm1_lb; cdnkm1 < cdnkm1_ub; ++cdnkm1) {
                                            EnergyT cdn_sc_km1 =
                                                (k - 1) / 3 == k / 3
                                                    ? sr.One()
                                                    : fld.idx_cdn_score[{k - 1, cdnkm1}];
                                            Base bip1 = fld.idx2base[{i + 1, cdnip1}];
                                            Base bjm1 = fld.idx2base[{j - 1, cdnjm1}];
                                            Base bkm1 = fld.idx2base[{k - 1, cdnkm1}];
                                            Base blp1 = fld.idx2base[{l + 1, cdnjm1}];
                                            EnergyT extra_codons = sr.One();
                                            // Only need to consider extra_codons on the
                                            // 5' side
                                            if ((i + 1) / 3 + 1 <= (k - 1) / 3 - 1)
                                                extra_codons = fld.unpaired_codons[{
                                                    (i + 1) / 3 + 1, (k - 1) / 3 - 1}];

                                            relax(sr.MultMany(fld.p_table[{k, l, cdnk, cdnl}],
                                                              extra_codons, cdn_sc_ip1, cdn_sc_jm1,
                                                              cdn_sc_k, cdn_sc_l, cdn_sc_km1,
                                                              fld.ctx.em.NormalInternalLoop(
                                                                  bp, bpkl, bip1, bjm1, bkm1, blp1,
                                                                  lup, rup)),
                                                  {{.table = DpTable::P,
                                                    .i = k,
                                                    .j = l,
                                                    .cdni = cdnk,
                                                    .cdnj = cdnl},
                                                   {.table = DpTable::UpCodons,
                                                    .i = (i + 1) / 3 + 1,
                                                    .j = (k - 1) / 3 - 1},
                                                   {.table = DpTable::SC_UpCodon,
                                                    .i = i + 1,
                                                    .cdni = cdnip1},
                                                   {.table = DpTable::SC_UpCodon,
                                                    .i = k - 1,
                                                    .cdni = cdnkm1},
                                                   {.table = DpTable::SC_UpCodon,
                                                    .i = j - 1,
                                                    .cdni = cdnjm1}});
                                        }
                                    } else {
                                        // No codons are shared
                                        // impossible to reach
                                        throw new std::runtime_error(
                                            "Impossible to reach no shared codons "
                                            "case");
                                    }
                                }
                            }
                        }
                    }
                } else {
                    // No codons shared between mismatches
                    EnergyT extra_codons = sr.One();
                    if ((i + 1) / 3 + 1 <= (k - 1) / 3 - 1)
                        extra_codons = fld.unpaired_codons[{(i + 1) / 3 + 1, (k - 1) / 3 - 1}];
                    if ((l + 1) / 3 + 1 <= (j - 1) / 3 - 1)
                        extra_codons = sr.MultMany(
                            extra_codons, fld.unpaired_codons[{(l + 1) / 3 + 1, (j - 1) / 3 - 1}]);
                    relax(
                        sr.MultMany(fld.omm_table[{k, l}], extra_codons, e_inner_mismatches,
                                    fld.ctx.em.ILInit(lup + rup), fld.ctx.em.ILAsymmetry(lup, rup)),
                        {{.table = DpTable::OMM, .i = k, .j = l},
                         {.table = DpTable::UpCodons, .i = (i + 1) / 3 + 1, .j = (k - 1) / 3 - 1},
                         {.table = DpTable::UpCodons, .i = (l + 1) / 3 + 1, .j = (j - 1) / 3 - 1},
                         {.table = DpTable::IMM, .i = i, .j = j, .cdni = cdni, .cdnj = cdnj}});
                }
            }
        }

        // Hairpin loops
        // To do hairpin loops we need to deal with all cases for i+1 and j-1 since
        // there can be dangles.
        // Some special cases are avoided since we assume min_hairpin_loop_size >= 3

        // Special logic for small hairpins (<=2 internal codons)
        // This handles special hairpins
        // ASSUMPTION: Special hairpins have length <= 6
        auto special_hairpins = fld.ctx.em.SpecialHairpins(j - i - 1);
        Primary pri;
        if (i / 3 + 1 == j / 3) {
            // 0 internal codons
            pri.assign(j - i + 1, Base::N);
            for (int k = i; k < i / 3 * 3 + 3; ++k) pri[k - i] = fld.idx2base[{k, cdni}];
            for (int k = i / 3 * 3 + 3; k <= j; ++k) pri[k - i] = fld.idx2base[{k, cdnj}];
            auto it = std::find_if(special_hairpins.begin(), special_hairpins.end(),
                                   [&pri](const auto& p) { return p.first == pri; });
            if (it != special_hairpins.end()) {
                relax(it->second,
                      {{.table = DpTable::SC_SpecialHairpin,
                        .i = i,
                        .j = j,
                        .cdni = static_cast<int>(std::distance(special_hairpins.begin(), it))}});
            } else {
                relax(fld.ctx.em.HairpinNormal(bp, fld.idx2base[{i + 1, cdni}],
                                               fld.idx2base[{j - 1, cdnj}], j - i - 1),
                      {});
            }
        } else if (i / 3 + 2 == j / 3) {
            // 1 internal codon
            pri.assign(j - i + 1, Base::N);
            for (int k = i; k < (i / 3) * 3 + 3; ++k) pri[k - i] = fld.idx2base[{k, cdni}];
            for (int k = j; k >= (j / 3) * 3; --k)
                pri[pri.size() - 1 - (j - k)] = fld.idx2base[{k, cdnj}];
            int cdnip1, cdnjm1;
            for (int cdn = 0; cdn < fld.num_codons[i / 3 * 3 + 3]; ++cdn) {
                EnergyT cdn_sc = fld.idx_cdn_score[{i / 3 * 3 + 3, cdn}];
                cdnip1 = cdni;
                cdnjm1 = cdnj;
                if ((i + 1) / 3 != i / 3) cdnip1 = cdn;
                if ((j - 1) / 3 != j / 3) cdnjm1 = cdn;
                for (int k = i / 3 * 3 + 3; k < i / 3 * 3 + 6; ++k)
                    pri[k - i] = fld.idx2base[{k, cdn}];
                auto it = std::find_if(special_hairpins.begin(), special_hairpins.end(),
                                       [&pri](const auto& p) { return p.first == pri; });
                if (it != special_hairpins.end()) {
                    relax(sr.MultMany(cdn_sc, it->second),
                          {{.table = DpTable::SC_SpecialHairpin,
                            .i = i,
                            .j = j,
                            .cdni = static_cast<int>(std::distance(special_hairpins.begin(), it))},
                           {.table = DpTable::SC_UpCodon, .i = i / 3 * 3 + 3, .cdni = cdn}});
                } else {
                    relax(sr.MultMany(cdn_sc, fld.ctx.em.HairpinNormal(
                                                  bp, fld.idx2base[{i + 1, cdnip1}],
                                                  fld.idx2base[{j - 1, cdnjm1}], j - i - 1)),
                          {{.table = DpTable::SC_UpCodon, .i = i / 3 * 3 + 3, .cdni = cdn}});
                }
            }
        } else if (i / 3 + 3 == j / 3) {
            // 2 internal codons
            pri.assign(j - i + 1, Base::N);
            for (int k = i; k < (i / 3) * 3 + 3; ++k) pri[k - i] = fld.idx2base[{k, cdni}];
            for (int k = j; k >= (j / 3) * 3; --k)
                pri[pri.size() - 1 - (j - k)] = fld.idx2base[{k, cdnj}];
            int cdnip1, cdnjm1;
            for (int cdnk = 0; cdnk < fld.num_codons[i / 3 * 3 + 3]; ++cdnk) {
                EnergyT cdn_sc_k = fld.idx_cdn_score[{i / 3 * 3 + 3, cdnk}];
                cdnip1 = cdni;
                if ((i + 1) / 3 != i / 3) cdnip1 = cdnk;
                for (int k = (i / 3) * 3 + 3; k < (i / 3) * 3 + 6; ++k)
                    pri[k - i] = fld.idx2base[{k, cdnk}];
                for (int cdnl = 0; cdnl < fld.num_codons[i / 3 * 3 + 6]; ++cdnl) {
                    EnergyT cdn_sc_l = fld.idx_cdn_score[{i / 3 * 3 + 6, cdnl}];
                    cdnjm1 = cdnj;
                    if ((j - 1) / 3 != j / 3) cdnjm1 = cdnl;
                    for (int l = (i / 3) * 3 + 6; l < (i / 3) * 3 + 9; ++l)
                        pri[l - i] = fld.idx2base[{l, cdnl}];
                    auto it = std::find_if(special_hairpins.begin(), special_hairpins.end(),
                                           [&pri](const auto& p) { return p.first == pri; });
                    if (it != special_hairpins.end()) {
                        relax(sr.MultMany(cdn_sc_k, cdn_sc_l, it->second),
                              {{.table = DpTable::SC_SpecialHairpin,
                                .i = i,
                                .j = j,
                                .cdni =
                                    static_cast<int>(std::distance(special_hairpins.begin(), it))},
                               {.table = DpTable::SC_UpCodon, .i = i / 3 * 3 + 3, .cdni = cdnk},
                               {.table = DpTable::SC_UpCodon, .i = i / 3 * 3 + 6, .cdni = cdnl}});
                    } else {
                        relax(sr.MultMany(cdn_sc_k, cdn_sc_l,
                                          fld.ctx.em.HairpinNormal(
                                              bp, fld.idx2base[{i + 1, cdnip1}],
                                              fld.idx2base[{j - 1, cdnjm1}], j - i - 1)),
                              {{.table = DpTable::SC_UpCodon, .i = i / 3 * 3 + 3, .cdni = cdnk},
                               {.table = DpTable::SC_UpCodon, .i = i / 3 * 3 + 6, .cdni = cdnl}});
                    }
                }
            }
        } else {
            // 3+ internal codons
            // Assumption: There are no special hairpins with 3+ internal codons
            EnergyT internal_combos = fld.unpaired_codons[{(i + 1) / 3 + 1, (j - 1) / 3 - 1}];

            if (i / 3 != (i + 1) / 3) {
                for (int cdnip1 = 0; cdnip1 < fld.num_codons[i + 1]; ++cdnip1) {
                    EnergyT cdn_sc_ip1 = fld.idx_cdn_score[{i + 1, cdnip1}];
                    if (j / 3 != (j - 1) / 3) {
                        // j-1 is in a different codon to j and i+1

                        for (int cdnjm1 = 0; cdnjm1 < fld.num_codons[j - 1]; ++cdnjm1) {
                            EnergyT cdn_sc_jm1 = fld.idx_cdn_score[{j - 1, cdnjm1}];
                            relax(sr.MultMany(internal_combos, cdn_sc_ip1, cdn_sc_jm1,
                                              fld.ctx.em.HairpinNormal(
                                                  bp, fld.idx2base[{i + 1, cdnip1}],
                                                  fld.idx2base[{j - 1, cdnjm1}], j - i - 1)),
                                  {{.table = DpTable::UpCodons,
                                    .i = (i + 1) / 3 + 1,
                                    .j = (j - 1) / 3 - 1},
                                   {.table = DpTable::SC_UpCodon, .i = i + 1, .cdni = cdnip1},
                                   {.table = DpTable::SC_UpCodon, .i = j - 1, .cdni = cdnjm1}});
                        }
                    } else if (j / 3 == (j - 1) / 3) {
                        // j-1 is in the same codon as j
                        relax(sr.MultMany(
                                  internal_combos, cdn_sc_ip1,
                                  fld.ctx.em.HairpinNormal(bp, fld.idx2base[{i + 1, cdnip1}],
                                                           fld.idx2base[{j - 1, cdnj}], j - i - 1)),
                              {{.table = DpTable::UpCodons,
                                .i = (i + 1) / 3 + 1,
                                .j = (j - 1) / 3 - 1},
                               {.table = DpTable::SC_UpCodon, .i = i + 1, .cdni = cdnip1}});
                    }
                }
            } else {
                if (j / 3 == (j - 1) / 3) {
                    // j-1 is in the same codon as j
                    relax(
                        sr.MultMany(internal_combos, fld.ctx.em.HairpinNormal(
                                                         bp, fld.idx2base[{i + 1, cdni}],
                                                         fld.idx2base[{j - 1, cdnj}], j - i - 1)),
                        {{.table = DpTable::UpCodons, .i = (i + 1) / 3 + 1, .j = (j - 1) / 3 - 1}});
                } else {
                    // j-1 is in a different codon to j and i+1
                    for (int cdnjm1 = 0; cdnjm1 < fld.num_codons[j - 1]; ++cdnjm1) {
                        EnergyT cdn_sc_jm1 = fld.idx_cdn_score[{j - 1, cdnjm1}];
                        relax(sr.MultMany(internal_combos, cdn_sc_jm1,
                                          fld.ctx.em.HairpinNormal(bp, fld.idx2base[{i + 1, cdni}],
                                                                   fld.idx2base[{j - 1, cdnjm1}],
                                                                   j - i - 1)),
                              {{.table = DpTable::UpCodons,
                                .i = (i + 1) / 3 + 1,
                                .j = (j - 1) / 3 - 1},
                               {.table = DpTable::SC_UpCodon, .i = j - 1, .cdni = cdnjm1}});
                    }
                }
            }
        }

        // Multi loops
        cdnip1_lb = 0, cdnip1_ub = fld.num_codons[i + 1];
        cdnjm1_lb = 0, cdnjm1_ub = fld.num_codons[j - 1];
        // (i+1)/3 can never equal (j-1)/3 since multiloops contain >= 2 branches and
        // min hairpin loop size is >= 3
        if ((i + 1) / 3 == i / 3) {
            cdnip1_lb = cdni;
            cdnip1_ub = cdni + 1;
        }
        if ((j - 1) / 3 == j / 3) {
            cdnjm1_lb = cdnj;
            cdnjm1_ub = cdnj + 1;
        }
        for (int cdnip1 = cdnip1_lb; cdnip1 < cdnip1_ub; ++cdnip1) {
            EnergyT cdn_sc_ip1 =
                (i + 1) / 3 == i / 3 ? sr.One() : fld.idx_cdn_score[{i + 1, cdnip1}];
            // We don't bother dealing with the case where i+1 and j-1 are in the same
            // codon since the ML table will be invalid
            for (int cdnjm1 = cdnjm1_lb; cdnjm1 < cdnjm1_ub; ++cdnjm1) {
                EnergyT cdn_sc_jm1 =
                    (j - 1) / 3 == j / 3 ? sr.One() : fld.idx_cdn_score[{j - 1, cdnjm1}];
                relax(sr.MultMany(fld.ml_table[{2, j - 1, i + 1, cdnip1, cdnjm1}], cdn_sc_ip1,
                                  cdn_sc_jm1, fld.ctx.em.MultiClosing(bp, Base::N, Base::N)),
                      {{.table = DpTable::ML,
                        .br = 2,
                        .i = i + 1,
                        .j = j - 1,
                        .cdni = cdnip1,
                        .cdnj = cdnjm1}});
            }
        }
    }

    void TraceOmm(const DpIdx& idx, const PartialTrace& pt) {
        int i = idx.i, j = idx.j;

        EnergyT base_e = fld.omm_table[{i, j}];

        auto relax_mismatch = [&](EnergyT e, std::initializer_list<DpIdx> idxs, int cdnim1,
                                  int cdnjp1) {
            auto new_pt = pt;
            AddCodon(new_pt.pri, i - 1, cdnim1);
            AddCodon(new_pt.pri, j + 1, cdnjp1);
            RelaxTrace(e, idxs, base_e, new_pt);
        };

        // Avoid edge cases that cannot be outer mismatches anyway
        if (i == 0 || j == fld.n - 1) return;
        for (int cdni = 0; cdni < fld.num_codons[i]; ++cdni) {
            EnergyT cdn_sc_i = fld.idx_cdn_score[{i, cdni}];
            for (int cdnj = 0; cdnj < fld.num_codons[j]; ++cdnj) {
                EnergyT cdn_sc_j = fld.idx_cdn_score[{j, cdnj}];
                BasePair bp = BasesToPair(fld.idx2base[{i, cdni}], fld.idx2base[{j, cdnj}]);
                if (bp == BasePair::NN) continue;
                int cdnim1_lb = 0, cdnim1_ub = fld.num_codons[i - 1];
                int cdnjp1_lb = 0, cdnjp1_ub = fld.num_codons[j + 1];
                if (i / 3 == (i - 1) / 3) {
                    cdnim1_lb = cdni;
                    cdnim1_ub = cdni + 1;
                }
                if (j / 3 == (j + 1) / 3) {
                    cdnjp1_lb = cdnj;
                    cdnjp1_ub = cdnj + 1;
                }

                for (int cdnim1 = cdnim1_lb; cdnim1 < cdnim1_ub; ++cdnim1) {
                    EnergyT cdn_sc_im1 =
                        (i - 1) / 3 == i / 3 ? sr.One() : fld.idx_cdn_score[{i - 1, cdnim1}];
                    for (int cdnjp1 = cdnjp1_lb; cdnjp1 < cdnjp1_ub; ++cdnjp1) {
                        EnergyT cdn_sc_jp1 =
                            (j + 1) / 3 == j / 3 ? sr.One() : fld.idx_cdn_score[{j + 1, cdnjp1}];
                        Base bim1 = fld.idx2base[{i - 1, cdnim1}];
                        Base bjp1 = fld.idx2base[{j + 1, cdnjp1}];
                        relax_mismatch(
                            sr.MultMany(fld.p_table[{i, j, cdni, cdnj}], cdn_sc_i, cdn_sc_im1,
                                        cdn_sc_j, cdn_sc_jp1,
                                        fld.ctx.em.ILOuterMismatch(bp, bim1, bjp1)),
                            {{.table = DpTable::P, .i = i, .j = j, .cdni = cdni, .cdnj = cdnj}},
                            cdnim1, cdnjp1);
                    }
                }
            }
        }
    }

    void TraceImm(const DpIdx& idx, const PartialTrace& pt) {
        int i = idx.i, j = idx.j, cdni = idx.cdni, cdnj = idx.cdnj;

        BasePair bp = BasesToPair(fld.idx2base[{i, cdni}], fld.idx2base[{j, cdnj}]);

        // Since we don't store this DP table, recompute it

        EnergyT base_e = sr.Zero();

        int cdnip1_lb = 0, cdnip1_ub = fld.num_codons[i + 1];
        int cdnjm1_lb = 0, cdnjm1_ub = fld.num_codons[j - 1];
        if (i / 3 == (i + 1) / 3) {
            cdnip1_lb = cdni;
            cdnip1_ub = cdni + 1;
        }
        if (j / 3 == (j - 1) / 3) {
            cdnjm1_lb = cdnj;
            cdnjm1_ub = cdnj + 1;
        }
        for (int cdnip1 = cdnip1_lb; cdnip1 < cdnip1_ub; ++cdnip1) {
            EnergyT cdn_sc_ip1 =
                (i + 1) / 3 == i / 3 ? sr.One() : fld.idx_cdn_score[{i + 1, cdnip1}];
            for (int cdnjm1 = cdnjm1_lb; cdnjm1 < cdnjm1_ub; ++cdnjm1) {
                EnergyT cdn_sc_jm1 =
                    (j - 1) / 3 == j / 3 ? sr.One() : fld.idx_cdn_score[{j - 1, cdnjm1}];
                Base bip1 = fld.idx2base[{i + 1, cdnip1}];
                Base bjm1 = fld.idx2base[{j - 1, cdnjm1}];
                base_e = sr.Add(base_e, sr.MultMany(cdn_sc_ip1, cdn_sc_jm1,
                                                    fld.ctx.em.ILInnerMismatch(bp, bip1, bjm1)));
            }
        }

        auto relax = [&](EnergyT e, std::initializer_list<DpIdx> idxs, int cdnip1, int cdnjm1) {
            auto new_pt = pt;
            AddCodon(new_pt.pri, i + 1, cdnip1);
            AddCodon(new_pt.pri, j - 1, cdnjm1);
            RelaxTrace(e, idxs, base_e, new_pt);
        };

        for (int cdnip1 = cdnip1_lb; cdnip1 < cdnip1_ub; ++cdnip1) {
            EnergyT cdn_sc_ip1 =
                (i + 1) / 3 == i / 3 ? sr.One() : fld.idx_cdn_score[{i + 1, cdnip1}];
            for (int cdnjm1 = cdnjm1_lb; cdnjm1 < cdnjm1_ub; ++cdnjm1) {
                EnergyT cdn_sc_jm1 =
                    (j - 1) / 3 == j / 3 ? sr.One() : fld.idx_cdn_score[{j - 1, cdnjm1}];
                Base bip1 = fld.idx2base[{i + 1, cdnip1}];
                Base bjm1 = fld.idx2base[{j - 1, cdnjm1}];
                base_e = sr.Add(base_e, sr.MultMany(cdn_sc_ip1, cdn_sc_jm1,
                                                    fld.ctx.em.ILInnerMismatch(bp, bip1, bjm1)));
                relax(
                    sr.MultMany(cdn_sc_ip1, cdn_sc_jm1, fld.ctx.em.ILInnerMismatch(bp, bip1, bjm1)),
                    {}, cdnip1, cdnjm1);
            }
        }
    }

    void TraceMl(const DpIdx& idx, PartialTrace& pt) {
        int i = idx.i, j = idx.j, br = idx.br, cdni = idx.cdni, cdnj = idx.cdnj;
        EnergyT base_e = fld.ml_table[{br, j, i, cdni, cdnj}];

        auto relax = [&](EnergyT e, std::initializer_list<DpIdx> idxs) {
            RelaxTrace(e, idxs, base_e, pt);
        };

        // If the codon is tbd, add all choices for the codon
        if (cdni == fld.tbd_cdn) {
            for (int cdn = 0; cdn < fld.num_codons[i]; ++cdn) {
                relax(
                    sr.MultMany(fld.idx_cdn_score[{i, cdn}], fld.ml_table[{br, j, i, cdn, cdnj}]),
                    {{.table = DpTable::ML, .br = br, .i = i, .j = j, .cdni = cdn, .cdnj = cdnj}});
            }
            return;
        }

        // Add the codons only if it's not tbd
        AddCodon(pt.pri, i, cdni);
        // This will redundantly add j multiple times. Could probability optimize this away
        AddCodon(pt.pri, j, cdnj);

        if (i + 1 <= j) {
            if ((i + 1) / 3 == i / 3)
                relax(sr.MultMany(fld.ctx.em.MultiUnpaired(),
                                  fld.ml_table[{br, j, i + 1, cdni, cdnj}]),
                      {{.table = DpTable::ML,
                        .br = br,
                        .i = i + 1,
                        .j = j,
                        .cdni = cdni,
                        .cdnj = cdnj}});
            else if ((i + 1) / 3 == j / 3)
                relax(sr.MultMany(fld.ctx.em.MultiUnpaired(),
                                  fld.ml_table[{br, j, i + 1, cdnj, cdnj}]),
                      {{.table = DpTable::ML,
                        .br = br,
                        .i = i + 1,
                        .j = j,
                        .cdni = cdnj,
                        .cdnj = cdnj}});
            else
                relax(sr.MultMany(fld.ctx.em.MultiUnpaired(),
                                  fld.ml_table[{br, j, i + 1, fld.tbd_cdn, cdnj}]),
                      {{.table = DpTable::ML,
                        .br = br,
                        .i = i + 1,
                        .j = j,
                        .cdni = fld.tbd_cdn,
                        .cdnj = cdnj}});
        } else if (br == 0) {
            relax(fld.ctx.em.MultiUnpaired(), {});
        }
        int nbr = std::max(0, br - 1);

        for (int k = i + fld.ctx.min_hairpin_loop_size + 1; k <= j; ++k) {
            int cdnk_lb = 0, cdnk_ub = fld.num_codons[k];

            // k/3 != i/3 because min_hairpin_loop_size >= 3
            if (k / 3 == j / 3) {
                cdnk_lb = cdnj;
                cdnk_ub = cdnj + 1;
            }
            for (int cdnk = cdnk_lb; cdnk < cdnk_ub; ++cdnk) {
                EnergyT cdn_sc_k = k / 3 == j / 3 ? sr.One() : fld.idx_cdn_score[{k, cdnk}];
                BasePair bp = fld.bp_table[{i, k, cdni, cdnk}];
                if (bp == BasePair::NN) continue;
                int cdnkp1 = fld.tbd_cdn;
                if ((k + 1) / 3 == k / 3)
                    cdnkp1 = cdnk;
                else if ((k + 1) / 3 == j / 3)
                    cdnkp1 = cdnj;
                if (k + 1 <= j) {
                    relax(sr.MultMany(fld.p_table[{i, k, cdni, cdnk}],
                                      fld.ml_table[{nbr, j, k + 1, cdnkp1, cdnj}], cdn_sc_k,
                                      fld.ctx.em.MultiBranch(Base::N, bp, Base::N)),
                          {{.table = DpTable::P, .i = i, .j = k, .cdni = cdni, .cdnj = cdnk},
                           {.table = DpTable::ML,
                            .br = nbr,
                            .i = k + 1,
                            .j = j,
                            .cdni = cdnkp1,
                            .cdnj = cdnj}});
                } else if (nbr == 0) {
                    relax(sr.MultMany(fld.p_table[{i, k, cdni, cdnk}], cdn_sc_k,
                                      fld.ctx.em.MultiBranch(Base::N, bp, Base::N)),
                          {{.table = DpTable::P, .i = i, .j = k, .cdni = cdni, .cdnj = cdnk}});
                }
            }
        }
    }

    void TraceEl(const DpIdx& idx, PartialTrace& pt) {
        int i = idx.i, cdni = idx.cdni;
        EnergyT base_e = fld.el_table[{i, cdni}];

        auto relax = [&](EnergyT e, std::initializer_list<DpIdx> idxs) {
            RelaxTrace(e, idxs, base_e, pt);
        };

        // If the codon is tbd, add all choices for the codon
        if (cdni == fld.tbd_cdn) {
            for (int cdn = 0; cdn < fld.num_codons[i]; ++cdn) {
                relax(sr.MultMany(fld.idx_cdn_score[{i, cdn}], fld.el_table[{i, cdn}]),
                      {{.table = DpTable::EL, .i = i, .cdni = cdn}});
            }
            return;
        }

        // Add the codons only if it's not tbd
        AddCodon(pt.pri, i, cdni);

        if (i + 1 < fld.n) {
            if ((i + 1) / 3 != i / 3)
                relax(fld.el_table[{i + 1, fld.tbd_cdn}],
                      {{.table = DpTable::EL, .i = i + 1, .cdni = fld.tbd_cdn}});
            else
                relax(fld.el_table[{i + 1, cdni}],
                      {{.table = DpTable::EL, .i = i + 1, .cdni = cdni}});
        } else {
            relax(sr.One(), {});
        }

        Base bi = fld.idx2base[{i, cdni}];
        for (int j = i + fld.ctx.min_hairpin_loop_size + 1; j < fld.n; ++j) {
            // Asumption: min hairpin size is >= 3, so cdnj is not dependent on cdni
            for (int cdnj = 0; cdnj < fld.num_codons[j]; ++cdnj) {
                EnergyT cdn_sc_j = fld.idx_cdn_score[{j, cdnj}];
                Base bj = fld.idx2base[{j, cdnj}];
                BasePair bp = BasesToPair(bi, bj);
                if (bp == BasePair::NN) continue;
                if (j + 1 < fld.n) {
                    if ((j + 1) / 3 == j / 3) {
                        relax(sr.MultMany(fld.ctx.em.ExtBranch(Base::N, bp, Base::N),
                                          fld.el_table[{j + 1, cdnj}], cdn_sc_j,
                                          fld.p_table[{i, j, cdni, cdnj}]),
                              {{.table = DpTable::P, .i = i, .j = j, .cdni = cdni, .cdnj = cdnj},
                               {.table = DpTable::EL, .i = j + 1, .cdni = cdnj}});
                    } else {
                        relax(sr.MultMany(fld.ctx.em.ExtBranch(Base::N, bp, Base::N),
                                          fld.el_table[{j + 1, fld.tbd_cdn}], cdn_sc_j,
                                          fld.p_table[{i, j, cdni, cdnj}]),
                              {{.table = DpTable::P, .i = i, .j = j, .cdni = cdni, .cdnj = cdnj},
                               {.table = DpTable::EL, .i = j + 1, .cdni = fld.tbd_cdn}});
                    }
                } else {
                    relax(sr.MultMany(fld.ctx.em.ExtBranch(Base::N, bp, Base::N), cdn_sc_j,
                                      fld.p_table[{i, j, cdni, cdnj}]),
                          {{.table = DpTable::P, .i = i, .j = j, .cdni = cdni, .cdnj = cdnj}});
                }
            }
        }
    }

    void TraceUnpairedCodons(const DpIdx& idx, const PartialTrace& pt) {
        int i = idx.i, j = idx.j;

        // Don't process empty sequences
        // Simply pass back onto the pq
        if (i > j) {
            RelaxTrace(pt.min_energy, {}, pt.min_energy, pt);
            return;
        }

        EnergyT base_e = fld.unpaired_codons[{i, j}];

        auto relax_j = [&](EnergyT e, std::initializer_list<DpIdx> idxs, int cdnj) {
            auto pt_new = pt;
            AddCodon(pt_new.pri, j * 3, cdnj);
            RelaxTrace(e, idxs, base_e, pt_new);
        };

        for (int cdnj = 0; cdnj < fld.num_codons[j * 3]; ++cdnj) {
            relax_j(sr.MultMany(fld.idx_cdn_score[{j * 3, cdnj}], fld.unpaired_codons[{i, j - 1}]),
                    {{.table = DpTable::UpCodons, .i = i, .j = j - 1}}, cdnj);
        }
    }

    void TraceScUpCodon(const DpIdx& idx, PartialTrace& pt) {
        int i = idx.i, cdni = idx.cdni;
        AddCodon(pt.pri, i, cdni);
        RelaxTrace(pt.min_energy, {}, pt.min_energy, pt);
    }

    void TraceScSpecialHairpin(const DpIdx& idx, PartialTrace& pt) {
        int i = idx.i, j = idx.j, cdni = idx.cdni;
        auto special_hairpins = fld.ctx.em.SpecialHairpins(j - i - 1);
        for (int k = i; k <= j; ++k) pt.pri[k] = special_hairpins[cdni].first[k - i];
        RelaxTrace(pt.min_energy, {}, pt.min_energy, pt);
    }
    double subopt_randomness;
    int max_pq_size;
    std::multiset<PartialTrace> set_pq;
    EnergyT trace_threshold;
    Fold<EnergyModelT> fld;
    SemiringT sr;
    std::default_random_engine rng;
};
}  // namespace mwmrna::fold_codons

#endif  // MWMRNA_LIB_TRACE_CODONS_HPP