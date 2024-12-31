#ifndef MWMRNA_LIB_CODON_FOLD_HPP
#define MWMRNA_LIB_CODON_FOLD_HPP

#include <stddef.h>

#include <algorithm>
#include <execution>
#include <optional>
#include <queue>
#include <random>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "energy_model.hpp"
#include "energy_semiring.hpp"
#include "md_array.hpp"
#include "protein.hpp"
#include "rna.hpp"
#include "structure.hpp"

using std::vector;

namespace mwmrna::fold_codons {

template <typename EnergyModelT>
struct FoldContext {
    using ViennaParamsT = typename EnergyModelT::params_type;
    using SemiringT = typename ViennaParamsT::semiring_type;
    using EnergyT = typename SemiringT::energy_type;
    EnergyModelT em;
    CodonTable codon_table;
    int max_two_loop_size = 30;
    int min_hairpin_loop_size = 3;
    AASeq aa_seq;
    vector<vector<EnergyT>> codon_score;
    bool parallel = false;
};

template <typename EnergyModelT>
class Fold {
   private:
    using EnergyT = typename FoldContext<EnergyModelT>::EnergyT;
    using SemiringT = typename FoldContext<EnergyModelT>::SemiringT;

   public:
    Fold(FoldContext<EnergyModelT> _ctx) : ctx(_ctx) { FillTables(); }
    Fold(FoldContext<EnergyModelT>&& _ctx) : ctx(_ctx) { FillTables(); }

    /// @brief Represents an undecided codon
    int tbd_cdn;
    /// @brief Number of bases in the sequence
    int n;

    /// @brief Converts (nt index, codon index) to a base
    MdArray<Base, 2> idx2base;
    /// @brief Converts an nt index to number of valid codons
    vector<int> num_codons;
    /// @brief Converts (nt index, codon index) to codon score
    MdArray<EnergyT, 2> idx_cdn_score;
    /// @brief Product (under semiring) of all codons between two amino acid indices
    MdArray<EnergyT, 2> unpaired_codons;
    /// @brief Indicates if (nt index i, nt index j, codon i, codon j) can form a base pair
    MdArray<BasePair, 4> bp_table;

    /// @brief DP table. (i,j,cdni,cdnj) -> sum (under semiring) of structues where i and j are
    /// paired and have codons cdni and cdnj
    MdArray<EnergyT, 4> p_table;
    /// @brief DP table for multiloops. (branch, j, i, cdni, cdnj) -> sum (under semiring) of all
    /// multiloops fragments containing at least branch branches between i and j with codons cdni
    /// and cdnj at i and j
    MdArray<EnergyT, 5> ml_table;

    /// @brief DP table for external loops (i, cdni) -> sum (under semiring) of all external loops
    /// for the suffix starting at i with codon cdni at i
    MdArray<EnergyT, 2> el_table;
    /// @brief DP table for outer mismatches. (i,j) -> sum (under semiring) of all outermismathed
    /// calls to p_table[i,j]. Useful for internal loop optimizations
    MdArray<EnergyT, 2> omm_table;
    SemiringT sr;

    /// @brief The context under which folding was computed
    FoldContext<EnergyModelT> ctx;

   protected:
    void FillTables() {
        if (ctx.min_hairpin_loop_size < 3)
            throw std::runtime_error("Hairpin loop size < 3 not supported");

        int max_codons = ctx.codon_table.MaxCodons();
        tbd_cdn = max_codons;

        n = ctx.aa_seq.size() * 3;
        // The ML table is indexed by branches, j, i, cdni, cdnj
        // j and i and flipped to allow for a more cache-friendly access pattern
        ml_table = MdArray<EnergyT, 5>(sr.Zero(), {3, n, n, max_codons + 1, max_codons});
        p_table = MdArray<EnergyT, 4>(sr.Zero(), {n, n, max_codons, max_codons});
        el_table = MdArray<EnergyT, 2>(sr.Zero(), {n, max_codons + 1});
        omm_table = MdArray<EnergyT, 2>(sr.Zero(), {n, n});

        idx2base = MdArray<Base, 2>(Base::N, {n, max_codons});
        num_codons = vector<int>(n);
        for (int i = 0; i < n; ++i) {
            auto codons = ctx.codon_table.GetCodons(ctx.aa_seq[i / 3]);
            num_codons[i] = static_cast<int>(codons.size());
            for (int j = 0; j < num_codons[i]; ++j) {
                idx2base[{i, j}] = CodonToPrimary(codons[j])[i % 3];
            }
        }

        int aa_len = ctx.aa_seq.size();

        // Cache codon scores in a 2D array, since nested vectors are slow
        // Store using nt indexes because it is faster (no /3) and convenient
        idx_cdn_score = MdArray<EnergyT, 2>(sr.One(), {n, max_codons});
        for (int i = 0; i < n; ++i) {
            if (static_cast<int>(ctx.codon_score[i / 3].size()) > max_codons)
                throw std::runtime_error(
                    "Codon score vector contains more than max_codons entries");
            for (int j = 0; j < static_cast<int>(ctx.codon_score[i / 3].size()); ++j) {
                idx_cdn_score[{i, j}] = ctx.codon_score[i / 3][j];
            }
        }

        unpaired_codons = MdArray<EnergyT, 2>(sr.One(), {aa_len, aa_len});
        for (int i = 0; i < aa_len; ++i) {
            EnergyT prod = sr.One();
            for (int j = i; j < aa_len; ++j) {
                EnergyT sum = sr.Zero();
                for (int cdnj = 0; cdnj < num_codons[j * 3]; ++cdnj) {
                    sum = sr.Add(sum, idx_cdn_score[{j * 3, cdnj}]);
                }
                prod = sr.Multiply(prod, sum);
                unpaired_codons[{i, j}] = prod;
            }
        }

        // Precompute base pair table
        // This is an optimization to avoid calling BasesToPair in the bottleneck loops
        bp_table = MdArray<BasePair, 4>(BasePair::NN, {n, n, max_codons, max_codons});
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                for (int cdni = 0; cdni < num_codons[i]; ++cdni) {
                    for (int cdnj = 0; cdnj < num_codons[j]; ++cdnj) {
                        bp_table[{i, j, cdni, cdnj}] =
                            BasesToPair(idx2base[{i, cdni}], idx2base[{j, cdnj}]);
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
                std::for_each(std::execution::par_unseq,
                              indexes.begin() + i + ctx.min_hairpin_loop_size + 1, indexes.end(),
                              [&](int j) { FillOmm(i, j); });
                std::for_each(std::execution::par_unseq, indexes.begin() + i, indexes.end(),
                              [&](int j) { FillMl(i, j); });
            } else {
                for (int j = i + ctx.min_hairpin_loop_size + 1; j < n; ++j) {
                    FillP(i, j);
                }
                for (int j = i + ctx.min_hairpin_loop_size + 1; j < n; ++j) {
                    FillOmm(i, j);
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

        for (int cdni = 0; cdni < num_codons[i]; ++cdni) {
            for (int cdnj = 0; cdnj < num_codons[j]; ++cdnj) {
                BasePair bp = BasesToPair(idx2base[{i, cdni}], idx2base[{j, cdnj}]);

                // Ignore invalid base pairs
                if (bp == BasePair::NN) continue;

                e = sr.Zero();

                // Stacks
                // We need to try all codons for i+1 and j-1 if they are in different codons
                // to i and j. Note that we can ignore cases where i+1=j-1, since these
                // are invalid for min_hairpin_loop_size >= 3 and the paired table lookup
                // will be zero.
                if (i / 3 != (i + 1) / 3) {
                    for (int cdnip1 = 0; cdnip1 < num_codons[i + 1]; ++cdnip1) {
                        if (j / 3 != (j - 1) / 3) {
                            for (int cdnjm1 = 0; cdnjm1 < num_codons[j - 1]; ++cdnjm1) {
                                e = sr.Add(e, sr.MultMany(
                                                  idx_cdn_score[{i + 1, cdnip1}],
                                                  idx_cdn_score[{j - 1, cdnjm1}],
                                                  p_table[{i + 1, j - 1, cdnip1, cdnjm1}],
                                                  ctx.em.Stack(
                                                      bp, BasesToPair(idx2base[{i + 1, cdnip1}],
                                                                      idx2base[{j - 1, cdnjm1}]))));
                            }
                        } else {
                            e = sr.Add(e,
                                       sr.MultMany(
                                           idx_cdn_score[{i + 1, cdnip1}],
                                           p_table[{i + 1, j - 1, cdnip1, cdnj}],
                                           ctx.em.Stack(bp, BasesToPair(idx2base[{i + 1, cdnip1}],
                                                                        idx2base[{j - 1, cdnj}]))));
                        }
                    }
                } else {
                    if (j / 3 != (j - 1) / 3) {
                        for (int cdnjm1 = 0; cdnjm1 < num_codons[j - 1]; ++cdnjm1) {
                            e = sr.Add(
                                e, sr.MultMany(
                                       idx_cdn_score[{j - 1, cdnjm1}],
                                       p_table[{i + 1, j - 1, cdni, cdnjm1}],
                                       ctx.em.Stack(bp, BasesToPair(idx2base[{i + 1, cdni}],
                                                                    idx2base[{j - 1, cdnjm1}]))));
                        }
                    } else {
                        e = sr.Add(
                            e, sr.MultMany(p_table[{i + 1, j - 1, cdni, cdnj}],
                                           ctx.em.Stack(bp, BasesToPair(idx2base[{i + 1, cdni}],
                                                                        idx2base[{j - 1, cdnj}]))));
                    }
                }

                // Bulges
                // 5' bulges
                for (int k = i + 2; (j - 1) - k - 1 >= ctx.min_hairpin_loop_size &&
                                    k - i - 1 <= ctx.max_two_loop_size;
                     ++k) {
                    int cdnjm1_lb = 0, cdnjm1_ub = num_codons[j - 1];
                    int cdnk_lb = 0, cdnk_ub = num_codons[k];
                    if (j / 3 == (j - 1) / 3) {
                        cdnjm1_lb = cdnj;
                        cdnjm1_ub = cdnj + 1;
                    }
                    if (k / 3 == i / 3) {
                        cdnk_lb = cdni;
                        cdnk_ub = cdni + 1;
                    }
                    EnergyT extra_codons = sr.One();
                    if (i / 3 + 1 <= k / 3 - 1)
                        extra_codons = unpaired_codons[{i / 3 + 1, k / 3 - 1}];
                    for (int cdnjm1 = cdnjm1_lb; cdnjm1 < cdnjm1_ub; ++cdnjm1) {
                        EnergyT cdn_sc_jm1 =
                            j / 3 == (j - 1) / 3 ? sr.One() : idx_cdn_score[{j - 1, cdnjm1}];
                        for (int cdnk = cdnk_lb; cdnk < cdnk_ub; ++cdnk) {
                            EnergyT cdn_sc_k = k / 3 == i / 3 ? sr.One() : idx_cdn_score[{k, cdnk}];
                            BasePair bpkl =
                                BasesToPair(idx2base[{k, cdnk}], idx2base[{j - 1, cdnjm1}]);
                            if (bpkl == BasePair::NN) continue;
                            e = sr.Add(e, sr.MultMany(extra_codons, cdn_sc_k, cdn_sc_jm1,
                                                      p_table[{k, j - 1, cdnk, cdnjm1}],
                                                      ctx.em.Bulge(bp, bpkl, k - i - 1)));
                        }
                    }
                }
                // 3' bulges
                for (int k = j - 2; j - k - 1 <= ctx.max_two_loop_size &&
                                    k - (i + 1) - 1 >= ctx.min_hairpin_loop_size;
                     --k) {
                    int cdnip1_lb = 0, cdnip1_ub = num_codons[i + 1];
                    int cdnk_lb = 0, cdnk_ub = num_codons[k];
                    if ((i + 1) / 3 == i / 3) {
                        cdnip1_lb = cdni;
                        cdnip1_ub = cdni + 1;
                    }
                    if (k / 3 == j / 3) {
                        cdnk_lb = cdnj;
                        cdnk_ub = cdnj + 1;
                    }
                    EnergyT extra_codons = sr.One();
                    if (k / 3 + 1 <= j / 3 - 1)
                        extra_codons = unpaired_codons[{k / 3 + 1, j / 3 - 1}];
                    for (int cdnip1 = cdnip1_lb; cdnip1 < cdnip1_ub; ++cdnip1) {
                        EnergyT cdn_sc_ip1 =
                            (i + 1) / 3 == i / 3 ? sr.One() : idx_cdn_score[{i + 1, cdnip1}];
                        for (int cdnk = cdnk_lb; cdnk < cdnk_ub; ++cdnk) {
                            EnergyT cdn_sc_k = k / 3 == j / 3 ? sr.One() : idx_cdn_score[{k, cdnk}];
                            BasePair bpkl =
                                BasesToPair(idx2base[{i + 1, cdnip1}], idx2base[{k, cdnk}]);
                            if (bpkl == BasePair::NN) continue;
                            e = sr.Add(e, sr.MultMany(extra_codons, cdn_sc_k, cdn_sc_ip1,
                                                      p_table[{i + 1, k, cdnip1, cdnk}],
                                                      ctx.em.Bulge(bp, bpkl, j - k - 1)));
                        }
                    }
                }

                // Internal loops

                // 1x1
                int cdnip2_lb = 0, cdnip2_ub = num_codons[i + 2];
                int cdnjm2_lb = 0, cdnjm2_ub = num_codons[j - 2];
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
                        (i + 2) / 3 == i / 3 ? sr.One() : idx_cdn_score[{i + 2, cdnip2}];
                    for (int cdnjm2 = cdnjm2_lb; cdnjm2 < cdnjm2_ub; ++cdnjm2) {
                        EnergyT cdn_sc_jm2 =
                            (j - 2) / 3 == j / 3 ? sr.One() : idx_cdn_score[{j - 2, cdnjm2}];
                        BasePair bpkl =
                            BasesToPair(idx2base[{i + 2, cdnip2}], idx2base[{j - 2, cdnjm2}]);
                        if (bpkl == BasePair::NN) continue;
                        Base bip1 = idx2base[{i + 1, (i + 1) / 3 == i / 3 ? cdni : cdnip2}];
                        Base bjm1 = idx2base[{j - 1, (j - 1) / 3 == j / 3 ? cdnj : cdnjm2}];
                        e = sr.Add(e,
                                   sr.MultMany(p_table[{i + 2, j - 2, cdnip2, cdnjm2}], cdn_sc_ip2,
                                               cdn_sc_jm2, ctx.em.IL11(bp, bpkl, bip1, bjm1)));
                    }
                }

                // 1x2
                for (int cdnip2 = cdnip2_lb; cdnip2 < cdnip2_ub; ++cdnip2) {
                    EnergyT cdn_sc_ip2 =
                        (i + 2) / 3 == i / 3 ? sr.One() : idx_cdn_score[{i + 2, cdnip2}];
                    for (int cdnjm3 = 0; cdnjm3 < num_codons[j - 3]; ++cdnjm3) {
                        EnergyT cdn_sc_jm3 = idx_cdn_score[{j - 3, cdnjm3}];
                        Base bip1 = idx2base[{i + 1, (i + 1) / 3 == i / 3 ? cdni : cdnip2}];
                        Base bjm1 = idx2base[{j - 1, (j - 1) / 3 == j / 3 ? cdnj : cdnjm3}];
                        Base bjm2 = idx2base[{j - 2, (j - 2) / 3 == j / 3 ? cdnj : cdnjm3}];
                        BasePair bpkl =
                            BasesToPair(idx2base[{i + 2, cdnip2}], idx2base[{j - 3, cdnjm3}]);
                        if (bpkl == BasePair::NN) continue;
                        e = sr.Add(
                            e, sr.MultMany(p_table[{i + 2, j - 3, cdnip2, cdnjm3}], cdn_sc_ip2,
                                           cdn_sc_jm3, ctx.em.IL12(bp, bpkl, bip1, bjm1, bjm2)));
                    }
                }

                // 2x1
                for (int cdnip3 = 0; cdnip3 < num_codons[i + 3]; ++cdnip3) {
                    EnergyT cdn_sc_ip3 = idx_cdn_score[{i + 3, cdnip3}];
                    Base bip1 = idx2base[{i + 1, (i + 1) / 3 == i / 3 ? cdni : cdnip3}];
                    Base bip2 = idx2base[{i + 2, (i + 2) / 3 == i / 3 ? cdni : cdnip3}];
                    for (int cdnjm2 = cdnjm2_lb; cdnjm2 < cdnjm2_ub; ++cdnjm2) {
                        EnergyT cdn_sc_jm2 =
                            (j - 2) / 3 == j / 3 ? sr.One() : idx_cdn_score[{j - 2, cdnjm2}];
                        Base bjm1 = idx2base[{j - 1, (j - 1) / 3 == j / 3 ? cdnj : cdnjm2}];
                        BasePair bpkl =
                            BasesToPair(idx2base[{i + 3, cdnip3}], idx2base[{j - 2, cdnjm2}]);
                        if (bpkl == BasePair::NN) continue;
                        e = sr.Add(
                            e, sr.MultMany(p_table[{i + 3, j - 2, cdnip3, cdnjm2}], cdn_sc_ip3,
                                           cdn_sc_jm2, ctx.em.IL21(bp, bpkl, bip1, bjm1, bip2)));
                    }
                }

                // 2x2
                for (int cdnip3 = 0; cdnip3 < num_codons[i + 3]; ++cdnip3) {
                    EnergyT cdn_sc_ip3 = idx_cdn_score[{i + 3, cdnip3}];
                    Base bip1 = idx2base[{i + 1, (i + 1) / 3 == i / 3 ? cdni : cdnip3}];
                    Base bip2 = idx2base[{i + 2, (i + 2) / 3 == i / 3 ? cdni : cdnip3}];

                    for (int cdnjm3 = 0; cdnjm3 < num_codons[j - 3]; ++cdnjm3) {
                        EnergyT cdn_sc_jm3 = idx_cdn_score[{j - 3, cdnjm3}];
                        Base bjm1 = idx2base[{j - 1, (j - 1) / 3 == j / 3 ? cdnj : cdnjm3}];
                        Base bjm2 = idx2base[{j - 2, (j - 2) / 3 == j / 3 ? cdnj : cdnjm3}];
                        BasePair bpkl =
                            BasesToPair(idx2base[{i + 3, cdnip3}], idx2base[{j - 3, cdnjm3}]);
                        if (bpkl == BasePair::NN) continue;
                        e = sr.Add(e, sr.MultMany(p_table[{i + 3, j - 3, cdnip3, cdnjm3}],
                                                  cdn_sc_ip3, cdn_sc_jm3,
                                                  ctx.em.IL22(bp, bpkl, bip1, bjm1, bip2, bjm2)));
                    }
                }

                // 2x3
                for (int cdnip3 = 0; cdnip3 < num_codons[i + 3]; ++cdnip3) {
                    EnergyT cdn_sc_ip3 = idx_cdn_score[{i + 3, cdnip3}];
                    Base bip1 = idx2base[{i + 1, (i + 1) / 3 == i / 3 ? cdni : cdnip3}];
                    Base bip2 = idx2base[{i + 2, (i + 2) / 3 == i / 3 ? cdni : cdnip3}];
                    for (int cdnjm4 = 0; cdnjm4 < num_codons[j - 4]; ++cdnjm4) {
                        EnergyT cdn_sc_jm4 = idx_cdn_score[{j - 4, cdnjm4}];
                        Base bjm1, bjm3;
                        BasePair bpkl =
                            BasesToPair(idx2base[{i + 3, cdnip3}], idx2base[{j - 4, cdnjm4}]);
                        if (bpkl == BasePair::NN) continue;
                        EnergyT p = p_table[{i + 3, j - 4, cdnip3, cdnjm4}];
                        // This can only happen when the 3nt unpaired region is a single codon
                        if ((j - 1) / 3 == (j - 3) / 3) {
                            for (int cdnjm1 = 0; cdnjm1 < num_codons[j - 1]; ++cdnjm1) {
                                EnergyT cdn_sc_jm1 = idx_cdn_score[{j - 1, cdnjm1}];
                                bjm1 = idx2base[{j - 1, cdnjm1}];
                                bjm3 = idx2base[{j - 3, cdnjm1}];
                                e = sr.Add(
                                    e, sr.MultMany(p, cdn_sc_ip3, cdn_sc_jm4, cdn_sc_jm1,
                                                   ctx.em.IL23(bp, bpkl, bip1, bjm1, bip2, bjm3)));
                            }
                        } else {
                            // TODO: There probably doesn't need to be a branch here
                            bjm1 = idx2base[{j - 1, (j - 1) / 3 == j / 3 ? cdnj : cdnjm4}];
                            bjm3 = idx2base[{j - 3, (j - 3) / 3 == j / 3 ? cdnj : cdnjm4}];

                            e = sr.Add(e,
                                       sr.MultMany(p, cdn_sc_ip3, cdn_sc_jm4,
                                                   ctx.em.IL23(bp, bpkl, bip1, bjm1, bip2, bjm3)));
                        }
                    }
                }

                // 3x2
                for (int cdnjm3 = 0; cdnjm3 < num_codons[j - 3]; ++cdnjm3) {
                    EnergyT cdn_sc_jm3 = idx_cdn_score[{j - 3, cdnjm3}];
                    Base bjm1 = idx2base[{j - 1, (j - 1) / 3 == j / 3 ? cdnj : cdnjm3}];
                    Base bjm2 = idx2base[{j - 2, (j - 2) / 3 == j / 3 ? cdnj : cdnjm3}];
                    for (int cdnip4 = 0; cdnip4 < num_codons[i + 4]; ++cdnip4) {
                        EnergyT cdn_sc_ip4 = idx_cdn_score[{i + 4, cdnip4}];
                        Base bip1, bip3;
                        BasePair bpkl =
                            BasesToPair(idx2base[{i + 4, cdnip4}], idx2base[{j - 3, cdnjm3}]);
                        if (bpkl == BasePair::NN) continue;
                        EnergyT p = p_table[{i + 4, j - 3, cdnip4, cdnjm3}];
                        // This can only happen when the 3nt unpaired region is a single codon
                        if ((i + 1) / 3 == (i + 3) / 3) {
                            for (int cdnip1 = 0; cdnip1 < num_codons[i + 1]; ++cdnip1) {
                                EnergyT cdn_sc_ip1 = idx_cdn_score[{i + 1, cdnip1}];
                                bip1 = idx2base[{i + 1, cdnip1}];
                                bip3 = idx2base[{i + 3, cdnip1}];
                                e = sr.Add(
                                    e, sr.MultMany(p, cdn_sc_ip1, cdn_sc_ip4, cdn_sc_jm3,
                                                   ctx.em.IL32(bp, bpkl, bip1, bjm1, bip3, bjm2)));
                            }
                        } else {
                            // TODO: There probably doesn't need to be a branch here
                            bip1 = idx2base[{i + 1, (i + 1) / 3 == i / 3 ? cdni : cdnip4}];
                            bip3 = idx2base[{i + 3, (i + 3) / 3 == i / 3 ? cdni : cdnip4}];
                            e = sr.Add(e,
                                       sr.MultMany(p, cdn_sc_ip4, cdn_sc_jm3,
                                                   ctx.em.IL32(bp, bpkl, bip1, bjm1, bip3, bjm2)));
                        }
                    }
                }

                // 1xn
                for (int k = j - 4;
                     j - k <= ctx.max_two_loop_size && k - (i + 2) - 1 >= ctx.min_hairpin_loop_size;
                     --k) {
                    // Note for j-4: 1x1 and 1x2 are handled above
                    // cdnk cannot be the same as cdnj since we know there are >2 nts
                    int cdnk_lb = 0, cdnk_ub = num_codons[k];
                    int cdnkp1_lb = 0, cdnkp1_ub = num_codons[k + 1];

                    for (int cdnip2 = cdnip2_lb; cdnip2 < cdnip2_ub; ++cdnip2) {
                        EnergyT cdn_sc_ip2 =
                            (i + 2) / 3 == i / 3 ? sr.One() : idx_cdn_score[{i + 2, cdnip2}];
                        int cdnip1 = i / 3 == (i + 1) / 3 ? cdni : cdnip2;
                        for (int cdnk = cdnk_lb; cdnk < cdnk_ub; ++cdnk) {
                            EnergyT cdn_sc_k = idx_cdn_score[{k, cdnk}];
                            BasePair bpkl =
                                BasesToPair(idx2base[{i + 2, cdnip2}], idx2base[{k, cdnk}]);
                            if (bpkl == BasePair::NN) continue;
                            if ((k + 1) / 3 == k / 3) {
                                cdnkp1_lb = cdnk;
                                cdnkp1_ub = cdnk + 1;
                            }
                            for (int cdnkp1 = cdnkp1_lb; cdnkp1 < cdnkp1_ub; ++cdnkp1) {
                                EnergyT cdn_sc_kp1 = (k + 1) / 3 == k / 3
                                                         ? sr.One()
                                                         : idx_cdn_score[{k + 1, cdnkp1}];
                                int cdnjm1_lb = 0, cdnjm1_ub = num_codons[j - 1];
                                if ((j - 1) / 3 == j / 3) {
                                    cdnjm1_lb = cdnj;
                                    cdnjm1_ub = cdnj + 1;
                                } else if ((j - 1) / 3 == (k + 1) / 3) {
                                    cdnjm1_lb = cdnkp1;
                                    cdnjm1_ub = cdnkp1 + 1;
                                }
                                for (int cdnjm1 = cdnjm1_lb; cdnjm1 < cdnjm1_ub; ++cdnjm1) {
                                    EnergyT cdn_sc_jm1 =
                                        (j - 1) / 3 == j / 3 || (j - 1) / 3 == (k + 1) / 3
                                            ? sr.One()
                                            : idx_cdn_score[{j - 1, cdnjm1}];
                                    Base bip1 = idx2base[{i + 1, cdnip1}];
                                    Base bjm1 = idx2base[{j - 1, cdnjm1}];
                                    Base bkp1 = idx2base[{k + 1, cdnkp1}];
                                    EnergyT extra_codons = sr.One();
                                    if ((k + 1) / 3 + 1 <= (j - 1) / 3 - 1)
                                        extra_codons =
                                            unpaired_codons[{(k + 1) / 3 + 1, (j - 1) / 3 - 1}];
                                    e = sr.Add(
                                        e, sr.MultMany(
                                               p_table[{i + 2, k, cdnip2, cdnk}], extra_codons,
                                               cdn_sc_ip2, cdn_sc_k, cdn_sc_kp1, cdn_sc_jm1,
                                               ctx.em.IL1N(bp, bpkl, bip1, bjm1, bkp1, j - k - 1)));
                                }
                            }
                        }
                    }
                }

                // nx1
                for (int k = i + 4;
                     (j - 2) - k - 1 >= ctx.min_hairpin_loop_size && k - i <= ctx.max_two_loop_size;
                     ++k) {
                    int cdnk_lb = 0, cdnk_ub = num_codons[k];
                    int cdnkm1_lb = 0, cdnkm1_ub = num_codons[k - 1];
                    for (int cdnjm2 = cdnjm2_lb; cdnjm2 < cdnjm2_ub; ++cdnjm2) {
                        EnergyT cdn_sc_jm2 =
                            (j - 2) / 3 == j / 3 ? sr.One() : idx_cdn_score[{j - 2, cdnjm2}];
                        int cdnjm1 = j / 3 == (j - 1) / 3 ? cdnj : cdnjm2;
                        for (int cdnk = cdnk_lb; cdnk < cdnk_ub; ++cdnk) {
                            EnergyT cdn_sc_k = idx_cdn_score[{k, cdnk}];
                            BasePair bpkl =
                                BasesToPair(idx2base[{k, cdnk}], idx2base[{j - 2, cdnjm2}]);
                            if (bpkl == BasePair::NN) continue;
                            if ((k - 1) / 3 == k / 3) {
                                cdnkm1_lb = cdnk;
                                cdnkm1_ub = cdnk + 1;
                            }
                            for (int cdnkm1 = cdnkm1_lb; cdnkm1 < cdnkm1_ub; ++cdnkm1) {
                                EnergyT cdn_sc_km1 = (k - 1) / 3 == k / 3
                                                         ? sr.One()
                                                         : idx_cdn_score[{k - 1, cdnkm1}];
                                int cdnip1_lb = 0, cdnip1_ub = num_codons[i + 1];
                                if ((i + 1) / 3 == i / 3) {
                                    cdnip1_lb = cdni;
                                    cdnip1_ub = cdni + 1;
                                } else if ((i + 1) / 3 == (k - 1) / 3) {
                                    cdnip1_lb = cdnkm1;
                                    cdnip1_ub = cdnkm1 + 1;
                                }
                                for (int cdnip1 = cdnip1_lb; cdnip1 < cdnip1_ub; ++cdnip1) {
                                    EnergyT cdn_sc_ip1 =
                                        (i + 1) / 3 == i / 3 || (i + 1) / 3 == (k - 1) / 3
                                            ? sr.One()
                                            : idx_cdn_score[{i + 1, cdnip1}];
                                    Base bip1 = idx2base[{i + 1, cdnip1}];
                                    Base bjm1 = idx2base[{j - 1, cdnjm1}];
                                    Base bkm1 = idx2base[{k - 1, cdnkm1}];
                                    EnergyT extra_codons = sr.One();
                                    if ((i + 1) / 3 + 1 <= (k - 1) / 3 - 1)
                                        extra_codons =
                                            unpaired_codons[{(i + 1) / 3 + 1, (k - 1) / 3 - 1}];
                                    e = sr.Add(
                                        e, sr.MultMany(
                                               p_table[{k, j - 2, cdnk, cdnjm2}], extra_codons,
                                               cdn_sc_k, cdn_sc_km1, cdn_sc_jm2, cdn_sc_ip1,
                                               ctx.em.ILN1(bp, bpkl, bip1, bjm1, bkm1, k - i - 1)));
                                }
                            }
                        }
                    }
                }

                // arbitrary internal loops

                // Precompute sum for inner mismatches
                // Used when inner mismatch codons are independent of outer mismatch codons
                int cdnip1_lb = 0, cdnip1_ub = num_codons[i + 1];
                int cdnjm1_lb = 0, cdnjm1_ub = num_codons[j - 1];
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
                        (i + 1) / 3 == i / 3 ? sr.One() : idx_cdn_score[{i + 1, cdnip1}];
                    for (int cdnjm1 = cdnjm1_lb; cdnjm1 < cdnjm1_ub; ++cdnjm1) {
                        EnergyT cdn_sc_jm1 =
                            (j - 1) / 3 == j / 3 ? sr.One() : idx_cdn_score[{j - 1, cdnjm1}];
                        Base bip1 = idx2base[{i + 1, cdnip1}];
                        Base bjm1 = idx2base[{j - 1, cdnjm1}];
                        e_inner_mismatches =
                            sr.Add(e_inner_mismatches,
                                   sr.MultMany(cdn_sc_ip1, cdn_sc_jm1,
                                               ctx.em.ILInnerMismatch(bp, bip1, bjm1)));
                    }
                }
                for (int k = i + 3; k - i - 1 <= ctx.max_two_loop_size && k < j - 3; ++k) {
                    for (int l = j - 3; j - l - 1 + k - i - 1 <= ctx.max_two_loop_size &&
                                        l - k - 1 >= ctx.min_hairpin_loop_size;
                         --l) {
                        int lup = k - i - 1, rup = j - l - 1;

                        // Ignore 2x2 and 2x3 internal loops as these are handled above
                        if ((lup == 3 && rup == 2) || (lup == 2 && rup == 3)) continue;
                        if (lup == 2 && rup == 2) continue;

                        int cdnip1_lb = 0, cdnip1_ub = num_codons[i + 1];
                        int cdnjm1_lb = 0, cdnjm1_ub = num_codons[j - 1];
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
                            int cdnk_lb = 0, cdnk_ub = num_codons[k];
                            int cdnl_lb = 0, cdnl_ub = num_codons[l];
                            for (int cdnk = cdnk_lb; cdnk < cdnk_ub; ++cdnk) {
                                EnergyT cdn_sc_k = idx_cdn_score[{k, cdnk}];
                                for (int cdnl = cdnl_lb; cdnl < cdnl_ub; ++cdnl) {
                                    EnergyT cdn_sc_l = idx_cdn_score[{l, cdnl}];
                                    BasePair bpkl =
                                        BasesToPair(idx2base[{k, cdnk}], idx2base[{l, cdnl}]);
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
                                        EnergyT cdn_sc_ip1 =
                                            (i + 1) / 3 == i / 3 || (i + 1) / 3 == k / 3
                                                ? sr.One()
                                                : idx_cdn_score[{i + 1, cdnip1}];
                                        for (int cdnjm1 = cdnjm1_lb; cdnjm1 < cdnjm1_ub; ++cdnjm1) {
                                            EnergyT cdn_sc_jm1 =
                                                (j - 1) / 3 == j / 3 || (j - 1) / 3 == l / 3
                                                    ? sr.One()
                                                    : idx_cdn_score[{j - 1, cdnjm1}];
                                            if ((k - 1) / 3 == (i + 1) / 3 &&
                                                (l + 1) / 3 == (j - 1) / 3) {
                                                // All mismatches share codons. Only happens for
                                                // 3x3s
                                                Base bip1 = idx2base[{i + 1, cdnip1}];
                                                Base bjm1 = idx2base[{j - 1, cdnjm1}];
                                                Base bkm1 = idx2base[{k - 1, cdnip1}];
                                                Base blp1 = idx2base[{l + 1, cdnjm1}];
                                                // No need for combos since this is only for 3x3
                                                e = sr.Add(
                                                    e, sr.MultMany(p_table[{k, l, cdnk, cdnl}],
                                                                   cdn_sc_ip1, cdn_sc_jm1, cdn_sc_k,
                                                                   cdn_sc_l,
                                                                   ctx.em.NormalInternalLoop(
                                                                       bp, bpkl, bip1, bjm1, bkm1,
                                                                       blp1, lup, rup)));
                                            } else if ((k - 1) / 3 == (i + 1) / 3) {
                                                // Only 5' mismatch shares codons
                                                int cdnlp1_lb = 0, cdnlp1_ub = num_codons[l + 1];
                                                if ((l + 1) / 3 == l / 3) {
                                                    cdnlp1_lb = cdnl;
                                                    cdnlp1_ub = cdnl + 1;
                                                }
                                                for (int cdnlp1 = cdnlp1_lb; cdnlp1 < cdnlp1_ub;
                                                     ++cdnlp1) {
                                                    EnergyT cdn_sc_lp1 =
                                                        (l + 1) / 3 == l / 3
                                                            ? sr.One()
                                                            : idx_cdn_score[{l + 1, cdnlp1}];
                                                    Base bip1 = idx2base[{i + 1, cdnip1}];
                                                    Base bjm1 = idx2base[{j - 1, cdnjm1}];
                                                    Base bkm1 = idx2base[{k - 1, cdnip1}];
                                                    Base blp1 = idx2base[{l + 1, cdnlp1}];
                                                    EnergyT extra_codons = sr.One();
                                                    // Only need to consider extra_codons on the
                                                    // 3' side
                                                    if ((l + 1) / 3 + 1 <= (j - 1) / 3 - 1)
                                                        extra_codons = unpaired_codons[{
                                                            (l + 1) / 3 + 1, (j - 1) / 3 - 1}];
                                                    e = sr.Add(
                                                        e, sr.MultMany(p_table[{k, l, cdnk, cdnl}],
                                                                       extra_codons, cdn_sc_ip1,
                                                                       cdn_sc_jm1, cdn_sc_k,
                                                                       cdn_sc_l, cdn_sc_lp1,
                                                                       ctx.em.NormalInternalLoop(
                                                                           bp, bpkl, bip1, bjm1,
                                                                           bkm1, blp1, lup, rup)));
                                                }
                                            } else if ((l + 1) / 3 == (j - 1) / 3) {
                                                // Only 3' mismatch shares codons
                                                int cdnkm1_lb = 0, cdnkm1_ub = num_codons[k - 1];
                                                if ((k - 1) / 3 == k / 3) {
                                                    cdnkm1_lb = cdnk;
                                                    cdnkm1_ub = cdnk + 1;
                                                }
                                                for (int cdnkm1 = cdnkm1_lb; cdnkm1 < cdnkm1_ub;
                                                     ++cdnkm1) {
                                                    EnergyT cdn_sc_km1 =
                                                        (k - 1) / 3 == k / 3
                                                            ? sr.One()
                                                            : idx_cdn_score[{k - 1, cdnkm1}];
                                                    Base bip1 = idx2base[{i + 1, cdnip1}];
                                                    Base bjm1 = idx2base[{j - 1, cdnjm1}];
                                                    Base bkm1 = idx2base[{k - 1, cdnkm1}];
                                                    Base blp1 = idx2base[{l + 1, cdnjm1}];
                                                    EnergyT extra_codons = sr.One();
                                                    // Only need to consider extra_codons on the
                                                    // 5' side
                                                    if ((i + 1) / 3 + 1 <= (k - 1) / 3 - 1)
                                                        extra_codons = unpaired_codons[{
                                                            (i + 1) / 3 + 1, (k - 1) / 3 - 1}];
                                                    e = sr.Add(
                                                        e, sr.MultMany(p_table[{k, l, cdnk, cdnl}],
                                                                       extra_codons, cdn_sc_ip1,
                                                                       cdn_sc_jm1, cdn_sc_k,
                                                                       cdn_sc_l, cdn_sc_km1,
                                                                       ctx.em.NormalInternalLoop(
                                                                           bp, bpkl, bip1, bjm1,
                                                                           bkm1, blp1, lup, rup)));
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
                                extra_codons = unpaired_codons[{(i + 1) / 3 + 1, (k - 1) / 3 - 1}];
                            if ((l + 1) / 3 + 1 <= (j - 1) / 3 - 1)
                                extra_codons = sr.MultMany(
                                    extra_codons,
                                    unpaired_codons[{(l + 1) / 3 + 1, (j - 1) / 3 - 1}]);
                            e = sr.Add(e, sr.MultMany(omm_table[{k, l}], extra_codons,
                                                      e_inner_mismatches, ctx.em.ILInit(lup + rup),
                                                      ctx.em.ILAsymmetry(lup, rup)));
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
                auto special_hairpins = ctx.em.SpecialHairpins(j - i - 1);
                Primary pri;
                if (i / 3 + 1 == j / 3) {
                    // 0 internal codons
                    pri.assign(j - i + 1, Base::N);
                    for (int k = i; k < i / 3 * 3 + 3; ++k) pri[k - i] = idx2base[{k, cdni}];
                    for (int k = i / 3 * 3 + 3; k <= j; ++k) pri[k - i] = idx2base[{k, cdnj}];
                    auto it = std::find_if(special_hairpins.begin(), special_hairpins.end(),
                                           [&pri](const auto& p) { return p.first == pri; });
                    if (it != special_hairpins.end()) {
                        e = sr.Add(e, it->second);
                    } else {
                        e = sr.Add(e, ctx.em.HairpinNormal(bp, idx2base[{i + 1, cdni}],
                                                           idx2base[{j - 1, cdnj}], j - i - 1));
                    }
                } else if (i / 3 + 2 == j / 3) {
                    // 1 internal codon
                    pri.assign(j - i + 1, Base::N);
                    for (int k = i; k < (i / 3) * 3 + 3; ++k) pri[k - i] = idx2base[{k, cdni}];
                    for (int k = j; k >= (j / 3) * 3; --k)
                        pri[pri.size() - 1 - (j - k)] = idx2base[{k, cdnj}];
                    int cdnip1, cdnjm1;
                    for (int cdn = 0; cdn < num_codons[i / 3 * 3 + 3]; ++cdn) {
                        EnergyT cdn_sc = idx_cdn_score[{i / 3 * 3 + 3, cdn}];
                        cdnip1 = cdni;
                        cdnjm1 = cdnj;
                        if ((i + 1) / 3 != i / 3) cdnip1 = cdn;
                        if ((j - 1) / 3 != j / 3) cdnjm1 = cdn;
                        for (int k = i / 3 * 3 + 3; k < i / 3 * 3 + 6; ++k)
                            pri[k - i] = idx2base[{k, cdn}];
                        auto it = std::find_if(special_hairpins.begin(), special_hairpins.end(),
                                               [&pri](const auto& p) { return p.first == pri; });
                        if (it != special_hairpins.end()) {
                            e = sr.Add(e, sr.MultMany(cdn_sc, it->second));
                        } else {
                            e = sr.Add(
                                e, sr.MultMany(cdn_sc, ctx.em.HairpinNormal(
                                                           bp, idx2base[{i + 1, cdnip1}],
                                                           idx2base[{j - 1, cdnjm1}], j - i - 1)));
                        }
                    }
                } else if (i / 3 + 3 == j / 3) {
                    // 2 internal codons
                    pri.assign(j - i + 1, Base::N);
                    for (int k = i; k < (i / 3) * 3 + 3; ++k) pri[k - i] = idx2base[{k, cdni}];
                    for (int k = j; k >= (j / 3) * 3; --k)
                        pri[pri.size() - 1 - (j - k)] = idx2base[{k, cdnj}];
                    int cdnip1, cdnjm1;
                    for (int cdnk = 0; cdnk < num_codons[i / 3 * 3 + 3]; ++cdnk) {
                        EnergyT cdn_sc_k = idx_cdn_score[{i / 3 * 3 + 3, cdnk}];
                        cdnip1 = cdni;
                        if ((i + 1) / 3 != i / 3) cdnip1 = cdnk;
                        for (int k = (i / 3) * 3 + 3; k < (i / 3) * 3 + 6; ++k)
                            pri[k - i] = idx2base[{k, cdnk}];
                        for (int cdnl = 0; cdnl < num_codons[i / 3 * 3 + 6]; ++cdnl) {
                            EnergyT cdn_sc_l = idx_cdn_score[{i / 3 * 3 + 6, cdnl}];
                            cdnjm1 = cdnj;
                            if ((j - 1) / 3 != j / 3) cdnjm1 = cdnl;
                            for (int l = (i / 3) * 3 + 6; l < (i / 3) * 3 + 9; ++l)
                                pri[l - i] = idx2base[{l, cdnl}];
                            auto it =
                                std::find_if(special_hairpins.begin(), special_hairpins.end(),
                                             [&pri](const auto& p) { return p.first == pri; });
                            if (it != special_hairpins.end()) {
                                e = sr.Add(e, sr.MultMany(cdn_sc_k, cdn_sc_l, it->second));
                            } else {
                                e = sr.Add(e,
                                           sr.MultMany(cdn_sc_k, cdn_sc_l,
                                                       ctx.em.HairpinNormal(
                                                           bp, idx2base[{i + 1, cdnip1}],
                                                           idx2base[{j - 1, cdnjm1}], j - i - 1)));
                            }
                        }
                    }
                } else {
                    // 3+ internal codons
                    // Assumption: There are no special hairpins with 3+ internal codons
                    EnergyT internal_combos = unpaired_codons[{(i + 1) / 3 + 1, (j - 1) / 3 - 1}];

                    if (i / 3 != (i + 1) / 3) {
                        for (int cdnip1 = 0; cdnip1 < num_codons[i + 1]; ++cdnip1) {
                            EnergyT cdn_sc_ip1 = idx_cdn_score[{i + 1, cdnip1}];
                            if (j / 3 != (j - 1) / 3) {
                                // j-1 is in a different codon to j and i+1

                                for (int cdnjm1 = 0; cdnjm1 < num_codons[j - 1]; ++cdnjm1) {
                                    EnergyT cdn_sc_jm1 = idx_cdn_score[{j - 1, cdnjm1}];
                                    e = sr.Add(
                                        e, sr.MultMany(internal_combos, cdn_sc_ip1, cdn_sc_jm1,
                                                       ctx.em.HairpinNormal(
                                                           bp, idx2base[{i + 1, cdnip1}],
                                                           idx2base[{j - 1, cdnjm1}], j - i - 1)));
                                }
                            } else if (j / 3 == (j - 1) / 3) {
                                // j-1 is in the same codon as j
                                e = sr.Add(e, sr.MultMany(internal_combos, cdn_sc_ip1,
                                                          ctx.em.HairpinNormal(
                                                              bp, idx2base[{i + 1, cdnip1}],
                                                              idx2base[{j - 1, cdnj}], j - i - 1)));
                            }
                        }
                    } else {
                        if (j / 3 == (j - 1) / 3) {
                            // j-1 is in the same codon as j
                            e = sr.Add(e, sr.MultMany(internal_combos,
                                                      ctx.em.HairpinNormal(
                                                          bp, idx2base[{i + 1, cdni}],
                                                          idx2base[{j - 1, cdnj}], j - i - 1)));
                        } else {
                            // j-1 is in a different codon to j and i+1
                            for (int cdnjm1 = 0; cdnjm1 < num_codons[j - 1]; ++cdnjm1) {
                                EnergyT cdn_sc_jm1 = idx_cdn_score[{j - 1, cdnjm1}];
                                e = sr.Add(
                                    e, sr.MultMany(internal_combos, cdn_sc_jm1,
                                                   ctx.em.HairpinNormal(bp, idx2base[{i + 1, cdni}],
                                                                        idx2base[{j - 1, cdnjm1}],
                                                                        j - i - 1)));
                            }
                        }
                    }
                }

                // Multi loops
                cdnip1_lb = 0, cdnip1_ub = num_codons[i + 1];
                cdnjm1_lb = 0, cdnjm1_ub = num_codons[j - 1];
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
                        (i + 1) / 3 == i / 3 ? sr.One() : idx_cdn_score[{i + 1, cdnip1}];
                    // We don't bother dealing with the case where i+1 and j-1 are in the same
                    // codon since the ML table will be invalid
                    for (int cdnjm1 = cdnjm1_lb; cdnjm1 < cdnjm1_ub; ++cdnjm1) {
                        EnergyT cdn_sc_jm1 =
                            (j - 1) / 3 == j / 3 ? sr.One() : idx_cdn_score[{j - 1, cdnjm1}];
                        e = sr.Add(
                            e, sr.MultMany(ml_table[{2, j - 1, i + 1, cdnip1, cdnjm1}], cdn_sc_ip1,
                                           cdn_sc_jm1, ctx.em.MultiClosing(bp, Base::N, Base::N)));
                    }
                }

                p_table[{i, j, cdni, cdnj}] = e;
            }
        }
    }

    void FillOmm(int i, int j) {
        // Avoid edge cases that cannot be outer mismatches anyway
        if (i == 0 || j == n - 1) return;
        EnergyT e = sr.Zero();
        for (int cdni = 0; cdni < num_codons[i]; ++cdni) {
            EnergyT cdn_sc_i = idx_cdn_score[{i, cdni}];
            for (int cdnj = 0; cdnj < num_codons[j]; ++cdnj) {
                EnergyT cdn_sc_j = idx_cdn_score[{j, cdnj}];
                BasePair bp = BasesToPair(idx2base[{i, cdni}], idx2base[{j, cdnj}]);
                if (bp == BasePair::NN) continue;
                int cdnim1_lb = 0, cdnim1_ub = num_codons[i - 1];
                int cdnjp1_lb = 0, cdnjp1_ub = num_codons[j + 1];
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
                        (i - 1) / 3 == i / 3 ? sr.One() : idx_cdn_score[{i - 1, cdnim1}];
                    for (int cdnjp1 = cdnjp1_lb; cdnjp1 < cdnjp1_ub; ++cdnjp1) {
                        EnergyT cdn_sc_jp1 =
                            (j + 1) / 3 == j / 3 ? sr.One() : idx_cdn_score[{j + 1, cdnjp1}];
                        Base bim1 = idx2base[{i - 1, cdnim1}];
                        Base bjp1 = idx2base[{j + 1, cdnjp1}];
                        e = sr.Add(e, sr.MultMany(p_table[{i, j, cdni, cdnj}], cdn_sc_i, cdn_sc_im1,
                                                  cdn_sc_j, cdn_sc_jp1,
                                                  ctx.em.ILOuterMismatch(bp, bim1, bjp1)));
                    }
                }
            }
        }
        omm_table[{i, j}] = e;
    }

    void FillMl(int i, int j) {
        EnergyT e;
        for (int br = 0; br <= 2; ++br) {
            for (int cdni = 0; cdni < num_codons[i]; ++cdni) {
                int cdnj_lb = 0, cdnj_ub = num_codons[j];
                // Only fill the table for valid cdni/cdnj pairs
                if (i / 3 == j / 3) {
                    cdnj_lb = cdni;
                    cdnj_ub = cdni + 1;
                }
                for (int cdnj = cdnj_lb; cdnj < cdnj_ub; ++cdnj) {
                    e = sr.Zero();
                    if (i + 1 <= j) {
                        if ((i + 1) / 3 == i / 3)
                            e = sr.MultMany(ctx.em.MultiUnpaired(),
                                            ml_table[{br, j, i + 1, cdni, cdnj}]);
                        else if ((i + 1) / 3 == j / 3)
                            e = sr.MultMany(ctx.em.MultiUnpaired(),
                                            ml_table[{br, j, i + 1, cdnj, cdnj}]);
                        else
                            e = sr.MultMany(ctx.em.MultiUnpaired(),
                                            ml_table[{br, j, i + 1, tbd_cdn, cdnj}]);
                    } else if (br == 0) {
                        e = ctx.em.MultiUnpaired();
                    }

                    int nbr = std::max(0, br - 1);
                    for (int k = i + ctx.min_hairpin_loop_size + 1; k <= j; ++k) {
                        int cdnk_lb = 0, cdnk_ub = num_codons[k];

                        // k/3 != i/3 because min_hairpin_loop_size >= 3
                        if (k / 3 == j / 3) {
                            cdnk_lb = cdnj;
                            cdnk_ub = cdnj + 1;
                        }
                        for (int cdnk = cdnk_lb; cdnk < cdnk_ub; ++cdnk) {
                            EnergyT cdn_sc_k = k / 3 == j / 3 ? sr.One() : idx_cdn_score[{k, cdnk}];
                            BasePair bp = bp_table[{i, k, cdni, cdnk}];
                            if (bp == BasePair::NN) continue;
                            int cdnkp1 = tbd_cdn;
                            if ((k + 1) / 3 == k / 3)
                                cdnkp1 = cdnk;
                            else if ((k + 1) / 3 == j / 3)
                                cdnkp1 = cdnj;
                            EnergyT nml;
                            if (nbr == 0) {
                                nml = sr.One();
                            } else {
                                nml = sr.Zero();
                            }
                            if (k + 1 <= j) nml = ml_table[{nbr, j, k + 1, cdnkp1, cdnj}];
                            e = sr.Add(e, sr.MultMany(p_table[{i, k, cdni, cdnk}], nml, cdn_sc_k,
                                                      ctx.em.MultiBranch(Base::N, bp, Base::N)));
                        }
                    }

                    ml_table[{br, j, i, cdni, cdnj}] = e;
                    ml_table[{br, j, i, tbd_cdn, cdnj}] =
                        sr.Add(ml_table[{br, j, i, tbd_cdn, cdnj}],
                               sr.MultMany(e, idx_cdn_score[{i, cdni}]));
                }
            }
        }
    }

    void FillEl(int i) {
        EnergyT e;
        for (int cdni = 0; cdni < num_codons[i]; ++cdni) {
            e = sr.One();
            if (i + 1 < n) {
                if ((i + 1) / 3 != i / 3)
                    e = el_table[{i + 1, tbd_cdn}];
                else
                    e = el_table[{i + 1, cdni}];
            }
            Base bi = idx2base[{i, cdni}];
            for (int j = i + ctx.min_hairpin_loop_size + 1; j < n; ++j) {
                // Asumption: min hairpin size is >= 3, so cdnj is not dependent on cdni
                for (int cdnj = 0; cdnj < num_codons[j]; ++cdnj) {
                    EnergyT cdn_sc_j = idx_cdn_score[{j, cdnj}];
                    Base bj = idx2base[{j, cdnj}];
                    BasePair bp = BasesToPair(bi, bj);
                    if (bp == BasePair::NN) continue;
                    EnergyT e_jp1 = sr.One();
                    if (j + 1 < n) {
                        if ((j + 1) / 3 == j / 3)
                            e_jp1 = el_table[{j + 1, cdnj}];
                        else
                            e_jp1 = el_table[{j + 1, tbd_cdn}];
                    }
                    e = sr.Add(e, sr.MultMany(ctx.em.ExtBranch(Base::N, bp, Base::N), e_jp1,
                                              cdn_sc_j, p_table[{i, j, cdni, cdnj}]));
                }
            }

            el_table[{i, cdni}] = e;
        }
        e = sr.Zero();
        for (int cdn = 0; cdn < num_codons[i]; ++cdn)
            e = sr.Add(e, sr.MultMany(idx_cdn_score[{i, cdn}], el_table[{i, cdn}]));
        el_table[{i, tbd_cdn}] = e;
    }
};

}  // namespace mwmrna::fold_codons
#endif  // MWMRNA_LIB_CODON_FOLD_HPP