#ifndef MWMRNA_LIB_FOLD_HPP
#define MWMRNA_LIB_FOLD_HPP

#include <algorithm>
#include <array>
#include <cstdarg>
#include <execution>
#include <vector>

#include "energy_model.hpp"
#include "energy_semiring.hpp"
#include "rna.hpp"

namespace mwmrna {

using std::array;
using std::vector;

enum class UpLoopType { Hairpin = 0, Bulge, Internal, Multi, External, NumLoopTypes };

template <typename EnergyModelT>
struct FoldConfig {
   private:
    using ParamsT = typename EnergyModelT::params_type;
    using SemiringT = typename ParamsT::semiring_type;
    using EnergyT = typename SemiringT::energy_type;

   public:
    EnergyModelT em;
    int max_two_loop_size = 30;
    int min_hairpin_loop_size = 3;
    EnergyT unpaired_factor = SemiringT().One();
    EnergyT nt_factor = SemiringT().One();
    EnergyT mismatch_up_factor = SemiringT().One();
    EnergyT hairpinloop_up_factor = SemiringT().One();
    EnergyT bulgeloop_up_factor = SemiringT().One();
    EnergyT internalloop_up_factor = SemiringT().One();
    EnergyT multiloop_up_factor = SemiringT().One();
    EnergyT externalloop_up_factor = SemiringT().One();

    bool parallel = false;
};

template <typename EnergyModelT>
class Fold {
   private:
    using ParamsT = typename EnergyModelT::params_type;
    using SemiringT = typename ParamsT::semiring_type;
    using EnergyT = typename SemiringT::energy_type;

   public:
    struct Result {
        EnergyT energy;
        vector<vector<EnergyT>> paired;
    };

    Fold() = delete;
    Fold(const Fold&) = default;
    Fold(Fold&&) = default;
    Fold& operator=(const Fold&) = default;

    Fold(const FoldConfig<EnergyModelT>& config)
        : em_(config.em),
          max_two_loop_size_(config.max_two_loop_size),
          min_hairpin_size_(config.min_hairpin_loop_size),
          unpaired_factor_(config.unpaired_factor),
          nt_factor_(config.nt_factor),
          mismatch_factor_(config.mismatch_up_factor),
          hairpinloop_up_factor_(config.hairpinloop_up_factor),
          bulgeloop_up_factor_(config.bulgeloop_up_factor),
          internalloop_up_factor_(config.internalloop_up_factor),
          multiloop_up_factor_(config.multiloop_up_factor),
          externalloop_up_factor_(config.externalloop_up_factor),
          parallel_(config.parallel) {}

    Result Inside(const Primary& pri) {
        pri_ = pri;
        const int N = pri_.size();
        SeqEnergyModel em(em_, pri_);

        // Allocate tables
        P = vector<vector<EnergyT>>(N, vector<EnergyT>(N, semiring_.Zero()));
        E = vector<EnergyT>(N + 1, semiring_.One());
        // N+1 x 3 x N might suit the iteration order better
        ML = vector<vector<vector<EnergyT>>>(
            3, vector<vector<EnergyT>>(N + 1, vector<EnergyT>(N, semiring_.Zero())));
        ML[0] = vector<vector<EnergyT>>(N + 1, vector<EnergyT>(N, semiring_.One()));

        // Precompute scaling factors
        unpaired_precomp = nt_factor_percomp = hairpin_precomp = bulge_precomp = internal_precomp =
            vector<EnergyT>(N + 1, semiring_.One());
        for (int i = 1; i <= N; ++i) {
            unpaired_precomp[i] = semiring_.MultMany(unpaired_precomp[i - 1], unpaired_factor_);
            nt_factor_percomp[i] = semiring_.MultMany(nt_factor_percomp[i - 1], nt_factor_);
            hairpin_precomp[i] = semiring_.MultMany(hairpin_precomp[i - 1], hairpinloop_up_factor_);
            bulge_precomp[i] = semiring_.MultMany(bulge_precomp[i - 1], bulgeloop_up_factor_);
            internal_precomp[i] =
                semiring_.MultMany(internal_precomp[i - 1], internalloop_up_factor_);
        }

        vector<int> indexes(N);
        std::iota(indexes.begin(), indexes.end(), 0);

        for (int i = N - 1; i >= 0; --i) {
            // External loop unapired case
            E[i] =
                semiring_.MultMany(E[i + 1], unpaired_factor_, nt_factor_, externalloop_up_factor_);

            if (parallel_) {
                // Parallelize the inner loops
                std::for_each(std::execution::par_unseq, begin(indexes) + i + min_hairpin_size_ + 1,
                              end(indexes), [&](int j) { ProcessP(em, i, j); });
                for (int b = 0; b < 3; ++b) {
                    // i>0 and j<N-1 because a multiloops are assumed to have a closing pair in the
                    // energy model
                    if (i > 0) {
                        std::for_each(std::execution::par_unseq, begin(indexes) + i,
                                      end(indexes) - 1, [&](int j) { ProcessMl(em, b, i, j); });
                    }
                }
            } else {
                for (int j = i + min_hairpin_size_ + 1; j < N; ++j) ProcessP(em, i, j);
                // Multiloop table
                // Must come after paired table because it depends on it
                for (int b = 0; b < 3; ++b) {
                    // i>0 and j<N-1 because a multiloops are assumed to have a closing pair in the
                    // energy model
                    if (i > 0) {
                        for (int j = i; j < N - 1; ++j) ProcessMl(em, b, i, j);
                    }
                }
            }
            // External loop branch
            // Must come after paired table because it depends on it
            for (int j = i + min_hairpin_size_ + 1; j < N; ++j) {
                E[i] =
                    semiring_.Add(E[i], semiring_.MultMany(P[i][j], em.ExtBranch(i, j), E[j + 1]));
            }
        }

        // Construct result
        Result res{E[0], P};
        // Deallocate tables
        P = {};
        ML = {};
        E = {};
        unpaired_precomp = {};
        nt_factor_percomp = {};

        return res;
    }

   private:
    inline void ProcessP(SeqEnergyModel<EnergyModelT>& em, int i, int j) {
        if (!IsValidPair(pri_[i], pri_[j])) return;
        // Paired table hairpins
        P[i][j] = semiring_.MultMany(unpaired_precomp[j - i - 1], nt_factor_percomp[j - i + 1],
                                     em.OneLoop(i, j), mismatch_factor_, mismatch_factor_,
                                     hairpin_precomp[j - i - 1]);
        // Paired table two loops (internal, bulges, and stacks)
        for (int k = 0; k <= max_two_loop_size_ && i + k + 1 < j - 1; ++k) {
            for (int l = 0;
                 l + k <= max_two_loop_size_ && (j - l - 1) - (i + k) >= min_hairpin_size_; ++l) {
                EnergyT two_loop_mismatch_factor = semiring_.One();
                for (int m = 0; m < std::min(k, 2) + std::min(l, 2); ++m) {
                    two_loop_mismatch_factor =
                        semiring_.Multiply(two_loop_mismatch_factor, mismatch_factor_);
                }

                EnergyT bulge_or_internal_loop_up_factor;
                if (k == 0 || l == 0) {
                    bulge_or_internal_loop_up_factor = bulge_precomp[k + l];
                } else {
                    bulge_or_internal_loop_up_factor = internal_precomp[k + l];
                }

                P[i][j] = semiring_.Add(
                    P[i][j], semiring_.MultMany(
                                 P[i + k + 1][j - l - 1], em.TwoLoop(i, j, i + k + 1, j - l - 1),
                                 unpaired_precomp[k + l], nt_factor_percomp[k + l + 2],
                                 two_loop_mismatch_factor, bulge_or_internal_loop_up_factor));
            }
        }
        // Paired table multi loops
        P[i][j] = semiring_.Add(
            P[i][j],
            semiring_.MultMany(ML[2][i + 1][j - 1], em.MultiClosing(i, j), nt_factor_percomp[2]));
    }

    inline void ProcessMl(SeqEnergyModel<EnergyModelT>& em, int b, int i, int j) {
        // Multiloop table unpaired case
        ML[b][i][j] = semiring_.MultMany(unpaired_factor_, nt_factor_, em.MultiUnpaired(),
                                         ML[b][i + 1][j], multiloop_up_factor_);
        // Multiloop tables branch cases
        // These depend on the paired table so must come after
        for (int k = i + min_hairpin_size_ + 1; k <= j; ++k) {
            ML[b][i][j] = semiring_.Add(
                ML[b][i][j], semiring_.MultMany(P[i][k], ML[std::max(b - 1, 0)][k + 1][j],
                                                em.MultiBranch(i, k)));
        }
    }

    SemiringT semiring_;
    EnergyModelT em_;
    int max_two_loop_size_;
    int min_hairpin_size_;
    EnergyT unpaired_factor_, nt_factor_, mismatch_factor_, hairpinloop_up_factor_,
        bulgeloop_up_factor_, internalloop_up_factor_, multiloop_up_factor_,
        externalloop_up_factor_;
    bool parallel_;
    vector<vector<EnergyT>> P;
    vector<EnergyT> E;
    vector<vector<vector<EnergyT>>> ML;
    vector<EnergyT> unpaired_precomp, nt_factor_percomp, hairpin_precomp, bulge_precomp,
        internal_precomp;
    Primary pri_;
};

}  // namespace mwmrna

#endif  // MWMRNA_LIB_FOLD_HPP