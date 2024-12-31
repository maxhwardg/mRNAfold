#ifndef MWMRNA_LIB_ENERGY_MODEL_HPP
#define MWMRNA_LIB_ENERGY_MODEL_HPP

#include <cstdlib>
#include <functional>
#include <random>
#include <type_traits>

#include "energy_semiring.hpp"
#include "rna.hpp"
#include "vienna_params.hpp"

namespace mwmrna {

/// @brief Nearest neighbor energy model
/// @tparam ParamsT ViennaRNA parameter type
/// @tparam Dangles Dangle model to use
/// @note Dangles d=0 and d=2 are supported
/// @note Dangles have the same meaning as in ViennaRNA
template <class ParamsT, int Dangles = 2>
class EnergyModel {
   private:
    using SemiringT = typename ParamsT::semiring_type;
    using EnergyT = typename SemiringT::energy_type;
    static_assert(EnergySemiring<SemiringT, EnergyT>, "SemiringT not an EnergySemiring");
    static_assert(Dangles == 0 || Dangles == 2, "Dangles must be 0 or 2");

   public:
    using params_type = ParamsT;
    using dangles_constant = std::integral_constant<int, Dangles>;
    using semiring_type = SemiringT;

    EnergyModel() = delete;
    EnergyModel(const EnergyModel&) = default;
    EnergyModel(EnergyModel&&) = default;
    EnergyModel& operator=(const EnergyModel&) = default;
    EnergyModel(ParamsT params) : params_(params) { Precomp(); }

    EnergyT HairpinNormal(BasePair bij, Base bip1, Base bjm1, int unpaired) const {
        EnergyT e;
        if (unpaired < static_cast<int>(params_.Data().hairpin_init.size()))
            e = params_.Data().hairpin_init[unpaired];
        else
            e = semiring_.Multiply(params_.Data().hairpin_init.back(),
                                   params_.LoopExtrapolation(unpaired));
        if (unpaired == 3) {
            if (IsNonGC(bij)) e = semiring_.Multiply(e, params_.Data().non_gc_closing_penalty);
        } else {
            e = semiring_.Multiply(e, params_.Data().hairpin_mismatch[bij][bip1][bjm1]);
        }
        return e;
    }

    EnergyT HairpinSpecial(int id, int len) const {
        switch (len) {
            case 4:
                return params_.Data().tetraloops[id].second;
            case 3:
                return params_.Data().triloops[id].second;
            case 6:
                return params_.Data().hexaloops[id].second;
            default:
                return ParamsT::kINF;
        }
    }

    EnergyT Stack(BasePair bij, BasePair bkl) const {
        // The bkl pair needs to be swapped due to vienna's convention
        return params_.Data().stack[bij][InvertPair(bkl)];
    }

    EnergyT ExtBranch(Base bim1, BasePair bij, Base bjp1) const {
        EnergyT e = semiring_.One();
        if constexpr (Dangles == 2) {
            // Special case for dangling ends where one of the bases is N
            // This is for branches at the ends of the sequence
            if (bim1 == Base::N) {
                if (bjp1 == Base::N) {
                    e = semiring_.One();
                } else {
                    e = params_.Data().dangle3[bij][bjp1];
                }
            } else if (bjp1 == Base::N) {
                e = params_.Data().dangle5[bij][bim1];
            } else {
                e = params_.Data().mismatch_exterior[bij][bim1][bjp1];
            }
        }
        if (IsNonGC(bij)) e = semiring_.Multiply(e, params_.Data().non_gc_closing_penalty);
        return e;
    }

    EnergyT MultiBranch(Base bim1, BasePair bij, Base bjp1) const {
        EnergyT e = params_.Data().ml_branch;
        // Always use mismatch since this can never go out of bounds
        if constexpr (Dangles == 2) {
            e = semiring_.Multiply(e, params_.Data().mismatch_multi[bij][bim1][bjp1]);
        }
        if (IsNonGC(bij)) e = semiring_.Multiply(e, params_.Data().non_gc_closing_penalty);
        return e;
    }

    EnergyT MultiUnpaired() const { return params_.Data().ml_unpaired; }

    EnergyT MultiClosing(BasePair bij, Base bip1, Base bjm1) const {
        // Note that 5' and 3' are swapped here
        return semiring_.Multiply(params_.Data().ml_init, MultiBranch(bjm1, InvertPair(bij), bip1));
    }

    EnergyT Bulge(BasePair bij, BasePair bkl, int unpaired) const {
        EnergyT e;
        if (unpaired < static_cast<int>(params_.Data().bulge_init.size()))
            e = params_.Data().bulge_init[unpaired];
        else
            e = semiring_.Multiply(params_.Data().bulge_init.back(),
                                   params_.LoopExtrapolation(unpaired));
        if (unpaired == 1) {
            return semiring_.Multiply(e, Stack(bij, bkl));
        }
        if (IsNonGC(bij)) e = semiring_.Multiply(e, params_.Data().non_gc_closing_penalty);
        if (IsNonGC(bkl)) e = semiring_.Multiply(e, params_.Data().non_gc_closing_penalty);
        return e;
    }

    const std::vector<std::pair<Primary, EnergyT>>& SpecialHairpins(int len) const {
        // This is needed since a reference must be returned
        static std::vector<std::pair<Primary, EnergyT>> empty_vec_;
        switch (len) {
            case 3:
                return params_.Data().triloops;
            case 4:
                return params_.Data().tetraloops;
            case 6:
                return params_.Data().hexaloops;
            default:
                return empty_vec_;
        }
    }

    // TODO: Test or use in IlNormal
    EnergyT ILInnerMismatch(BasePair bij, Base bip1, Base bjm1) const {
        return params_.Data().mismatch_interior[bij][bip1][bjm1];
    }

    // TODO: Test or use in IlNormal
    EnergyT ILOuterMismatch(BasePair bij, Base bim1, Base bjp1) const {
        // The bij pair needs to be swapped due to vienna's convention
        return params_.Data().mismatch_interior[InvertPair(bij)][bjp1][bim1];
    }

    EnergyT ILInit(int sz) const {
        if (sz < static_cast<int>(params_.Data().interior_init.size()))
            return params_.Data().interior_init[sz];
        else
            return semiring_.Multiply(params_.Data().interior_init.back(),
                                      params_.LoopExtrapolation(sz));
    }

    EnergyT ILAsymmetry(int asym) const {
        if (asym < static_cast<int>(asymmetry_precomp_.size())) return asymmetry_precomp_[asym];

        EnergyT e = asymmetry_precomp_.back();
        // Check if we can early exit due to threshold
        if (semiring_.LessThan(e, params_.Data().il_asym_thrshold)) {
            return params_.Data().il_asym_thrshold;
        }
        // Otherwise, use the slow power operation
        e = semiring_.Power(params_.Data().il_asym, asym);
        if (semiring_.LessThan(e, params_.Data().il_asym_thrshold)) {
            return params_.Data().il_asym_thrshold;
        }
        return e;
    }

    EnergyT ILAsymmetry(int lup, int rup) const { return ILAsymmetry(std::abs(lup - rup)); }

    virtual EnergyT IL11(BasePair bij, BasePair bkl, Base bip1, Base bjm1) const {
        // The bkl pair needs to be swapped due to vienna's convention
        return params_.Data().int11[bij][InvertPair(bkl)][bip1][bjm1];
    }

    virtual EnergyT IL21(BasePair bij, BasePair bkl, Base bip1, Base bjm1, Base bkm1) const {
        // The bkl pair needs to be swapped due to vienna's convention
        return params_.Data().int21[InvertPair(bkl)][bij][bjm1][bip1][bkm1];
    }

    virtual EnergyT IL12(BasePair bij, BasePair bkl, Base bip1, Base bjm1, Base blp1) const {
        // The bkl pair needs to be swapped due to vienna's convention
        return params_.Data().int21[bij][InvertPair(bkl)][bip1][blp1][bjm1];
    }

    virtual EnergyT IL22(BasePair bij, BasePair bkl, Base bip1, Base bjm1, Base bkm1,
                         Base blp1) const {
        // The bkl pair needs to be swapped due to vienna's convention
        return params_.Data().int22[bij][InvertPair(bkl)][bip1][bkm1][blp1][bjm1];
    }

    virtual EnergyT IL23(BasePair bij, BasePair bkl, Base bip1, Base bjm1, Base bkm1,
                         Base blp1) const {
        // The bkl pair needs to be swapped due to vienna's convention
        EnergyT e =
            semiring_.Multiply(params_.Data().mismatch_interior_23[bij][bip1][bjm1],
                               params_.Data().mismatch_interior_23[InvertPair(bkl)][blp1][bkm1]);
        e = semiring_.Multiply(e, ILInit(2 + 3));
        e = semiring_.Multiply(e, ILAsymmetry(2, 3));
        return e;
    }

    virtual EnergyT IL32(BasePair bij, BasePair bkl, Base bip1, Base bjm1, Base bkm1,
                         Base blp1) const {
        return IL23(bij, bkl, bip1, bjm1, bkm1, blp1);
    }

    virtual EnergyT IL1N(BasePair bij, BasePair bkl, Base bip1, Base bjm1, Base blp1, int n) const {
        // The bkl pair needs to be swapped due to vienna's convention
        EnergyT e =
            semiring_.Multiply(params_.Data().mismatch_interior_1n[bij][bip1][bjm1],
                               params_.Data().mismatch_interior_1n[InvertPair(bkl)][blp1][bip1]);
        e = semiring_.Multiply(e, ILInit(1 + n));
        e = semiring_.Multiply(e, ILAsymmetry(1, n));
        return e;
    }

    virtual EnergyT ILN1(BasePair bij, BasePair bkl, Base bip1, Base bjm1, Base bkm1, int n) const {
        // The bkl pair needs to be swapped due to vienna's convention
        EnergyT e =
            semiring_.Multiply(params_.Data().mismatch_interior_1n[bij][bip1][bjm1],
                               params_.Data().mismatch_interior_1n[InvertPair(bkl)][bjm1][bkm1]);
        e = semiring_.Multiply(e, ILInit(1 + n));
        e = semiring_.Multiply(e, ILAsymmetry(1, n));
        return e;
    }

    /// @brief Computes the energy of an internal loop assuming it is not a special case
    /// @note Will *not* infer the type of loop
    virtual EnergyT NormalInternalLoop(BasePair bij, BasePair bkl, Base bip1, Base bjm1, Base bkm1,
                                       Base blp1, int lup, int rup) const {
        EnergyT e =
            semiring_.Multiply(ILInnerMismatch(bij, bip1, bjm1), ILOuterMismatch(bkl, bkm1, blp1));
        e = semiring_.Multiply(e, ILInit(lup + rup));
        e = semiring_.Multiply(e, ILAsymmetry(lup, rup));
        return e;
    }

    /// @brief Computes the energy of an internal loop
    /// @note Will infer the type of loop
    EnergyT InternalLoop(BasePair bij, BasePair bkl, Base bip1, Base bjm1, Base bkm1, Base blp1,
                         int lup, int rup) const {
        switch (lup) {
            case 1:
                switch (rup) {
                    case 1:
                        return IL11(bij, bkl, bip1, bjm1);

                    case 2:
                        return IL12(bij, bkl, bip1, bjm1, blp1);

                    default:
                        return IL1N(bij, bkl, bip1, bjm1, blp1, rup);
                }
            case 2:
                switch (rup) {
                    case 1:
                        return IL21(bij, bkl, bip1, bjm1, bkm1);

                    case 2:
                        return IL22(bij, bkl, bip1, bjm1, bkm1, blp1);

                    case 3:
                        return IL23(bij, bkl, bip1, bjm1, bkm1, blp1);

                    default:
                        return NormalInternalLoop(bij, bkl, bip1, bjm1, bkm1, blp1, lup, rup);
                }

            case 3:
                switch (rup) {
                    case 1:
                        return ILN1(bij, bkl, bip1, bjm1, bkm1, 3);

                    case 2:
                        return IL32(bij, bkl, bip1, bjm1, bkm1, blp1);

                    default:
                        return NormalInternalLoop(bij, bkl, bip1, bjm1, bkm1, blp1, lup, rup);
                }

            default:
                switch (rup) {
                    case 1:
                        return ILN1(bij, bkl, bip1, bjm1, bkm1, lup);

                    default:
                        return NormalInternalLoop(bij, bkl, bip1, bjm1, bkm1, blp1, lup, rup);
                }
                return NormalInternalLoop(bij, bkl, bip1, bjm1, bkm1, blp1, lup, rup);
        }
    }

   private:
    /// @brief Maximum precomp for asymmetry penalty
    /// Aymmetry uses the Power operation for the Boltzmann semiring
    /// Since this is slow, precomp is used to speed it up
    static constexpr int kMaxAsymPrecomp = 100;

    void Precomp() {
        asymmetry_precomp_.resize(kMaxAsymPrecomp + 1);
        asymmetry_precomp_[0] = semiring_.One();
        for (int i = 1; i <= kMaxAsymPrecomp; ++i) {
            asymmetry_precomp_[i] =
                semiring_.Multiply(asymmetry_precomp_[i - 1], params_.Data().il_asym);
            // Clamp to threshold
            if (semiring_.LessThan(asymmetry_precomp_[i], params_.Data().il_asym_thrshold)) {
                asymmetry_precomp_[i] = params_.Data().il_asym_thrshold;
            }
        }
    }

    ParamsT params_;
    SemiringT semiring_;
    std::vector<EnergyT> asymmetry_precomp_;
};

enum class Motif {
    Hairpin = 0,
    Bulge,
    Internal,
    Stack,
    Multi,
    ExtBranch,
    SpecialHairpin,
    Il11,
    Il12,
    Il22,
    Il23,
    Il1N,
    IlNormal,
    NUM_MOTIFS
};

template <class ParamsT, int Dangles>
class RandomEnergyModel : public EnergyModel<ParamsT, Dangles> {
   private:
    using SemiringT = typename ParamsT::semiring_type;
    using EnergyT = typename SemiringT::energy_type;
    static_assert(EnergySemiring<SemiringT, EnergyT>, "SemiringT not an EnergySemiring");
    static_assert(Dangles == 0 || Dangles == 2, "Dangles must be 0 or 2");

   public:
    using params_type = ParamsT;
    using dangles_constant = std::integral_constant<int, Dangles>;
    using semiring_type = SemiringT;

    RandomEnergyModel(ParamsT params, unsigned seed) : EnergyModel<ParamsT, Dangles>(params) {
        std::default_random_engine generator(seed);
        generator.seed(seed);
        // Warm up the generator
        for (int i = 0; i < 1007; ++i) generator();
        base_ = std::uniform_int_distribution<int>(0, MOD - 1)(generator);
    }

    void SetAllowSparsity(bool allow) { allow_sparsity_ = allow; }

    EnergyT HairpinNormal(BasePair bij, Base bip1, Base bjm1, int unpaired) const {
        return RandomEnergy({static_cast<int>(Motif::Hairpin), static_cast<int>(bij),
                             static_cast<int>(bip1), static_cast<int>(bjm1), unpaired});
    }

    EnergyT HairpinSpecial(int id, int len) const {
        return RandomEnergy({static_cast<int>(Motif::SpecialHairpin), id, len});
    }

    EnergyT Stack(BasePair bij, BasePair bkl) const {
        return RandomEnergy(
            {static_cast<int>(Motif::Stack), static_cast<int>(bij), static_cast<int>(bkl)});
    }

    EnergyT ExtBranch(Base bim1, BasePair bij, Base bjp1) const {
        if constexpr (Dangles == 0) {
            return RandomEnergy({static_cast<int>(Motif::ExtBranch), static_cast<int>(bij)});
        }
        return RandomEnergy({static_cast<int>(Motif::ExtBranch), static_cast<int>(bim1),
                             static_cast<int>(bij), static_cast<int>(bjp1)});
    }

    EnergyT MultiBranch(Base bim1, BasePair bij, Base bjp1) const {
        if constexpr (Dangles == 0) {
            return RandomEnergy({
                static_cast<int>(Motif::Multi),
                1,
                static_cast<int>(bij),
            });
        }
        return RandomEnergy({static_cast<int>(Motif::Multi), 1, static_cast<int>(bim1),
                             static_cast<int>(bij), static_cast<int>(bjp1)});
    }

    EnergyT MultiUnpaired() const {
        if (allow_sparsity_) return sr_.One();
        return RandomEnergy({static_cast<int>(Motif::Multi), 2});
    }

    EnergyT MultiClosing(BasePair bij, Base bip1, Base bjm1) const {
        if constexpr (Dangles == 0) {
            return RandomEnergy({static_cast<int>(Motif::Multi), 3, static_cast<int>(bij)});
        }
        return RandomEnergy({static_cast<int>(Motif::Multi), 3, static_cast<int>(bij),
                             static_cast<int>(bip1), static_cast<int>(bjm1)});
    }

    EnergyT Bulge(BasePair bij, BasePair bkl, int unpaired) const {
        return RandomEnergy({static_cast<int>(Motif::Bulge), static_cast<int>(bij),
                             static_cast<int>(bkl), unpaired});
    }

    EnergyT IL11(BasePair bij, BasePair bkl, Base bip1, Base bjm1) const override {
        return RandomEnergy({static_cast<int>(Motif::Il11), static_cast<int>(bij),
                             static_cast<int>(bkl), static_cast<int>(bip1),
                             static_cast<int>(bjm1)});
    }

    EnergyT IL21(BasePair bij, BasePair bkl, Base bip1, Base bjm1, Base bkm1) const override {
        return RandomEnergy({static_cast<int>(Motif::Il12), 1, static_cast<int>(bij),
                             static_cast<int>(bkl), static_cast<int>(bip1), static_cast<int>(bjm1),
                             static_cast<int>(bkm1)});
    }

    EnergyT IL12(BasePair bij, BasePair bkl, Base bip1, Base bjm1, Base blp1) const override {
        return RandomEnergy({static_cast<int>(Motif::Il12), 2, static_cast<int>(bij),
                             static_cast<int>(bkl), static_cast<int>(bip1), static_cast<int>(bjm1),
                             static_cast<int>(blp1)});
    }

    EnergyT IL22(BasePair bij, BasePair bkl, Base bip1, Base bjm1, Base bkm1,
                 Base blp1) const override {
        return RandomEnergy({static_cast<int>(Motif::Il22), static_cast<int>(bij),
                             static_cast<int>(bkl), static_cast<int>(bip1), static_cast<int>(bjm1),
                             static_cast<int>(bkm1), static_cast<int>(blp1)});
    }

    EnergyT IL23(BasePair bij, BasePair bkl, Base bip1, Base bjm1, Base bkm1,
                 Base blp1) const override {
        return RandomEnergy({static_cast<int>(Motif::Il23), 1, static_cast<int>(bij),
                             static_cast<int>(bkl), static_cast<int>(bip1), static_cast<int>(bjm1),
                             static_cast<int>(bkm1), static_cast<int>(blp1)});
    }

    EnergyT IL32(BasePair bij, BasePair bkl, Base bip1, Base bjm1, Base bkm1,
                 Base blp1) const override {
        return RandomEnergy({static_cast<int>(Motif::Il23), 2, static_cast<int>(bij),
                             static_cast<int>(bkl), static_cast<int>(bip1), static_cast<int>(bjm1),
                             static_cast<int>(bkm1), static_cast<int>(blp1)});
    }

    EnergyT IL1N(BasePair bij, BasePair bkl, Base bip1, Base bjm1, Base blp1,
                 int n) const override {
        return RandomEnergy({static_cast<int>(Motif::Il1N), 1, static_cast<int>(bij),
                             static_cast<int>(bkl), static_cast<int>(bip1), static_cast<int>(bjm1),
                             static_cast<int>(blp1), n});
    }

    EnergyT ILN1(BasePair bij, BasePair bkl, Base bip1, Base bjm1, Base bkm1,
                 int n) const override {
        return RandomEnergy({static_cast<int>(Motif::Il1N), 2, static_cast<int>(bij),
                             static_cast<int>(bkl), static_cast<int>(bip1), static_cast<int>(bjm1),
                             static_cast<int>(bkm1), n});
    }

    EnergyT ILInnerMismatch(BasePair bij, Base bip1, Base bjm1) const {
        return RandomEnergy({static_cast<int>(Motif::IlNormal), 1, static_cast<int>(bij),
                             static_cast<int>(bip1), static_cast<int>(bjm1)});
    }

    EnergyT ILOuterMismatch(BasePair bij, Base bim1, Base bjp1) const {
        return RandomEnergy({static_cast<int>(Motif::IlNormal), 2, static_cast<int>(bij),
                             static_cast<int>(bim1), static_cast<int>(bjp1)});
    }

    EnergyT ILInit(int sz) const {
        return RandomEnergy({static_cast<int>(Motif::IlNormal), 3, sz});
    }

    EnergyT ILAsymmetry(int asym) const {
        return RandomEnergy({static_cast<int>(Motif::IlNormal), 4, asym});
    }

    EnergyT ILAsymmetry(int lup, int rup) const { return this->ILAsymmetry(std::abs(lup - rup)); }

    EnergyT NormalInternalLoop(BasePair bij, BasePair bkl, Base bip1, Base bjm1, Base bkm1,
                               Base blp1, int lup, int rup) const override {
        return sr_.MultMany(ILInnerMismatch(bij, bip1, bjm1), ILOuterMismatch(bkl, bkm1, blp1),
                            ILInit(lup + rup), ILAsymmetry(lup, rup));
    }

   private:
    const int MOD = 6131, MUL = 5881;
    int base_;
    SemiringT sr_;
    // If true, enforces the triangle inequality by making MultiUnpaired always sr_.One() (which is
    // 0 for MFE folding)
    bool allow_sparsity_ = false;

    EnergyT RandomEnergy(std::initializer_list<int> args) const {
        int num = base_;
        for (auto arg : args) {
            num = (num * MUL + arg) % MOD;
        }
        double d = static_cast<double>(num) / MOD;
        // Boltzmann needs to be >= 0
        if constexpr (std::is_base_of<BoltzSemiring<EnergyT>, SemiringT>::value)
            return static_cast<EnergyT>(d * 2);
        return static_cast<EnergyT>(d * 2 - 1);
    }
};

/// @brief A version of EnergyModel that allows for injecting changes such as banned motifs
/// Usefull for testing
template <class ParamsT, int Dangles>
class InjectableEnergyModel : public EnergyModel<ParamsT, Dangles> {
   private:
    using SemiringT = typename ParamsT::semiring_type;
    using EnergyT = typename SemiringT::energy_type;
    using SuperT = EnergyModel<ParamsT, Dangles>;
    static_assert(EnergySemiring<SemiringT, EnergyT>, "SemiringT not an EnergySemiring");
    static_assert(Dangles == 0 || Dangles == 2, "Dangles must be 0 or 2");

   public:
    using params_type = ParamsT;
    using dangles_constant = std::integral_constant<int, Dangles>;
    using semiring_type = SemiringT;

    InjectableEnergyModel(ParamsT params) : SuperT(params) {
        allowed_loops_.assign(static_cast<int>(Motif::NUM_MOTIFS), true);
    }

    InjectableEnergyModel(ParamsT params, std::vector<Motif> banned_motifs) : SuperT(params) {
        allowed_loops_.assign(static_cast<int>(Motif::NUM_MOTIFS), true);
        for (auto motif : banned_motifs) {
            DisallowMotif(motif);
        }
    }

    void DisallowMotif(Motif type) {
        if (type == Motif::NUM_MOTIFS) return;
        allowed_loops_[static_cast<int>(type)] = false;
    }

    EnergyT HairpinNormal(BasePair bij, Base bip1, Base bjm1, int unpaired) const {
        if (!IsAllowed(Motif::Hairpin)) return semiring_.Zero();
        return SuperT::HairpinNormal(bij, bip1, bjm1, unpaired);
    }

    EnergyT HairpinSpecial(int id, int len) const {
        if (!IsAllowed(Motif::Hairpin)) return semiring_.Zero();
        if (!IsAllowed(Motif::SpecialHairpin)) return semiring_.Zero();
        return SuperT::HairpinSpecial(id, len);
    }

    EnergyT Stack(BasePair bij, BasePair bkl) const {
        if (!IsAllowed(Motif::Stack)) return semiring_.Zero();
        return SuperT::Stack(bij, bkl);
    }

    EnergyT ExtBranch(Base bim1, BasePair bij, Base bjp1) const {
        if (!IsAllowed(Motif::ExtBranch)) return semiring_.Zero();
        return SuperT::ExtBranch(bim1, bij, bjp1);
    }

    EnergyT MultiBranch(Base bim1, BasePair bij, Base bjp1) const {
        if (!IsAllowed(Motif::Multi)) return semiring_.Zero();
        return SuperT::MultiBranch(bim1, bij, bjp1);
    }

    EnergyT MultiUnpaired() const {
        if (!IsAllowed(Motif::Multi)) return semiring_.Zero();
        return SuperT::MultiUnpaired();
    }

    EnergyT MultiClosing(BasePair bij, Base bip1, Base bjm1) const {
        if (!IsAllowed(Motif::Multi)) return semiring_.Zero();
        return SuperT::MultiClosing(bij, bip1, bjm1);
    }

    EnergyT Bulge(BasePair bij, BasePair bkl, int unpaired) const {
        if (!IsAllowed(Motif::Bulge)) return semiring_.Zero();
        return SuperT::Bulge(bij, bkl, unpaired);
    }

    EnergyT InternalLoop(BasePair bij, BasePair bkl, Base bip1, Base bjm1, Base bkm1, Base blp1,
                         int lup, int rup) const {
        if (!IsAllowed(Motif::Internal)) return semiring_.Zero();
        switch (lup) {
            case 1:
                switch (rup) {
                    case 1:
                        if (!IsAllowed(Motif::Il11)) return semiring_.Zero();
                        return this->IL11(bij, bkl, bip1, bjm1);

                    case 2:
                        if (!IsAllowed(Motif::Il12)) return semiring_.Zero();
                        return this->IL12(bij, bkl, bip1, bjm1, blp1);

                    default:
                        if (!IsAllowed(Motif::Il1N)) return semiring_.Zero();
                        return this->IL1N(bij, bkl, bip1, bjm1, blp1, rup);
                }
            case 2:
                switch (rup) {
                    case 1:
                        if (!IsAllowed(Motif::Il12)) return semiring_.Zero();
                        return this->IL21(bij, bkl, bip1, bjm1, bkm1);

                    case 2:
                        if (!IsAllowed(Motif::Il22)) return semiring_.Zero();
                        return this->IL22(bij, bkl, bip1, bjm1, bkm1, blp1);

                    case 3:
                        if (!IsAllowed(Motif::Il23)) return semiring_.Zero();
                        return this->IL23(bij, bkl, bip1, bjm1, bkm1, blp1);

                    default:
                        if (!IsAllowed(Motif::IlNormal)) return semiring_.Zero();
                        return this->NormalInternalLoop(bij, bkl, bip1, bjm1, bkm1, blp1, lup, rup);
                }

            case 3:
                switch (rup) {
                    case 1:
                        if (!IsAllowed(Motif::Il1N)) return semiring_.Zero();
                        return this->ILN1(bij, bkl, bip1, bjm1, bkm1, 3);

                    case 2:
                        if (!IsAllowed(Motif::Il23)) return semiring_.Zero();
                        return this->IL32(bij, bkl, bip1, bjm1, bkm1, blp1);

                    default:
                        if (!IsAllowed(Motif::IlNormal)) return semiring_.Zero();
                        return this->NormalInternalLoop(bij, bkl, bip1, bjm1, bkm1, blp1, lup, rup);
                }

            default:
                switch (rup) {
                    case 1:
                        if (!IsAllowed(Motif::Il1N)) return semiring_.Zero();
                        return this->ILN1(bij, bkl, bip1, bjm1, bkm1, lup);

                    default:
                        if (!IsAllowed(Motif::IlNormal)) return semiring_.Zero();
                        return this->NormalInternalLoop(bij, bkl, bip1, bjm1, bkm1, blp1, lup, rup);
                }
                if (!IsAllowed(Motif::IlNormal)) return semiring_.Zero();
                return this->NormalInternalLoop(bij, bkl, bip1, bjm1, bkm1, blp1, lup, rup);
        }
    }

   private:
    bool IsAllowed(Motif type) const { return allowed_loops_[static_cast<int>(type)]; }
    std::vector<bool> allowed_loops_;
    SemiringT semiring_;
};

template <class EnergyModelT>
class SeqEnergyModel {
   private:
    using SemiringT = typename EnergyModelT::semiring_type;
    using EnergyT = typename SemiringT::energy_type;
    using ParamsT = typename EnergyModelT::params_type;

   public:
    using semiring_type = SemiringT;
    using params_type = ParamsT;

    SeqEnergyModel(EnergyModelT em, Primary pri) : em_(em), pri_(pri) {}

    EnergyT OneLoop(int i, int j) const {
        int len = j - i - 1;
        Primary pri_sub(pri_.begin() + i, pri_.begin() + j + 1);
        auto special_loops = em_.SpecialHairpins(len);
        // TODO: Consider optimizing this
        // Since hairpins are not usually in the inner loop, it is a low priority
        auto it = std::find_if(special_loops.begin(), special_loops.end(),
                               [&pri_sub](const auto& p) { return p.first == pri_sub; });
        if (it == special_loops.end()) {
            return em_.HairpinNormal(BasesToPair(pri_[i], pri_[j]), pri_[i + 1], pri_[j - 1], len);
        } else {
            return it->second;
        }
    }

    EnergyT TwoLoop(int i, int j, int k, int l) const {
        if (k - i - 1 == 0 && j - l - 1 == 0)
            return em_.Stack(BasesToPair(pri_[i], pri_[j]), BasesToPair(pri_[k], pri_[l]));
        if (k - i - 1 == 0 || j - l - 1 == 0)
            return em_.Bulge(BasesToPair(pri_[i], pri_[j]), BasesToPair(pri_[k], pri_[l]),
                             k - i - 1 + j - l - 1);
        return em_.InternalLoop(BasesToPair(pri_[i], pri_[j]), BasesToPair(pri_[k], pri_[l]),
                                pri_[i + 1], pri_[j - 1], pri_[k - 1], pri_[l + 1], k - i - 1,
                                j - l - 1);
    }

    EnergyT MultiUnpaired() const { return em_.MultiUnpaired(); }

    EnergyT MultiClosing(int i, int j) const {
        return em_.MultiClosing(BasesToPair(pri_[i], pri_[j]), pri_[i + 1], pri_[j - 1]);
    }

    EnergyT MultiBranch(int i, int j) const {
        // Assumption: i-1 and j+1 are always valid since a Multiloop must have a closing pair
        return em_.MultiBranch(pri_[i - 1], BasesToPair(pri_[i], pri_[j]), pri_[j + 1]);
    }

    EnergyT ExtBranch(int i, int j) const {
        return em_.ExtBranch(i == 0 ? N : pri_[i - 1], BasesToPair(pri_[i], pri_[j]),
                             j == static_cast<int>(pri_.size()) - 1 ? N : pri_[j + 1]);
    }

    int NumNucleotides() const { return pri_.size(); }

   protected:
    EnergyModelT em_;
    Primary pri_;
    SemiringT sr_;
};

}  // namespace mwmrna

#endif  // MWMRNA_LIB_ENERGY_MODEL_HPP