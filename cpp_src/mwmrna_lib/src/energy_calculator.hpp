#ifndef MWMRNA_LIB_ENERGY_CALCULATOR_HPP
#define MWMRNA_LIB_ENERGY_CALCULATOR_HPP

#include <sstream>
#include <tuple>

#include "energy_model.hpp"
#include "rna.hpp"
#include "structure.hpp"

namespace mwmrna {
template <class EnergyModelT>
class EnergyCalculator {
   private:
    using SemiringT = typename EnergyModelT::semiring_type;
    using EnergyT = typename SemiringT::energy_type;
    using ParamsT = typename EnergyModelT::params_type;

   public:
    EnergyCalculator(SeqEnergyModel<EnergyModelT> sem) : sem_(sem) { Precomp(); }

    EnergyCalculator(EnergyModelT em, Primary rna) : sem_(em, rna) { Precomp(); }

    void SetModel(EnergyModelT em) {
        sem_ = SeqEnergyModel<EnergyModelT>(em);
        Precomp();
    }

    EnergyT Calculate(const PairTree& pt) const {
        EnergyT energy = semiring_.One();
        int rt = pt.Root();
        for (int c = 0; c < pt.Children(rt); ++c) {
            int c_id = pt.Child(rt, c);
            auto [i, j] = pt.Pair(c_id);
            energy = semiring_.Multiply(energy, sem_.ExtBranch(i, j));
            energy = semiring_.Multiply(energy, RecursiveHelper(pt, c_id));
        }
        return energy;
    }

    std::string Describe(const PairTree& pt) const {
        std::stringstream ss;
        EnergyT energy = semiring_.One();
        int rt = pt.Root();
        for (int c = 0; c < pt.Children(rt); ++c) {
            int c_id = pt.Child(rt, c);
            auto [i, j] = pt.Pair(c_id);
            energy = semiring_.Multiply(energy, sem_.ExtBranch(i, j));
            EnergyT rec_energy = RecursiveHelper(pt, c_id);
            ss << "ExtBranch(" << i << ", " << j << "): " << sem_.ExtBranch(i, j) << ", "
               << rec_energy << std::endl;
            energy = semiring_.Multiply(energy, RecursiveHelperDescribe(pt, c_id, ss));
        }
        ss << "Total: " << energy << std::endl;
        return ss.str();
    }

   private:
    EnergyT RecursiveHelper(const PairTree& pt, int at) const {
        auto [i, j] = pt.Pair(at);
        // Declare ahead of time because of switch
        int k, l, c, up, prevl;
        EnergyT energy;
        switch (pt.Children(at)) {
            case 0:
                return sem_.OneLoop(i, j);
            case 1:
                c = pt.Child(at, 0);
                std::tie(k, l) = pt.Pair(c);
                return semiring_.Multiply(sem_.TwoLoop(i, j, k, l), RecursiveHelper(pt, c));
            default:  // multiloops
                energy = sem_.MultiClosing(i, j);
                up = 0;
                prevl = i;
                for (int c = 0; c < pt.Children(at); ++c) {
                    std::tie(k, l) = pt.Pair(pt.Child(at, c));
                    up += k - prevl - 1;
                    energy = semiring_.Multiply(energy, sem_.MultiBranch(k, l));
                    energy = semiring_.Multiply(energy, RecursiveHelper(pt, pt.Child(at, c)));
                    prevl = l;
                }
                up += j - prevl - 1;
                return semiring_.Multiply(energy, unpaired_precomp_[up]);
        }
    }

    EnergyT RecursiveHelperDescribe(const PairTree& pt, int at, std::stringstream& ss,
                                    int depth = 1) const {
        auto [i, j] = pt.Pair(at);
        // Declare ahead of time because of switch
        int k, l, c, up, prevl;
        EnergyT energy;
        EnergyT rec_energy;
        switch (pt.Children(at)) {
            case 0:
                ss << std::string(depth, ' ') << "OneLoop(" << i << ", " << j
                   << "): " << sem_.OneLoop(i, j) << std::endl;
                return sem_.OneLoop(i, j);
            case 1:
                c = pt.Child(at, 0);
                std::tie(k, l) = pt.Pair(c);
                rec_energy = RecursiveHelper(pt, c);
                ss << std::string(depth, ' ') << "TwoLoop(" << i << ", " << j << ", " << k << ", "
                   << l << "): " << sem_.TwoLoop(i, j, k, l) << ", " << rec_energy << std::endl;
                return semiring_.Multiply(sem_.TwoLoop(i, j, k, l),
                                          RecursiveHelperDescribe(pt, c, ss, depth + 1));
            default:  // multiloops
                energy = sem_.MultiClosing(i, j);
                ss << std::string(depth, ' ') << "MultiClosing(" << i << ", " << j
                   << "): " << sem_.MultiClosing(i, j) << std::endl;
                up = 0;
                prevl = i;
                for (int c = 0; c < pt.Children(at); ++c) {
                    std::tie(k, l) = pt.Pair(pt.Child(at, c));
                    up += k - prevl - 1;
                    energy = semiring_.Multiply(energy, sem_.MultiBranch(k, l));
                    rec_energy = RecursiveHelper(pt, pt.Child(at, c));
                    ss << std::string(depth, ' ') << "MultiBranch(" << k << ", " << l
                       << "): " << sem_.MultiBranch(k, l) << ", " << rec_energy << std::endl;
                    energy = semiring_.Multiply(
                        energy, RecursiveHelperDescribe(pt, pt.Child(at, c), ss, depth + 1));
                    prevl = l;
                }
                up += j - prevl - 1;
                ss << std::string(depth, ' ') << "Multiloop unpaired(" << up
                   << "): " << unpaired_precomp_[up] << std::endl;
                ss << std::string(depth, ' ') << "Multiloop(" << i << ", " << j
                   << "): " << semiring_.Multiply(energy, unpaired_precomp_[up]) << std::endl;
                return semiring_.Multiply(energy, unpaired_precomp_[up]);
        }
    }

    void Precomp() {
        unpaired_precomp_.resize(sem_.NumNucleotides());
        unpaired_precomp_[0] = semiring_.One();
        for (int i = 1; i < sem_.NumNucleotides(); ++i) {
            unpaired_precomp_[i] =
                semiring_.Multiply(unpaired_precomp_[i - 1], sem_.MultiUnpaired());
        }
    }
    SeqEnergyModel<EnergyModelT> sem_;
    SemiringT semiring_;
    std::vector<EnergyT> unpaired_precomp_;
};
}  // namespace mwmrna

#endif  // MWMRNA_LIB_ENERGY_CALCULATOR_HPP