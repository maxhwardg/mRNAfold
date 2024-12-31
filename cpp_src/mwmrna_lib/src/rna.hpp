#ifndef MWMRNA_LIB_RNA_HPP
#define MWMRNA_LIB_RNA_HPP

#include <array>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace mwmrna {

enum Base { A, U, C, G, N, NUM_BASES };

enum BasePair { AU, UA, CG, GC, GU, UG, NN, NUM_ORDERED_PAIRS };

typedef std::vector<Base> Primary;

typedef std::vector<std::array<double, Base::N>> ContinuousPrimary;

// The following functions are inline so they can be optimized across translation units

inline BasePair BasesToPair(Base b1, Base b2) {
    static const BasePair pairs[NUM_BASES + 1][NUM_BASES + 1] = {
        {NN, AU, NN, NN, NN, NN},  // A
        {UA, NN, NN, UG, NN, NN},  // U
        {NN, NN, NN, CG, NN, NN},  // C
        {NN, GU, GC, NN, NN, NN},  // G
        {NN, NN, NN, NN, NN, NN},  // N
        {NN, NN, NN, NN, NN, NN},  // NUM_BASES
    };
    return pairs[b1][b2];
}

inline BasePair InvertPair(BasePair op) {
    switch (op) {
        case GC:
            return CG;
        case CG:
            return GC;
        case GU:
            return UG;
        case UG:
            return GU;
        case AU:
            return UA;
        case UA:
            return AU;
        default:
            return NN;
    }
}

inline bool IsValidPair(Base b1, Base b2) { return BasesToPair(b1, b2) != NN; }

inline bool IsNonGC(BasePair op) {
    switch (op) {
        case GC:
        case CG:
            return false;
        default:
            return true;
    }
}

inline bool IsNonGC(Base b1, Base b2) { return IsNonGC(BasesToPair(b1, b2)); }

inline std::pair<Base, Base> ExpandPair(BasePair op) {
    switch (op) {
        case GC:
            return {G, C};
        case CG:
            return {C, G};
        case GU:
            return {G, U};
        case UG:
            return {U, G};
        case AU:
            return {A, U};
        case UA:
            return {U, A};
        default:
            return {N, N};
    }
}

Base CharToBase(char c);

char BaseToChar(Base b);

std::string PrimaryToString(const Primary &seq);

Primary StringToPrimary(const std::string &str);

void EnumeratePrimaries(int N, std::function<void(const Primary &)> callback);

ContinuousPrimary PrimaryToContinuous(const Primary &pri);

}  // namespace mwmrna

#endif  // MWMRNA_LIB_RNA_HPP