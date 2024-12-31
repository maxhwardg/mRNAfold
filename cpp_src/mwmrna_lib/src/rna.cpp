#include "rna.hpp"

namespace mwmrna {

Base CharToBase(char c) {
    switch (c) {
        case 'G':
            return G;
        case 'C':
            return C;
        case 'U':
            return U;
        case 'A':
            return A;
        default:
            return N;
    }
}

char BaseToChar(Base b) {
    switch (b) {
        case G:
            return 'G';
        case C:
            return 'C';
        case U:
            return 'U';
        case A:
            return 'A';
        default:
            return 'N';
    }
}

std::string PrimaryToString(const Primary& seq) {
    std::string str;
    str.reserve(seq.size());
    for (auto b : seq) str += BaseToChar(b);
    return str;
}

Primary StringToPrimary(const std::string& str) {
    Primary seq(str.size());
    for (size_t i = 0; i < str.size(); ++i) seq[i] = CharToBase(str[i]);
    return seq;
}

void EnumerateHelper(int N, int i, Primary& seq, std::function<void(const Primary&)> callback) {
    if (i == N) {
        callback(seq);
        return;
    }
    for (int b = 0; b < Base::N; ++b) {
        seq[i] = static_cast<Base>(b);
        EnumerateHelper(N, i + 1, seq, callback);
    }
}

void EnumeratePrimaries(int N, std::function<void(const Primary&)> callback) {
    Primary seq(N);
    EnumerateHelper(N, 0, seq, callback);
}

ContinuousPrimary PrimaryToContinuous(const Primary& pri) {
    ContinuousPrimary cpri(pri.size());
    for (size_t i = 0; i < pri.size(); ++i) {
        for (int b = 0; b < Base::N; ++b) {
            cpri[i][b] = pri[i] == b ? 1 : 0;
        }
    }
    return cpri;
}

}  // namespace mwmrna