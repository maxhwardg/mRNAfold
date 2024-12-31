#ifndef MWMRNA_LIB_STRUCTURE_HPP
#define MWMRNA_LIB_STRUCTURE_HPP

#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "rna.hpp"

namespace mwmrna {

class Matching : public std::vector<int> {
   public:
    static constexpr int kUnmatched = -1;
    using std::vector<int>::vector;
};

Matching DbToMatching(const std::string &str);

std::string MatchingToDb(const Matching &m);

class PairTree {
   public:
    PairTree(const Matching &m);

    int Children(int i) const;
    int Child(int i, int j) const;
    std::pair<int, int> Pair(int i) const;
    int Root() const { return root_; }

   private:
    void Build(int i);
    Matching m_;
    int root_;
    std::vector<std::vector<int>> children_;
    std::vector<std::pair<int, int>> pairs_;
};

void EnumerateStructures(const Primary &pri, const std::function<void(const Matching &)> &callback,
                         int hairpin_min = 3);

int NumUnapired(const Matching &m);

}  // namespace mwmrna

#endif  // MWMRNA_LIB_STRUCTURE_HPP