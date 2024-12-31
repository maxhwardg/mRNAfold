#include "structure.hpp"

#include <functional>
#include <stack>

namespace mwmrna {
using namespace std;

Matching DbToMatching(const std::string& str) {
    Matching match(str.size(), Matching::kUnmatched);
    stack<int> stk;
    for (int i = 0; i < static_cast<int>(str.size()); ++i) {
        if (str[i] == '(') {
            stk.push(i);
        } else if (str[i] == ')') {
            match[i] = stk.top();
            match[stk.top()] = i;
            stk.pop();
        }
    }
    return match;
}

string MatchingToDb(const Matching& m) {
    string str(m.size(), '.');
    for (int i = 0; i < static_cast<int>(m.size()); ++i) {
        if (m[i] != Matching::kUnmatched && m[i] > i) {
            str[i] = '(';
            str[m[i]] = ')';
        }
    }
    return str;
}

namespace {
void EnumerateRec(const Primary& pri, const std::function<void(const Matching&)>& callback,
                  int hairpin_min, Matching& m, vector<int>& stk, int i) {
    if (i == static_cast<int>(pri.size())) {
        if (stk.empty()) callback(m);
        return;
    }
    EnumerateRec(pri, callback, hairpin_min, m, stk, i + 1);
    if (pri.size() - i - 1 >= stk.size() + 1) {
        stk.push_back(i);
        EnumerateRec(pri, callback, hairpin_min, m, stk, i + 1);
        stk.pop_back();
    }
    if (!stk.empty() && IsValidPair(pri[i], pri[stk.back()]) && i - stk.back() - 1 >= hairpin_min) {
        int j = stk.back();
        stk.pop_back();
        m[i] = j;
        m[j] = i;
        EnumerateRec(pri, callback, hairpin_min, m, stk, i + 1);
        m[i] = m[j] = Matching::kUnmatched;
        stk.push_back(j);
    }
}
}  // namespace

void EnumerateStructures(const Primary& pri, const std::function<void(const Matching&)>& callback,
                         int hairpin_min) {
    Matching m(pri.size(), Matching::kUnmatched);
    vector<int> stk;
    EnumerateRec(pri, callback, hairpin_min, m, stk, 0);
}

int NumUnapired(const Matching& m) {
    int num = 0;
    for (int i = 0; i < static_cast<int>(m.size()); ++i) {
        if (m[i] == Matching::kUnmatched) ++num;
    }
    return num;
}

PairTree::PairTree(const Matching& m) : m_(m) {
    children_.resize(1);
    pairs_.resize(1);
    root_ = 0;
    pairs_[root_] = {-1, m_.size()};
    Build(root_);
}

int PairTree::Children(int i) const { return children_[i].size(); }

int PairTree::Child(int i, int j) const { return children_[i][j]; }

std::pair<int, int> PairTree::Pair(int i) const { return pairs_[i]; }

void PairTree::Build(int i) {
    int at = pairs_[i].first + 1;
    while (at < pairs_[i].second) {
        if (m_[at] == Matching::kUnmatched) {
            ++at;
        } else {
            pairs_.emplace_back(at, m_[at]);
            int ch_i = children_.size();
            children_[i].push_back(ch_i);
            children_.emplace_back();
            Build(ch_i);
            at = m_[at] + 1;
        }
    }
}

}  // namespace mwmrna