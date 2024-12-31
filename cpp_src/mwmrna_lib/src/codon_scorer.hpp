#ifndef MWMRNA_LIB_CODON_SCORER_HPP
#define MWMRNA_LIB_CODON_SCORER_HPP

#include <cmath>
#include <vector>

#include "md_array.hpp"
#include "protein.hpp"

using std::vector;

namespace mwmrna::fold_codons {
class CodonScorer {
   public:
    CodonScorer(const CodonTable& codon_table, const AASeq& aa_seq, double lamdba = 1.0)
        : codon_table_(codon_table), aa_seq_(aa_seq) {
        n_ = aa_seq.size();
        max_codons_ = codon_table.MaxCodons();
        codons_.resize(aa_seq.size());
        for (int i = 0; i < n_; ++i) {
            codons_[i] = codon_table.GetCodons(aa_seq[i]);
        }
        scores_ = MdArray<double, 2>(0, {static_cast<int>(aa_seq.size()), max_codons_});
        for (int i = 0; i < static_cast<int>(aa_seq.size()); ++i) {
            for (int j = 0; j < static_cast<int>(codons_[i].size()); ++j) {
                scores_[{i, j}] =
                    -std::log(codon_table.CodonAdaptationWeight(codons_[i][j])) * lamdba;
            }
        }
    }

    vector<Codon> GetCodons(int idx) const { return codons_[idx]; }

    double GetScore(int idx, int cdn) { return scores_[{idx, cdn}]; }

    void BanCodon(int idx, Codon codon) {
        auto it = std::find(codons_[idx].begin(), codons_[idx].end(), codon);
        if (it == codons_[idx].end()) {
            throw std::invalid_argument("Banned codon (" + CodonToString(codon) +
                                        ") not found for amino acid " +
                                        AASeqToString(aa_seq_).substr(idx, 1));
        }
        int cdnid = it - codons_[idx].begin();
        codons_[idx].erase(it);
        for (int i = cdnid; i + 1 < static_cast<int>(codons_[idx].size()); ++i) {
            scores_[{idx, i}] = scores_[{idx, i + 1}];
        }
    }

    void ForceCodon(int idx, Codon codon) {
        codons_[idx] = {codon};
        for (int i = 0; i < max_codons_; ++i) {
            scores_[{idx, i}] = 0;
        }
    }

   private:
    int max_codons_, n_;
    CodonTable codon_table_;
    AASeq aa_seq_;
    MdArray<double, 2> scores_;
    vector<vector<Codon>> codons_;
};
};  // namespace mwmrna::fold_codons

#endif  // MWMRNA_LIB_CODON_SCORER_HPP