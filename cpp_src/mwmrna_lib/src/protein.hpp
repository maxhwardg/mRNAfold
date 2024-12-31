#ifndef MWMRNA_LIB_PROTEIN_HPP
#define MWMRNA_LIB_PROTEIN_HPP

#include <array>
#include <istream>
#include <string>
#include <vector>

#include "rna.hpp"

namespace mwmrna {

enum AminoAcid {
    ALA = 0,
    ARG = 1,
    ASN = 2,
    ASP = 3,
    CYS = 4,
    GLN = 5,
    GLU = 6,
    GLY = 7,
    HIS = 8,
    ILE = 9,
    LEU = 10,
    LYS = 11,
    MET = 12,
    PHE = 13,
    PRO = 14,
    SER = 15,
    THR = 16,
    TRP = 17,
    TYR = 18,
    VAL = 19,
    END = 20,
    NUM_AAS = 21
};

constexpr char kAALetters[] = "ARNDCQEGHILKMFPSTWYV*";

enum Codon {
    AAA = 0,
    AAU = 1,
    AAC = 2,
    AAG = 3,
    AUA = 4,
    AUU = 5,
    AUC = 6,
    AUG = 7,
    ACA = 8,
    ACU = 9,
    ACC = 10,
    ACG = 11,
    AGA = 12,
    AGU = 13,
    AGC = 14,
    AGG = 15,
    UAA = 16,
    UAU = 17,
    UAC = 18,
    UAG = 19,
    UUA = 20,
    UUU = 21,
    UUC = 22,
    UUG = 23,
    UCA = 24,
    UCU = 25,
    UCC = 26,
    UCG = 27,
    UGA = 28,
    UGU = 29,
    UGC = 30,
    UGG = 31,
    CAA = 32,
    CAU = 33,
    CAC = 34,
    CAG = 35,
    CUA = 36,
    CUU = 37,
    CUC = 38,
    CUG = 39,
    CCA = 40,
    CCU = 41,
    CCC = 42,
    CCG = 43,
    CGA = 44,
    CGU = 45,
    CGC = 46,
    CGG = 47,
    GAA = 48,
    GAU = 49,
    GAC = 50,
    GAG = 51,
    GUA = 52,
    GUU = 53,
    GUC = 54,
    GUG = 55,
    GCA = 56,
    GCU = 57,
    GCC = 58,
    GCG = 59,
    GGA = 60,
    GGU = 61,
    GGC = 62,
    GGG = 63,
    NUM_CODONS = 64,
};

/// @brief Stands for CoDon Sequence or CoDing Sequence
typedef std::vector<Codon> CDS;

/// @brief An amino acid sequence
typedef std::vector<AminoAcid> AASeq;

AASeq StringToAASeq(const std::string &aa_seq);
std::string AASeqToString(const AASeq &aa_seq);
AminoAcid ShortStringToAA(const std::string &aa);
std::string ConvertTtoU(const std::string &rna);
Codon StringToCodon(const std::string &codon);
Codon PrimaryToCodon(const Primary &primary);
Primary CodonToPrimary(Codon codon);
std::string CodonToString(mwmrna::Codon codon);
bool IsAALetter(char c);
CDS StringToCDS(const std::string &cds_str);
CDS PrimaryToCDS(const Primary &primary);
Primary CDSToPrimary(const CDS &cds);
std::string CDSToString(const CDS &cds);

class CodonTable {
   public:
    CodonTable(const std::string &filename);
    CodonTable(std::istream &file);
    CodonTable(std::istream &&file);

    CodonTable() = delete;
    CodonTable(const CodonTable &) = default;
    CodonTable(CodonTable &&) = default;
    CodonTable &operator=(const CodonTable &) = default;
    CodonTable &operator=(CodonTable &&) = default;

    double GetCodonFreq(Codon codon) const;
    double GetAAMaxFreq(AminoAcid aa) const;
    const std::vector<Codon> &GetCodons(AminoAcid aa) const;
    int MaxCodons() const;
    AminoAcid GetAA(Codon codon) const;
    double CodonAdaptationWeight(Codon codon) const;
    double CodonAdaptationIndex(const CDS &cds) const;
    double LogCodonAdaptationIdex(const CDS &cds) const;

   private:
    void Init(std::istream &file);
    std::array<AminoAcid, static_cast<size_t>(Codon::NUM_CODONS)> codon_to_aa_;
    std::array<std::vector<Codon>, static_cast<size_t>(AminoAcid::NUM_AAS)> aa_to_codons_;
    std::array<double, static_cast<size_t>(AminoAcid::NUM_AAS)> aa_max_freq_;
    std::array<double, static_cast<size_t>(Codon::NUM_CODONS)> codon_freq_;
    int max_codons_;
};

void EnumerateCds(const AASeq &aa_seq, const CodonTable &codon_table,
                  std::function<void(const CDS &)> callback);

}  // namespace mwmrna

#endif  // MWMRNA_LIB_PROTEIN_HPP