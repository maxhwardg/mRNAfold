#include "protein.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

using namespace std;

namespace mwmrna {
AASeq StringToAASeq(const std::string& aa_seq) {
    AASeq seq;
    for (auto c : aa_seq) {
        auto it = find(kAALetters, kAALetters + AminoAcid::NUM_AAS, c);
        if (it == kAALetters + AminoAcid::NUM_AAS) return AASeq();
        seq.push_back(static_cast<AminoAcid>(distance(kAALetters, it)));
    }
    return seq;
}

std::string AASeqToString(const AASeq& aa_seq) {
    string seq;
    for (auto aa : aa_seq) {
        if (aa == AminoAcid::NUM_AAS) return "";
        seq += kAALetters[static_cast<int>(aa)];
    }
    return seq;
}

AminoAcid ShortStringToAA(const std::string& aa) {
    static const map<string, AminoAcid> kMap = {
        {string("Ala"), AminoAcid::ALA}, {string("Arg"), AminoAcid::ARG},
        {string("Asn"), AminoAcid::ASN}, {string("Asp"), AminoAcid::ASP},
        {string("Cys"), AminoAcid::CYS}, {string("Gln"), AminoAcid::GLN},
        {string("Glu"), AminoAcid::GLU}, {string("Gly"), AminoAcid::GLY},
        {string("His"), AminoAcid::HIS}, {string("Ile"), AminoAcid::ILE},
        {string("Leu"), AminoAcid::LEU}, {string("Lys"), AminoAcid::LYS},
        {string("Met"), AminoAcid::MET}, {string("Phe"), AminoAcid::PHE},
        {string("Pro"), AminoAcid::PRO}, {string("Ser"), AminoAcid::SER},
        {string("Thr"), AminoAcid::THR}, {string("Trp"), AminoAcid::TRP},
        {string("Tyr"), AminoAcid::TYR}, {string("Val"), AminoAcid::VAL},
        {string("End"), AminoAcid::END}};

    auto it = kMap.find(aa);
    if (it == kMap.end()) return AminoAcid::NUM_AAS;
    return it->second;
}

std::string ConvertTtoU(const std::string& str) {
    string codon_str = str;
    for (auto& c : codon_str)
        if (c == 'T') c = 'U';
    return codon_str;
}

Codon StringToCodon(const std::string& codon) {
    return PrimaryToCodon(StringToPrimary(ConvertTtoU(codon)));
}

Codon PrimaryToCodon(const Primary& primary) {
    if (primary.size() != 3) return Codon::NUM_CODONS;
    int idx = 0;
    for (int i = 0; i < 3; ++i) {
        if (primary[i] == Base::N || primary[i] == Base::NUM_BASES) return Codon::NUM_CODONS;
        idx *= 4;
        idx += primary[i];
    }
    return static_cast<Codon>(idx);
}

Primary CodonToPrimary(Codon codon) {
    if (codon == Codon::NUM_CODONS) return Primary();
    Primary primary(3);
    int idx = static_cast<int>(codon);
    for (int i = 2; i >= 0; --i) {
        primary[i] = static_cast<Base>(idx % 4);
        idx /= 4;
    }
    return primary;
}

string CodonToString(mwmrna::Codon codon) {
    return mwmrna::PrimaryToString(mwmrna::CodonToPrimary(codon));
}

bool IsAALetter(char c) {
    return find(kAALetters, kAALetters + AminoAcid::NUM_AAS, c) != kAALetters + AminoAcid::NUM_AAS;
}

CDS StringToCDS(const std::string& cds_str) {
    return PrimaryToCDS(StringToPrimary(ConvertTtoU(cds_str)));
}

CDS PrimaryToCDS(const Primary& primary) {
    CDS cds;
    for (size_t i = 0; i < primary.size(); i += 3) {
        Codon codon = PrimaryToCodon(Primary(primary.begin() + i, primary.begin() + i + 3));
        if (codon == Codon::NUM_CODONS) return CDS();
        cds.push_back(codon);
    }
    return cds;
}

Primary CDSToPrimary(const CDS& cds) {
    Primary primary;
    for (auto codon : cds) {
        Primary codon_primary = CodonToPrimary(codon);
        if (codon_primary.size() != 3) return Primary();
        primary.insert(primary.end(), codon_primary.begin(), codon_primary.end());
    }
    return primary;
}

std::string CDSToString(const CDS &cds) {
    return PrimaryToString(CDSToPrimary(cds));
}

void EnumerateCdsRec(const AASeq& aa_seq, const CodonTable& codon_table, int aa_idx, CDS& cds,
                     std::function<void(const CDS&)> callback) {
    if (aa_idx == static_cast<int>(aa_seq.size())) {
        callback(cds);
        return;
    }
    for (auto codon : codon_table.GetCodons(aa_seq[aa_idx])) {
        cds.push_back(codon);
        EnumerateCdsRec(aa_seq, codon_table, aa_idx + 1, cds, callback);
        cds.pop_back();
    }
}

void EnumerateCds(const AASeq& aa_seq, const CodonTable& codon_table,
                  std::function<void(const CDS&)> callback) {
    CDS cds;
    EnumerateCdsRec(aa_seq, codon_table, 0, cds, callback);
}

CodonTable::CodonTable(const std::string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Could not open codon table file " + filename);
    }
    Init(file);
    file.close();
}

CodonTable::CodonTable(std::istream& file) { Init(file); }

CodonTable::CodonTable(std::istream&& file) { Init(file); }

double CodonTable::GetCodonFreq(Codon codon) const {
    return codon_freq_[static_cast<size_t>(codon)];
}

double CodonTable::GetAAMaxFreq(AminoAcid aa) const {
    return aa_max_freq_[static_cast<size_t>(aa)];
}

const std::vector<Codon>& CodonTable::GetCodons(AminoAcid aa) const {
    return aa_to_codons_[static_cast<size_t>(aa)];
}

int CodonTable::MaxCodons() const { return max_codons_; }

AminoAcid CodonTable::GetAA(Codon codon) const { return codon_to_aa_[static_cast<size_t>(codon)]; }

double CodonTable::CodonAdaptationWeight(Codon codon) const {
    return codon_freq_[static_cast<size_t>(codon)] /
           static_cast<double>(aa_max_freq_[static_cast<size_t>(GetAA(codon))]);
}

double CodonTable::CodonAdaptationIndex(const CDS& cds) const {
    double cai = 1;
    for (auto codon : cds) cai *= CodonAdaptationWeight(codon);
    return std::pow(cai, 1.0 / cds.size());
}

double CodonTable::LogCodonAdaptationIdex(const CDS& cds) const {
    double cai = 0;
    for (auto codon : cds) cai += std::log(CodonAdaptationWeight(codon));
    return cai / cds.size();
}

void CodonTable::Init(std::istream& file) {
    // Init members
    codon_freq_.fill(0);
    codon_to_aa_.fill(AminoAcid::NUM_AAS);
    aa_to_codons_.fill(vector<Codon>());
    aa_max_freq_.fill(0);
    max_codons_ = 0;

    string line;
    while (getline(file, line)) {
        // Strip leading whitespace
        line.erase(line.begin(),
                   find_if(line.begin(), line.end(), [](int ch) { return !isspace(ch); }));
        if (line.empty()) continue;
        stringstream ss(line);
        string aa_short_str, codon_str;
        double freq;
        ss >> aa_short_str >> codon_str >> freq;
        Codon codon = StringToCodon(codon_str);
        if (codon == Codon::NUM_CODONS) {
            throw runtime_error("Invalid codon in codon table: " + codon_str +
                                "\n Line was: " + line);
        }
        AminoAcid aa = ShortStringToAA(aa_short_str);
        if (aa == AminoAcid::NUM_AAS) {
            throw runtime_error("Invalid amino acid in codon table: " + aa_short_str);
        }
        int aa_idx = static_cast<int>(aa);
        int codon_idx = static_cast<int>(codon);
        codon_freq_[codon_idx] = freq;
        codon_to_aa_[codon_idx] = aa;
        aa_to_codons_[aa_idx].push_back(codon);
        aa_max_freq_[aa_idx] = max(aa_max_freq_[aa_idx], freq);
        max_codons_ = max(max_codons_, static_cast<int>(aa_to_codons_[aa_idx].size()));
    }
}

}  // namespace mwmrna