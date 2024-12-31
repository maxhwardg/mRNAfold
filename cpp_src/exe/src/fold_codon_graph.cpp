#include <algorithm>
#include <cmath>
#include <execution>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "codon_graph.hpp"
#include "codon_graph_fold.hpp"
#include "codon_graph_trace.hpp"
#include "codon_table_strings.hpp"
#include "energy_model.hpp"
#include "energy_semiring.hpp"
#include "fold.hpp"
#include "logging.hpp"
#include "protein.hpp"
#include "rna.hpp"
#include "structure.hpp"
#include "vienna_param_strings.hpp"
#include "vienna_params.hpp"

using namespace std;
using ParamT = mwmrna::ViennaParamsScaled;
using EnergyModelT = mwmrna::EnergyModel<ParamT, 0>;
using SemiringT = typename ParamT::semiring_type;
using EnergyT = typename SemiringT::energy_type;

const map<string, string> kOrganismCodonFreqStr = {{"human", mwmrna::kHomosapiensCodonFreq},
                                                   {"mouse", mwmrna::kMouseCodonFreq}};
const double kBannedPenalty = 100;
const double kPfnNtFactor = 0.37;

struct Config {
    mwmrna::AASeq aa_seq;
    unique_ptr<mwmrna::CodonTable> codon_table;
    vector<vector<EnergyT>> codon_score;
    mwmrna::Primary utr_5p, utr_3p;
    bool parallel = false;
    double subopt_randomness = 0.0;
    double lambda = 1.0;
    int num_subopt_traces = 100;
    double keep_chance = 1.0;
    vector<pair<int, int>> encourage_unpaired;
    int num_subopt_rounds = 1;
    int rng_seed = 0;
};

void AddCodonScore(Config& config, const vector<pair<int, mwmrna::Codon>>& banned_codons,
                   const vector<pair<int, mwmrna::Codon>>& forced_codons) {
    config.codon_score.assign(config.aa_seq.size(),
                              vector<EnergyT>(config.codon_table->MaxCodons()));
    for (int i = 0; i < static_cast<int>(config.aa_seq.size()); ++i) {
        auto codons = config.codon_table->GetCodons(config.aa_seq[i]);
        for (int j = 0; j < static_cast<int>(codons.size()); ++j) {
            config.codon_score[i][j] =
                -std::log(config.codon_table->CodonAdaptationWeight(codons[j])) * config.lambda;
        }
    }
    for (auto [pos, codon] : banned_codons) {
        auto codons = config.codon_table->GetCodons(config.aa_seq[pos]);
        auto it = std::find(codons.begin(), codons.end(), codon);
        if (it == codons.end()) {
            throw std::runtime_error("Banned codon not found: " + std::to_string(pos) + " " +
                                     mwmrna::CodonToString(codon));
        }
        config.codon_score[pos][it - codons.begin()] = kBannedPenalty;
    }
    for (auto [pos, codon] : forced_codons) {
        auto codons = config.codon_table->GetCodons(config.aa_seq[pos]);
        auto it = std::find(codons.begin(), codons.end(), codon);
        if (it == codons.end()) {
            throw std::runtime_error("Forced codon not found: " + std::to_string(pos) + " " +
                                     mwmrna::CodonToString(codon));
        }
        for (int i = 0; i < static_cast<int>(codons.size()); ++i) {
            if (i != it - codons.begin()) {
                config.codon_score[pos][i] = kBannedPenalty;
            }
        }
    }
}

Config ParseConfigFile(const string& filename) {
    std::ifstream ifs(filename);
    if (!ifs) {
        throw std::runtime_error("Could not open file " + filename);
    }
    string name, value;
    Config config;
    vector<pair<int, mwmrna::Codon>> banned_codons, forced_codons;
    // Init codon freq table
    std::stringstream ss(mwmrna::kHomosapiensCodonFreq);
    config.codon_table = make_unique<mwmrna::CodonTable>(ss);
    for (string line; std::getline(ifs, line);) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        iss >> name >> value;
        if (name == "aa_seq") {
            config.aa_seq = mwmrna::StringToAASeq(value);
        } else if (name == "forced_codon" || name == "banned_codon") {
            string cdn;
            iss >> cdn;
            auto codon = mwmrna::StringToCodon(cdn);

            if (codon == mwmrna::Codon::NUM_CODONS) {
                throw std::runtime_error("Couldn't parse forced codon");
            }
            int aa_pos = stoi(value);
            if (aa_pos >= static_cast<int>(config.aa_seq.size()) || aa_pos < 0) {
                throw std::runtime_error("Forced codon position out of bounds");
            }
            if (name == "forced_codon") {
                forced_codons.emplace_back(aa_pos, codon);
            } else {
                banned_codons.emplace_back(aa_pos, codon);
            }
        } else if (name == "parallel") {
            config.parallel = value == "true";
        } else if (name == "subopt_randomness") {
            config.subopt_randomness = stod(value);
        } else if (name == "lambda") {
            config.lambda = stod(value);
        } else if (name == "num_subopt_traces") {
            config.num_subopt_traces = stoi(value);
        } else if (name == "utr_5p") {
            config.utr_5p = mwmrna::StringToPrimary(value);
        } else if (name == "utr_3p") {
            config.utr_3p = mwmrna::StringToPrimary(value);
        } else if (name == "keep_chance") {
            config.keep_chance = stod(value);
        } else if (name == "encourage_unpaired") {
            // split by comma
            auto it = value.find(',');
            int i = stoi(value.substr(0, it));
            int j = stoi(value.substr(it + 1));
            config.encourage_unpaired.push_back({i, j});
        } else if (name == "num_subopt_rounds") {
            config.num_subopt_rounds = stoi(value);
        } else if (name == "rng_seed") {
            config.rng_seed = stoi(value);
        } else if (name == "organism") {
            // Load codon table
            auto it = kOrganismCodonFreqStr.find(value);
            if (it == kOrganismCodonFreqStr.end()) {
                throw std::runtime_error("No codon table for organism: " + value);
            }
            std::clog << "Using codon table for organism: " << value << std::endl;
            config.codon_table = make_unique<mwmrna::CodonTable>(std::istringstream(it->second));
        } else if (name == "codon_freq_path") {
            std::ifstream ifs(value);
            std::clog << "Using custom codon table: " << value << std::endl;
            if (!ifs) {
                throw std::runtime_error("Could not open file " + value);
            }
            config.codon_table = make_unique<mwmrna::CodonTable>(ifs);
        } else {
            throw std::runtime_error("Unknown config option: " + name);
        }
    }

    // Compute codon scores
    AddCodonScore(config, banned_codons, forced_codons);

    return config;
}

int FindMotif(const mwmrna::Primary& primary, const mwmrna::Primary& motif) {
    for (int i = 0; i < static_cast<int>(primary.size()) - static_cast<int>(motif.size()) + 1;
         ++i) {
        bool found = true;
        for (int j = 0; j < static_cast<int>(motif.size()); ++j) {
            if (primary[i + j] != motif[j]) {
                found = false;
                break;
            }
        }
        if (found) {
            return i;
        }
    }
    return -1;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file_path>" << std::endl;
        return 1;
    }
    mwmrna::allow_logging = false;
    Config config = ParseConfigFile(argv[1]);

    auto& codon_table = *config.codon_table;
    std::stringstream ss(mwmrna::kViennaParams2004);
    ParamT params(ss);
    EnergyModelT em(params);
    std::stringstream ssb(mwmrna::kViennaParams2004);
    mwmrna::ViennaParamsBoltz boltz_params(ssb);
    mwmrna::EnergyModel<mwmrna::ViennaParamsBoltz, 0> boltz_em(boltz_params);

    while (true) {
        mwmrna::fold_codons::CodonGraph graph(mwmrna::fold_codons::CodonGraphConfig<SemiringT>{
            .aa_seq = config.aa_seq,
            .codon_table = codon_table,
            .codon_score =
                [&config](int aa_pos, mwmrna::Codon codon) {
                    auto codons = config.codon_table->GetCodons(config.aa_seq[aa_pos]);
                    auto it = std::find(codons.begin(), codons.end(), codon);
                    if (it == codons.end()) {
                        return kBannedPenalty;
                    }
                    return config.codon_score[aa_pos][static_cast<int>(it - codons.begin())];
                },
            .utr_5p = config.utr_5p,
            .utr_3p = config.utr_3p,
        });

        vector<vector<EnergyT>> pair_bonus(graph.Nts(), vector<EnergyT>(graph.Nts(), 0));

        for (auto [i, j] : config.encourage_unpaired) {
            for (int x = i; x <= j; ++x) {
                for (int y = 0; y < graph.Nts(); ++y) {
                    pair_bonus[x][y] = kBannedPenalty;
                    pair_bonus[y][x] = kBannedPenalty;
                }
            }
        }

        mwmrna::fold_codons::GraphFoldContext<EnergyModelT> ctx{
            .em = em,
            .graph = graph,
            .parallel = config.parallel,
            .pair_bonus = pair_bonus,
        };

        mwmrna::fold_codons::GraphFold fold(ctx);

        std::clog << "DP table filled. Suboptimal trace generation starting..." << std::endl;

        std::cout << std::fixed << std::setprecision(4);


        mwmrna::FoldConfig<decltype(boltz_em)> boltz_cfg{
            .em = boltz_em,
            .nt_factor = kPfnNtFactor,
            .parallel = ctx.parallel,
        };

        auto do_traces = [&](int seed) {
            mwmrna::fold_codons::GraphTraceIterator tracer(
                fold, mwmrna::fold_codons::GraphTraceContext<EnergyModelT>{
                          .max_pq_size = config.num_subopt_traces,
                          .subopt_randomness = config.subopt_randomness,
                          .rng_seed = static_cast<int>(seed),
                          .keep_chance = config.keep_chance,
                      });
            vector<tuple<double, double, mwmrna::fold_codons::CodonGraphTrace<EnergyT>>> results;
            for (int i = 0; i < config.num_subopt_traces && tracer.HasNext(); ++i) {
                auto trace = tracer.Next();

                auto pfn = mwmrna::Fold(boltz_cfg).Inside(trace.pri).energy;
                double efe = (log(pfn) - log(boltz_cfg.nt_factor) * trace.pri.size()) *
                             (-1 / mwmrna::ViennaParamsBoltz::kBeta);
                auto cai_pri_str = mwmrna::PrimaryToString(trace.pri)
                                       .substr(config.utr_5p.size())
                                       .substr(0, config.aa_seq.size() * 3);
                double cai = codon_table.CodonAdaptationIndex(mwmrna::StringToCDS(cai_pri_str));
                results.emplace_back(efe, cai, trace);

                if (i % 10 == 0) {
                    string message = "Subopt round with rng seed=" + to_string(seed) + ": " +
                                     to_string(i) + " traces generated...";
                    std::clog << message << std::endl;
                }
            }
            return results;
        };

        default_random_engine rng(config.rng_seed);
        vector<int> seeds(config.num_subopt_rounds);
        std::generate(seeds.begin(), seeds.end(), [&rng]() { return rng(); });

        vector<vector<tuple<double, double, mwmrna::fold_codons::CodonGraphTrace<EnergyT>>>>
            all_results(config.num_subopt_rounds);

        if (config.parallel) {
            std::transform(std::execution::par_unseq, seeds.begin(), seeds.end(),
                           all_results.begin(), [&](int seed) { return do_traces(seed); });
        } else {
            for (int i = 0; i < config.num_subopt_rounds; ++i) {
                all_results[i] = do_traces(seeds[i]);
            }
        }

        for (int i = 0; i < config.num_subopt_rounds; ++i) {
            for (const auto& [efe, cai, trace] : all_results[i]) {
                std::cout << mwmrna::PrimaryToString(trace.pri) << std::endl;
                std::cout << mwmrna::MatchingToDb(trace.match) << std::endl;
                std::cout << "Score: " << trace.energy << std::endl;
                std::cout << "EFE: " << efe << std::endl << "CAI: " << cai << std::endl;
            }
        }

        break;
    }
}
