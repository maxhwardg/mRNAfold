#ifndef MWMRNA_LIB_VIENNA_PARAMS_HPP
#define MWMRNA_LIB_VIENNA_PARAMS_HPP

#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <istream>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "energy_semiring.hpp"
#include "logging.hpp"
#include "rna.hpp"

namespace mwmrna {

typedef std::vector<std::vector<std::string>> ViennaDataChunk;

/// @brief A class for parsing Vienna parameter files
template <typename EnergyT>
class ViennaParams {
   public:
    /// @brief Contains all the data extracted from the parameter file
    struct ViennaData {
        EnergyT hairpin_mismatch[NUM_ORDERED_PAIRS][NUM_BASES][NUM_BASES];
        EnergyT stack[NUM_ORDERED_PAIRS][NUM_ORDERED_PAIRS];
        EnergyT int11[NUM_ORDERED_PAIRS][NUM_ORDERED_PAIRS][NUM_BASES][NUM_BASES];
        EnergyT int21[NUM_ORDERED_PAIRS][NUM_ORDERED_PAIRS][NUM_BASES][NUM_BASES][NUM_BASES];
        EnergyT int22[NUM_ORDERED_PAIRS][NUM_ORDERED_PAIRS][NUM_BASES][NUM_BASES][NUM_BASES]
                     [NUM_BASES];
        EnergyT mismatch_interior[NUM_ORDERED_PAIRS][NUM_BASES][NUM_BASES];
        EnergyT mismatch_interior_1n[NUM_ORDERED_PAIRS][NUM_BASES][NUM_BASES];
        EnergyT mismatch_interior_23[NUM_ORDERED_PAIRS][NUM_BASES][NUM_BASES];
        EnergyT mismatch_multi[NUM_ORDERED_PAIRS][NUM_BASES][NUM_BASES];
        EnergyT mismatch_exterior[NUM_ORDERED_PAIRS][NUM_BASES][NUM_BASES];
        EnergyT dangle5[NUM_ORDERED_PAIRS][NUM_BASES];
        EnergyT dangle3[NUM_ORDERED_PAIRS][NUM_BASES];

        EnergyT il_asym;
        EnergyT il_asym_thrshold;
        EnergyT non_gc_closing_penalty;
        EnergyT duplex_init;
        EnergyT ml_unpaired;
        EnergyT ml_init;
        EnergyT ml_branch;

        std::vector<std::pair<Primary, EnergyT>> tetraloops;
        std::vector<std::pair<Primary, EnergyT>> triloops;
        std::vector<std::pair<Primary, EnergyT>> hexaloops;

        std::vector<EnergyT> bulge_init;
        std::vector<EnergyT> hairpin_init;
        std::vector<EnergyT> interior_init;
    };
    static constexpr BasePair kPairOrder[] = {
        CG, GC, GU, UG, AU, UA, NN,
    };
    static constexpr Base kBaseOrder[] = {
        N, A, C, G, U,
    };
    static constexpr double kLoopExtrapolationFactor = 1.07856;
    static constexpr int kMaxLoopSize = 30;
    static constexpr int kDgDenom = 100;
    static constexpr EnergyT kINF = 1000000;

    /// @brief Constructs from a .par file
    /// @param par_stream Path to a .par file
    /// @note Will throw an exception if the file cannot be opened
    /// @param f A function to apply to each parameter such as the boltzmann transform
    /// @note The function is applied to each parameter before it is stored
    ViennaParams(const std::string &par_file_path, std::function<EnergyT(EnergyT)> f) {
        std::ifstream par_stream(par_file_path);
        if (!par_stream.is_open()) {
            throw std::runtime_error("Could not open Vienna parameter file: " + par_file_path);
        }
        f_ = f;
        ProcessParStream(par_stream);
    }

    /// @brief Constructs from a stream
    /// @param par_stream Input stream with the contents of a .par file
    /// @note The stream must be open and readable
    /// @note The stream is consumed by this function
    /// @note The stream is not closed by this function
    /// @param f A function to apply to each parameter such as the boltzmann transform
    /// @note The function is applied to each parameter before it is stored
    ViennaParams(std::istream &par_stream, std::function<EnergyT(EnergyT)> f) {
        f_ = f;
        ProcessParStream(par_stream);
    }

    /// @brief Getter for the data struct
    /// @return An immutable reference to the data
    const ViennaData &Data() const { return data_; }

    ViennaParams &ExtrapolateLoopInits(int to_size) {
        ExtrapolateLoop(data_.bulge_init, to_size);
        ExtrapolateLoop(data_.hairpin_init, to_size);
        ExtrapolateLoop(data_.interior_init, to_size);
        return *this;
    }

    EnergyT LoopExtrapolation(int size) const {
        // Round would be more correct, but floor is used to reproduce Vienna's buggy behavior
        return f_(std::floor(kLoopExtrapolationFactor *
                             std::log(static_cast<double>(size) / kMaxLoopSize) * kDgDenom));
    }

   private:
    ViennaData data_;
    std::function<EnergyT(EnergyT)> f_;

    /// @brief Processess a stream containing a .par file
    /// @param par_stream the input stream
    void ProcessParStream(std::istream &par_stream) {
        std::string line;
        ViennaDataChunk chunk;
        std::string chunk_name;
        while (getline(par_stream, line)) {
            std::stringstream ss(line);

            // Deal with end of file indicators
            if (line.substr(0, 4) == "#END") break;
            if (line.substr(0, 5) == "# END") break;

            std::string first;

            // Ignore empty lines
            if (!(ss >> first)) continue;

            // Ignore line comments
            if (first == "##") continue;

            // Catch chunk names
            if (first == "#") {
                if (!chunk.empty()) {
                    ProcessChunk(chunk_name, chunk);
                    chunk.clear();
                }
                ss >> chunk_name;
                continue;
            }

            // Ignore comments
            if (first.substr(0, 2) == "/*") continue;

            std::string dg;
            chunk.push_back({first});
            while (ss >> dg) {
                if (dg.substr(0, 2) == "/*") break;
                chunk.back().push_back(dg);
            }
        }

        if (!chunk.empty()) ProcessChunk(chunk_name, chunk);
    }

    void ProcessChunk(const std::string &name, const ViennaDataChunk &chunk) {
        if (name == "stack") {
            ProcessStack(chunk);
        } else if (name == "bulge") {
            ProcessBulge(chunk);
        } else if (name == "hairpin") {
            ProcessHairpin(chunk);
        } else if (name == "mismatch_hairpin") {
            ProcessMismatchHairpin(chunk);
        } else if (name == "interior") {
            ProcessInterior(chunk);
        } else if (name == "NINIO") {
            ProcessNinio(chunk);
        } else if (name == "Misc") {
            ProcessMisc(chunk);
        } else if (name == "ML_params") {
            ProcessMLParams(chunk);
        } else if (name == "int11") {
            ProcessInt11(chunk);
        } else if (name == "int21") {
            ProcessInt21(chunk);
        } else if (name == "int22") {
            ProcessInt22(chunk);
        } else if (name == "mismatch_interior") {
            ProcessMismatchInterior(chunk);
        } else if (name == "mismatch_interior_1n") {
            ProcessMismatchInterior1n(chunk);
        } else if (name == "mismatch_interior_23") {
            ProcessMismatchInterior23(chunk);
        } else if (name == "mismatch_multi") {
            ProcessMismatchMulti(chunk);
        } else if (name == "mismatch_exterior") {
            ProcessMismatchExterior(chunk);
        } else if (name == "dangle5") {
            ProcessDangle5(chunk);
        } else if (name == "dangle3") {
            ProcessDangle3(chunk);
        } else if (name == "Tetraloops") {
            ProcessTetraloop(chunk);
        } else if (name == "Triloops") {
            ProcessTriloop(chunk);
        } else if (name == "Hexaloops") {
            ProcessHexaloop(chunk);
        } else {
            mwmrna::log("Ignoring chunk: " + name);
        }
    }

    void ProcessStack(const ViennaDataChunk &chunk) {
        for (int i = 0; i < NUM_ORDERED_PAIRS; ++i) {
            auto nums = ConvertToNums(chunk[i]);
            for (int j = 0; j < NUM_ORDERED_PAIRS; ++j) {
                data_.stack[kPairOrder[i]][kPairOrder[j]] = nums[j];
            }
        }
    }

    void ProcessBulge(const ViennaDataChunk &chunk) { FillInits(data_.bulge_init, chunk); }

    void ProcessHairpin(const ViennaDataChunk &chunk) { FillInits(data_.hairpin_init, chunk); }

    void ProcessMismatchHairpin(const ViennaDataChunk &chunk) {
        auto nums = ConvertToNums(chunk);
        for (int i = 0; i < NUM_ORDERED_PAIRS; ++i) {
            for (int j = 0; j < NUM_BASES; ++j) {
                for (int k = 0; k < NUM_BASES; ++k) {
                    data_.hairpin_mismatch[kPairOrder[i]][kBaseOrder[j]][kBaseOrder[k]] =
                        nums[i * NUM_BASES + j][k];
                }
            }
        }
    }

    void ProcessInterior(const ViennaDataChunk &chunk) { FillInits(data_.interior_init, chunk); }

    void ProcessNinio(const ViennaDataChunk &chunk) {
        auto nums = ConvertToNums(chunk);
        data_.il_asym = nums[0][0];
        data_.il_asym_thrshold = nums[0][2];
    }

    void ProcessMisc(const ViennaDataChunk &chunk) {
        auto nums = ConvertToNums(chunk);
        data_.duplex_init = nums[0][0];
        data_.non_gc_closing_penalty = nums[0][2];
    }

    void ProcessMLParams(const ViennaDataChunk &chunk) {
        auto nums = ConvertToNums(chunk);
        data_.ml_unpaired = nums[0][0];
        data_.ml_init = nums[0][2];
        data_.ml_branch = nums[0][4];
    }

    void ProcessInt11(const ViennaDataChunk &chunk) {
        auto nums = ConvertToNums(chunk);
        for (int i = 0; i < NUM_ORDERED_PAIRS; ++i) {
            for (int j = 0; j < NUM_ORDERED_PAIRS; ++j) {
                for (int k = 0; k < NUM_BASES; ++k) {
                    for (int l = 0; l < NUM_BASES; ++l) {
                        int line_idx = i * NUM_ORDERED_PAIRS * NUM_BASES + j * NUM_BASES + k;
                        data_.int11[kPairOrder[i]][kPairOrder[j]][kBaseOrder[k]][kBaseOrder[l]] =
                            nums[line_idx][l];
                    }
                }
            }
        }
    }

    void ProcessInt21(const ViennaDataChunk &chunk) {
        auto nums = ConvertToNums(chunk);
        for (int i = 0; i < NUM_ORDERED_PAIRS; ++i) {
            for (int j = 0; j < NUM_ORDERED_PAIRS; ++j) {
                for (int k = 0; k < NUM_BASES; ++k) {
                    for (int l = 0; l < NUM_BASES; ++l) {
                        for (int m = 0; m < NUM_BASES; ++m) {
                            int line_idx = i * NUM_ORDERED_PAIRS * NUM_BASES * NUM_BASES +
                                           j * NUM_BASES * NUM_BASES + k * NUM_BASES + l;
                            data_.int21[kPairOrder[i]][kPairOrder[j]][kBaseOrder[k]][kBaseOrder[l]]
                                       [kBaseOrder[m]] = nums[line_idx][m];
                        }
                    }
                }
            }
        }
    }

    void ProcessInt22(const ViennaDataChunk &chunk) {
        auto nums = ConvertToNums(chunk);
        // There are no entries for NN pairs or N bases
        // So the loops go up to -1
        int num_pairs = NUM_ORDERED_PAIRS - 1;
        int num_bases = NUM_BASES - 1;
        for (int i = 0; i < num_pairs; ++i) {
            for (int j = 0; j < num_pairs; ++j) {
                for (int k = 0; k < num_bases; ++k) {
                    for (int l = 0; l < num_bases; ++l) {
                        for (int m = 0; m < num_bases; ++m) {
                            for (int n = 0; n < num_bases; ++n) {
                                int line_idx = i * num_pairs * num_bases * num_bases * num_bases +
                                               j * num_bases * num_bases * num_bases +
                                               k * num_bases * num_bases + l * num_bases + m;
                                // Base indexes are shifted by 1
                                // Since kBaseOrder[0] == N
                                data_.int22[kPairOrder[i]][kPairOrder[j]][kBaseOrder[k + 1]]
                                           [kBaseOrder[l + 1]][kBaseOrder[m + 1]]
                                           [kBaseOrder[n + 1]] = nums[line_idx][n];
                            }
                        }
                    }
                }
            }
        }
    }

    void ProcessMismatchInterior(const ViennaDataChunk &chunk) {
        FillMismatch(data_.mismatch_interior, chunk);
    }

    void ProcessMismatchInterior1n(const ViennaDataChunk &chunk) {
        FillMismatch(data_.mismatch_interior_1n, chunk);
    }

    void ProcessMismatchInterior23(const ViennaDataChunk &chunk) {
        FillMismatch(data_.mismatch_interior_23, chunk);
    }

    void ProcessMismatchMulti(const ViennaDataChunk &chunk) {
        FillMismatch(data_.mismatch_multi, chunk);
    }

    void ProcessMismatchExterior(const ViennaDataChunk &chunk) {
        FillMismatch(data_.mismatch_exterior, chunk);
    }

    void ProcessDangle5(const ViennaDataChunk &chunk) { FillDangle(data_.dangle5, chunk); }

    void ProcessDangle3(const ViennaDataChunk &chunk) { FillDangle(data_.dangle3, chunk); }

    void ProcessTetraloop(const ViennaDataChunk &chunk) {
        FillSpecialHairpin(data_.tetraloops, chunk);
    }

    void ProcessTriloop(const ViennaDataChunk &chunk) { FillSpecialHairpin(data_.triloops, chunk); }

    void ProcessHexaloop(const ViennaDataChunk &chunk) {
        FillSpecialHairpin(data_.hexaloops, chunk);
    }

    EnergyT ParseViennaNumber(const std::string &s) {
        if (s == "INF") return kINF;
        return f_(static_cast<EnergyT>(stoi(s)));
    }

    std::vector<EnergyT> ConvertToNums(const std::vector<std::string> &toks) {
        std::vector<EnergyT> nums;
        std::transform(toks.begin(), toks.end(), back_inserter(nums),
                       [this](const std::string &tok) { return this->ParseViennaNumber(tok); });
        return nums;
    }

    std::vector<std::vector<EnergyT>> ConvertToNums(const ViennaDataChunk &chunk) {
        std::vector<std::vector<EnergyT>> nums;
        std::transform(
            chunk.begin(), chunk.end(), back_inserter(nums),
            [this](const std::vector<std::string> &toks) { return this->ConvertToNums(toks); });
        return nums;
    }

    void FillMismatch(EnergyT (&table)[NUM_ORDERED_PAIRS][NUM_BASES][NUM_BASES],
                      const ViennaDataChunk &chunk) {
        auto nums = ConvertToNums(chunk);
        for (int i = 0; i < NUM_ORDERED_PAIRS; ++i) {
            for (int j = 0; j < NUM_BASES; ++j) {
                for (int k = 0; k < NUM_BASES; ++k) {
                    table[kPairOrder[i]][kBaseOrder[j]][kBaseOrder[k]] = nums[i * NUM_BASES + j][k];
                }
            }
        }
    }

    void FillDangle(EnergyT (&table)[NUM_ORDERED_PAIRS][NUM_BASES], const ViennaDataChunk &chunk) {
        auto nums = ConvertToNums(chunk);
        for (int i = 0; i < NUM_ORDERED_PAIRS; ++i) {
            for (int j = 0; j < NUM_BASES; ++j) {
                table[kPairOrder[i]][kBaseOrder[j]] = nums[i][j];
            }
        }
    }

    void FillSpecialHairpin(std::vector<std::pair<Primary, EnergyT>> &loops,
                            const ViennaDataChunk &chunk) {
        for (auto &toks : chunk)
            loops.emplace_back(StringToPrimary(toks[0]), ParseViennaNumber(toks[1]));
    }

    void ExtrapolateLoop(std::vector<EnergyT> &loop, int to_size) {
        while (loop.size() < to_size + 1) {
            loop.push_back(loop[kMaxLoopSize] + LoopExtrapolation(loop.size()));
        }
    }

    void FillInits(std::vector<EnergyT> &inits, const ViennaDataChunk &chunk) {
        for (auto &toks : chunk) {
            auto nums = ConvertToNums(toks);
            inits.insert(inits.end(), nums.begin(), nums.end());
        }
    }
};

/// @brief A class for parsing Vienna parameter files using the raw units from the files
class ViennaParamsRaw : public ViennaParams<int> {
   public:
    using semiring_type = MFESemiring<int>;
    ViennaParamsRaw(const std::string &par_file_path) : ViennaParams(par_file_path, Identity) {}
    ViennaParamsRaw(std::istream &par_stream) : ViennaParams(par_stream, Identity) {}

   private:
    static int Identity(int x) { return x; }
};

/// @brief A class for parsing Vienna parameter files with units scaled to kcal/mol
class ViennaParamsScaled : public ViennaParams<double> {
   public:
    using semiring_type = MFESemiring<double>;
    ViennaParamsScaled(const std::string &par_file_path) : ViennaParams(par_file_path, Scale) {}
    ViennaParamsScaled(std::istream &par_stream) : ViennaParams(par_stream, Scale) {}
    static double Scale(double x) { return x / ViennaParams::kDgDenom; }
};

class ViennaParamsBoltz : public ViennaParams<double> {
   public:
    static constexpr double kCellTempK = 310.15;
    static constexpr double kGasConstant = 1.98717e-3;
    static constexpr double kBeta = 1.0 / (kGasConstant * kCellTempK);
    using semiring_type = BoltzSemiring<double>;

    ViennaParamsBoltz(const std::string &par_file_path)
        : ViennaParams(par_file_path, ScaledBoltz) {}
    ViennaParamsBoltz(std::istream &par_stream) : ViennaParams(par_stream, ScaledBoltz) {}
    static double ScaledBoltz(double x) { return Boltz(ViennaParamsScaled::Scale(x)); }
    static double Boltz(double x) { return std::exp(-kBeta * x); }
    static double InverseBoltz(double x) { return -kGasConstant * kCellTempK * std::log(x); }
};

// TODO: Make a ViennaParramsBoltzmann class

}  // namespace mwmrna

#endif  // MWMRNA_LIB_VIENNA_PARAMS_HPP
