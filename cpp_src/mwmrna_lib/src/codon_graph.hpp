#ifndef MWMRNA_LIB_CODON_GRAPH_HPP
#define MWMRNA_LIB_CODON_GRAPH_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <tuple>

#include "md_array.hpp"
#include "protein.hpp"
#include "rna.hpp"

using std::function;
using std::pair;
using std::tuple;
using std::vector;

namespace mwmrna::fold_codons {

template <typename SemiringT>
struct CodonGraphConfig {
    using EnergyT = typename SemiringT::energy_type;
    AASeq aa_seq;
    CodonTable& codon_table;
    // (codon index, codon) -> score, usually CAI
    function<EnergyT(int, Codon)> codon_score;
    Primary utr_5p{};
    Primary utr_3p{};
    // TODO: Banned and forced codons
};

template <typename EnergyT>
struct CodonGraphEdge {
    int to;
    Base base;
    EnergyT weight;
};

template <typename T>
bool operator==(const CodonGraphEdge<T>& lhs, const CodonGraphEdge<T>& rhs) {
    return std::tie(lhs.to, lhs.base, lhs.weight) == std::tie(rhs.to, rhs.base, rhs.weight);
}

template <typename SemiringT>
class CodonGraph {
   private:
    using EnergyT = typename SemiringT::energy_type;
    vector<CodonGraphEdge<EnergyT>> edges;
    vector<std::pair<int, int>> edges_l2r;
    vector<std::pair<int, int>> edges_r2l;
    vector<vector<int>> nt_nodes;
    vector<Base> node_base;
    int num_nodes;
    CodonGraphConfig<SemiringT> config;

    /// @brief Builds a codon graph assuming a variant of the standard codon table is used. In
    /// particular, a heuristic is used where the 2nt 5' prefix of each codon is assumed to have at
    /// most 2 different combinations. Also, assumes CAI is used to score codons.
    /// @return The codon graph.
    void BuildStandardCodonGraph() {
        using EnergyT = typename SemiringT::energy_type;
        SemiringT sr;
        vector<vector<CodonGraphEdge<EnergyT>>> adj;

        this->nt_nodes.resize(config.aa_seq.size() * 3 + config.utr_5p.size() +
                              config.utr_3p.size());

        auto add_node = [&](Base base, int nt) {
            int id = adj.size();
            this->nt_nodes[nt].push_back({static_cast<int>(adj.size())});
            adj.push_back({});
            this->node_base.push_back(base);
            return id;
        };

        vector<int> prev_nodes;
        for (int i = 0; i < static_cast<int>(config.utr_5p.size()); ++i) {
            int id = add_node(config.utr_5p[i], i);
            for (int nid : prev_nodes) {
                adj[nid].push_back({id, config.utr_5p[i], sr.One()});
            }
            prev_nodes = {id};
        }

        for (int i = 0; i < static_cast<int>(config.aa_seq.size()); ++i) {
            auto codons = config.codon_table.GetCodons(config.aa_seq[i]);
            std::map<Primary, std::set<Base>> pref_map;
            for (const auto& codon : codons) {
                auto prim = CodonToPrimary(codon);
                pref_map[{prim[0], prim[1]}].insert(prim[2]);
            }

            if (pref_map.size() > 2) {
                std::stringstream ss;
                for (const auto& codon : codons) {
                    ss << CodonToString(codon) << ", ";
                }
                throw std::invalid_argument(
                    "BuildStandardCodonGraph assumes that there are at most 2 codon prefixes for "
                    "each "
                    "amino acid. Too many codon prefixes for amino acid " +
                    AASeqToString(config.aa_seq).substr(i, 1) + ". Codons: " + ss.str());
            }

            int base_adj_idx = adj.size();

            // Add prefix nodes and preceding transitions
            int pref_idx = 0;
            vector<Base> all_bases;
            for (auto& [pref, bases] : pref_map) {
                int id0 = add_node(pref[0], config.utr_5p.size() + i * 3);
                int id1 = add_node(pref[1], config.utr_5p.size() + i * 3 + 1);
                for (int nid : prev_nodes) {
                    adj[nid].push_back({id0, pref[0], sr.One()});
                }
                adj[id0].push_back({id1, pref[1], sr.One()});
                for (Base b : bases) all_bases.push_back(b);
                pref_idx++;
            }

            // Add 3rd base nodes and preceding transitions

            // Remove duplicates from all_bases
            std::sort(all_bases.begin(), all_bases.end());
            all_bases.erase(std::unique(all_bases.begin(), all_bases.end()), all_bases.end());
            prev_nodes.clear();
            //  Add 3rd base nodes
            for (Base b : all_bases) {
                prev_nodes.push_back(add_node(b, config.utr_5p.size() + i * 3 + 2));
            }
            // Add transitions from 2nd base to 3rd base
            pref_idx = 0;
            for (auto& [pref, bases] : pref_map) {
                for (Base b : bases) {
                    int bid = std::distance(all_bases.begin(),
                                            std::find(all_bases.begin(), all_bases.end(), b));
                    auto codon = PrimaryToCodon({pref[0], pref[1], b});
                    EnergyT weight = config.codon_score(i, codon);
                    adj[base_adj_idx + pref_idx * 2 + 1].push_back(
                        {base_adj_idx + static_cast<int>(pref_map.size()) * 2 + bid, b, weight});
                }
                ++pref_idx;
            }
        }

        for (int i = 0; i < static_cast<int>(config.utr_3p.size()); ++i) {
            if (config.utr_3p[i] == N) {
                // Special case for 3' UTR. Allow N to represent any base.
                vector<int> ids;
                for (Base b : {A, U, G, C}) {
                    int id = add_node(b, config.utr_5p.size() + config.aa_seq.size() * 3 + i);
                    ids.push_back(id);
                    for (int nid : prev_nodes) {
                        adj[nid].push_back({id, b, sr.One()});
                    }
                }
                prev_nodes = ids;
            } else {
                int id =
                    add_node(config.utr_3p[i], config.utr_5p.size() + config.aa_seq.size() * 3 + i);
                for (int nid : prev_nodes) {
                    adj[nid].push_back({id, config.utr_3p[i], sr.One()});
                }
                prev_nodes = {id};
            }
        }

        this->num_nodes = adj.size();

        vector<int> vis(this->num_nodes, 0);
        this->edges_l2r = this->edges_r2l = vector<std::pair<int, int>>(this->num_nodes);

        vector<vector<CodonGraphEdge<EnergyT>>> adj_rev(this->num_nodes);

        // DFS to build reverse adjacency list and populate l2r edges
        auto dfs = [&](auto dfs, int at) {
            if (vis[at]) return;
            vis[at] = 1;
            int st = static_cast<int>(this->edges.size());
            for (const auto& e : adj[at]) {
                this->edges.push_back(e);
                adj_rev[e.to].push_back({at, this->node_base[at], e.weight});
            }
            this->edges_l2r[at] = {st, static_cast<int>(this->edges.size())};
            for (const auto& e : adj[at]) {
                dfs(dfs, e.to);
            }
        };

        for (int nid : this->nt_nodes[0]) {
            dfs(dfs, nid);
        }

        for (int i = 0; i < this->num_nodes; ++i) {
            int st = static_cast<int>(this->edges.size());
            for (const auto& e : adj_rev[i]) {
                this->edges.push_back(e);
            }
            this->edges_r2l[i] = {st, static_cast<int>(this->edges.size())};
        }
    }

    struct NodeEdges {
        const CodonGraph& graph;
        int st, en;

        auto begin() const { return graph.edges.begin() + st; }
        auto end() const { return graph.edges.begin() + en; }
        auto size() const { return en - st; }
        auto operator[](int idx) const { return graph.edges[st + idx]; }
    };

   public:
    CodonGraph(CodonGraphConfig<SemiringT> config) : config(config) { BuildStandardCodonGraph(); }

    NodeEdges EdgesL2R(int node) const {
        return {*this, this->edges_l2r[node].first, this->edges_l2r[node].second};
    }
    NodeEdges EdgesR2L(int node) const {
        return {*this, this->edges_r2l[node].first, this->edges_r2l[node].second};
    }

    int NumNodes() const { return num_nodes; }

    const vector<int>& NtNodes(int nt) const { return nt_nodes[nt]; }

    Base NodeBase(int node) const { return node_base[node]; }

    int Nts() const { return nt_nodes.size(); }

    const CodonGraphConfig<SemiringT>& Config() const { return config; }
};

}  // namespace mwmrna::fold_codons

#endif  // MWMRNA_LIB_CODON_GRAPH_HPP