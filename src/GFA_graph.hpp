#ifndef GFA_GRAPH_HPP
#define GFA_GRAPH_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <filesystem>
#include <stdexcept>
#include <iostream>
#include <algorithm>

namespace fs = std::filesystem;

class GFAGraph {
public:
    enum class Orientation : char {
        FORWARD = '+',
        REVERSE = '-',
        UNKNOWN = '?'
    };

    enum class EdgeStats : char {
        CIRCULAR = 'c',
        LINEAR = 'l',
        COMPLEX = '-'
    };

    struct Edge {
        std::string node_id;
        Orientation node_orientation;
        size_t overlap;

        Edge(std::string id, Orientation o, size_t ol = 0);
    };

    struct EdgeView {
        const std::vector<Edge>& out_forward;
        const std::vector<Edge>& out_reverse;
        const std::vector<Edge>& in_forward;
        const std::vector<Edge>& in_reverse;
    };

    struct Node {
        std::string id;
        std::string sequence;
        EdgeStats status = EdgeStats::COMPLEX;

        std::vector<Edge> out_forward;
        std::vector<Edge> out_reverse;
        std::vector<Edge> in_forward;
        std::vector<Edge> in_reverse;

        const std::vector<Edge>& outgoing(Orientation o) const;
        const std::vector<Edge>& incoming(Orientation o) const;

        Node(const std::string& id, const std::string& seq) :id(id), sequence(seq) {};
    };

    GFAGraph() = default;
    explicit GFAGraph(const fs::path& gfa_path);

    void add_edge(const std::string& from_id, Orientation from_orient,
        const std::string& to_id, Orientation to_orient,
        size_t overlap = 0);

    const std::vector<Edge>& outgoing_edges(const std::string& node_id,
        Orientation source_orient) const;
    const std::vector<Edge>& incoming_edges(const std::string& node_id,
        Orientation target_orient) const;
    EdgeView all_edges(const std::string& node_id) const;

    void write_fasta(fs::path& output);

    // assign edge as linear, circular, or in complex component
    void assign_edge_type();

    inline int get_node_size() { return nodes.size(); }
    Node& get_node(const std::string& id);

    static bool has_node_in_edges(const std::vector<Edge>& edges,
        const std::string& node_id,
        Orientation orientation = Orientation::UNKNOWN) {
        return std::any_of(edges.begin(), edges.end(),
            [&](const Edge& e) {
                return e.node_id == node_id &&
                    (orientation == Orientation::UNKNOWN ||
                        e.node_orientation == orientation);
            });
    }

    // Non-const version
    static bool has_node_in_edges(std::vector<Edge>& edges,
        const std::string& node_id,
        Orientation orientation = Orientation::UNKNOWN) {
        return has_node_in_edges(const_cast<const std::vector<Edge>&>(edges),
            node_id, orientation);
    }

    inline static void set_verbose(bool verbose) { verbose_ = verbose; }
    inline static bool get_verbose() { return verbose_; }

private:

    std::unordered_map<std::string, Node> nodes;

    void parse_gfa_line(const std::string& line);
    void add_edge_internal(const std::string& from_id, Orientation from_orient,
        const std::string& to_id, Orientation to_orient,
        size_t overlap);

    static Orientation complement(Orientation o);

    static inline bool verbose_ = false;

    void log_verbose(const std::string& message) const {
        if (verbose_) {
            std::cout << "[GFA] " << message << '\n';
        }
    }
};

#endif // GFA_GRAPH_HPP