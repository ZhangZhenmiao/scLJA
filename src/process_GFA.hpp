#ifndef GFA_READER_H
#define GFA_READER_H

#include <string>
#include <unordered_map>
#include <vector>
#include <memory>
#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;

namespace gfa {

    // Orientation of a node in an edge
    enum class Orientation {
        FORWARD = '+',  // +
        REVERSE = '-'   // -
    };

    enum class Circular {
        CIRCULAR = 'c',
        LINEAR = 'l',
        COMPLEX = '-'
    };

    // Edge structure
    struct Edge {
        size_t id;
        std::string from_id;
        std::string to_id;
        Orientation from_orient;
        Orientation to_orient;
        std::string overlap;
    };

    struct EdgeInNode {
        std::shared_ptr<Edge> edge;
        bool is_forward = true;
        EdgeInNode(std::shared_ptr<Edge>& edge, bool is_forward) : edge(edge), is_forward(is_forward) {}
    };

    // Node structure with separated plus/minus outgoing edges
    struct Node {
        std::string id;
        std::string sequence;
        Circular status;

        // Outgoing edges separated by orientation
        std::unordered_map<std::string, std::vector<EdgeInNode>> outgoing_plus_plus;
        std::unordered_map<std::string, std::vector<EdgeInNode>> outgoing_plus_minus;
        std::unordered_map<std::string, std::vector<EdgeInNode>> outgoing_minus_plus;
        std::unordered_map<std::string, std::vector<EdgeInNode>> outgoing_minus_minus;

        // Incoming edges separated by orientation
        std::unordered_map<std::string, std::vector<EdgeInNode>> incoming_plus_minus;
        std::unordered_map<std::string, std::vector<EdgeInNode>> incoming_minus_minus;
        std::unordered_map<std::string, std::vector<EdgeInNode>> incoming_plus_plus;
        std::unordered_map<std::string, std::vector<EdgeInNode>> incoming_minus_plus;

        std::unordered_map<std::string, std::vector<EdgeInNode>>& getOutEdges(Orientation& from_o, Orientation& to_o);
        std::unordered_map<std::string, std::vector<EdgeInNode>>& getInEdges(Orientation& from_o, Orientation& to_o);
        bool nodeInOutEdge(Orientation from_o, Orientation to_o, std::string to_id);
        inline size_t getPlusOutSize() { return outgoing_plus_plus.size() + outgoing_plus_minus.size(); };
        inline size_t getPlusInSize() { return incoming_plus_plus.size() + incoming_plus_minus.size(); };
        inline size_t getMinusOutSize() { return outgoing_minus_plus.size() + outgoing_minus_minus.size(); };
        inline size_t getMinusInSize() { return incoming_minus_plus.size() + incoming_minus_minus.size(); };
    };

    // Graph structure
    class GFAGraph {
    public:
        GFAGraph() = default;

        // Add a node to the graph
        void addNode(const std::string& id, const std::string& sequence);

        // Add an edge to the graph
        void addEdge(
            const std::string& from_id, Orientation from_orient,
            const std::string& to_id, Orientation to_orient,
            const std::string& overlap);

        // Get node by ID
        std::shared_ptr<Node> getNode(const std::string& id) const;

        // Get all nodes
        std::unordered_map<std::string, std::shared_ptr<Node>>& getNodes();

        // Get all edges (key: edge ID)
        const std::unordered_map<size_t, std::shared_ptr<Edge>>& getEdges() const;

        // Parse GFA file
        void parseGFA(const std::string& filename);

        void assign_node_type();

        void collapse_short_edges(int short_length = 50000);

        inline Orientation complement(Orientation o) {
            return (o == Orientation::FORWARD) ? Orientation::REVERSE : Orientation::FORWARD;
        }

        inline Orientation complement(char o) {
            auto ori = charToOrientation(o);
            return (ori == Orientation::FORWARD) ? Orientation::REVERSE : Orientation::FORWARD;
        }

        inline static void set_verbose(bool verbose) { verbose_ = verbose; }

        inline static bool get_verbose() { return verbose_; }

        GFAGraph(fs::path& filename) { parseGFA(filename); };

        void write_fasta(fs::path& output);

        void move_to_contracted(std::unordered_map<std::string, std::vector<EdgeInNode>>& contract_edges, std::unordered_map<std::string, std::vector<EdgeInNode>>& ori_edges, std::string& contract_name, std::string& ori_name);
        void rename_connected_node(std::unordered_map<std::string, std::vector<EdgeInNode>>& edges, std::string& contract_name, std::string& ori_name);

        void write_gfa(fs::path output);

    private:
        std::unordered_map<std::string, std::shared_ptr<Node>> nodes;
        std::unordered_map<size_t, std::shared_ptr<Edge>> edges;
        size_t next_edge_id = 0;
        static inline bool verbose_ = false;

        // Helper functions
        void parseSegmentLine(const std::string& line);
        void parseLinkLine(const std::string& line);
        Orientation charToOrientation(char c) const;

        // Find existing edge
        std::shared_ptr<Edge> findEdge(
            const std::string& from_id, Orientation from_orient,
            const std::string& to_id, Orientation to_orient) const;

        void log_verbose(const std::string& message) const {
            if (verbose_) {
                std::cout << "[GFA] " << message << std::endl;
            }
        }
    };
}

#endif // GFA_READER_H