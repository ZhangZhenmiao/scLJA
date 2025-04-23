#include "GFA_graph.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>

GFAGraph::Edge::Edge(std::string id, Orientation o, size_t ol)
    : node_id(std::move(id)), node_orientation(o), overlap(ol) {
}

GFAGraph::GFAGraph(const fs::path& gfa_path) {
    std::ifstream file(gfa_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open GFA file: " + gfa_path.string());
    }

    std::string line;
    while (std::getline(file, line)) {
        parse_gfa_line(line);
    }
    log_verbose("Loaded " + std::to_string(nodes.size()) + " nodes");
}

GFAGraph::Orientation GFAGraph::complement(Orientation o) {
    return (o == Orientation::FORWARD) ? Orientation::REVERSE :
        (o == Orientation::REVERSE) ? Orientation::FORWARD :
        Orientation::UNKNOWN;
}

const std::vector<GFAGraph::Edge>& GFAGraph::Node::outgoing(Orientation o) const {
    return (o == Orientation::FORWARD) ? out_forward : out_reverse;
}

const std::vector<GFAGraph::Edge>& GFAGraph::Node::incoming(Orientation o) const {
    return (o == Orientation::FORWARD) ? in_forward : in_reverse;
}

GFAGraph::Node& GFAGraph::get_node(const std::string& id) {
    auto it = nodes.find(id);
    if (it == nodes.end()) throw std::runtime_error("Node not found: " + id);
    return it->second;
}

void GFAGraph::parse_gfa_line(const std::string& line) {
    if (line.empty()) return;

    std::istringstream iss(line);
    char record_type;
    iss >> record_type;

    try {
        switch (record_type) {
        case 'S': {
            std::string id, seq;
            iss >> id >> seq;
            nodes.try_emplace(id, Node{ id, seq });
            log_verbose("Loaded segment " + id + ", len " + std::to_string(seq.size()));
            break;
        }
        case 'L': {
            std::string from, to;
            char from_orient, to_orient;
            size_t overlap = 0;
            iss >> from >> from_orient >> to >> to_orient;
            try { iss >> overlap; }
            catch (...) {}

            add_edge(
                from, static_cast<Orientation>(from_orient),
                to, static_cast<Orientation>(to_orient),
                overlap
            );
            break;
        }
        }
    }
    catch (...) {
        // Skip malformed lines
    }
}

void GFAGraph::add_edge_internal(const std::string& from_id, Orientation from_orient,
    const std::string& to_id, Orientation to_orient,
    size_t overlap) {
    Node& from_node = get_node(from_id);
    Node& to_node = get_node(to_id);

    auto& out_edges = (from_orient == Orientation::FORWARD)
        ? from_node.out_forward
        : from_node.out_reverse;
    out_edges.emplace_back(to_id, to_orient, overlap);

    auto& in_edges = (to_orient == Orientation::FORWARD)
        ? to_node.in_forward
        : to_node.in_reverse;
    in_edges.emplace_back(from_id, from_orient, overlap);

    log_verbose("Loaded edge " + from_id + static_cast<char>(from_orient) + " -> " + to_id + static_cast<char>(to_orient) + ", overlap " + std::to_string(overlap));
}

void GFAGraph::add_edge(const std::string& from_id, Orientation from_orient,
    const std::string& to_id, Orientation to_orient,
    size_t overlap) {
    add_edge_internal(from_id, from_orient, to_id, to_orient, overlap);

    //in case of LJA, add RC edge
    if (from_id.find('s') == std::string::npos) {
        add_edge_internal(to_id, complement(to_orient), from_id, complement(from_orient), overlap);
    }
}

const std::vector<GFAGraph::Edge>& GFAGraph::outgoing_edges(const std::string& node_id,
    Orientation source_orient) const {
    static const std::vector<Edge> empty;
    auto it = nodes.find(node_id);
    return (it != nodes.end()) ? it->second.outgoing(source_orient) : empty;
}

const std::vector<GFAGraph::Edge>& GFAGraph::incoming_edges(const std::string& node_id,
    Orientation target_orient) const {
    static const std::vector<Edge> empty;
    auto it = nodes.find(node_id);
    return (it != nodes.end()) ? it->second.incoming(target_orient) : empty;
}

GFAGraph::EdgeView GFAGraph::all_edges(const std::string& node_id) const {
    static const std::vector<Edge> empty;
    auto it = nodes.find(node_id);
    if (it == nodes.end()) return { empty, empty, empty, empty };
    return {
        it->second.out_forward,
        it->second.out_reverse,
        it->second.in_forward,
        it->second.in_reverse
    };
}

void GFAGraph::write_fasta(fs::path& output) {
    std::ofstream file(output);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to write fasta file: " + output.string());
    }

    for (auto&& n : nodes) {
        file << ">" << n.second.id << "\n" << n.second.sequence << "\n";
    }

    file.close();
}
void GFAGraph::assign_edge_type() {
    int cnt_circular = 0, cnt_linear = 0, cnt_complex = 0;
    for (auto&& node : nodes) {
        assert(node.second.in_forward.size() == node.second.out_reverse.size());
        assert(node.second.in_reverse.size() == node.second.out_forward.size());
        if (node.second.in_forward.size() == 0 && node.second.out_forward.size() == 0) {
            node.second.status = EdgeStats::LINEAR;
            cnt_linear++;
        }
        else if (node.second.in_forward.size() == 1 && node.second.out_forward.size() == 1 && has_node_in_edges(node.second.in_forward, node.first, Orientation::FORWARD)) {
            node.second.status = EdgeStats::CIRCULAR;
            cnt_circular++;
        }
        else {
            node.second.status = EdgeStats::COMPLEX;
            cnt_complex++;
        }
    }

    log_verbose("Linear contigs: " + std::to_string(cnt_linear) + " , circular contigs: " + std::to_string(cnt_circular) + ", contigs in complex components: " + std::to_string(cnt_complex));
}