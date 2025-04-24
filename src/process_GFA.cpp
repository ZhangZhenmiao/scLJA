#include "process_GFA.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cassert>

using namespace gfa;

std::unordered_map <std::string, std::vector<EdgeInNode>>& Node::getOutEdges(Orientation& from_orient, Orientation& to_orient) {
    if (from_orient == Orientation::FORWARD && to_orient == Orientation::FORWARD) {
        return outgoing_plus_plus;
    }
    else if (from_orient == Orientation::FORWARD && to_orient == Orientation::REVERSE) {
        return outgoing_plus_minus;
    }
    else if (from_orient == Orientation::REVERSE && to_orient == Orientation::FORWARD) {
        return outgoing_minus_plus;
    }
    else {
        return outgoing_minus_minus;
    }
}

std::unordered_map <std::string, std::vector<EdgeInNode>>& Node::getInEdges(Orientation& from_orient, Orientation& to_orient) {
    if (from_orient == Orientation::FORWARD && to_orient == Orientation::FORWARD) {
        return incoming_plus_plus;
    }
    else if (from_orient == Orientation::FORWARD && to_orient == Orientation::REVERSE) {
        return incoming_plus_minus;
    }
    else if (from_orient == Orientation::REVERSE && to_orient == Orientation::FORWARD) {
        return incoming_minus_plus;
    }
    else {
        return incoming_minus_minus;
    }
}

bool Node::nodeInOutEdge(Orientation from_o, Orientation to_o, std::string to_id) {
    auto& edges = getOutEdges(from_o, to_o);
    return edges.find(to_id) != edges.end();
}

void GFAGraph::addNode(const std::string& id, const std::string& sequence) {
    if (nodes.find(id) != nodes.end()) {
        throw std::runtime_error("Node with id " + id + " already exists");
    }

    auto node = std::make_shared<Node>();
    node->id = id;
    node->sequence = sequence;
    nodes[id] = node;
}

std::shared_ptr<Edge> GFAGraph::findEdge(
    const std::string& from_id, Orientation from_orient,
    const std::string& to_id, Orientation to_orient) const {

    for (const auto& [id, edge] : edges) {
        if (edge->from_id == from_id && edge->from_orient == from_orient &&
            edge->to_id == to_id && edge->to_orient == to_orient) {
            return edge;
        }
    }
    return nullptr;
}

void GFAGraph::addEdge(
    const std::string& from_id, Orientation from_orient,
    const std::string& to_id, Orientation to_orient,
    const std::string& overlap) {

    auto from_node = getNode(from_id);
    auto to_node = getNode(to_id);

    if (!from_node || !to_node) {
        throw std::runtime_error("One or both nodes not found for edge creation");
    }

    // Check if edge already exists (in either direction)
    auto existing_edge = findEdge(from_id, from_orient, to_id, to_orient);
    if (existing_edge) {
        // Edge exists, just add connections to nodes
        bool is_forward = (existing_edge->from_id == from_id &&
            existing_edge->from_orient == from_orient);

        from_node->getOutEdges(from_orient, to_orient)[to_id].emplace_back(EdgeInNode(existing_edge, is_forward));
        to_node->getInEdges(to_orient, from_orient)[from_id].emplace_back(EdgeInNode(existing_edge, is_forward));
        assert(from_node->getOutEdges(from_orient, to_orient)[to_id].size() == 1);
        return;
    }

    // Create new edge
    auto edge = std::make_shared<Edge>();
    edge->id = next_edge_id++;
    edge->from_id = from_id;
    edge->to_id = to_id;
    edge->from_orient = from_orient;
    edge->to_orient = to_orient;
    edge->overlap = overlap;

    // Add to edges list
    edges[edge->id] = edge;

    from_node->getOutEdges(from_orient, to_orient)[to_id].emplace_back(EdgeInNode(edge, true));
    to_node->getInEdges(to_orient, from_orient)[from_id].emplace_back(EdgeInNode(edge, true));
    assert(to_node->getInEdges(to_orient, from_orient)[from_id].size() == 1);
    log_verbose("Loaded edge " + from_id + static_cast<char>(from_orient) + " -> " + to_id + static_cast<char>(to_orient) + ", overlap " + overlap);
}

std::shared_ptr<Node> GFAGraph::getNode(const std::string& id) const {
    auto it = nodes.find(id);
    if (it != nodes.end()) {
        return it->second;
    }
    return nullptr;
}

std::unordered_map<std::string, std::shared_ptr<Node>>& GFAGraph::getNodes() {
    return nodes;
}

const std::unordered_map<size_t, std::shared_ptr<Edge>>& GFAGraph::getEdges() const {
    return edges;
}

void GFAGraph::parseGFA(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    while (getline(file, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') {
            continue;
        }

        // Split line into tokens
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;
        while (iss >> token) {
            tokens.push_back(token);
        }

        if (tokens.empty()) {
            continue;
        }

        // Handle different record types
        if (tokens[0] == "S") {  // Segment line
            parseSegmentLine(line);
        }
        else if (tokens[0] == "L") {  // Link line
            parseLinkLine(line);
        }
    }
    file.close();
    log_verbose("Loaded " + std::to_string(nodes.size()) + " nodes, " + std::to_string(edges.size()) + " edges");
}

void GFAGraph::parseSegmentLine(const std::string& line) {
    std::istringstream iss(line);
    std::string record_type, id, sequence;

    if (!(iss >> record_type >> id >> sequence)) {
        throw std::runtime_error("Invalid segment line: " + line);
    }

    addNode(id, sequence);
}

void GFAGraph::parseLinkLine(const std::string& line) {
    std::istringstream iss(line);
    std::string record_type, from_id, to_id;
    char from_orient_char, to_orient_char;
    std::string overlap;

    if (!(iss >> record_type >> from_id >> from_orient_char >> to_id >> to_orient_char)) {
        throw std::runtime_error("Invalid link line: " + line);
    }

    iss >> overlap;

    auto from_ori = charToOrientation(from_orient_char), to_ori = charToOrientation(to_orient_char);

    addEdge(from_id, from_ori, to_id, to_ori, overlap);
    if (from_id.find('s') == std::string::npos) {
        addEdge(to_id, complement(to_ori), from_id, complement(from_ori), overlap);
    }
}

Orientation GFAGraph::charToOrientation(char c) const {
    c = toupper(c);
    if (c == '+') {
        return Orientation::FORWARD;
    }
    else if (c == '-') {
        return Orientation::REVERSE;
    }
    throw std::runtime_error("Invalid orientation character: " + std::string(1, c));
}

void GFAGraph::assign_node_type() {
    int cnt_circular = 0, cnt_linear = 0, cnt_complex = 0;
    for (auto&& node : nodes) {
        assert(node.second->getPlusOutSize() == node.second->getMinusInSize());
        assert(node.second->getMinusOutSize() == node.second->getPlusInSize());

        if (node.second->getPlusInSize() == 0 && node.second->getPlusOutSize() == 0) {
            node.second->status = Circular::LINEAR;
            cnt_linear++;
        }
        else if (node.second->getPlusInSize() == 1 && node.second->getPlusOutSize() == 1 && node.second->nodeInOutEdge(Orientation::FORWARD, Orientation::FORWARD, node.first)) {
            node.second->status = Circular::CIRCULAR;
            cnt_circular++;
        }
        else {
            node.second->status = Circular::COMPLEX;
            cnt_complex++;
        }
    }

    log_verbose("Linear contigs: " + std::to_string(cnt_linear) + " , circular contigs: " + std::to_string(cnt_circular) + ", contigs in complex components: " + std::to_string(cnt_complex));
}

void GFAGraph::write_fasta(fs::path& output) {
    std::ofstream file(output);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to write fasta file: " + output.string());
    }

    for (auto&& n : nodes) {
        file << ">" << n.first << "\n" << n.second->sequence << "\n";
    }

    file.close();
}

void GFAGraph::move_to_contracted(std::unordered_map<std::string, std::vector<EdgeInNode>>& contract_edges, std::unordered_map<std::string, std::vector<EdgeInNode>>& ori_edges, std::string& contract_name, std::string& ori_name) {
    for (auto&& e : ori_edges) {
        if (e.first == contract_name)
            continue;
        if (e.first == ori_name)
            continue;
        if (contract_edges.find(e.first) != contract_edges.end())
            contract_edges[e.first].insert(contract_edges[e.first].end(), e.second.begin(), e.second.end());
        else
            contract_edges[e.first] = e.second;
    }
}

void GFAGraph::rename_connected_node(std::unordered_map<std::string, std::vector<EdgeInNode>>& edges, std::string& contract_name, std::string& ori_name) {
    if (edges.find(ori_name) != edges.end()) {
        edges[contract_name].insert(edges[contract_name].end(), edges[ori_name].begin(), edges[ori_name].end());
        edges.erase(ori_name);
    }
}

void GFAGraph::collapse_short_edges(int short_length) {
    std::cout << "Contract nodes" << std::endl;
    // collect all nodes to be contracted
    std::vector<std::string> segments_to_contract;
    for (auto&& s : nodes) {
        if (s.second->status == Circular::CIRCULAR || s.second->status == Circular::LINEAR)
            continue;
        if (s.second->sequence.size() <= size_t(short_length))
            segments_to_contract.push_back(s.first);
    }

    // contract nodes one by one
    int contracted_id = 1;
    std::unordered_map<std::string, std::string> segment_to_contracted_name;
    for (auto&& s : segments_to_contract) {
        auto node = std::make_shared<Node>();

        if (segment_to_contracted_name.find(s) != segment_to_contracted_name.end()) {
            assert(nodes.find(segment_to_contracted_name[s]) != nodes.end());
            node = nodes[segment_to_contracted_name[s]];
        }
        else {
            node->id = "contracted_" + std::to_string(contracted_id++);
            node->sequence = "NNNNNNNNNNNNNNNNNNNN";
        }

        std::cout << "Move links from " << s << " to " << node->id << std::endl;
        move_to_contracted(node->incoming_minus_minus, nodes[s]->incoming_minus_minus, node->id, s);
        move_to_contracted(node->incoming_minus_plus, nodes[s]->incoming_minus_plus, node->id, s);
        move_to_contracted(node->incoming_plus_minus, nodes[s]->incoming_plus_minus, node->id, s);
        move_to_contracted(node->incoming_plus_plus, nodes[s]->incoming_plus_plus, node->id, s);
        move_to_contracted(node->outgoing_minus_minus, nodes[s]->outgoing_minus_minus, node->id, s);
        move_to_contracted(node->outgoing_minus_plus, nodes[s]->outgoing_minus_plus, node->id, s);
        move_to_contracted(node->outgoing_plus_minus, nodes[s]->outgoing_plus_minus, node->id, s);
        move_to_contracted(node->outgoing_plus_plus, nodes[s]->outgoing_plus_plus, node->id, s);

        node->incoming_minus_minus.erase(s);
        node->incoming_minus_plus.erase(s);
        node->incoming_plus_minus.erase(s);
        node->incoming_plus_plus.erase(s);
        node->outgoing_minus_minus.erase(s);
        node->outgoing_minus_plus.erase(s);
        node->outgoing_plus_minus.erase(s);
        node->outgoing_plus_plus.erase(s);
        for (auto&& neighbor : node->incoming_minus_minus) {
            rename_connected_node(nodes[neighbor.first]->outgoing_minus_minus, node->id, s);
            segment_to_contracted_name[neighbor.first] = node->id;
        }
        for (auto&& neighbor : node->incoming_minus_plus) {
            rename_connected_node(nodes[neighbor.first]->outgoing_plus_minus, node->id, s);
            segment_to_contracted_name[neighbor.first] = node->id;
        }
        for (auto&& neighbor : node->incoming_plus_minus) {
            rename_connected_node(nodes[neighbor.first]->outgoing_minus_plus, node->id, s);
            segment_to_contracted_name[neighbor.first] = node->id;
        }
        for (auto&& neighbor : node->incoming_plus_plus) {
            rename_connected_node(nodes[neighbor.first]->outgoing_plus_plus, node->id, s);
            segment_to_contracted_name[neighbor.first] = node->id;
        }
        for (auto&& neighbor : node->outgoing_plus_minus) {
            rename_connected_node(nodes[neighbor.first]->incoming_minus_plus, node->id, s);
            segment_to_contracted_name[neighbor.first] = node->id;
        }
        for (auto&& neighbor : node->outgoing_plus_plus) {
            rename_connected_node(nodes[neighbor.first]->incoming_plus_plus, node->id, s);
            segment_to_contracted_name[neighbor.first] = node->id;
        }
        for (auto&& neighbor : node->outgoing_minus_minus) {
            rename_connected_node(nodes[neighbor.first]->incoming_minus_minus, node->id, s);
            segment_to_contracted_name[neighbor.first] = node->id;
        }
        for (auto&& neighbor : node->outgoing_minus_plus) {
            rename_connected_node(nodes[neighbor.first]->incoming_plus_minus, node->id, s);
            segment_to_contracted_name[neighbor.first] = node->id;
        }

        nodes.erase(s);
        nodes[node->id] = node;
    }
    std::cout << "Contract nodes further" << std::endl;
    std::vector<std::string> contracted_to_contract;
    for (auto&& s : nodes) {
        if (s.first.find("contracted") == std::string::npos)
            continue;
        bool flag = false;
        for (auto&& e : s.second->outgoing_minus_minus) {
            if (e.first.find("contracted") != std::string::npos)
                flag = true;
        }
        for (auto&& e : s.second->outgoing_minus_plus) {
            if (e.first.find("contracted") != std::string::npos)
                flag = true;
        }
        for (auto&& e : s.second->outgoing_plus_minus) {
            if (e.first.find("contracted") != std::string::npos)
                flag = true;
        }
        for (auto&& e : s.second->outgoing_plus_plus) {
            if (e.first.find("contracted") != std::string::npos)
                flag = true;
        }
        if (flag)
            contracted_to_contract.push_back(s.first);
    }

    for (auto&& s : contracted_to_contract) {
        std::string id;
        for (auto&& e : nodes[s]->outgoing_minus_minus) {
            if (e.first.find("contracted") != std::string::npos)
                id = e.first;
        }
        for (auto&& e : nodes[s]->outgoing_minus_plus) {
            if (e.first.find("contracted") != std::string::npos)
                id = e.first;
        }
        for (auto&& e : nodes[s]->outgoing_plus_minus) {
            if (e.first.find("contracted") != std::string::npos)
                id = e.first;
        }
        for (auto&& e : nodes[s]->outgoing_plus_plus) {
            if (e.first.find("contracted") != std::string::npos)
                id = e.first;
        }

        if (id.empty())
            continue;
        assert(nodes.find(id) != nodes.end());
        auto node = nodes[id];
        assert(id != s);
        std::cout << "Move links from " << s << " to " << node->id << std::endl;
        move_to_contracted(node->incoming_minus_minus, nodes[s]->incoming_minus_minus, node->id, s);
        move_to_contracted(node->incoming_minus_plus, nodes[s]->incoming_minus_plus, node->id, s);
        move_to_contracted(node->incoming_plus_minus, nodes[s]->incoming_plus_minus, node->id, s);
        move_to_contracted(node->incoming_plus_plus, nodes[s]->incoming_plus_plus, node->id, s);
        move_to_contracted(node->outgoing_minus_minus, nodes[s]->outgoing_minus_minus, node->id, s);
        move_to_contracted(node->outgoing_minus_plus, nodes[s]->outgoing_minus_plus, node->id, s);
        move_to_contracted(node->outgoing_plus_minus, nodes[s]->outgoing_plus_minus, node->id, s);
        move_to_contracted(node->outgoing_plus_plus, nodes[s]->outgoing_plus_plus, node->id, s);

        node->incoming_minus_minus.erase(s);
        node->incoming_minus_plus.erase(s);
        node->incoming_plus_minus.erase(s);
        node->incoming_plus_plus.erase(s);
        node->outgoing_minus_minus.erase(s);
        node->outgoing_minus_plus.erase(s);
        node->outgoing_plus_minus.erase(s);
        node->outgoing_plus_plus.erase(s);
        for (auto&& neighbor : node->incoming_minus_minus) {
            rename_connected_node(nodes[neighbor.first]->outgoing_minus_minus, node->id, s);
            segment_to_contracted_name[neighbor.first] = node->id;
        }
        for (auto&& neighbor : node->incoming_minus_plus) {
            rename_connected_node(nodes[neighbor.first]->outgoing_plus_minus, node->id, s);
            segment_to_contracted_name[neighbor.first] = node->id;
        }
        for (auto&& neighbor : node->incoming_plus_minus) {
            rename_connected_node(nodes[neighbor.first]->outgoing_minus_plus, node->id, s);
            segment_to_contracted_name[neighbor.first] = node->id;
        }
        for (auto&& neighbor : node->incoming_plus_plus) {
            rename_connected_node(nodes[neighbor.first]->outgoing_plus_plus, node->id, s);
            segment_to_contracted_name[neighbor.first] = node->id;
        }
        for (auto&& neighbor : node->outgoing_plus_minus) {
            rename_connected_node(nodes[neighbor.first]->incoming_minus_plus, node->id, s);
            segment_to_contracted_name[neighbor.first] = node->id;
        }
        for (auto&& neighbor : node->outgoing_plus_plus) {
            rename_connected_node(nodes[neighbor.first]->incoming_plus_plus, node->id, s);
            segment_to_contracted_name[neighbor.first] = node->id;
        }
        for (auto&& neighbor : node->outgoing_minus_minus) {
            rename_connected_node(nodes[neighbor.first]->incoming_minus_minus, node->id, s);
            segment_to_contracted_name[neighbor.first] = node->id;
        }
        for (auto&& neighbor : node->outgoing_minus_plus) {
            rename_connected_node(nodes[neighbor.first]->incoming_plus_minus, node->id, s);
            segment_to_contracted_name[neighbor.first] = node->id;
        }

        nodes.erase(s);
        nodes[node->id] = node;
    }
    std::cout << "Contract nodes done" << std::endl;
}

void GFAGraph::write_gfa(fs::path output) {
    std::ofstream file(output);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to write fasta file: " + output.string());
    }
    file << "H\tVN:Z:1.0\n";

    for (auto&& s : nodes) {
        file << "S\t" << s.first << "\t" << s.second->sequence << "\n";
    }

    for (auto&& s : nodes) {
        for (auto&& n : s.second->outgoing_minus_minus) {
            for (auto&& e : n.second) {
                file << "L\t" << s.first << "\t-\t" << n.first << "\t-\t" << e.edge->overlap << "\n";
            }
        }
        for (auto&& n : s.second->outgoing_minus_plus) {
            for (auto&& e : n.second) {
                file << "L\t" << s.first << "\t-\t" << n.first << "\t+\t" << e.edge->overlap << "\n";
            }
        }
        for (auto&& n : s.second->outgoing_plus_plus) {
            for (auto&& e : n.second) {
                file << "L\t" << s.first << "\t+\t" << n.first << "\t+\t" << e.edge->overlap << "\n";
            }
        }
        for (auto&& n : s.second->outgoing_plus_minus) {
            for (auto&& e : n.second) {
                file << "L\t" << s.first << "\t+\t" << n.first << "\t-\t" << e.edge->overlap << "\n";
            }
        }
    }

    file.close();
}