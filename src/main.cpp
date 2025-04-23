#include <iostream>
#include <string>
#include <filesystem>
#include <map>
#include <algorithm>
#include "align_SAGs.hpp"
#include "GFA_graph.hpp"

namespace fs = std::filesystem;

struct ProgramOptions {
    fs::path gfa_path;
    fs::path sag_path;
    fs::path output_path;
    bool verbose = false;
    bool help_requested = false;
};

void print_help(const char* program_name) {
    std::cout << "Usage: " << program_name << " [OPTIONS]\n\n"
        << "Required options:\n"
        << "  --assembly <path>    Path to assembly file (FASTA format)\n"
        << "  --sag <path>         Path to SAG dir (FASTQ format)\n"
        << "  --output <path>      Path to output dir\n\n"
        << "Optional options:\n"
        << "  --verbose            Enable verbose output\n"
        << "  --help               Show this help message\n\n"
        << "Example:\n"
        << "  " << program_name << " --assembly genome.fasta --sag reads.fastq --verbose\n";
}

ProgramOptions parse_arguments(int argc, char* argv[]) {
    ProgramOptions options;
    std::map<std::string, fs::path*> arg_map = {
        {"--assembly", &options.gfa_path},
        {"--sag", &options.sag_path},
        {"--output", &options.output_path}
    };

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--help") {
            options.help_requested = true;
            return options;
        }
        else if (arg == "--verbose") {
            options.verbose = true;
        }
        else if (arg_map.find(arg) != arg_map.end()) {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value for option: " + arg);
            }
            *arg_map[arg] = argv[++i];
        }
        else {
            throw std::runtime_error("Unknown option: " + arg);
        }
    }

    // Validate required options
    if (options.gfa_path.empty() || options.sag_path.empty() || options.output_path.empty()) {
        throw std::runtime_error("--assembly, --sag, and --output options are required");
    }

    return options;
}

void validate_paths(const ProgramOptions& options) {
    auto check_path = [](const fs::path& p, const std::string& name) {
        if (!fs::exists(p)) {
            throw std::runtime_error(name + " path does not exist: " + p.string());
        }
        };

    check_path(options.gfa_path, "Assembly");
    check_path(options.sag_path, "SAG");
}

int main(int argc, char* argv[]) {
    ProgramOptions options = parse_arguments(argc, argv);

    if (options.help_requested || argc == 1) {
        print_help(argv[0]);
        return 0;
    }

    validate_paths(options);

    if (options.verbose) {
        std::cout << "=== Processing Parameters ===\n"
            << "GFA file: " << fs::absolute(options.gfa_path) << "\n"
            << "SAG dir:     " << fs::absolute(options.sag_path) << "\n"
            << "Output dir:   " << fs::absolute(options.output_path) << "\n"
            << "Verbose mode: " << (options.verbose ? "ON" : "OFF") << "\n\n";
    }

    if (options.verbose) {
        std::cout << "Processing completed successfully\n";
    }

    // ==============================================
    // Program start
    // ==============================================

    if (options.verbose) {
        std::cout << "=== Program start ===\n";
    }
    GFAGraph::set_verbose(options.verbose);
    GFAGraph graph(options.gfa_path);

    graph.assign_edge_type();

    fs::create_directory(options.output_path);
    fs::path assembly_path = options.output_path / "contigs.fasta";
    if (!fs::exists(assembly_path))
        graph.write_fasta(assembly_path);

    SAGAligner::set_verbose(options.verbose);
    SAGAligner aligner(assembly_path, options.sag_path, options.output_path);
    aligner.countReads(graph);

    return 0;
}