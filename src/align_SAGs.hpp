#include <string>
#include <filesystem>
#include <map>
#include <unordered_map>
#include "process_GFA.hpp"

namespace fs = std::filesystem;

class SAGAligner {
public:
    std::string assembly_path, SAGs_dir, output_dir;
    std::map<std::string, fs::path> bam_files;
    std::unordered_map<std::string, std::map<std::string, int>>sag_contig_reads;
    int thread = 100;

    SAGAligner(const fs::path& assembly_path, const std::string& SAGs_dir, const std::string& output_dir, int thread = 100);

    std::map<std::string, int> getReadCounts(const fs::path& bam_path, const fs::path& output_file, gfa::GFAGraph& graph, float min_identity = 0.99f, int min_mapq = 40);
    void countReads(gfa::GFAGraph& graph);

    inline static void set_verbose(bool verbose) { verbose_ = verbose; }
    inline static bool get_verbose() { return verbose_; }

private:
    void executeBwa(const std::string& ref, const std::string& reads1, const std::string& reads2, const std::string& output, bool remove_duplicate = false);
    std::string getBarcode(const std::string& path);
    void sortFastqByName(const std::filesystem::path& input_path, const std::filesystem::path& output_path);

    static inline bool verbose_ = false;
    void log_verbose(const std::string& message) const {
        if (verbose_) {
            std::cout << "[SAG] " << message << '\n';
        }
    }
};