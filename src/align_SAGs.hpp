#include <string>

class SAGAligner {
public:
    std::string assembly_path, SAGs_dir, output_dir;
    int thread = 100;

    SAGAligner(std::string assembly_path, std::string SAGs_dir, std::string output_dir, int thread = 100);

private:
    void execute_bwa(const std::string& ref, const std::string& reads1, const std::string& reads2, const std::string& output);
    std::string get_barcode(const std::string& path);
    void sort_fastq_by_name(const std::filesystem::path& input_path, const std::filesystem::path& output_path);
};