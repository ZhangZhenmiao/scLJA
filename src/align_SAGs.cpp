#include <string>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <omp.h>
#include <zlib.h>
#include "align_SAGs.hpp"
#include "utils.hpp"

#define BUFFER_SIZE 4096

namespace fs = std::filesystem;

void SAGAligner::sort_fastq_by_name(const fs::path& input_path, const fs::path& output_path) {
    struct FastqRecord {
        std::string header;
        std::string sequence;
        std::string quality;
        std::string separator; // The '+' line

        // For sorting
        bool operator<(const FastqRecord& other) const {
            return header < other.header;
        }
    };

    // Check input file exists
    if (!fs::exists(input_path)) {
        throw std::runtime_error("Input file does not exist: " + input_path.string());
    }

    // Create output directory if needed
    fs::create_directories(output_path.parent_path());

    // Read all records
    std::vector<FastqRecord> records;

    gzFile file_in = gzopen(input_path.c_str(), "rb");
    if (!file_in) {
        throw std::runtime_error("Error opening file: " + input_path.string());
    }

    char buffer[BUFFER_SIZE];

    while (true) {
        FastqRecord record;

        if (!gzgets(file_in, buffer, BUFFER_SIZE)) break; // Header
        record.header = buffer;

        if (!gzgets(file_in, buffer, BUFFER_SIZE)) break; // Sequence
        record.sequence = buffer;

        if (!gzgets(file_in, buffer, BUFFER_SIZE)) break; // Separator ('+')
        record.separator = buffer;

        if (!gzgets(file_in, buffer, BUFFER_SIZE)) break; // Quality score
        record.quality = buffer;

        // Remove trailing newline characters
        record.header.erase(record.header.find_last_not_of("\r\n") + 1);
        record.sequence.erase(record.sequence.find_last_not_of("\r\n") + 1);
        record.separator.erase(record.separator.find_last_not_of("\r\n") + 1);
        record.quality.erase(record.quality.find_last_not_of("\r\n") + 1);

        records.push_back(std::move(record));
    }

    gzclose(file_in);

    std::sort(records.begin(), records.end());

    // Write sorted output
    gzFile file_out = gzopen(output_path.c_str(), "wb");
    if (!file_out) {
        throw std::runtime_error("Error opening file for writing: " + output_path.string());
    }

    for (size_t i = 0; i < records.size(); ++i) {
        gzprintf(file_out, "%s\n%s\n%s\n%s\n",
            records[i].header.c_str(),
            records[i].sequence.c_str(),
            records[i].separator.c_str(),
            records[i].quality.c_str());
    }

    gzclose(file_out);
}

void SAGAligner::execute_bwa(const std::string& ref, const std::string& reads1, const std::string& reads2, const std::string& output) {
    if (fs::exists(output)) {
        std::cout << "Skip bwa for " << output << ": the bam already exists." << std::endl;
        return;
    }
    std::string cmd = "bwa mem " + ref + " " + reads1 + " " + reads2 + " -t " + std::to_string(thread) + " | samtools sort -@ " + std::to_string(thread) + " -o " + output;
    std::cout << "Run BWA: " + cmd << std::endl;
    if (execute_command(cmd, output + ".log") != 0) {
        throw std::runtime_error("Run BWA: " + cmd + " failed.");
    }
}

std::string SAGAligner::get_barcode(const std::string& path) {
    size_t last_dash = path.find_last_of("-");
    if (last_dash == std::string::npos)
        return "";
    return path.substr(last_dash + 1, path.find("_", last_dash + 1) - last_dash - 1);
}

SAGAligner::SAGAligner(std::string assembly_path, std::string SAGs_dir, std::string output_dir, int thread) {
    thread = thread;
    omp_set_num_threads(thread);
    fs::path dir = output_dir;
    fs::create_directory(dir);

    fs::path target_asm = dir / "contigs.fasta";
    if (!fs::exists(target_asm))
        fs::copy_file(assembly_path, target_asm);

    if (!fs::exists(target_asm.string() + ".bwt")) {
        std::cout << "Creating bwa index for assembly ..." << std::endl;
        if (execute_command("bwa index " + target_asm.string()) != 0) {
            throw std::runtime_error("Run BWA index: bwa index " + target_asm.string() + " failed.");
        }
    }

    std::vector <std::string> SAG_reads;
    for (const auto& entry : fs::directory_iterator(SAGs_dir)) {
        if (entry.is_regular_file() && entry.path().extension() == ".gz") {
            SAG_reads.push_back(fs::absolute(entry.path()));
        }
    }

    std::sort(SAG_reads.begin(), SAG_reads.end());
#pragma omp parallel for
    for (size_t i = 0; i < SAG_reads.size() / 2; ++i) {
        std::string barcode1 = get_barcode(SAG_reads[2 * i]);
        std::string barcode2 = get_barcode(SAG_reads[2 * i + 1]);
        assert(barcode1 == barcode2);
        sort_fastq_by_name(SAG_reads[2 * i], dir / (barcode1 + "_R1.fq.gz"));
        sort_fastq_by_name(SAG_reads[2 * i + 1], dir / (barcode2 + "_R2.fq.gz"));
        fs::path output_bam = dir / (barcode1 + ".bam");
        execute_bwa(target_asm, dir / (barcode1 + "_R1.fq.gz"), dir / (barcode2 + "_R2.fq.gz"), output_bam);
    }
}