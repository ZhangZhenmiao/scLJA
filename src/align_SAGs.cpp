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
#include <htslib/sam.h>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include "process_GFA.hpp"

#define BUFFER_SIZE 4096

namespace fs = std::filesystem;

void SAGAligner::sortFastqByName(const fs::path& input_path, const fs::path& output_path) {
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

void SAGAligner::executeBwa(const std::string& ref, const std::string& reads1, const std::string& reads2, const std::string& output, bool remove_duplicate) {
    if (fs::exists(output)) {
#pragma omp critical
        {
            log_verbose("Skip bwa for " + output + ": the bam already exists");
        }
        return;
    }

    std::string bc = reads1.substr(reads1.find_last_of('/') + 1, reads1.find_last_of('_') - reads1.find_last_of('/') - 1);
    if (remove_duplicate) {
        std::string cmd = "bwa mem " + ref + " " + reads1 + " " + reads2 + " -R \"@RG\\tID:" + bc + "\\tSM:" + bc + "\\tPL:ILLUMINA\"" + " -t " + std::to_string(thread) + " | samtools sort -@ " + std::to_string(thread) + " -o " + output + ".tmp.bam";
#pragma omp critical
        {
            log_verbose("Run BWA for " + reads1 + " and " + reads2 + " ...");
        }
        if (execute_command(cmd, output + ".tmp.bam.log") != 0) {
            throw std::runtime_error("Run BWA: " + cmd + " failed.");
        }

        cmd = "picard MarkDuplicates I=" + output + ".tmp.bam O=" + output + " M=" + output + ".dedup_metrics.txt REMOVE_DUPLICATES=true";
        if (execute_command(cmd, output + ".log") != 0) {
            throw std::runtime_error("Run Picard: " + cmd + " failed.");
        }

        if (execute_command("rm " + output + ".tmp.bam") != 0) {
            throw std::runtime_error("Run rm " + output + ".tmp.bam failed.");
        }

        if (execute_command("rm " + output + ".tmp.bam.log") != 0) {
            throw std::runtime_error("Run rm " + output + ".tmp.bam.log failed.");
        }
    }
    else {
        std::string cmd = "bwa mem " + ref + " " + reads1 + " " + reads2 + " -t " + std::to_string(thread) + " | samtools sort -@ " + std::to_string(thread) + " -o " + output;
#pragma omp critical
        {
            log_verbose("Run BWA for " + reads1 + " and " + reads2 + " ...");
        }
        if (execute_command(cmd, output + ".log") != 0) {
            throw std::runtime_error("Run BWA: " + cmd + " failed.");
        }
    }
}

std::string SAGAligner::getBarcode(const std::string& path) {
    size_t last_dash = path.find_last_of("-");
    if (last_dash == std::string::npos)
        return "";
    return path.substr(last_dash + 1, path.find("_", last_dash + 1) - last_dash - 1);
}

SAGAligner::SAGAligner(const fs::path& assembly_path, const std::string& SAGs_dir, const std::string& output_dir, int thread) {
    thread = thread;
    omp_set_num_threads(thread);
    fs::path dir = output_dir;
    fs::create_directory(dir);

    if (!fs::exists(assembly_path.string() + ".bwt")) {
        log_verbose("Creating bwa index for assembly ...");
        if (execute_command("bwa index " + assembly_path.string()) != 0) {
            throw std::runtime_error("Run BWA index: bwa index " + assembly_path.string() + " failed.");
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
        std::string barcode1 = getBarcode(SAG_reads[2 * i]);
        std::string barcode2 = getBarcode(SAG_reads[2 * i + 1]);
        assert(barcode1 == barcode2);
        if (!fs::exists(dir / (barcode1 + "_R1.fq.gz")) || !fs::exists(dir / (barcode1 + "_R2.fq.gz"))) {
            sortFastqByName(SAG_reads[2 * i], dir / (barcode1 + "_R1.fq.gz"));
            sortFastqByName(SAG_reads[2 * i + 1], dir / (barcode2 + "_R2.fq.gz"));
        }
        fs::path output_bam = dir / (barcode1 + ".bam");
        executeBwa(assembly_path, dir / (barcode1 + "_R1.fq.gz"), dir / (barcode2 + "_R2.fq.gz"), output_bam);
#pragma omp critical
        {
            bam_files[barcode1] = output_bam;
        }
    }
}


std::map<std::string, int> SAGAligner::getReadCounts(const fs::path& bam_path, const fs::path& output_file, gfa::GFAGraph& graph, float min_identity, int min_mapq) {
    // Open BAM file
    samFile* bam_file = sam_open(bam_path.c_str(), "r");
    if (!bam_file) {
        throw std::runtime_error("Failed to open BAM file: " + std::string(bam_path));
    }

    // Load the header
    bam_hdr_t* header = sam_hdr_read(bam_file);
    if (!header) {
        sam_close(bam_file);
        throw std::runtime_error("Failed to read BAM header");
    }

    // Initialize BAM record
    bam1_t* read = bam_init1();

    // Data structures
    std::map<std::string, int> high_confidence_counts;
    std::unordered_map<std::string, std::string> read_to_contig;
    std::unordered_set<std::string> ambiguous_reads;

    // Initialize contig counts
    for (int i = 0; i < header->n_targets; ++i) {
        high_confidence_counts[header->target_name[i]] = 0;
    }

    // Process each read
    while (sam_read1(bam_file, header, read) >= 0) {
        // Skip unmapped, secondary, and supplementary alignments
        if ((read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) ||
            read->core.tid < 0) {
            continue;
        }

        // Check mapping quality
        if (read->core.qual < min_mapq) {
            continue;
        }

        std::string read_name = bam_get_qname(read);
        std::string contig_name = header->target_name[read->core.tid];

        // Calculate alignment identity from CIGAR
        uint32_t* cigar = bam_get_cigar(read);
        int matches = 0;
        int total_bases = read->core.l_qseq;

        for (uint32_t k = 0; k < read->core.n_cigar; ++k) {
            int op = bam_cigar_op(cigar[k]);
            int len = bam_cigar_oplen(cigar[k]);

            if (op == BAM_CMATCH || op == BAM_CEQUAL) {
                matches += len;
            }
        }

        // Check identity threshold
        if (total_bases == 0 || static_cast<float>(matches) / total_bases < min_identity) {
            continue;
        }

        // Check for ambiguous mappings
        auto it = read_to_contig.find(read_name);
        if (it != read_to_contig.end()) {
            if (it->second != contig_name) {
                ambiguous_reads.insert(read_name);
            }
        }
        else {
            read_to_contig[read_name] = contig_name;
        }
    }

    // Count high-confidence reads (non-ambiguous)
    for (const auto& [read_name, contig_name] : read_to_contig) {
        if (ambiguous_reads.count(read_name) == 0) {
            high_confidence_counts[contig_name]++;
        }
    }

    // Write results to file if requested
    std::ofstream out(output_file);
    if (!out) {
        throw std::runtime_error("Failed to open output file: " + std::string(output_file));
    }

    int cnt_0_1 = 0;
    out << "contig_name\tlen\tstatus\tread_count\tdensity\n";
    for (const auto& [contig, count] : high_confidence_counts) {
        auto status = graph.getNode(contig)->status == gfa::Circular::CIRCULAR ? "circular" : (graph.getNode(contig)->status == gfa::Circular::LINEAR ? "linear" : "complex");
        out << contig << "\t" << graph.getNode(contig)->sequence.size() << "\t" << status << "\t" << count << "\t" << std::fixed << static_cast<float>(count) * 150 / graph.getNode(contig)->sequence.size() << "\n";
        if (static_cast<float>(count) * 150 / graph.getNode(contig)->sequence.size() > 0.1)
            cnt_0_1++;
    }

#pragma omp critical
    {
        log_verbose("Density > 0.1 for " + bam_path.string() + ": " + std::to_string(cnt_0_1));
    }

    // Cleanup
    bam_destroy1(read);
    bam_hdr_destroy(header);
    sam_close(bam_file);

    return high_confidence_counts;
}

void SAGAligner::countReads(gfa::GFAGraph& graph) {
#pragma omp parallel for
    for (size_t i = 0; i < bam_files.size(); ++i) {
        auto bam_it = bam_files.begin();
        advance(bam_it, i);
#pragma omp critical
        {
            log_verbose("Process bam " + bam_it->second.string());
        }
        std::map<std::string, int> contig_to_reads = getReadCounts(bam_it->second, bam_it->second.string() + ".stats", graph);
        sag_contig_reads[bam_it->first] = std::move(contig_to_reads);
    }
}