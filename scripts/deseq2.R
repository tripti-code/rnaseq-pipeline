library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)
counts_dir <- args[1]
output_dir <- args[2]
metadata_file <- args[3]

sample_info <- read.csv(metadata_file, row.names = 1)

count_files <- list.files(counts_dir, pattern = "\\.txt$", full.names = TRUE)

count_list <- list()

for (file in count_files) {
    df <- read.table(file, skip=2, header=TRUE)
    counts <- df[, ncol(df)]
    names(counts) <- df[, 1]
    sample_name <- gsub("_Aligned.sortedByCoord.out.txt", "", basename(file))
    count_list[[sample_name]] <- counts
    }

counts_matrix <- do.call(cbind, count_list)

dds <- DESeqDataSetFromMatrix(
    countData = counts_matrix,
    colData = sample_info,
    design = ~ condition
)

dds <- DESeq(dds)

results <- results(dds)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(results, file.path(output_dir, "results.csv"))

write.csv(counts(dds, normalized=TRUE), file.path(output_dir, "normalized_counts.csv"))

print("DESeq2 analysis is completed.")