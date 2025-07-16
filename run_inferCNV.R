if (!requireNamespace("infercnv", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("infercnv")
}

install.packages("rjags")
library(rjags)
library(infercnv)

# Set working directories
base_dir <- "/data/groups/gomesr/REU_2025/Warren/scGPT_files"
output_base <- "scGPT_files/infercnv_results"
gene_order_file <- file.path(base_dir, "gene_ordering_file.txt")


# Define non-epithelial reference cell types
reference_types <- c("T/NK Cells", "Endothelial", "Myeloid", "B/plasma", "Mesenchyme", "Mast")

# List GSM sample folders
samples <- list.dirs(file.path(base_dir, "infercnv_input"), full.names = TRUE, recursive = FALSE)
samples <- samples[!grepl("/\\.", samples)] # Exclude hidden directories

# Loop over each GSM sample
for (sample_dir in samples) {
  sample_name <- basename(sample_dir)
  message("Processing sample: ", sample_name)
  
  # File paths
  counts_file <- file.path(sample_dir, "counts.txt")
  annot_file  <- file.path(sample_dir, "cell_annotations.txt")
  output_dir  <- file.path(output_base, sample_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load Cell Counts Matrix
  counts_matrix <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t")
  counts_matrix <- as.matrix(counts_matrix)
  
  # Load and Review Cell Annotations
  cell_annotations <- read.table(annot_file, header = TRUE, sep="\t")
  colnames(cell_annotations) <- c("CellID", "ClusterLabel")
  
  # Identify Reference Cells
  ref_cells <- cell_annotations$CellID[cell_annotations$ClusterLabel %in% reference_types]
  
  # Downsample Reference Cells
  set.seed(34)
  if (length(ref_cells) > 1000) {
    ref_cells <- sample(ref_cells, 1000)
  }
  
  # Add ~10% of Reference Cells to Observations
  n_smoothing <- round(0.1 * length(ref_cells))
  smoothing_cells <- sample(ref_cells, n_smoothing)
  obs_cells <- cell_annotations$CellID[cell_annotations$ClusterLabel == "epithelial"]
  obs_cell <- unique(c(obs_cells, smoothing_cells)) 
  
  final_annotations <- data.frame(
    CellID = c(ref_cells, obs_cells),
    Group = c(rep("reference", length(ref_cells)), rep("observation", length(obs_cells)))
  )
  final_annot_file <- file.path(output_dir, "infercnv_annotations.txt")
  write.table(final_annotations, final_annot_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  infercnv_object <- CreateInfercnvObject(
    raw_counts_matrix = counts_matrix,
    annotations_file = final_annot_file,
    delim = "\t", 
    gene_order_file = gene_order_file, 
    ref_group_names = reference_types
  )
  
  infercnv_object <- infercnv::run(
    infercnv_object, 
    cutoff = 1,
    out_dir = output_dir, 
    cluster_by_groups = TRUE,
    HMM = TRUE,
    sliding_window = 101,
    analysis_mode = "subsample"
  )
}
