library(tidyverse)

###############################################################################
# 1. Settings
###############################################################################

genome_dir <- "/Users/Mac/Documents/INRAE_PaleoLab/05_Vertebrae/00_genomes"
output_dir <- "A273_vs_Salmo_trutta"

dir.create(output_dir, showWarnings = FALSE)

window_size <- 50
minimum_detected_per_window <- 10

copy_levels <- c("0 detected", "1 copy", "2 copies", "3 copies", "4+ copies")

copy_colours <- c(
  "0 detected" = "#D9D9D9",
  "1 copy"     = "#0072B2",
  "2 copies"   = "#E69F00",
  "3 copies"   = "#CC79A7",
  "4+ copies"  = "#D55E00"
)

###############################################################################
# 2. Read chromosome tables and gene-to-HOG assignments
###############################################################################

read_table <- function(species, filename) {
  read_delim(
    file.path(genome_dir, species, filename),
    delim = ";",
    show_col_types = FALSE
  )
}

avk_chr <- read_table("A273", "A273.chr") %>%
  arrange(ChrNum)

trout_chr <- read_table("Salmo_trutta", "Salmo_trutta.chr") %>%
  arrange(ChrNum)

avk <- read_table("A273", "A273.orthogroups.csv.gz") %>%
  transmute(
    HOG,
    AVK_gene = GeneName,
    AVK_chr = ChrName,
    AVK_position = as.numeric(chromStart)
  ) %>%
  distinct() %>%
  filter(!is.na(HOG), !is.na(AVK_gene))

trout_raw <- read_table(
  "Salmo_trutta",
  "Salmo_trutta.orthogroups.csv.gz"
)

###############################################################################
# 3. Define modern gene identifiers
#
# Count a gene once per chromosome label, even if multiple proteins are listed.
# Protein identifiers are used only when a gene name is unavailable.
# Unplaced genes remain included in copy-number counts.
###############################################################################

trout <- trout_raw %>%
  mutate(
    chromosome_key = case_when(
      !is.na(ChrNum) ~ paste0("ChrNum_", ChrNum),
      !is.na(ChrName) ~ ChrName,
      TRUE ~ "unplaced"
    ),
    gene_identifier = case_when(
      !is.na(GeneName) & GeneName != "" ~ GeneName,
      !is.na(ProtName) & ProtName != "" ~ ProtName,
      TRUE ~ NA_character_
    ),
    Trout_gene = paste(chromosome_key, gene_identifier, sep = "::")
  ) %>%
  filter(!is.na(HOG), !is.na(gene_identifier)) %>%
  transmute(
    HOG,
    Trout_gene,
    Trout_chr_num = ChrNum,
    Trout_start = as.numeric(chromStart),
    Trout_end = as.numeric(chromEnd)
  ) %>%
  distinct()

# A modern gene must not be assigned to several different HOGs.
conflicting_genes <- trout %>%
  distinct(Trout_gene, HOG) %>%
  count(Trout_gene) %>%
  filter(n > 1)

if (nrow(conflicting_genes) > 0) {
  write_tsv(
    semi_join(trout, conflicting_genes, by = "Trout_gene"),
    file.path(output_dir, "conflicting_modern_gene_assignments.tsv")
  )
  stop("Some modern genes have multiple HOG assignments. See the exported table.")
}

# Keep one record per modern gene and HOG.
# If several transcript coordinates exist, retain their total span.
trout <- trout %>%
  group_by(HOG, Trout_gene, Trout_chr_num) %>%
  summarise(
    Trout_start = if (all(is.na(Trout_start))) {
      NA_real_
    } else {
      min(Trout_start, na.rm = TRUE)
    },
    Trout_end = if (all(is.na(Trout_end))) {
      NA_real_
    } else {
      max(Trout_end, na.rm = TRUE)
    },
    .groups = "drop"
  )

###############################################################################
# 4. Restrict gene-level analyses to HOGs occurring once in A273
###############################################################################

ancestral_hog_counts <- avk %>%
  distinct(HOG, AVK_gene) %>%
  count(HOG, name = "n_AVK_genes")

ambiguous_hogs <- ancestral_hog_counts %>%
  filter(n_AVK_genes > 1)

write_tsv(
  avk %>% inner_join(ambiguous_hogs, by = "HOG"),
  file.path(output_dir, "ambiguous_ancestral_HOGs.tsv")
)

avk_unique <- avk %>%
  semi_join(
    ancestral_hog_counts %>% filter(n_AVK_genes == 1),
    by = "HOG"
  )

if (nrow(avk_unique) == 0) {
  stop("No HOG occurs exactly once in A273.")
}

message("Ancestral genes retained for analysis: ", nrow(avk_unique))
message("Ambiguous ancestral HOGs excluded: ", nrow(ambiguous_hogs))

###############################################################################
# 5. Count modern copies, including unplaced genes
###############################################################################

modern_counts <- trout %>%
  distinct(HOG, Trout_gene) %>%
  count(HOG, name = "n_modern")

retention <- avk_unique %>%
  left_join(modern_counts, by = "HOG") %>%
  mutate(
    n_modern = replace_na(n_modern, 0L),
    Copy_class = case_when(
      n_modern == 0 ~ "0 detected",
      n_modern == 1 ~ "1 copy",
      n_modern == 2 ~ "2 copies",
      n_modern == 3 ~ "3 copies",
      TRUE ~ "4+ copies"
    ),
    Copy_class = factor(Copy_class, levels = copy_levels),
    AVK_chr = factor(AVK_chr, levels = avk_chr$ChrName)
  )

write_tsv(
  retention,
  file.path(output_dir, "ancestral_gene_copy_retention.tsv")
)

###############################################################################
# 6. Build chromosome offsets for the global dotplot
#
# X: ancestral gene-order units.
# Y: modern physical coordinates in megabases.
###############################################################################

avk_layout <- avk_chr %>%
  transmute(
    AVK_chr = ChrName,
    length = as.numeric(ChrSize),
    offset = lag(cumsum(length), default = 0),
    centre = offset + length / 2
  )

trout_layout <- trout_chr %>%
  transmute(
    Trout_chr_num = ChrNum,
    Trout_chr = ChrName,
    length = as.numeric(ChrSize),
    offset = lag(cumsum(length), default = 0),
    centre = (offset + length / 2) / 1e6
  )

anchors <- avk_unique %>%
  inner_join(trout, by = "HOG") %>%
  inner_join(
    trout_layout,
    by = "Trout_chr_num"
  ) %>%
  filter(is.finite(Trout_start), is.finite(Trout_end)) %>%
  inner_join(
    avk_layout %>% select(AVK_chr, AVK_offset = offset),
    by = "AVK_chr"
  ) %>%
  left_join(
    retention %>% select(HOG, Copy_class),
    by = "HOG"
  ) %>%
  mutate(
    x = AVK_offset + AVK_position,
    y = (offset + (Trout_start + Trout_end) / 2) / 1e6
  )

if (nrow(anchors) == 0) {
  stop("No chromosome-positioned anchors were found.")
}

write_tsv(
  anchors,
  file.path(output_dir, "positioned_anchors.tsv")
)

###############################################################################
# 7. Plot theme and export function
###############################################################################

theme_article <- function() {
  theme_classic(base_size = 10) +
    theme(
      axis.text = element_text(colour = "grey20"),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

save_plot <- function(p, name, width, height) {
  ggsave(
    file.path(output_dir, paste0(name, ".pdf")),
    p, width = width, height = height,
    device = grDevices::pdf, useDingbats = FALSE
  )
  ggsave(
    file.path(output_dir, paste0(name, ".png")),
    p, width = width, height = height,
    dpi = 400, bg = "white"
  )
}

###############################################################################
# 8. Figure 1: whole-genome dotplot
###############################################################################

p1 <- ggplot(anchors, aes(x, y)) +
  geom_vline(
    xintercept = avk_layout$offset,
    colour = "grey88", linewidth = 0.25
  ) +
  geom_hline(
    yintercept = trout_layout$offset / 1e6,
    colour = "grey88", linewidth = 0.25
  ) +
  geom_point(
    aes(colour = Copy_class),
    size = 0.35, alpha = 0.65
  ) +
  scale_colour_manual(
    values = copy_colours,
    drop = FALSE,
    name = "Modern copies"
  ) +
  scale_x_continuous(
    breaks = avk_layout$centre,
    labels = sub("^A273_", "", avk_layout$AVK_chr),
    expand = expansion(mult = 0.005)
  ) +
  scale_y_continuous(
    breaks = trout_layout$centre,
    labels = trout_layout$Trout_chr_num,
    expand = expansion(mult = 0.005)
  ) +
  labs(
    title = "Ancestral vertebrate karyotype versus brown trout",
    subtitle = "Unambiguous ancestral HOGs; chromosome-positioned modern genes",
    x = "A273 chromosome — ancestral gene order",
    y = "Salmo trutta chromosome — physical position"
  ) +
  theme_article() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

save_plot(p1, "01_AVK_trout_dotplot", 13, 11)

###############################################################################
# 9. Figure 2: modern copy-number composition by ancestral chromosome
###############################################################################

copy_summary <- retention %>%
  count(AVK_chr, Copy_class, .drop = FALSE) %>%
  group_by(AVK_chr) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup()

p2 <- ggplot(
  copy_summary,
  aes(x = AVK_chr, y = Proportion, fill = Copy_class)
) +
  geom_col(width = 0.8) +
  scale_fill_manual(
    values = copy_colours,
    drop = FALSE,
    name = "Modern copies"
  ) +
  scale_y_continuous(
    labels = scales::label_percent(),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "Modern copy retention by ancestral chromosome",
    subtitle = "Denominator: all analysed ancestral genes, including undetected genes",
    x = "A273 chromosome",
    y = "Proportion of ancestral genes"
  ) +
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

save_plot(p2, "02_copy_retention_by_AVK_chromosome", 12, 5)

###############################################################################
# 10. Figure 3: singleton profiles along ancestral chromosomes
#
# Use non-overlapping windows of ancestral gene-order positions.
# Singleton proportion is calculated among genes with >= 1 detected copy.
###############################################################################

window_summary <- retention %>%
  mutate(
    Window = floor((AVK_position - 1) / window_size)
  ) %>%
  group_by(AVK_chr, Window) %>%
  summarise(
    Position = (first(Window) + 0.5) * window_size,
    n_analysed = n(),
    n_detected = sum(n_modern > 0),
    n_singleton = sum(n_modern == 1),
    .groups = "drop"
  ) %>%
  mutate(
    Singleton_fraction = if_else(
      n_detected >= minimum_detected_per_window,
      n_singleton / pmax(n_detected, 1),
      NA_real_
    )
  )

p3 <- ggplot(
  window_summary,
  aes(x = Position, y = Singleton_fraction)
) +
  geom_line(colour = "#0072B2", linewidth = 0.5, na.rm = TRUE) +
  geom_point(
    aes(size = n_detected),
    colour = "#0072B2", alpha = 0.8,
    na.rm = TRUE
  ) +
  facet_wrap(vars(AVK_chr), scales = "free_x", ncol = 4) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::label_percent()
  ) +
  scale_size_continuous(range = c(1, 3), name = "Detected genes") +
  labs(
    title = "Local singleton retention along the ancestral genome",
    subtitle = paste0(
      window_size, "-position windows; at least ",
      minimum_detected_per_window, " detected genes per window"
    ),
    x = "Ancestral gene-order position",
    y = "Singletons among detected ancestral genes"
  ) +
  theme_article()

save_plot(p3, "03_local_singleton_retention", 13, 10)

###############################################################################
# 11. Figure 4: chromosome correspondence heatmap
#
# Each cell counts distinct ancestral genes with a modern match.
# Values are not normalized and do not establish collinearity.
###############################################################################

chromosome_links <- anchors %>%
  distinct(AVK_chr, Trout_chr_num, HOG) %>%
  count(AVK_chr, Trout_chr_num, name = "n_ancestral_genes") %>%
  complete(
    AVK_chr = avk_layout$AVK_chr,
    Trout_chr_num = trout_layout$Trout_chr_num,
    fill = list(n_ancestral_genes = 0)
  ) %>%
  mutate(
    AVK_chr = factor(AVK_chr, levels = avk_layout$AVK_chr),
    Trout_chr_num = factor(
      Trout_chr_num,
      levels = rev(trout_layout$Trout_chr_num)
    )
  )

p4 <- ggplot(
  chromosome_links,
  aes(x = AVK_chr, y = Trout_chr_num, fill = n_ancestral_genes)
) +
  geom_tile(colour = "white", linewidth = 0.15) +
  scale_fill_viridis_c(
    trans = "sqrt",
    name = "Ancestral genes"
  ) +
  labs(
    title = "Chromosome-level ancestral correspondence",
    x = "A273 chromosome",
    y = "Salmo trutta chromosome"
  ) +
  theme_article() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

save_plot(p4, "04_chromosome_correspondence", 12, 9)

###############################################################################
# 12. Export numerical summaries and mapping coverage
###############################################################################

chromosome_summary <- retention %>%
  group_by(AVK_chr) %>%
  summarise(
    Analysed_genes = n(),
    Undetected_genes = sum(n_modern == 0),
    Detected_genes = sum(n_modern > 0),
    Singleton_genes = sum(n_modern == 1),
    Multicopy_genes = sum(n_modern >= 2),
    Singleton_fraction_among_detected =
      Singleton_genes / pmax(Detected_genes, 1),
    .groups = "drop"
  )

write_tsv(
  chromosome_summary,
  file.path(output_dir, "retention_by_ancestral_chromosome.tsv")
)

write_tsv(
  window_summary,
  file.path(output_dir, "local_singleton_retention.tsv")
)

write_tsv(
  chromosome_links,
  file.path(output_dir, "chromosome_correspondence.tsv")
)

mapped_genes <- anchors %>%
  distinct(Trout_gene)

message("Modern genes counted: ", n_distinct(trout$Trout_gene))
message("Modern genes plotted: ", nrow(mapped_genes))
message("Plots and tables exported to: ", output_dir)

print(p1)
print(p2)
print(p3)
print(p4)






library(tidyverse)

###############################################################################
# 1. Settings
###############################################################################

input_dir <- "A273_vs_Salmo_trutta"
output_dir <- file.path(input_dir, "evolutionary_analysis")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Exploratory chain-detection parameters.
min_chain_genes <- 6
max_ancestral_gap <- 25
max_modern_gap_bp <- 2e6

# Minimum number of shared chain anchors between candidate regions.
min_shared_anchors <- 4

# Number of candidate region pairs illustrated in detail.
n_pairs_to_plot <- 6

###############################################################################
# 2. Read previous outputs
###############################################################################

retention <- read_tsv(
  file.path(input_dir, "ancestral_gene_copy_retention.tsv"),
  show_col_types = FALSE
)

anchors <- read_tsv(
  file.path(input_dir, "positioned_anchors.tsv"),
  show_col_types = FALSE
) %>%
  distinct(
    HOG, AVK_chr, AVK_position,
    Trout_gene, Trout_chr_num, Trout_start, Trout_end,
    .keep_all = TRUE
  )

required_retention <- c("HOG", "AVK_chr", "AVK_position", "n_modern")
required_anchors <- c(
  "HOG", "AVK_chr", "AVK_position", "Trout_gene",
  "Trout_chr_num", "Trout_start", "Trout_end"
)

stopifnot(
  all(required_retention %in% names(retention)),
  all(required_anchors %in% names(anchors)),
  !anyDuplicated(retention$HOG)
)

chr_order <- retention %>%
  distinct(AVK_chr) %>%
  arrange(as.numeric(sub("^A273_", "", AVK_chr))) %>%
  pull(AVK_chr)

anchors <- anchors %>%
  mutate(Modern_position = (Trout_start + Trout_end) / 2)

###############################################################################
# 3. Plot theme and export function
###############################################################################

theme_article <- function() {
  theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      axis.text = element_text(colour = "grey20"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

save_plot <- function(plot, filename, width, height) {
  
  ggsave(
    file.path(output_dir, paste0(filename, ".pdf")),
    plot, width = width, height = height,
    device = grDevices::pdf, useDingbats = FALSE
  )
  
  ggsave(
    file.path(output_dir, paste0(filename, ".png")),
    plot, width = width, height = height,
    dpi = 400, bg = "white"
  )
}

###############################################################################
# 4. Singleton retention with Wilson intervals
#
# Intervals describe binomial sampling uncertainty.
# They do not account for dependence among neighbouring genes.
###############################################################################

singleton_summary <- retention %>%
  group_by(AVK_chr) %>%
  summarise(
    n_analysed = n(),
    n_detected = sum(n_modern > 0),
    n_singleton = sum(n_modern == 1),
    n_multicopy = sum(n_modern >= 2),
    n_undetected = sum(n_modern == 0),
    .groups = "drop"
  ) %>%
  mutate(
    proportion = if_else(
      n_detected > 0,
      n_singleton / pmax(n_detected, 1),
      NA_real_
    ),
    z = qnorm(0.975),
    denominator = 1 + z^2 / pmax(n_detected, 1),
    centre = (
      proportion + z^2 / (2 * pmax(n_detected, 1))
    ) / denominator,
    half_width = z * sqrt(
      proportion * (1 - proportion) / pmax(n_detected, 1) +
        z^2 / (4 * pmax(n_detected, 1)^2)
    ) / denominator,
    lower = pmax(0, centre - half_width),
    upper = pmin(1, centre + half_width)
  ) %>%
  select(-z, -denominator, -centre, -half_width)

write_tsv(
  singleton_summary,
  file.path(output_dir, "01_singleton_summary_with_intervals.tsv")
)

global_singleton_fraction <- with(
  singleton_summary,
  sum(n_singleton) / sum(n_detected)
)

p1 <- singleton_summary %>%
  mutate(AVK_chr = fct_reorder(AVK_chr, proportion)) %>%
  ggplot(aes(x = proportion, y = AVK_chr)) +
  geom_vline(
    xintercept = global_singleton_fraction,
    linetype = 2, colour = "grey55"
  ) +
  geom_segment(
    aes(x = lower, xend = upper, yend = AVK_chr),
    colour = "#7DAFC3", linewidth = 0.5
  ) +
  geom_point(aes(size = n_detected), colour = "#0072B2") +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_size_continuous(range = c(1.5, 4), name = "Detected genes") +
  labs(
    title = "Singleton retention across ancestral units",
    subtitle = "95% Wilson intervals; dashed line: pooled singleton fraction",
    x = "Singletons among detected ancestral genes",
    y = NULL
  ) +
  theme_article()

save_plot(p1, "01_singleton_retention_intervals", 8, 11)

###############################################################################
# 5. Mapping coverage and copy-number sensitivity
#
# Compare copy counts using all annotated genes with counts using only
# chromosome-positioned genes. A singleton can otherwise be misclassified.
###############################################################################

positioned_counts <- anchors %>%
  distinct(HOG, Trout_gene) %>%
  count(HOG, name = "n_positioned")

mapping_qc <- retention %>%
  left_join(positioned_counts, by = "HOG") %>%
  mutate(
    n_positioned = replace_na(n_positioned, 0L),
    n_not_plotted = n_modern - n_positioned,
    singleton_all_annotations = n_modern == 1,
    singleton_positioned_only = n_positioned == 1,
    apparent_singleton_from_mapping_filter =
      n_positioned == 1 & n_modern > 1
  )

if (any(mapping_qc$n_not_plotted < 0)) {
  stop("Positioned counts exceed total counts: check input consistency.")
}

write_tsv(
  mapping_qc,
  file.path(output_dir, "02_mapping_and_singleton_sensitivity.tsv")
)

qc_summary <- mapping_qc %>%
  summarise(
    analysed_ancestral_genes = n(),
    detected_ancestral_genes = sum(n_modern > 0),
    modern_gene_matches = sum(n_modern),
    positioned_modern_gene_matches = sum(n_positioned),
    true_observed_singletons = sum(singleton_all_annotations),
    apparent_singletons_if_unpositioned_genes_ignored =
      sum(apparent_singleton_from_mapping_filter)
  )

write_tsv(
  qc_summary,
  file.path(output_dir, "02_mapping_QC_summary.tsv")
)

###############################################################################
# 6. Normalized chromosome correspondence
#
# Each cell is the fraction of analysed ancestral genes with at least one
# positioned descendant on the indicated modern chromosome.
# A gene may contribute to several modern chromosomes.
###############################################################################

chromosome_links <- anchors %>%
  distinct(AVK_chr, Trout_chr_num, HOG) %>%
  count(AVK_chr, Trout_chr_num, name = "n_genes") %>%
  complete(
    AVK_chr = chr_order,
    Trout_chr_num = sort(unique(anchors$Trout_chr_num)),
    fill = list(n_genes = 0)
  ) %>%
  left_join(
    singleton_summary %>% select(AVK_chr, n_analysed),
    by = "AVK_chr"
  ) %>%
  mutate(
    fraction = n_genes / n_analysed,
    AVK_chr = factor(AVK_chr, levels = chr_order),
    Trout_chr_num = factor(
      Trout_chr_num,
      levels = rev(sort(unique(anchors$Trout_chr_num)))
    )
  )

write_tsv(
  chromosome_links,
  file.path(output_dir, "03_normalized_chromosome_correspondence.tsv")
)

p2 <- ggplot(
  chromosome_links,
  aes(AVK_chr, Trout_chr_num, fill = fraction)
) +
  geom_tile(colour = "white", linewidth = 0.1) +
  scale_fill_viridis_c(
    labels = scales::label_percent(),
    name = "Ancestral genes\nrepresented"
  ) +
  labs(
    title = "Normalized ancestral chromosome correspondence",
    subtitle = "Fractions are not additive across modern chromosomes",
    x = "Ancestral unit", y = "Modern chromosome index"
  ) +
  theme_article() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

save_plot(p2, "03_normalized_correspondence", 13, 9)

###############################################################################
# 7. Exploratory detection of ordered anchor chains
#
# Find increasing and decreasing chains with bounded gaps.
# Retain the longest chain, remove its modern genes, and repeat.
#
# These are candidate chains, not validated WGD-derived blocks.
# Ancestral positions must represent a reconstructed gene order.
###############################################################################

longest_chain <- function(d, orientation) {
  
  d <- d %>% arrange(AVK_position, Modern_position)
  
  n <- nrow(d)
  score <- rep(1L, n)
  previous <- integer(n)
  
  for (i in seq_len(n)) {
    
    if (i == 1L) next
    
    j <- seq_len(i - 1L)
    
    ancestral_gap <- d$AVK_position[i] - d$AVK_position[j]
    modern_gap <- orientation *
      (d$Modern_position[i] - d$Modern_position[j])
    
    eligible <- j[
      ancestral_gap > 0 &
        ancestral_gap <= max_ancestral_gap &
        modern_gap > 0 &
        modern_gap <= max_modern_gap_bp
    ]
    
    if (length(eligible) > 0) {
      best <- eligible[which.max(score[eligible])]
      score[i] <- score[best] + 1L
      previous[i] <- best
    }
  }
  
  endpoint <- which.max(score)
  indices <- endpoint
  
  while (previous[endpoint] > 0L) {
    endpoint <- previous[endpoint]
    indices <- c(endpoint, indices)
  }
  
  d[indices, ] %>%
    mutate(Orientation = if_else(orientation == 1, "+", "-"))
}

detect_chains <- function(d, key) {
  
  remaining <- d
  results <- list()
  block_number <- 0L
  
  while (nrow(remaining) >= min_chain_genes) {
    
    forward <- longest_chain(remaining, 1)
    reverse <- longest_chain(remaining, -1)
    
    best <- if (nrow(forward) >= nrow(reverse)) forward else reverse
    
    if (nrow(best) < min_chain_genes) break
    
    block_number <- block_number + 1L
    
    best$Block <- paste(
      key$AVK_chr,
      paste0("T", key$Trout_chr_num),
      paste0("B", block_number),
      sep = "_"
    )
    
    results[[block_number]] <- best
    
    remaining <- remaining %>%
      filter(!Trout_gene %in% best$Trout_gene)
  }
  
  if (length(results) == 0) {
    return(tibble())
  }
  
  bind_rows(results)
}

chain_anchors <- anchors %>%
  group_by(AVK_chr, Trout_chr_num) %>%
  group_modify(~ detect_chains(.x, .y)) %>%
  ungroup()

if (nrow(chain_anchors) == 0) {
  stop(
    "No candidate chains found. Earlier QC outputs were saved. ",
    "Inspect gene order and gap thresholds."
  )
}

blocks <- chain_anchors %>%
  group_by(AVK_chr, Trout_chr_num, Block, Orientation) %>%
  summarise(
    n_anchors = n_distinct(HOG),
    AVK_start = min(AVK_position),
    AVK_end = max(AVK_position),
    Modern_start = min(Modern_position),
    Modern_end = max(Modern_position),
    .groups = "drop"
  )

write_tsv(
  chain_anchors,
  file.path(output_dir, "04_candidate_chain_anchors.tsv")
)

write_tsv(
  blocks,
  file.path(output_dir, "04_candidate_chains.tsv")
)

###############################################################################
# 8. Compare candidate regions covering the same ancestral interval
#
# Require different modern chromosomes and shared chain anchors.
# Regions on the same chromosome are not tested in this first-pass analysis.
###############################################################################

candidate_rows <- list()
counter <- 0L

for (ancestral_unit in unique(blocks$AVK_chr)) {
  
  b <- blocks %>% filter(AVK_chr == ancestral_unit)
  
  if (nrow(b) < 2L) next
  
  combinations <- combn(seq_len(nrow(b)), 2)
  
  for (k in seq_len(ncol(combinations))) {
    
    a <- b[combinations[1, k], ]
    z <- b[combinations[2, k], ]
    
    if (a$Trout_chr_num == z$Trout_chr_num) next
    
    overlap_start <- max(a$AVK_start, z$AVK_start)
    overlap_end <- min(a$AVK_end, z$AVK_end)
    
    if (overlap_start >= overlap_end) next
    
    h1 <- chain_anchors %>%
      filter(
        Block == a$Block,
        between(AVK_position, overlap_start, overlap_end)
      ) %>%
      pull(HOG)
    
    h2 <- chain_anchors %>%
      filter(
        Block == z$Block,
        between(AVK_position, overlap_start, overlap_end)
      ) %>%
      pull(HOG)
    
    shared <- length(intersect(h1, h2))
    
    if (shared < min_shared_anchors) next
    
    counter <- counter + 1L
    
    candidate_rows[[counter]] <- tibble(
      AVK_chr = ancestral_unit,
      Block_A = a$Block,
      Block_B = z$Block,
      Chr_A = a$Trout_chr_num,
      Chr_B = z$Trout_chr_num,
      AVK_start = overlap_start,
      AVK_end = overlap_end,
      A_start = a$Modern_start,
      A_end = a$Modern_end,
      B_start = z$Modern_start,
      B_end = z$Modern_end,
      Shared_chain_anchors = shared
    )
  }
}

if (length(candidate_rows) == 0) {
  stop(
    "Candidate chains were exported, but no region pair passed ",
    "the shared-anchor threshold."
  )
}

pairs <- bind_rows(candidate_rows) %>%
  arrange(desc(Shared_chain_anchors)) %>%
  mutate(Pair = paste0("Pair_", row_number()), .before = 1)

write_tsv(
  pairs,
  file.path(output_dir, "05_candidate_region_pairs.tsv")
)

###############################################################################
# 9. Retention within candidate region pairs
#
# Presence is evaluated inside modern region boundaries, not across the
# entire modern chromosome.
#
# "Neither region" does NOT mean genome-wide gene loss.
# Copy counts and matches elsewhere are exported separately.
###############################################################################

pair_gene_tables <- vector("list", nrow(pairs))

for (i in seq_len(nrow(pairs))) {
  
  p <- pairs[i, ]
  
  ancestral_genes <- retention %>%
    filter(
      AVK_chr == p$AVK_chr,
      between(AVK_position, p$AVK_start, p$AVK_end)
    )
  
  in_a <- anchors %>%
    filter(
      AVK_chr == p$AVK_chr,
      Trout_chr_num == p$Chr_A,
      between(Modern_position, p$A_start, p$A_end)
    ) %>%
    pull(HOG) %>%
    unique()
  
  in_b <- anchors %>%
    filter(
      AVK_chr == p$AVK_chr,
      Trout_chr_num == p$Chr_B,
      between(Modern_position, p$B_start, p$B_end)
    ) %>%
    pull(HOG) %>%
    unique()
  
  pair_gene_tables[[i]] <- ancestral_genes %>%
    mutate(
      Pair = p$Pair,
      Present_A = HOG %in% in_a,
      Present_B = HOG %in% in_b,
      Retention = case_when(
        Present_A & Present_B ~ "Both regions",
        Present_A ~ "A only",
        Present_B ~ "B only",
        TRUE ~ "Neither region"
      ),
      Found_elsewhere_or_unpositioned =
        !Present_A & !Present_B & n_modern > 0
    )
}

pair_genes <- bind_rows(pair_gene_tables)

pair_summary <- pair_genes %>%
  group_by(Pair) %>%
  summarise(
    n_ancestral_genes = n(),
    Both = sum(Present_A & Present_B),
    A_only = sum(Present_A & !Present_B),
    B_only = sum(!Present_A & Present_B),
    Neither = sum(!Present_A & !Present_B),
    Neither_but_detected_elsewhere =
      sum(Found_elsewhere_or_unpositioned),
    .groups = "drop"
  ) %>%
  mutate(
    # Positive values indicate more genes retained only in region A.
    Asymmetry = if_else(
      A_only + B_only > 0,
      (A_only - B_only) / pmax(A_only + B_only, 1),
      NA_real_
    )
  ) %>%
  left_join(pairs, by = "Pair")

write_tsv(
  pair_genes,
  file.path(output_dir, "06_gene_retention_in_candidate_pairs.tsv")
)

write_tsv(
  pair_summary,
  file.path(output_dir, "06_candidate_pair_retention_summary.tsv")
)

###############################################################################
# 10. Figure: retention composition of the strongest candidate pairs
###############################################################################

selected_pairs <- pairs %>%
  slice_head(n = n_pairs_to_plot) %>%
  pull(Pair)

pair_labels <- pairs %>%
  mutate(
    Label = paste0(
      Pair, " | ", AVK_chr,
      " | modern ", Chr_A, " vs ", Chr_B
    )
  ) %>%
  select(Pair, Label)

composition <- pair_genes %>%
  filter(Pair %in% selected_pairs) %>%
  count(Pair, Retention) %>%
  group_by(Pair) %>%
  mutate(Fraction = n / sum(n)) %>%
  ungroup() %>%
  left_join(pair_labels, by = "Pair")

p3 <- ggplot(composition, aes(Label, Fraction, fill = Retention)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Both regions" = "#009E73",
      "A only" = "#0072B2",
      "B only" = "#E69F00",
      "Neither region" = "#D9D9D9"
    )
  ) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(
    title = "Gene retention in candidate duplicated regions",
    subtitle = "Neither region does not imply absence from the whole genome",
    x = NULL, y = "Ancestral genes in the shared interval",
    fill = NULL
  ) +
  theme_article()

save_plot(p3, "05_candidate_pair_retention", 11, 6)

###############################################################################
# 11. Figure: gene-by-gene retention tracks
###############################################################################

track_data <- pair_genes %>%
  filter(Pair %in% selected_pairs) %>%
  select(Pair, AVK_position, Present_A, Present_B) %>%
  pivot_longer(
    cols = c(Present_A, Present_B),
    names_to = "Region",
    values_to = "Present"
  ) %>%
  mutate(
    Region = recode(Present_A = "A", Present_B = "B", .x = Region)
  )

p4 <- ggplot(
  track_data,
  aes(AVK_position, Region, colour = Present)
) +
  geom_point(shape = 15, size = 1.5) +
  facet_wrap(vars(Pair), scales = "free_x", ncol = 2) +
  scale_colour_manual(
    values = c("FALSE" = "#D9D9D9", "TRUE" = "#0072B2"),
    labels = c("Not detected in region", "Detected in region"),
    name = NULL
  ) +
  labs(
    title = "Complementary and shared gene retention",
    x = "Ancestral gene-order position",
    y = "Candidate region"
  ) +
  theme_article()

save_plot(p4, "06_gene_retention_tracks", 12, 7)

###############################################################################
# 12. Figure: local dotplots for candidate region pairs
###############################################################################

zoom_tables <- list()

for (i in seq_along(selected_pairs)) {
  
  p <- pairs %>% filter(Pair == selected_pairs[i])
  
  local_anchors <- anchors %>%
    filter(
      AVK_chr == p$AVK_chr,
      between(AVK_position, p$AVK_start, p$AVK_end)
    )
  
  zoom_tables[[2 * i - 1]] <- local_anchors %>%
    filter(
      Trout_chr_num == p$Chr_A,
      between(Modern_position, p$A_start, p$A_end)
    ) %>%
    mutate(
      Pair = p$Pair,
      Region = paste0("A: chromosome ", p$Chr_A),
      Relative_Mb = (Modern_position - p$A_start) / 1e6
    )
  
  zoom_tables[[2 * i]] <- local_anchors %>%
    filter(
      Trout_chr_num == p$Chr_B,
      between(Modern_position, p$B_start, p$B_end)
    ) %>%
    mutate(
      Pair = p$Pair,
      Region = paste0("B: chromosome ", p$Chr_B),
      Relative_Mb = (Modern_position - p$B_start) / 1e6
    )
}

zoom_data <- bind_rows(zoom_tables)

p5 <- ggplot(
  zoom_data,
  aes(AVK_position, Relative_Mb, colour = Region)
) +
  geom_point(size = 1, alpha = 0.8) +
  facet_wrap(vars(Pair), scales = "free", ncol = 2) +
  labs(
    title = "Local correspondence of candidate region pairs",
    x = "Ancestral gene-order position",
    y = "Position from candidate region start (Mb)"
  ) +
  theme_article() +
  theme(legend.position = "right")

save_plot(p5, "07_candidate_pair_dotplots", 13, 9)

###############################################################################
# 13. Save parameters and software information
###############################################################################

write_tsv(
  tibble(
    Parameter = c(
      "min_chain_genes",
      "max_ancestral_gap",
      "max_modern_gap_bp",
      "min_shared_anchors"
    ),
    Value = c(
      min_chain_genes,
      max_ancestral_gap,
      max_modern_gap_bp,
      min_shared_anchors
    )
  ),
  file.path(output_dir, "analysis_parameters.tsv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

print(qc_summary)

message("Candidate ordered chains: ", nrow(blocks))
message("Candidate region pairs: ", nrow(pairs))
message("Results exported to: ", output_dir)









library(tidyverse)

###############################################################################
# 1. Settings
###############################################################################

base_dir <- "A273_vs_Salmo_trutta"
previous_dir <- file.path(base_dir, "evolutionary_analysis")
output_dir <- file.path(previous_dir, "validation")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pairs_to_plot <- paste0("Pair_", 1:6)

# Expand each modern region by these distances on both sides.
boundary_margins <- c(0, 250000, 500000, 1000000)

# Conditional permutations for descriptive order checks.
# These are not selection-corrected significance tests.
n_permutations <- 499
set.seed(123)

###############################################################################
# 2. Read input tables
###############################################################################

retention <- read_tsv(
  file.path(base_dir, "ancestral_gene_copy_retention.tsv"),
  show_col_types = FALSE
)

anchors <- read_tsv(
  file.path(base_dir, "positioned_anchors.tsv"),
  show_col_types = FALSE
) %>%
  mutate(Modern_position = (Trout_start + Trout_end) / 2) %>%
  distinct(
    HOG, AVK_chr, AVK_position,
    Trout_gene, Trout_chr_num, Modern_position,
    .keep_all = TRUE
  )

pairs <- read_tsv(
  file.path(previous_dir, "05_candidate_region_pairs.tsv"),
  show_col_types = FALSE
)

stopifnot(
  !anyDuplicated(retention$HOG),
  !anyDuplicated(pairs$Pair)
)

positioned_counts <- anchors %>%
  distinct(HOG, Trout_gene) %>%
  count(HOG, name = "n_positioned")

retention <- retention %>%
  left_join(positioned_counts, by = "HOG") %>%
  mutate(n_positioned = replace_na(n_positioned, 0L))

if (any(retention$n_positioned > retention$n_modern)) {
  stop("Inconsistent inputs: positioned copies exceed total modern copies.")
}

###############################################################################
# 3. Plot settings
###############################################################################

theme_article <- function() {
  theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

save_plot <- function(p, filename, width, height) {
  
  ggsave(
    file.path(output_dir, paste0(filename, ".pdf")),
    p, width = width, height = height,
    device = grDevices::pdf, useDingbats = FALSE
  )
  
  ggsave(
    file.path(output_dir, paste0(filename, ".png")),
    p, width = width, height = height,
    dpi = 400, bg = "white"
  )
}

###############################################################################
# 4. Reclassify gene retention within each candidate pair
#
# Keep the ancestral interval fixed when expanding modern boundaries.
# Distinguish "outside these regions" from genome-wide non-detection.
###############################################################################

classify_pair <- function(p, margin = 0) {
  
  genes <- retention %>%
    filter(
      AVK_chr == p$AVK_chr,
      between(AVK_position, p$AVK_start, p$AVK_end)
    )
  
  region_a <- anchors %>%
    filter(
      AVK_chr == p$AVK_chr,
      Trout_chr_num == p$Chr_A,
      between(
        Modern_position,
        max(0, p$A_start - margin),
        p$A_end + margin
      ),
      HOG %in% genes$HOG
    ) %>%
    distinct(HOG, Trout_gene)
  
  region_b <- anchors %>%
    filter(
      AVK_chr == p$AVK_chr,
      Trout_chr_num == p$Chr_B,
      between(
        Modern_position,
        max(0, p$B_start - margin),
        p$B_end + margin
      ),
      HOG %in% genes$HOG
    ) %>%
    distinct(HOG, Trout_gene)
  
  genes %>%
    left_join(
      count(region_a, HOG, name = "n_A"),
      by = "HOG"
    ) %>%
    left_join(
      count(region_b, HOG, name = "n_B"),
      by = "HOG"
    ) %>%
    mutate(
      n_A = replace_na(n_A, 0L),
      n_B = replace_na(n_B, 0L),
      Pair = p$Pair,
      Margin_bp = margin,
      
      Status = case_when(
        n_A > 0 & n_B > 0 ~ "Both regions",
        n_A > 0 ~ "A only",
        n_B > 0 ~ "B only",
        n_positioned > 0 ~ "Elsewhere on listed chromosomes",
        n_modern > 0 ~ "Detected but not plotted",
        TRUE ~ "No modern descendant detected"
      ),
      
      # A region-specific gene is not necessarily a genomic singleton.
      Genomic_singleton = n_modern == 1,
      
      Exclusive_status = case_when(
        Status == "A only" & Genomic_singleton ~
          "A only: genomic singleton",
        Status == "B only" & Genomic_singleton ~
          "B only: genomic singleton",
        Status == "A only" ~
          "A only: additional copies elsewhere",
        Status == "B only" ~
          "B only: additional copies elsewhere",
        TRUE ~ NA_character_
      )
    )
}

baseline <- map_dfr(
  seq_len(nrow(pairs)),
  function(i) classify_pair(pairs[i, ], margin = 0)
)

write_tsv(
  baseline,
  file.path(output_dir, "01_reclassified_gene_retention.tsv")
)

###############################################################################
# 5. Summarize shared, exclusive and externally retained descendants
###############################################################################

pair_summary <- baseline %>%
  group_by(Pair) %>%
  summarise(
    n_ancestral = n(),
    Both = sum(Status == "Both regions"),
    A_only = sum(Status == "A only"),
    B_only = sum(Status == "B only"),
    Elsewhere_positioned =
      sum(Status == "Elsewhere on listed chromosomes"),
    Detected_not_plotted =
      sum(Status == "Detected but not plotted"),
    No_descendant_detected =
      sum(Status == "No modern descendant detected"),
    
    A_only_genomic_singletons =
      sum(Status == "A only" & Genomic_singleton),
    B_only_genomic_singletons =
      sum(Status == "B only" & Genomic_singleton),
    
    .groups = "drop"
  ) %>%
  mutate(
    n_in_either_region = Both + A_only + B_only,
    Shared_fraction_among_region_retained =
      Both / pmax(n_in_either_region, 1),
    Exclusive_gene_count = A_only + B_only,
    Exclusive_retention_asymmetry = if_else(
      Exclusive_gene_count > 0,
      (A_only - B_only) / pmax(Exclusive_gene_count, 1),
      NA_real_
    )
  ) %>%
  left_join(
    pairs %>% select(Pair, AVK_chr, Chr_A, Chr_B),
    by = "Pair"
  )

write_tsv(
  pair_summary,
  file.path(output_dir, "02_pair_validation_summary.tsv")
)

###############################################################################
# 6. Export modern positions for exclusive and missing-from-region genes
#
# Multiple rows per ancestral gene are intentional: all positioned descendants
# are listed. Genes without positioned descendants retain an empty position.
###############################################################################

genes_to_audit <- baseline %>%
  filter(Status != "Both regions")

audit <- genes_to_audit %>%
  left_join(
    anchors %>%
      select(
        HOG, Trout_gene, Trout_chr_num,
        Trout_start, Trout_end, Modern_position
      ) %>%
      distinct(),
    by = "HOG",
    relationship = "many-to-many"
  )

write_tsv(
  audit,
  file.path(output_dir, "03_exclusive_and_external_gene_locations.tsv")
)

###############################################################################
# 7. Figure: replace the ambiguous grey category
###############################################################################

status_colours <- c(
  "Both regions" = "#009E73",
  "A only" = "#0072B2",
  "B only" = "#E69F00",
  "Elsewhere on listed chromosomes" = "#CC79A7",
  "Detected but not plotted" = "#999999",
  "No modern descendant detected" = "#E5E5E5"
)

composition <- baseline %>%
  filter(Pair %in% pairs_to_plot) %>%
  count(Pair, Status) %>%
  group_by(Pair) %>%
  mutate(Fraction = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    Pair = factor(Pair, levels = rev(pairs_to_plot)),
    Status = factor(Status, levels = names(status_colours))
  )

p1 <- ggplot(composition, aes(Pair, Fraction, fill = Status)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = status_colours, drop = FALSE) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(
    title = "Where are the ancestral gene descendants?",
    subtitle = "Absence from a candidate pair is separated from genome-wide non-detection",
    x = NULL,
    y = "Ancestral genes in the shared interval",
    fill = NULL
  ) +
  theme_article() +
  guides(fill = guide_legend(ncol = 2))

save_plot(p1, "01_resolved_retention_categories", 11, 6)

###############################################################################
# 8. Direct modern-versus-modern gene-order comparison
#
# Retain HOGs with exactly one positioned gene in each candidate region.
# Do not choose the copy giving the best correlation.
# Ancestral gene order is not used in the correlation calculation.
###############################################################################

modern_matches <- list()
order_results <- list()

for (i in seq_len(nrow(pairs))) {
  
  p <- pairs[i, ]
  
  genes <- baseline %>%
    filter(Pair == p$Pair) %>%
    pull(HOG)
  
  a <- anchors %>%
    filter(
      HOG %in% genes,
      Trout_chr_num == p$Chr_A,
      between(Modern_position, p$A_start, p$A_end)
    ) %>%
    distinct(HOG, Trout_gene, Modern_position) %>%
    add_count(HOG, name = "n") %>%
    filter(n == 1) %>%
    transmute(HOG, Gene_A = Trout_gene, Position_A = Modern_position)
  
  b <- anchors %>%
    filter(
      HOG %in% genes,
      Trout_chr_num == p$Chr_B,
      between(Modern_position, p$B_start, p$B_end)
    ) %>%
    distinct(HOG, Trout_gene, Modern_position) %>%
    add_count(HOG, name = "n") %>%
    filter(n == 1) %>%
    transmute(HOG, Gene_B = Trout_gene, Position_B = Modern_position)
  
  matched <- inner_join(a, b, by = "HOG") %>%
    mutate(
      Pair = p$Pair,
      Local_A_Mb = (Position_A - p$A_start) / 1e6,
      Local_B_Mb = (Position_B - p$B_start) / 1e6
    )
  
  modern_matches[[i]] <- matched
  
  rho <- NA_real_
  null_q95 <- NA_real_
  conditional_p <- NA_real_
  
  if (
    nrow(matched) >= 6 &&
    n_distinct(matched$Position_A) > 1 &&
    n_distinct(matched$Position_B) > 1
  ) {
    
    rank_a <- rank(matched$Position_A)
    rank_b <- rank(matched$Position_B)
    
    rho <- cor(rank_a, rank_b)
    
    null_scores <- replicate(
      n_permutations,
      abs(cor(rank_a, sample(rank_b)))
    )
    
    null_q95 <- unname(quantile(null_scores, 0.95))
    
    conditional_p <- (
      1 + sum(null_scores >= abs(rho))
    ) / (n_permutations + 1)
  }
  
  order_results[[i]] <- tibble(
    Pair = p$Pair,
    n_unambiguous_shared_HOGs = nrow(matched),
    Spearman_rho = rho,
    Absolute_rho = abs(rho),
    Conditional_shuffle_q95 = null_q95,
    Conditional_shuffle_p = conditional_p
  )
}

modern_matches <- bind_rows(modern_matches)
order_results <- bind_rows(order_results)

write_tsv(
  modern_matches,
  file.path(output_dir, "04_direct_modern_gene_matches.tsv")
)

write_tsv(
  order_results,
  file.path(output_dir, "04_direct_modern_order_check.tsv")
)

# The conditional shuffle p-values are descriptive only.
# Region selection and chain detection were performed on the same data.
# No genome-wide significance or WGD assignment is inferred.

p2 <- modern_matches %>%
  filter(Pair %in% pairs_to_plot) %>%
  ggplot(aes(Local_A_Mb, Local_B_Mb)) +
  geom_point(size = 1.7, alpha = 0.8, colour = "#0072B2") +
  facet_wrap(vars(Pair), scales = "free", ncol = 2) +
  labs(
    title = "Direct correspondence between modern regions",
    subtitle = "One gene per HOG in each region; ancestral order not used on either axis",
    x = "Position from region A start (Mb)",
    y = "Position from region B start (Mb)"
  ) +
  theme_article()

save_plot(p2, "02_direct_modern_vs_modern_dotplots", 10, 9)

###############################################################################
# 9. Boundary sensitivity analysis
#
# Expand modern boundaries while keeping the ancestral gene set fixed.
# Large changes suggest sensitivity to the initial region delimitation.
###############################################################################

boundary_summary <- map_dfr(
  boundary_margins,
  function(margin) {
    
    map_dfr(
      seq_len(nrow(pairs)),
      function(i) classify_pair(pairs[i, ], margin)
    ) %>%
      group_by(Pair, Margin_bp) %>%
      summarise(
        n_ancestral = n(),
        Both = sum(n_A > 0 & n_B > 0),
        A_only = sum(n_A > 0 & n_B == 0),
        B_only = sum(n_A == 0 & n_B > 0),
        Neither = sum(n_A == 0 & n_B == 0),
        .groups = "drop"
      ) %>%
      mutate(
        Both_fraction = Both / n_ancestral,
        Neither_fraction = Neither / n_ancestral
      )
  }
)

boundary_change <- boundary_summary %>%
  group_by(Pair) %>%
  arrange(Margin_bp, .by_group = TRUE) %>%
  summarise(
    Baseline_both_fraction = first(Both_fraction),
    Expanded_both_fraction = last(Both_fraction),
    Change_in_both_fraction =
      Expanded_both_fraction - Baseline_both_fraction,
    Baseline_neither_fraction = first(Neither_fraction),
    Expanded_neither_fraction = last(Neither_fraction),
    .groups = "drop"
  )

write_tsv(
  boundary_summary,
  file.path(output_dir, "05_boundary_sensitivity.tsv")
)

write_tsv(
  boundary_change,
  file.path(output_dir, "05_boundary_sensitivity_summary.tsv")
)

p3 <- boundary_summary %>%
  filter(Pair %in% pairs_to_plot) %>%
  select(Pair, Margin_bp, Both_fraction, Neither_fraction) %>%
  pivot_longer(
    c(Both_fraction, Neither_fraction),
    names_to = "Category",
    values_to = "Fraction"
  ) %>%
  mutate(
    Category = recode(
      Category,
      Both_fraction = "Both regions",
      Neither_fraction = "Neither region"
    )
  ) %>%
  ggplot(aes(Margin_bp / 1e6, Fraction, colour = Category)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  facet_wrap(vars(Pair), ncol = 2) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::label_percent()
  ) +
  scale_colour_manual(
    values = c(
      "Both regions" = "#009E73",
      "Neither region" = "#777777"
    )
  ) +
  labs(
    title = "Sensitivity to modern region boundaries",
    subtitle = "The ancestral gene interval remains fixed",
    x = "Boundary expansion on each side (Mb)",
    y = "Proportion of ancestral genes",
    colour = NULL
  ) +
  theme_article()

save_plot(p3, "03_boundary_sensitivity", 10, 8)

###############################################################################
# 10. Integrated validation table
###############################################################################

validation <- pair_summary %>%
  left_join(order_results, by = "Pair") %>%
  left_join(boundary_change, by = "Pair") %>%
  mutate(
    Exclusive_genomic_singletons =
      A_only_genomic_singletons + B_only_genomic_singletons,
    
    Exclusive_genes_with_additional_copies =
      A_only + B_only - Exclusive_genomic_singletons
  )

write_tsv(
  validation,
  file.path(output_dir, "06_integrated_validation.tsv")
)

print(
  validation %>%
    filter(Pair %in% pairs_to_plot) %>%
    select(
      Pair, AVK_chr, Chr_A, Chr_B,
      Both, A_only, B_only,
      Exclusive_genomic_singletons,
      Elsewhere_positioned,
      Detected_not_plotted,
      No_descendant_detected,
      n_unambiguous_shared_HOGs,
      Spearman_rho,
      Change_in_both_fraction
    )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

message("Validation tables and figures exported to: ", output_dir)

