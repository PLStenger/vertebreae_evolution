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
