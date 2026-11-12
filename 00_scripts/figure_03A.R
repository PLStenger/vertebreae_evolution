library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

###############################################################################
# 1. Read data and build the ultra-HOG matrix
###############################################################################

orthogroup <- read.delim("OG_count_counts.tsv", check.names = FALSE)

species <- names(orthogroup)[
  match("Acanthochromis_polyacanthus", names(orthogroup)):
    match("Zootoca_vivipara", names(orthogroup))
]

ucg <- orthogroup %>% filter(Genome == 271)

uCG.mat <- as.matrix(ucg[, species])
rownames(uCG.mat) <- ucg$HOG
uCG.mat[uCG.mat > 9] <- 10

###############################################################################
# 2. Assign taxonomic groups using genera from the supplied species list
###############################################################################

split_genera <- function(x) strsplit(x, "\\s+")[[1]]

group_genera <- list(
  Mammalia = split_genera(paste(
    "Acinonyx Acomys Ailuropoda Antechinus Apodemus Arvicanthis Arvicola",
    "Balaenoptera Bos Bubalus Budorcas Callithrix Camelus Canis Capra",
    "Cervus Chionomys Choloepus Cricetulus Dasypus Desmodus Diceros",
    "Dromiciops Elephas Eptesicus Equus Felis Gorilla Gracilinanus",
    "Hippopotamus Homo Jaculus Lemur Leopardus Lutra Lynx Macaca Manis",
    "Meles Microcebus Microtus Mirounga Monodelphis Mus Mustela Neofelis",
    "Neogale Neomonachus Nomascus Nycticebus Ochotona Onychomys Orcinus",
    "Ornithorhynchus Oryctolagus Ovis Pan Panthera Papio Perognathus",
    "Peromyscus Phacochoerus Phocoena Phyllostomus Physeter Piliocolobus",
    "Pongo Prionailurus Rattus Rhinolophus Rhinopithecus Sarcophilus",
    "Sciurus Sorex Suncus Suricata Sus Symphalangus Tachyglossus",
    "Theropithecus Trichosurus Tursiops Vulpes Zalophus"
  )),
  
  Aves = split_genera(paste(
    "Accipiter Agelaius Ammospiza Anas Apus Aquila Aythya Calypte",
    "Camarhynchus Catharus Chiroxiphia Corvus Coturnix Cuculus Cygnus",
    "Dryobates Falco Ficedula Gallus Grus Gymnogyps Harpia Hirundo",
    "Indicator Lagopus Lonchura Meleagris Melopsittacus Melospiza",
    "Molothrus Motacilla Numida Oenanthe Oxyura Parus Poecile Rissa",
    "Serinus Strigops Taeniopygia Vidua"
  )),
  
  Chondrichthyes = split_genera(paste(
    "Amblyraja Carcharodon Chiloscyllium Leucoraja Pristis",
    "Rhincodon Scyliorhinus Stegostoma"
  )),
  
  Cladistia = split_genera("Erpetoichthys Polypterus"),
  
  Hyperoartia = "Petromyzon",
  
  Lepidosauria = split_genera(paste(
    "Ahaetulla Anolis Eublepharis Euleptes Hemicordylus Lacerta",
    "Podarcis Sceloporus Sphaerodactylus Thamnophis Zootoca"
  )),
  
  Amphibia = split_genera(paste(
    "Bombina Bufo Geotrypetes Hyla Microcaecilia Rana",
    "Rhinatrema Spea Xenopus"
  )),
  
  Reptilia = split_genera(
    "Caretta Chelonia Chrysemys Dermochelys Gopherus Malaclemys Mauremys Trachemys"
  ),
  
  Sarcopterygii = "Protopterus"
)

# All remaining genera in this specific dataset belong to Actinopteri.
genus <- sub("_.*", "", colnames(uCG.mat))
species_group <- rep("Actinopteri", length(genus))

for (group in names(group_genera)) {
  species_group[genus %in% group_genera[[group]]] <- group
}

###############################################################################
# 3. Define group colours and display order
###############################################################################

group_colours <- c(
  Actinopteri    = "#CC33AA",
  Mammalia       = "#3333CC",
  Chondrichthyes = "#66A61E",
  Cladistia      = "#E6A6CF",
  Hyperoartia    = "#E6B3DD",
  Aves           = "#E64B64",
  Lepidosauria   = "#E68632",
  Amphibia       = "#368093",
  Reptilia       = "#A65D64",
  Sarcopterygii  = "#C58CB4"
)

species_group <- factor(species_group, levels = names(group_colours))
label_colours <- unname(group_colours[as.character(species_group)])

# Export assignments for inspection.
write_tsv(
  tibble(Species = colnames(uCG.mat), Class = species_group),
  "species_taxonomic_groups.tsv"
)

print(table(species_group))

###############################################################################
# 4. Define the copy-number palette
###############################################################################

col_fun <- colorRamp2(
  1:10,
  c(
    "#C6DBEF", "#4292C6", "#ADDD8E", "#FFEDA0", "#FEB24C",
    "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#67000D"
  )
)

###############################################################################
# 5. Shared heatmap function
###############################################################################

make_heatmap <- function(k = 1) {
  Heatmap(
    uCG.mat,
    name = "CNV",
    col = col_fun,
    show_row_names = FALSE,
    row_km = k,
    clustering_method_rows = "ward.D2",
    row_dend_width = unit(25, "mm"),
    
    # Group species by class and cluster within each class.
    column_split = species_group,
    cluster_column_slices = FALSE,
    clustering_method_columns = "ward.D2",
    column_dend_height = unit(20, "mm"),
    column_gap = unit(2, "mm"),
    column_title = NULL,
    
    column_labels = gsub("_", " ", colnames(uCG.mat)),
    column_names_gp = gpar(
      col = label_colours,
      fontsize = 5,
      fontface = "italic"
    ),
    column_names_rot = 90,
    column_names_max_height = unit(75, "mm"),
    
    top_annotation = HeatmapAnnotation(
      Class = species_group,
      col = list(Class = group_colours),
      annotation_legend_param = list(
        Class = list(title = "Class", at = names(group_colours))
      )
    ),
    
    use_raster = TRUE,
    heatmap_legend_param = list(
      title = "CNV",
      at = 1:10,
      labels = c(as.character(1:9), "10+")
    )
  )
}

###############################################################################
# 6. Export both heatmaps
###############################################################################

pdf("ultra_HOGs_coloured_labels.pdf", width = 32, height = 18)
ht_simple <- draw(
  make_heatmap(),
  column_title = "Ultra-HOGs",
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

set.seed(123)

pdf("ultra_HOGs_coloured_labels_km11.pdf", width = 16, height = 7)
ht_drawn <- draw(
  make_heatmap(k = 11),
 # column_title = "Ultra-HOGs - 11 row clusters",
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

###############################################################################
# 7. Export HOGs from cluster 4
###############################################################################

cluster_rows <- row_order(ht_drawn)
genes_cluster4 <- rownames(uCG.mat)[cluster_rows[["4"]]]

writeLines(genes_cluster4, "genes_for_annotation.txt")

message("Both heatmaps exported.")
message("Cluster 4: ", length(genes_cluster4), " HOGs exported.")

