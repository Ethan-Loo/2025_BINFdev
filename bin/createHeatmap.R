################################################
## LOAD LIBRARIES                             ##
################################################
################################################

# Check if the required libraries are installed, and install them if not
if (!require("optparse")) {
    install.packages("optparse", repos="http://cran.us.r-project.org")
}
if (!require("ggplot2")) {
    install.packages("ggplot2", repos="http://cran.us.r-project.org")
}
if (!require("RColorBrewer")) {
    install.packages("RColorBrewer", repos="http://cran.us.r-project.org")
}
if (!require("pheatmap")) {
    install.packages("pheatmap", repos="http://cran.us.r-project.org")
}

library(optparse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################
option_list <- list(
  make_option(c("-i", "--input_file"), type="character", default=NULL, metavar="path", help="Input sample file"),
  make_option(c("-g", "--geneFunctions_file"), type="character", default=NULL, metavar="path", help="Gene Functions file."),
  make_option(c("-a", "--annoData_file"), type="character", default=NULL, metavar="path", help="Annotation Data file."),
  make_option(c("-p", "--outprefix"), type="character", default='projectID', metavar="string", help="Output prefix.")
)


opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

sampleInput=opt$input_file
geneInput=opt$geneFunctions_file
annoInput=opt$annoData_file
outprefix=opt$outprefix

testing="Y"
if (testing == "Y"){
  sampleInput="./data/sampleData.csv"
  geneInput="./data/geneFunctions.csv"
  annoInput="./data/annoData.csv"
  outprefix="test"
}


if (is.null(sampleInput)){
  print_help(opt_parser)
  stop("Please provide an input file.", call.=FALSE)
}

################################################
################################################
## READ IN FILES##
################################################
################################################
sampleData=read.csv(sampleInput,row.names=1)
annoData=read.csv(annoInput,row.names=1)
geneFunctions=read.csv(geneInput,row.names=1)

################################################
################################################
## Set colors##
################################################
################################################
annoColors <- list(
  gene_functions = c("Oxidative_phosphorylation" = "#F46D43",
                     "Cell_cycle" = "#708238",
                     "Immune_regulation" = "#9E0142",
                     "Signal_transduction" = "beige", 
                     "Transcription" = "violet"), 
  Group = c("Disease" = "darkgreen",
            "Control" = "blueviolet"),
  Lymphocyte_count = brewer.pal(5, 'PuBu')
)

################################################
################################################
## Create a basic heatmap##
################################################
################################################

# Convert matrix into a long format data frame
data_matrix <- as.matrix(sampleData)
data_long <- data.frame(
    gene = rep(rownames(data_matrix), ncol(data_matrix)),
    sample = rep(colnames(data_matrix), each = nrow(data_matrix)),
    expression = as.vector(data_matrix)
)

# Add a sorting column to maintain gene order from clustering
gene_clustering <- hclust(dist(data_matrix, method = "euclidean"), method = "ward.D")
gene_order <- rownames(data_matrix)[gene_clustering$order]
data_long$gene <- factor(data_long$gene, levels = gene_order)


ggplot(data_long, aes(x = sample, y = gene)) +
    geom_tile(aes(fill = expression), color = "white") +
    scale_fill_gradientn(
        colors = c("#3A90C4", "#FEF5E7", "#EE4444"), #BWR color scheme
        name = "Expression\nZ-score",
    ) +
    

    theme_minimal() +
    theme(
        plot.title = element_text(size = 12, face = "plain", hjust = 0.5),
      
        axis.title.x = element_text(size = 10, margin = margin(t = 10)),
        axis.title.y = element_text(size = 10, margin = margin(r = 10)),
        
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 7),
        
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.key.height = unit(1, "cm"),
        legend.margin = margin(l = 5, r = 5),
        
        panel.grid = element_blank(),
        aspect.ratio = 3
    ) +
    
    # Add labels
    labs(
        title = paste0("basic_heatmap_", outprefix),
        x = "Samples",
        y = "Genes"
    )

# Save the plot with adjusted dimensions
ggsave(paste0("./computational_output/basic_heatmap_", outprefix, ".pdf"), 
       width = 5, height = 8,
       device = cairo_pdf)

################################################
################################################
## Create a complex heatmap##
################################################
################################################

# Define quantile breaks for our data
data_matrix <- as.matrix(sampleData)
quantile_breaks <- quantile(data_matrix, probs = c(0, 1/3, 2/3, 1))

# Create display matrix (keeping original data for clustering but using categories for display)
discrete_matrix <- matrix(0, nrow=nrow(data_matrix), ncol=ncol(data_matrix))
rownames(discrete_matrix) <- rownames(data_matrix)
colnames(discrete_matrix) <- colnames(data_matrix)

# Assign categories based on the quantile breaks
discrete_matrix[data_matrix <= quantile_breaks[2]] <- 1        # Low expression
discrete_matrix[data_matrix > quantile_breaks[2] & 
               data_matrix <= quantile_breaks[3]] <- 2         # Medium expression
discrete_matrix[data_matrix > quantile_breaks[3]] <- 3        # High expression

my_colors <- c("#2E1F3E", "#8B7B9E", "#E8E6EE") 
names(my_colors) <- c("1", "2", "3")

pheatmap(
    mat = discrete_matrix,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D",
    annotation_col = annoData,
    annotation_row = geneFunctions,
    annotation_colors = annoColors,
    annotation_names_row = FALSE,
    annotation_names_col = FALSE,
    color = my_colors,
    breaks = c(0.5, 1.5, 2.5, 3.5),
    legend_breaks = c(1, 2, 3),
    legend_labels = c("Low", "Medium", "High"),
    main = paste0("Gene Expression Heatmap - ", outprefix),
    fontsize = 8,
    filename = paste0("./computational_output/complex_heatmap_", outprefix, ".pdf")
)
  