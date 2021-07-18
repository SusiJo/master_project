#!/usr/bin/env Rscript



suppressWarnings(suppressMessages(library(tidyverse, quietly=TRUE)))
suppressWarnings(suppressMessages(library(sva, quietly=TRUE)))
suppressWarnings(suppressMessages(library(tibble, quietly=TRUE)))
suppressWarnings(suppressMessages(library(ggplot2, quietly=TRUE)))
suppressWarnings(suppressMessages(library(grid, quietly=TRUE)))
suppressWarnings(suppressMessages(library(gridExtra, quietly=TRUE)))
suppressWarnings(suppressMessages(library(ggpubr, quietly=TRUE)))
suppressWarnings(suppressMessages(library(optparse, quietly=TRUE)))
suppressWarnings(suppressMessages(library(data.table, quietly=TRUE)))
suppressWarnings(suppressMessages(library(limma, quietly=TRUE)))


option_list = list(
    make_option(c("-c", "--counts"), type="character", default=NULL, help="path to raw count table", metavar="character"),
    make_option(c("-m", "--metadata"), type="character", default=NULL, help="path to metadata table", metavar="character"),
    make_option(c("-l", "--lengths"), type="character", default=NULL, help="path to gene_lengths table", metavar="character"),
    make_option(c("-t", "--title"), type="character", default=NULL, help="tissue source", metavar="character"),
    make_option(c("-o", "--outpath"), type="character", default=NULL, help="path for output", metavar="character")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# Validate and read input
if (is.null(opt$counts)){
    print_help(opt_parser)
    stop("Counts table needs to be provided!")
} else {
    path_count_table = opt$counts
}
if (is.null(opt$metadata)){
    print_help(opt_parser)
    stop("Metadata table needs to be provided!")
} else {
    path_metadata = opt$metadata
}
if (is.null(opt$lengths)){
    print_help(opt_parser)
    stop("Gene_lengths table needs to be provided!")
} else {
    path_gene_lengths = opt$lengths
}
if (is.null(opt$title)){
    print_help(opt_parser)
    stop("Tissue source needs to be provided!")
} else {
    tissue = opt$title
}
if (is.null(opt$outpath)){
    print_help(opt_parser)
    stop("Outpath needs to be provided!")
} else {
    path_outpath = opt$outpath
}

# create directories needed
ifelse(!dir.exists(paste(path_outpath, "batch_correction", sep='')), dir.create(paste(path_outpath, "batch_correction", sep=''), showWarnings = FALSE), FALSE)
dir.create(paste(path_outpath, "batch_correction/plots/", sep=''), showWarnings = FALSE)

# new output paths
new_outpath = paste(path_outpath, 'batch_correction/', tissue, sep='')
plot_path = paste(path_outpath, 'batch_correction/plots/', tissue, sep='')

##################################### FUNCTIONS #####################################

# Functions for reading input data

parseInput <- function(gene_df){
    
    df <- as_tibble(gene_df)
    gene_ids <- df$GeneID
    gene_names <- df$GeneName
    # expression matrix
    expr_df <- df %>% select(3:ncol(df))    
    output = list("expr_data"= expr_df, "ids" = gene_ids, "names"= gene_names )
    return(output)    
}

formatSamples <- function(expr_df){
    # adjust sample_names format
    sample_names <- sub("X", "", colnames(expr_df))
    sample_names <- sample_names %>% str_replace_all("\\.", "\\-") 
    return(sample_names)
}

# TPM Function  
# source: https://gist.github.com/slowkow/counts_to_tpm.R
tpm <- function(counts, lengths) {
    rate <- counts / lengths
    rate / sum(rate) * 1e6
}

# mimics scikit-learn.preprocessing.MinMaxScaler()
# mi, ma = range(mi, ma)
MinMaxScaler <- function(x, mi, ma){
    X_sc <- (x- min(x)) / (max(x)-min(x))
    X_sc * (ma-mi) + mi
}

# conversion of counts to TPM
TPM_conversion_df <- function(df, fileName, gene_id, gene_name, table_header, gene_lengths) {
    
    # expression matrix with TPM
    expr_TPM <- apply(df, 2, tpm, lengths=gene_lengths)

    # convert matrix to dataframe    
    TPM_df <- data.frame(expr_TPM)
    TPM_df_out <- data.frame(expr_TPM)
    
    # add columns GeneId, GeneName to dataframe 
    TPM_df_out <- add_column(TPM_df_out, gene_id, .before=1)
    TPM_df_out <- add_column(TPM_df_out, gene_name, .before=2)
    colnames(TPM_df_out) <- table_header
    
    # save dataframe
    write.table(TPM_df_out, fileName, sep="\t", row.names = FALSE, col.names =  TRUE, quote=FALSE)   
    
    return(TPM_df)
    
}

# Function for transposing and scaling 
processing_func <- function(df, scaler, sample_names, gene_id){
    
    df_t <- transpose(df)
    row.names(df_t) <- sample_names
    colnames(df_t) <- gene_id

    # Scaling
    df_scaled <- apply(df_t, 2, scaler, -1,1)

    # replace infinity and is.na with 0
    is.na(df_scaled) <- sapply(df_scaled, is.infinite)
    df_scaled[is.na(df_scaled)] <- 0
    
    return(df_scaled)    
}


# Function for generating PCA plots
# on is either "condition" or "project"

plot_pca <- function(df, metadata, on){
    
    # compute pca
    df_pca <- prcomp(df)
    # compute explained variance on PC1, PC2
    ev <- summary(df_pca)
    pc1_ev <- round((ev$importance[2,1] *100),1)
    pc2_ev <- round((ev$importance[2,2] *100),1)
    
    # color according to condition, or project
    df_out <- as.data.frame(df_pca$x)
    if(on == "condition"){
        df_out$condition <- metadata$SampleType
        head(df_out)
        pca <- ggplot(df_out,aes(x=PC1,y=PC2,color=condition))
        pca <- pca + scale_color_manual(values=c("#636EFA", "#EF553B")) # blue, red
    } 
    else if (on == "project"){
        df_out$project <- metadata$Project
        pca <- ggplot(df_out,aes(x=PC1,y=PC2,color=project))
    }    
    
    # plot pca
    pca <- pca + geom_point() + 
                xlab(paste("PC1 - " , pc1_ev, " %")) +
                ylab(paste("PC2 - " , pc2_ev, " %")) + 
                theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
    
    return(pca)
}

####### LOADING AND PROCESSING COUNT, METADATA AND GENE-LENGTHS TABLE ###############

count_table = read.table(path_count_table, sep="\t", header=TRUE)

metadata_table = read.table(path_metadata, sep=";", header=TRUE)

gene_lengths = read.table(path_gene_lengths, sep=',', header=TRUE)

output <- parseInput(count_table)
expr_df <- output$expr_data
gene_ids = output$ids
gene_names = output$names
sample_names <- formatSamples(expr_df)


# sort sample names of metadata according to order in expr_df
sorted_meta <- metadata_table[order(match(metadata_table[,1],sample_names)),]
# sort geneIds of featureCounts_geneLenghts according to expr_df
sorted_gene_lengths <- gene_lengths[order(match(gene_lengths[,1], count_table$GeneID)),]

# vector of batches
batches = sorted_meta$Batch
table_header <- append(c("GeneID", "GeneName"), sample_names, after=2)



####### PREPROCESSING UNCORRECTED FEATURE COUNTS DATA ###############



# save uncorrected TPM from featureCounts
uncorr_tpm <- TPM_conversion_df(expr_df, paste(new_outpath, "_uncorrected_fcounts_tpm_reduced_geneset.txt", sep=''), 
                gene_ids, gene_names, table_header, sorted_gene_lengths$GeneLength) #sorted_gene_lengths['GeneLength']

# scale uncorrected TPM for PCA
df_uncorr_scaled <- processing_func(uncorr_tpm, MinMaxScaler, sample_names, gene_ids)



####### COMBAT BATCH CORRECTION ###############

# RUN COMBAT
# NOTE: expr_df not working on denbi --> convert to matrix: as_matrix(expr_df) 
print('Running ComBat correction...', quote=FALSE)

combat_edata = ComBat(dat=as.matrix(expr_df), batch=batches, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

combat_df <- TPM_conversion_df(combat_edata, paste(new_outpath, "_combat_corrected_fcounts_tpm_reduced_geneset.txt", sep=''),
                gene_ids, gene_names, table_header, sorted_gene_lengths$GeneLength)

df_combat_scaled <- processing_func(combat_df, MinMaxScaler, sample_names, gene_ids)



####### COMBAT-SEQ BATCH CORRECTION ###############

# combat-seq takes as input a count matrix
expr_mat <- as.matrix(expr_df)

# perform batch correction with batches specified
print('Running ComBat-Seq correction...', quote=FALSE)

combatseq_counts <- ComBat_seq(expr_mat, batch=batches, group=NULL, full_mod=FALSE)

combatseq_df <- TPM_conversion_df(combatseq_counts, paste(new_outpath, "_combatseq_corrected_fcounts_tpm_reduced_geneset.txt", sep=''),
                    gene_ids, gene_names, table_header, sorted_gene_lengths$GeneLength)

df_combatseq_scaled <- processing_func(combatseq_df, MinMaxScaler, sample_names, gene_ids)


####### REMOVE BATCH EFFECT BATCH CORRECTION ###############


# 1) RemoveBatchEffect without model
print('Running removeBatchEffect()...', quote=FALSE)

rbe <- removeBatchEffect(expr_df, batches) 

df_rbe <- TPM_conversion_df(rbe, paste(new_outpath, "_rbe_corrected_fcounts_tpm_reduced_geneset.txt", sep=''),
                gene_ids, gene_names, table_header, sorted_gene_lengths$GeneLength)

df_rbe_scaled <- processing_func(df_rbe, MinMaxScaler, sample_names, gene_ids)




########################################## PCA ANALYSIS ################################################

print('Performing PCA analysis...', quote=FALSE)

############ UNCORRECTED ###################

uc = plot_pca(df_uncorr_scaled, sorted_meta, "condition")
ggsave(uc, filename = paste(plot_path, "_pca_uncorrected_scaled_condition.pdf", sep=''))
ggsave(uc, dpi=300, width=20, height=10, units= 'cm', filename = paste(plot_path, "_pca_uncorrected_scaled_condition.png", sep=''))


up = plot_pca(df_uncorr_scaled, sorted_meta, "project")
ggsave(up, filename = paste(plot_path, "_pca_uncorrected_scaled_project.pdf", sep=''))
ggsave(up, dpi=300, width=20, height=10, units= 'cm', filename = paste(plot_path, "_pca_uncorrected_scaled_project.png", sep=''))



############ COMBAT ###################


cc = plot_pca(df_combat_scaled, sorted_meta, "condition")
ggsave(cc, filename = paste(plot_path, "_pca_combat_scaled_condition.pdf", sep=''))
ggsave(cc, dpi=300, width=20, height=10, units= 'cm', filename = paste(plot_path, "_pca_combat_scaled_condition.png", sep=''))

cp = plot_pca(df_combat_scaled, sorted_meta, "project") 
ggsave(cp, filename = paste(plot_path, "_pca_combat_scaled_project.pdf", sep=''))
ggsave(cp, dpi=300, width=20, height=10, units= 'cm', filename = paste(plot_path, "_pca_combat_scaled_project.png", sep=''))


############ COMBAT-SEQ ###################


csc = plot_pca(df_combatseq_scaled, sorted_meta, "condition") 
ggsave(csc, filename = paste(plot_path, "_pca_combatseq_scaled_condition.pdf", sep=''))
ggsave(csc, dpi=300, width=20, height=10, units= 'cm', filename = paste(plot_path, "_pca_combatseq_scaled_condition.png", sep=''))


csp = plot_pca(df_combatseq_scaled, sorted_meta, "project") 
ggsave(csp, filename = paste(plot_path, "_pca_combatseq_scaled_project.pdf", sep=''))
ggsave(csp, dpi=300, width=20, height=10, units= 'cm', filename = paste(plot_path, "_pca_combatseq_scaled_project.png", sep=''))



############ LIMMA ###################

lc = plot_pca(df_rbe_scaled, sorted_meta, "condition")
ggsave(lc, filename = paste(plot_path, "_pca_limma_scaled_condition.pdf", sep=''))
ggsave(lc, dpi=300, width=20, height=10, units= 'cm', filename = paste(plot_path, "_pca_limma_scaled_condition.png", sep=''))

lp = plot_pca(df_rbe_scaled, sorted_meta, "project")
ggsave(lp, filename = paste(plot_path, "_pca_limma_scaled_project.pdf", sep=''))
ggsave(lp, dpi=300, width=20, height=10, units= 'cm', filename = paste(plot_path, "_pca_limma_scaled_project.png", sep=''))



####### DRAW FIGURE AND EXPORT ###############

# draw all in one figure
figure <- ggarrange(uc, up, cc, cp, csc, csp, lc, lp,
            labels = c("U", "U", 
                        "C", "C",
                        "CS", "CS",
                        "L", "L"),
                        vjust = 0.2,
                        nrow = 4, ncol = 2,
                        legend="none")


figure <- annotate_figure(figure, top = text_grob("", color = "black", face = "bold", size = 12))

ggexport(figure, filename = paste(plot_path, "_all_pcas.pdf", sep=''))
ggexport(figure, res=300, width=3000, height=2000, units= 'cm', filename = paste(plot_path, "_all_pcas.png", sep=''))


conditions <- ggarrange(uc, cc, csc, lc, labels=c("U", "C", "CS", "L"), vjust = 0.2, nrow=4, ncol=1, common.legend=TRUE)
ggexport(conditions, filename=paste(plot_path, "_pca_bc_conditions.pdf", sep=""))
projects <- ggarrange(up, cp, csp, lp, labels=c("U", "C", "CS", "L"), vjust = 0.2, nrow=4, ncol=1, common.legend=TRUE)
ggexport(projects, filename=paste(plot_path, "_pca_bc_projects.pdf", sep=""))