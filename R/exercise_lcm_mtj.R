# library loading-----------------------------------------------------------------
source(here::here("R/library.R"))

# meta data loading -------------------------------------------------------
metadata <- read_excel("C:/Users/carol/OneDrive/Carolina Jacob/Doutorado/CBMR/MS Prep/Final protocol/MetaData_10Jan.xlsx")

metadata <- metadata %>%
    dplyr::select(sample_id, tissue, intervention, time, sampleprep_date) %>%
    mutate_at(1:5, as.factor) %>%
    arrange(tissue)

usethis::use_data(metadata, overwrite = T)
# data loading -----------------------------------------------

originaldata <- read.delim("C:/Users/carol/OneDrive/Carolina Jacob/Doutorado/CBMR/MS Prep/Final protocol/Caro_Ex_LCM_MTJ_Spectronaut/20240109_093728_20231215_ProteinQuantReport.tsv")
View(X20240109_093728_20231215_ProteinQuantReport)

caroldata <- read.delim("C:/Users/carol/OneDrive/Carolina Jacob/Doutorado/CBMR/MS Prep/Final protocol/Caro_Ex_LCM_MTJ_Spectronaut/20240109_093728_20231215_ProteinQuantReport.tsv", header = T, check.names = T, sep = "\t") %>%
    dplyr::select("PG.ProteinGroups", "PG.Genes", 6:113)


new_colnames <- sapply(strsplit(colnames(caroldata), "_"), '[',10)
colnames(caroldata) <- new_colnames

colnames(caroldata)[1] <- "protein_names"
colnames(caroldata)[2] <- "genes"
colnames(caroldata)[96] <- "S109"

counts <- table(caroldata$genes) #unique gene IDs

caroldata <- caroldata %>%
    mutate(genes = if_else(genes == "", protein_names, genes))  %>%
    mutate_at(3:110, as.numeric) %>%
    mutate_all(~ifelse(is.nan(.), NA, .)) %>%
    mutate(genes = make.names(genes, unique = T), genes) %>%
    column_to_rownames(var = "genes") %>%
    dplyr::select(-c("protein_names")) %>%
    log2()

names(caroldata) <- sub("^S", "", names(caroldata))

sample_order <- metadata$sample_id[order(metadata$tissue)] #
sample_order <- trimws(sample_order)
caroldata <- caroldata[ ,sample_order]

usethis::use_data(caroldata, overwrite = T)

#summary reports ----------
summary(caroldata)
str(caroldata)

#save your file


# NA visualization --------------------------------------------------------
#load your metadata and df again

vis_dat(caroldata)

vis_miss(caroldata, sort_miss = T)

#visualization of number of NA vs number of valid values per sample -----

counts_na <- sapply(caroldata, function (x) sum(is.na(x)))
counts_valid <- sapply(caroldata, function (x) sum(!is.na(x)))

df_na <- tibble(counts_na, counts_valid)
df_na$tissue <- metadata$tissue

df_na <- df_na %>%
    arrange(tissue)

df_na$sampleid <- metadata$sample_id

df_na <- df_na %>% mutate(sampleid = factor(sampleid, levels = sampleid))


ggplot(df_na, aes(x = sampleid)) +
    geom_bar(aes(y = counts_valid, fill = tissue), stat = "identity", position = "dodge", alpha = 0.5) +
    geom_bar(aes(y = counts_na, fill = tissue), stat = "identity", position = "dodge", alpha = 0.9) +
    labs(title = "counts_valid vs counts_na per sample spect", y = "counts") +
    theme_minimal()

#Scurve ------------------------------------------

df_scurve <- caroldata %>% #use the original df
    mutate(genes = ifelse(genes == "", protein_names, genes)) %>%
    mutate_at(3:47, as.numeric) %>%
    mutate_all(~ifelse(is.nan(.), NA, .)) %>%
    mutate(genes = make.names(genes, unique = T), genes) %>%
    column_to_rownames(var = "genes") %>%
    dplyr::select(-c("protein_names")) %>%
    log10()

df_scurve$median = apply(df_scurve, 1, median, na.rm = TRUE)
df_scurve <- df_scurve %>% arrange(desc(median))
df_scurve$rank = (1:2318)
ggplot(df_scurve, aes(rank, median)) +
    geom_line() +
    labs(title = "S-curve_allsamples") +
    ylab("Log10 protein intensities (median)") +
    xlab("Abundance rank") +
    theme_minimal()

#the 10 most abundant proteins ----------
row.names(df_scurve[1:10,])

#the 10 least abundant proteins ----------
row.names(df_scurve[2309:2318,])

#plot densities -----------
plotDensities(caroldata)
caroldata <- caroldata[, -c(91)]
plotDensities(caroldata)
use_data(caroldata, overwrite = T)

metadata <- metadata[-c(19),]

#filtering 70%-------------

load(here::here("data/caroldata.rda"))
caroldata70f <- selectOverallPercent(caroldata, 0.7)

saveRDS(caroldata70f, file = "C:\\Users\\carol\\OneDrive\\Área de Trabalho\\MTJ_CarolinaJacob\\caroldata70f.rds")

#Load the dataset
caroldata70f <- readRDS("C:\\Users\\carol\\OneDrive\\Área de Trabalho\\MTJ_CarolinaJacob\\caroldata70f.rds")

use_data(caroldata70f, overwrite = T)

caroldata <- normalizeQuantiles(caroldata70f)

plotDensities(caroldata)

help("selectOverallPercent")

#######filter
caroldata70f <- selectOverallPercent(caroldata, 0.7)

############normalize
caroldata_70f_qn <- normalizeQuantiles(caroldata70f)

#data visualization
metadata <- metadata %>% arrange(tissue) #mtj, muscle, tendon order
metadata.grp.col <- c("#F8766D", "#00BA38", "#619CFF")
grp.col <- metadata.grp.col[match(metadata$tissue, unique(metadata$tissue))]

boxplot(caroldata_70f_qn, main = "caro_spect_1sampleremoved_NA70f_qn", ylab = "
log2 intensities", col = grp.col)

unique_tissues <- unique(metadata$tissue)
tissue_colors <- grp.col[match(unique_tissues, metadata$tissue)]
legend("topleft", legend = unique_tissues, fill = tissue_colors, horiz = T )


#histogram
par(mfrow=c(2, 2))

for (i in 1:ncol(caroldata_70f_qn)) {
    hist(caroldata_70f_qn[[i]],
         main=paste("Histogram for Column", i),
         xlab="Value",
         ylab="Frequency",
         col="lightblue",
         border="black",
         breaks=20)
}

#filtering 60%-------------
caroldata60f <- selectOverallPercent(caroldata, 0.6)
saveRDS(caroldata60f, file = "C:\\Users\\carol\\OneDrive\\Área de Trabalho\\MTJ_CarolinaJacob\\caroldata60f.rds")

#Load the dataset
caroldata60f <- readRDS("C:\\Users\\carol\\OneDrive\\Área de Trabalho\\MTJ_CarolinaJacob\\caroldata60f.rds")

use_data(caroldata60f, overwrite = T)

caroldata <- normalizeQuantiles(caroldata60f)

plotDensities(caroldata)

help("selectOverallPercent")

#######filter
caroldata60f <- selectOverallPercent(caroldata, 0.6)

############normalize
caroldata_60f_qn <- normalizeQuantiles(caroldata60f)

#data visualization
metadata <- metadata %>% arrange(tissue) #mtj, muscle, tendon order
metadata.grp.col <- c("#F8766D", "#00BA38", "#619CFF")
grp.col <- metadata.grp.col[match(metadata$tissue, unique(metadata$tissue))]

boxplot(caroldata_60f_qn, main = "caro_spect_1sampleremoved_NA60f_qn", ylab = "
log2 intensities", col = grp.col)

unique_tissues <- unique(metadata$tissue)
tissue_colors <- grp.col[match(unique_tissues, metadata$tissue)]
legend("topleft", legend = unique_tissues, fill = tissue_colors, horiz = T )

#histogram
par(mfrow=c(2, 2))

for (i in 1:ncol(caroldata_60f_qn)) {
    hist(caroldata_60f_qn[[i]],
         main=paste("Histogram for Column", i),
         xlab="Value",
         ylab="Frequency",
         col="lightyellow",
         border="black",
         breaks=20)
}

#filtering 50%---------------
caroldata50f <- selectOverallPercent(caroldata, 0.5)
saveRDS(caroldata50f, file = "C:\\Users\\carol\\OneDrive\\Área de Trabalho\\MTJ_CarolinaJacob\\caroldata50f.rds")

#Load the dataset
caroldata50f <- readRDS("C:\\Users\\carol\\OneDrive\\Área de Trabalho\\MTJ_CarolinaJacob\\caroldata50f.rds")

use_data(caroldata50f, overwrite = T)

caroldata <- normalizeQuantiles(caroldata50f)

plotDensities(caroldata)

help("selectOverallPercent")

#######filter
caroldata50f <- selectOverallPercent(caroldata, 0.5)

############normalize
caroldata_50f_qn <- normalizeQuantiles(caroldata50f)

#data visualization
metadata <- metadata %>% arrange(tissue) #mtj, muscle, tendon order
metadata.grp.col <- c("#F8766D", "#00BA38", "#619CFF")
grp.col <- metadata.grp.col[match(metadata$tissue, unique(metadata$tissue))]

boxplot(caroldata_50f_qn, main = "caro_spect_1sampleremoved_NA50f_qn", ylab = "
log2 intensities", col = grp.col)

unique_tissues <- unique(metadata$tissue)
tissue_colors <- grp.col[match(unique_tissues, metadata$tissue)]
legend("topleft", legend = unique_tissues, fill = tissue_colors, horiz = T )

#histogram
par(mfrow=c(2, 2))

for (i in 1:ncol(caroldata_50f_qn)) {
    hist(caroldata_50f_qn[[i]],
         main=paste("Histogram for Column", i),
         xlab="Value",
         ylab="Frequency",
         col="lightgreen",
         border="black",
         breaks=20)
}

#filtering 40%---------------
caroldata40f <- selectOverallPercent(caroldata, 0.4)
saveRDS(caroldata40f, file = "C:\\Users\\carol\\OneDrive\\Área de Trabalho\\MTJ_CarolinaJacob\\caroldata40f.rds")

#Load the dataset
caroldata40f <- readRDS("C:\\Users\\carol\\OneDrive\\Área de Trabalho\\MTJ_CarolinaJacob\\caroldata40f.rds")

use_data(caroldata40f, overwrite = T)

caroldata <- normalizeQuantiles(caroldata40f)

plotDensities(caroldata)

help("selectOverallPercent")

#######filter
caroldata40f <- selectOverallPercent(caroldata, 0.4)

############normalize
caroldata_40f_qn <- normalizeQuantiles(caroldata40f)

#data visualization
metadata <- metadata %>% arrange(tissue) #mtj, muscle, tendon order
metadata.grp.col <- c("#F8766D", "#00BA38", "#619CFF")
grp.col <- metadata.grp.col[match(metadata$tissue, unique(metadata$tissue))]

boxplot(caroldata_40f_qn, main = "caro_spect_1sampleremoved_NA40f_qn", ylab = "
log2 intensities", col = grp.col)

unique_tissues <- unique(metadata$tissue)
tissue_colors <- grp.col[match(unique_tissues, metadata$tissue)]
legend("topleft", legend = unique_tissues, fill = tissue_colors, horiz = T )

#histogram
par(mfrow=c(2, 2))

for (i in 1:ncol(caroldata_40f_qn)) {
    hist(caroldata_40f_qn[[i]],
         main=paste("Histogram for Column", i),
         xlab="Value",
         ylab="Frequency",
         col="lightpink",
         border="black",
         breaks=20)
}

#100%filtered data (na.omit)---------------------
df_nona <- na.omit(caroldata)
df_nona_qn <- normalizeQuantiles(df_nona)

colSums(is.na(df_nona_qn))

pca_nona_qn <- prcomp(t(df_nona_qn), scale = T)
summary(pca_nona_qn)

fviz_pca_ind(pca_nona_qn) #quick view

plot(pca_nona_qn$x[,1], pca_nona_qn$x[,2]) #alternative

pca_var <- pca_nona_qn$sdev^2
pca_var_perc <- round(pca_var/sum(pca_var)*100, digits = 1)

fviz_eig(pca_nona_qn, addlabels = T)

pca_results <- as.data.frame(pca_nona_qn$x)
pca_pc1_pc2 <- pca_results %>%
    dplyr::select(PC1, PC2)
metadata$expgrp <- as.factor(paste(metadata$tissue, metadata$intervention, metadata
                                   $time, sep = "-"))


pca_pc1_pc2 <- pca_pc1_pc2 %>%
    mutate(
        sample = rownames(pca_results),
        exp_group = metadata$expgrp[1:107],
        tissue = metadata$tissue[1:107],
        intervention = metadata$intervention[1:107],
        duration = metadata$time[1:107],
        cryo_date = metadata$cryostat_date[1:107],
        lcm_date = metadata$lcm_date[1:107],
        sampleprep_occ = metadata$sampleprep_occ[1:107],
        sampleprep_date = metadata$sampleprep_date[1:107]
    )

#by tissue ---------
#code1
ggplot(pca_pc1_pc2, aes(x = PC1 , y = PC2, colour = tissue, label = sample)) +
    geom_point(size = 3) +
    geom_text(fontface = "bold", show.legend = FALSE, vjust = 0.5, hjust = -0.5, size = 3)

#hazals code
ggplot(pca_pc1_pc2, aes(x = PC1 , y = PC2, colour = tissue, label = sample)) +
    geom_point(size = 3) +
    #geom_text(fontface = "bold", show.legend = F, vjust = 0.5, hjust = -0.5, size = 3)
    xlab(paste("PC1 - ", pca_var_perc[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca_var_perc[2], "%", sep = "")) +
    #xlim(-50, 50) + ylim(-50, 50) +
    ggtitle("by_tissue_PCA_noNA_qn") +
    theme_minimal()

#by intervention and area -----------
ggplot(pca_pc1_pc2, aes(x = PC1 , y = PC2, colour = tissue, shape = intervention, label = sample)) +
    geom_point(size = 3) +
    scale_shape_manual(values = c("cc" = 1, "ecc" = 16, "sed" = 13)) +
    #geom_text(fontface = "bold", show.legend = F, vjust = 0.5, hjust = -0.5, size = 3) +
    xlab(paste("PC1 - ", pca_var_perc[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca_var_perc[2], "%", sep = "")) +
    #xlim(-50, 50) + ylim(-50, 50) +
    ggtitle("by_intervention_PCA_noNA_qn") +
    theme_minimal()

#by duration and area -----------
# first, convert duration to a factor
pca_pc1_pc2$duration <- factor(pca_pc1_pc2$duration)

ggplot(pca_pc1_pc2, aes(x = PC1 , y = PC2, colour = tissue, shape = duration, label = sample)) +
    geom_point(size = 3) +
    scale_shape_manual(values = c("0" = 0, "5" = 15, "10" = 12)) +
    #geom_text(fontface = "bold", show.legend = F, vjust = 0.5, hjust = -0.5, size = 3) +
    xlab(paste("PC1 - ", pca_var_perc[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca_var_perc[2], "%", sep = "")) +
    #xlim(-50, 50) + ylim(-50, 50) +
    ggtitle("by_time_point_noNA_qn") +
    theme_minimal()







#Creating a binary matrix indicating missing values---------
missing_matrix <- is.na(caroldata)

# Creating the heatmap
pheatmap(missing_matrix, cex_col = 0.5, cex_row = 0.5)


# Create a heatmap of missing values using pheatmap
missing_matrix <- is.na(caroldata)
pheatmap(missing_matrix, cex_col = 0.5, cex_row = 0.5)



# Visualize missing data using the aggr function
aggr(caroldata, col = c('navajowhite3', 'red'), numbers = TRUE, sortVars = TRUE, labels = names(caroldata), cex.axis = 0.7, gap = 3)



# Create a matrix indicating missing values (TRUE for missing, FALSE for present)
missing_matrix <- is.na(caroldata)


# Use heatmap to visualize missing data
heatmap(missing_matrix, cexCol = 0.5, cexRow = 0.5, col = c("white", "red"))


#by expgrp ------------
ggplot(pca_pc1_pc2, aes(x = PC1 , y = PC2, colour = exp_group, label = sample)) +
    geom_point(size = 3) +
    #geom_text(fontface = "bold", show.legend = F, vjust = 0.5, hjust = -0.5, size = 3) +
    xlab(paste("PC1 - ", pca_var_perc[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca_var_perc[2], "%", sep = "")) +
    #xlim(-50, 50) + ylim(-50, 50) +
    ggtitle("by_exp_group_PCA_noNA_qn") +
    theme_minimal()

#by tissue and sampleprep_date ----------
# Convert sampleprep_date to a factor
pca_pc1_pc2$sampleprep_date <- factor(pca_pc1_pc2$sampleprep_date)

ggplot(pca_pc1_pc2, aes(x = PC1 , y = PC2, shape = sampleprep_date, colour = tissue, label = sample)) +
    geom_point(size = 3) +
    scale_shape_manual(values = c("1" = 2, "2" = 17)) +
    #geom_text(fontface = "bold", show.legend = F, vjust = 0.5, hjust = -0.5, size = 3) +
    xlab(paste("PC1 - ", pca_var_perc[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca_var_perc[2], "%", sep = "")) +
    #xlim(-50, 50) + ylim(-50, 50) +
    ggtitle("by_sampleprep_date_PCA_noNA_qn") +
    theme_minimal()

#by area and sampleprep_occ --------------
pca_pc1_pc2$sampleprep_occ <- factor(pca_pc1_pc2$sampleprep_date)

ggplot(pca_pc1_pc2, aes(x = PC1 , y = PC2, shape = sampleprep_occ, colour = tissue, label = sample)) +
    geom_point(size = 3) +
    scale_shape_manual(values = c("1" = 1, "2" = 16, "3" = 10)) +
    #geom_text(fontface = "bold", show.legend = F, vjust = 0.5, hjust = -0.5, size = 3) +
    xlab(paste("PC1 - ", pca_var_perc[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca_var_perc[2], "%", sep = "")) +
    #xlim(-50, 50) + ylim(-50, 50) +
    ggtitle("by_occ_PCA_noNA_qn") +
    theme_minimal()

#by cryo_date ----------------
pca_pc1_pc2$cryo_date <- factor(pca_pc1_pc2$sampleprep_date)

ggplot(pca_pc1_pc2, aes(x = PC1 , y = PC2, colour = cryo_date, label = sample)) +
    geom_point(size = 3) +
    #geom_text(fontface = "bold", show.legend = F, vjust = 0.5, hjust = -0.5, size = 3) +
    xlab(paste("PC1 - ", pca_var_perc[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca_var_perc[2], "%", sep = "")) +
    #xlim(-50, 50) + ylim(-50, 50) +
    ggtitle("cryo_date_PCA_noNA_qn") +
    theme_minimal()

#contributions ----------
contr <- pca_nona_qn$rotation
contr_absolute <- abs(contr)

# 10 highest contributions on diml
dim1_10highest=row.names(contr_absolute[order(contr_absolute[,1], decreasing = T),])[1:10]
dim1_10highest

# 10 highest on dim2
dim2_10highest=row.names(contr_absolute[order(contr_absolute[,2], decreasing = T),])[1:10]
dim2_10highest

fviz_pca_biplot(pca_nona_qn,
                col.var = "contrib", # Color by contributions to the PC
                select.var = list(name=c(dim1_10highest, dim2_10highest)),
                #habillage = metadata$area,
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = T, # Avoid text overlapping
                geom.ind = c("point"),
                pointsize = 2,
                addEllipses = F,
                #xlim = c(-50, 50),
                #ylim = c(-50, 50),
)

#MTJ ------------
metadata_mtj <- metadata[metadata$tissue == "MTJ", ]
mtj_order = metadata_mtj$sample_id[order(metadata_mtj$intervention)]
mtj_order = trimws(mtj_order)
df_nona_qn_mtj = df_nona_qn[ ,mtj_order]
pca_nona_qn_mtj <- prcomp(t(df_nona_qn_mtj), scale = T)
summary(pca_nona_qn_mtj)

fviz_pca_ind(pca_nona_qn_mtj)

plot(pca_nona_qn_mtj$x[,1], pca_nona_qn_mtj$x[,2])

pca_var <- pca_nona_qn_mtj$sdev^2
pca_var_perc <- round(pca_var/sum(pca_var)*100, digits = 1) #percentage of each component
fviz_eig(pca_nona_qn_mtj, addlabels = T)

pca_results <- as.data.frame(pca_nona_qn_mtj$x)
pca_mtj_pc1_pc2 <- pca_results %>%
    dplyr::select(PC1, PC2)

pca_mtj_pc1_pc2 <- pca_mtj_pc1_pc2 %>%
    add_column(sample = rownames(pca_results),
               exp_group = metadata_mtj$expgrp,
               tissue = metadata_mtj$tissue,
               intervention = metadata_mtj$intervention,
               duration = metadata_mtj$time,
               cryo_date = metadata_mtj$cryostat_date,
               lcm_date = metadata_mtj$lcm_date,
               sampleprep_occ = metadata_mtj$sampleprep_occ,
               sampleprep_date = metadata_mtj$sampleprep_date)

#by intervention
ggplot(pca_mtj_pc1_pc2, aes(x = PC1 , y = PC2, colour = intervention, label = sample)) +
    geom_point(size = 3) +
    #geom_text(fontface = "bold", show.legend = F, vjust = 0.5, hjust = -0.5, size = 3) +
    xlab(paste("PC1 - ", pca_var_perc[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca_var_perc[2], "%", sep = "")) +
    #xlim(-50, 50) + ylim(-50, 50) +
    ggtitle("intervention_PCA_noNA_qn_MTJ") +
    theme_minimal()


#by duration and intervention
pca_mtj_pc1_pc2$intervention <- as.factor(pca_mtj_pc1_pc2$intervention)
pca_mtj_pc1_pc2$duration <- as.factor(pca_mtj_pc1_pc2$duration)

unique(pca_mtj_pc1_pc2$intervention)
unique(pca_mtj_pc1_pc2$duration)

ggplot(pca_mtj_pc1_pc2, aes(x = PC1, y = PC2, colour = intervention, shape = duration, label = sample)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("red", "green", "blue")) +
    scale_shape_manual(values = c(0, 15, 12)) +
    xlab(paste("PC1 - ", pca_var_perc[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca_var_perc[2], "%", sep = "")) +
    ggtitle("caro_lcm_spect_PCA_noNA_qn_MTJ") +
    theme_minimal()

#by expgrp
# Plotting by exp_group
ggplot(pca_mtj_pc1_pc2, aes(x = PC1, y = PC2, colour = exp_group, shape = duration, label = sample)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("cc" = "red", "ecc" = "green", "sed" = "blue")) +
    scale_shape_manual(values = c("0" = 0, "5" = 15, "10" = 12)) +
    xlab(paste("PC1 - ", pca_var_perc[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca_var_perc[2], "%", sep = "")) +
    ggtitle("caro_lcm_spect_PCA_noNA_qn_MTJ") +
    theme_minimal()

# Plotting by expgrp (assuming you meant exp_group)
ggplot(pca_mtj_pc1_pc2, aes(x = PC1, y = PC2, colour = exp_group, label = sample)) +
    geom_point(size = 3) +
    ggtitle("caro_lcm_spect_PCA_noNA_qn_MTJ by exp_group") +
    theme_minimal()


# Specify shapes manually based on the number of levels in exp_group
custom_shapes <- seq(0, 15, length.out = length(unique(pca_mtj_pc1_pc2$exp_group)))

# Plotting with manual shape specification
ggplot(pca_mtj_pc1_pc2, aes(x = PC1, y = PC2, shape = exp_group, label = sample)) +
    geom_point(size = 3) +
    scale_shape_manual(values = custom_shapes) +
    ggtitle("caro_lcm_spect_PCA_noNA_qn_MTJ by exp_group") +
    theme_minimal()

#by tissue and sampleprep_date
class(pca_mtj_pc1_pc2$sampleprep_date)
unique(pca_mtj_pc1_pc2$sampleprep_date)
pca_mtj_pc1_pc2_numeric_date <- pca_mtj_pc1_pc2[!is.na(as.numeric(pca_mtj_pc1_pc2$sampleprep_date)), ]

ggplot(pca_mtj_pc1_pc2, aes(x = PC1, y = PC2, shape = factor(sampleprep_date), colour = intervention, label = sample)) +
    geom_point(size = 3) +
    scale_shape_manual(values = c("1" = 2, "2" = 17)) +
    xlab(paste("PC1 - ", pca_var_perc[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca_var_perc[2], "%", sep = "")) +
    ggtitle("caro_lcm_spect_PCA_noNA_qn_MTJ") +
    theme_minimal()

