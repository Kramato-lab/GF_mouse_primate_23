library(DESeq2)
library(RITAN)
library(RITANdata)
library(ggplot2)
library(stringr)
library(biomaRt)

# Set working directory
setwd("/Users/derek/Desktop/Katie/")

#read in raw count data and metadata with columns of first matching rows of second
countdata<-read.csv("Liver_gene_expression_rawcounts.csv", row.names = "ensgene", check.names = FALSE)
metadata<-read.csv("Liver_Gene_expression_metadata.csv", row.names = 1)

# Map ENSEMBL IDs to Gene Symbols for Downstream Analyses
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
gene_names <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=rownames(countdata),mart= mart)
rownames(gene_names) <- gene_names$ensembl_gene_id


#species only
#set up data set
dds<-DESeqDataSetFromMatrix(countData = countdata, colData = metadata, design = ~ Species)

#filter out low count genes (rows with at least 10 reads total)
keep<-rowSums(counts(dds)) >= 10
dds<-dds[keep,]

#run DESeq2
dds<-DESeq(dds)

#get results
res<-results (dds, contrast=c("Species", "Macaque", "Squirrel_monkey"))
res2<-results (dds, contrast=c("Species", "Macaque", "Human"))
res3<-results (dds, contrast=c("Species", "Squirrel_monkey", "Human"))
res.sig<-subset(res, padj<0.05)
res.sig2<-subset(res2, padj<0.05)
res.sig3<-subset(res3, padj<0.05)
write.csv(res.sig, file="Macaque_SQM.csv")
write.csv(res.sig2, file="Macaque_Human.csv")
write.csv(res.sig3, file="SQM_Human.csv")


#output transformed data for plots
vsd <- vst(dds, blind=FALSE)

write.csv(as.data.frame(assay(vsd)), file="transformed_data.csv")




# Perform functional analysis of DE genes
resources <- c("ReactomePathways")

# FC threshold
up_sites <- res.sig[res.sig$log2FoldChange > 1,]
down_sites <- res.sig[res.sig$log2FoldChange < 1,]

# Order based on p-value
up_sites <- up_sites[order(up_sites$pvalue, decreasing = F),]
down_sites <- down_sites[order(down_sites$pvalue, decreasing = F),]

# Convert to gene names
gene_names <- gene_names[is.na(gene_names$mgi_symbol)==F,]

up_sites <- up_sites[intersect(rownames(up_sites), rownames(gene_names)),]
down_sites <- down_sites[intersect(rownames(down_sites), rownames(gene_names)),]

rownames(up_sites) <- gene_names[rownames(up_sites),]$mgi_symbol
rownames(down_sites) <- gene_names[rownames(down_sites),]$mgi_symbol


# Enrichment
enr_up <- term_enrichment(str_to_upper(rownames(up_sites)), resources = resources, all_symbols = cached_coding_genes)
enr_down <- term_enrichment(str_to_upper(rownames(down_sites)), resources = resources, all_symbols = cached_coding_genes)

# Trim pathway names
enr_up$name <- str_to_title(tolower(gsub("_", " ", substr(enr_up$name, 18, 100))))
enr_up$name <- factor(enr_up$name, levels=enr_up[order(enr_up$q, decreasing = T), "name"])

enr_down$name <- str_to_title(tolower(gsub("_", " ", substr(enr_down$name, 18, 100))))
enr_down$name <- factor(enr_down$name, levels=enr_down[order(enr_down$q, decreasing = T), "name"])

png("Images/HumanAdult_vs_HumanInfant/up_in_infant.png")
ggplot(head(enr_up[order(enr_up$q),], n=20), aes(x = name, y = -1 * log10(q))) +
  geom_bar(stat="identity") +
  coord_flip() +
  xlab("Annotation") +
  ylab("Negative Log Adj. P-value") + 
  ggtitle("Sites Up-Regulated in Infant")
dev.off()

png("Images/HumanAdult_vs_HumanInfant/down_in_infant.png")
ggplot(head(enr_down[order(enr_down$q),], n=20), aes(x = name, y = -1 * log10(q))) +
  geom_bar(stat="identity") +
  coord_flip() +
  xlab("Annotation") +
  ylab("Negative Log Adj. P-value") + 
  ggtitle("Sites Down-Regulated in Infant")
dev.off()

#MiMeNet Analysis required for below
'''
c0 <- read.csv("cluster_0.csv")$X0
c1 <- read.csv("cluster_1.csv")$X0
c2 <- read.csv("cluster_2.csv")$X0
c3 <- read.csv("cluster_3.csv")$X0
c4 <- read.csv("cluster_4.csv")$X0
c5 <- read.csv("cluster_5.csv")$X0

c0_enrichment <-  term_enrichment(str_to_upper(c0), resources = resources, all_symbols = cached_coding_genes)
c0_enrichment$name <- str_to_title(tolower(gsub("_", " ", substr(c0_enrichment$name, 18, 100))))
c0_enrichment$name <- factor(c0_enrichment$name, levels=c0_enrichment[order(c0_enrichment$q, decreasing = T), "name"])

ggplot(head(c0_enrichment[order(c0_enrichment$q),], n=20), aes(x = name, y = -1 * log10(q))) +
  geom_bar(stat="identity") +
  coord_flip() +
  xlab("Annotation") +
  ylab("Negative Log Adj. P-value") + 
  ggtitle("Sites Up-Regulated in C0")


c1_enrichment <-  term_enrichment(str_to_upper(c1), resources = resources, all_symbols = cached_coding_genes)
c1_enrichment$name <- str_to_title(tolower(gsub("_", " ", substr(c1_enrichment$name, 18, 100))))
c1_enrichment$name <- factor(c1_enrichment$name, levels=c1_enrichment[order(c1_enrichment$q, decreasing = T), "name"])

ggplot(head(c1_enrichment[order(c1_enrichment$q),], n=20), aes(x = name, y = -1 * log10(q))) +
  geom_bar(stat="identity") +
  coord_flip() +
  xlab("Annotation") +
  ylab("Negative Log Adj. P-value") + 
  ggtitle("Sites Up-Regulated in C1")



c2_enrichment <-  term_enrichment(str_to_upper(c2), resources = resources, all_symbols = cached_coding_genes)
c2_enrichment$name <- str_to_title(tolower(gsub("_", " ", substr(c2_enrichment$name, 18, 100))))
c2_enrichment$name <- factor(c2_enrichment$name, levels=c2_enrichment[order(c2_enrichment$q, decreasing = T), "name"])

ggplot(head(c2_enrichment[order(c2_enrichment$q),], n=20), aes(x = name, y = -1 * log10(q))) +
  geom_bar(stat="identity") +
  coord_flip() +
  xlab("Annotation") +
  ylab("Negative Log Adj. P-value") + 
  ggtitle("Sites Up-Regulated in C2")



c3_enrichment <-  term_enrichment(str_to_upper(c3), resources = resources, all_symbols = cached_coding_genes)
c3_enrichment$name <- str_to_title(tolower(gsub("_", " ", substr(c3_enrichment$name, 18, 100))))
c3_enrichment$name <- factor(c3_enrichment$name, levels=c3_enrichment[order(c3_enrichment$q, decreasing = T), "name"])

ggplot(head(c3_enrichment[order(c3_enrichment$q),], n=20), aes(x = name, y = -1 * log10(q))) +
  geom_bar(stat="identity") +
  coord_flip() +
  xlab("Annotation") +
  ylab("Negative Log Adj. P-value") + 
  ggtitle("Sites Up-Regulated in C3")





c4_enrichment <-  term_enrichment(str_to_upper(c4), resources = resources, all_symbols = cached_coding_genes)
c4_enrichment$name <- str_to_title(tolower(gsub("_", " ", substr(c4_enrichment$name, 18, 100))))
c4_enrichment$name <- factor(c4_enrichment$name, levels=c4_enrichment[order(c4_enrichment$q, decreasing = T), "name"])

ggplot(head(c4_enrichment[order(c4_enrichment$q),], n=20), aes(x = name, y = -1 * log10(q))) +
  geom_bar(stat="identity") +
  coord_flip() +
  xlab("Annotation") +
  ylab("Negative Log Adj. P-value") + 
  ggtitle("Sites Up-Regulated in C4")



c5_enrichment <-  term_enrichment(str_to_upper(c5), resources = resources, all_symbols = cached_coding_genes)
c5_enrichment$name <- str_to_title(tolower(gsub("_", " ", substr(c5_enrichment$name, 18, 100))))
c5_enrichment$name <- factor(c5_enrichment$name, levels=c5_enrichment[order(c5_enrichment$q, decreasing = T), "name"])

ggplot(head(c5_enrichment[order(c5_enrichment$q),], n=20), aes(x = name, y = -1 * log10(q))) +
  geom_bar(stat="identity") +
  coord_flip() +
  xlab("Annotation") +
  ylab("Negative Log Adj. P-value") + 
  ggtitle("Sites Up-Regulated in C0")

cluster_enrichment <-term_enrichment_by_subset(list("Cluster 0/Cluster 3"=str_to_upper(c(c0,c3)),"Cluster 1"=str_to_upper(c1),
                                                    "Cluster 2"=str_to_upper(c2), "Cluster 4"=str_to_upper(c4),
                                                    "Cluseter 5"=str_to_upper(c5)),
                                               resources = resources, all_symbols = cached_coding_genes)
cluster_enrichment <- cluster_enrichment[cluster_enrichment$name != "ReactomePathways.Metabolism", ]

png("MiMeNet_Reactome.png")
plot(cluster_enrichment, label_size_y = 7, label_size_x = 7 )
dev.off()
'''

# Subset for up regulated in pairwise compairsons
up_in_m_v_s <- res.sig[res.sig$log2FoldChange >= 1,]
up_in_s_v_m <- res.sig[res.sig$log2FoldChange <= 1,]
up_in_m_v_h <- res.sig2[res.sig2$log2FoldChange >= 1,]
up_in_h_v_m <- res.sig2[res.sig2$log2FoldChange <= 1,]
up_in_s_v_h <- res.sig3[res.sig3$log2FoldChange >= 1,]
up_in_h_v_s <- res.sig3[res.sig3$log2FoldChange <= 1,]

# Joint enrichment analysis of all gene sets

up_in_m_v_s <- up_in_m_v_s[order(up_in_m_v_s$pvalue, decreasing = F),]
up_in_m_v_s <- up_in_m_v_s[intersect(rownames(up_in_m_v_s), rownames(gene_names)),]
rownames(up_in_m_v_s) <- gene_names[rownames(up_in_m_v_s),]$mgi_symbol

up_in_s_v_m <- up_in_s_v_m[order(up_in_s_v_m$pvalue, decreasing = F),]
up_in_s_v_m <- up_in_s_v_m[intersect(rownames(up_in_s_v_m), rownames(gene_names)),]
rownames(up_in_s_v_m) <- gene_names[rownames(up_in_s_v_m),]$mgi_symbol

up_in_m_v_h <- up_in_m_v_h[order(up_in_m_v_h$pvalue, decreasing = F),]
up_in_m_v_h <- up_in_m_v_h[intersect(rownames(up_in_m_v_h), rownames(gene_names)),]
rownames(up_in_m_v_h) <- gene_names[rownames(up_in_m_v_h),]$mgi_symbol

up_in_h_v_m <- up_in_h_v_m[order(up_in_h_v_m$pvalue, decreasing = F),]
up_in_h_v_m <- up_in_h_v_m[intersect(rownames(up_in_h_v_m), rownames(gene_names)),]
rownames(up_in_h_v_m) <- gene_names[rownames(up_in_h_v_m),]$mgi_symbol

up_in_s_v_h <- up_in_s_v_h[order(up_in_s_v_h$pvalue, decreasing = F),]
up_in_s_v_h <- up_in_s_v_h[intersect(rownames(up_in_s_v_h), rownames(gene_names)),]
rownames(up_in_s_v_h) <- gene_names[rownames(up_in_s_v_h),]$mgi_symbol

up_in_h_v_s <- up_in_h_v_s[order(up_in_h_v_s$pvalue, decreasing = F),]
up_in_h_v_s <- up_in_h_v_s[intersect(rownames(up_in_h_v_s), rownames(gene_names)),]
rownames(up_in_h_v_s) <- gene_names[rownames(up_in_h_v_s),]$mgi_symbol

cluster_enrichment <-term_enrichment_by_subset(list("Up in Macaque vs SM"=str_to_upper(rownames(up_in_m_v_s)),
                                                    "Up in Macaque vs Human"=str_to_upper(rownames(up_in_m_v_h)),
                                                    "Up in SM vs Macaue"=str_to_upper(rownames(up_in_s_v_m)),
                                                    "Up in SM vs Human"=str_to_upper(rownames(up_in_s_v_h)),
                                                    "Up in Human vs Macaque"=str_to_upper(rownames(up_in_h_v_m)),
                                                    "Up in Human vs SM"=str_to_upper(rownames(up_in_h_v_s))), 
                                               q_value_threshold = 1e-4, resources = resources, all_symbols = cached_coding_genes)
cluster_enrichment <- cluster_enrichment[cluster_enrichment$name != "ReactomePathways.Metabolism", ]

png("Reactome_DEGenes.png", res = 300, width = 10, height = 10, units = "in")
par(mar = c(0, 0, 0, 0))
plot(cluster_enrichment, label_size_y = 12, label_size_x = 12, grid_line_color = 'black')
dev.off()

