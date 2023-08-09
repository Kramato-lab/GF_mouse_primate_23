library(DESeq2)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/GF_Mouse_Human_v_Primate/GF_Mouse_Human/First_experiment/Final_analysis")

#read in raw count data and metadata with columns of first matching rows of second
countdata<-read.csv("Liver_gene_expression_rawcounts.csv", row.names = "ensgene", check.names = FALSE)
metadata<-read.csv("Gene_expression_metadata.csv", row.names = 1)

#set up data set
dds<-DESeqDataSetFromMatrix(countData = countdata, colData = metadata, design = ~ Species)

#set reference
dds$Species <- relevel(dds$Species, ref = "Human")

#filter out low count genes (rows with at least 10 reads total)
keep<-rowSums(counts(dds)) >= 10
dds<-dds[keep,]

#run DESeq2
dds<-DESeq(dds)

#get results
res<-results (dds, contrast=c("Species", "Squirrel_monkey", "Human"))
res.sig<-subset(res, padj<0.05)
write.csv(res.sig, file="Adult_SQM_Human.csv")

res<-results (dds, contrast=c("Species", "Macaque", "Human"))
res.sig<-subset(res, padj<0.05)
write.csv(res.sig, file="Adult_Macaque_Human.csv")

#large vs small brain
#read in raw count data and metadata with columns of first matching rows of second
countdata<-read.csv("Liver_gene_expression_rawcounts.csv", row.names = "ensgene", check.names = FALSE)
metadata<-read.csv("Liver_Gene_expression_metadata.csv", row.names = 1)


#set up data set
dds<-DESeqDataSetFromMatrix(countData = countdata, colData = metadata, design = ~ Group)

#set reference
dds$Group <- relevel(dds$Group, ref = "Small")

#filter out low count genes (rows with at least 10 reads total)
keep<-rowSums(counts(dds)) >= 10
dds<-dds[keep,]

#run DESeq2
dds<-DESeq(dds)

#get results
res<-results (dds)
res.sig<-subset(res, padj<0.05)
write.csv(res.sig, file="Large_v_small.csv")


##not sure i need the below code
#output transformed data for plots (for one gene)
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="Species", 
                returnData=TRUE)
write.csv(d, file="transformed_data_to_plot_species.csv")

#output transformed data for plots
vsd <- vst(dds, blind=FALSE)

write.csv(as.data.frame(assay(vsd)), file="transformed_data_to_plot2.csv")
