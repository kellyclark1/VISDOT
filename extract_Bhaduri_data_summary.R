library(Seurat)
library(dplyr)
library(ggplot2)
library(gprofiler2)
library(gridExtra)


print("reading RDS")

#read in R object, replace with path to file
#this step takes several minutes
#R object from Bhaduri et al can be accessed at:
#     https://cellxgene.cziscience.com/collections/c8565c6a-01a1-435b-a549-f11b452a83a8
d <- readRDS("local.rds")

print("done")


#assign levels to order time points
f <- factor(d$development_stage, levels =c("14th week post-fertilization human stage",
"16th week post-fertilization human stage",
"17th week post-fertilization human stage", 
"18th week post-fertilization human stage",
"19th week post-fertilization human stage",
"20th week post-fertilization human stage",
"22nd week post-fertilization human stage",
"25th week post-fertilization human stage"
))
d$development_stage <- f


#vector of cell types
cells <-  c("cerebral cortex GABAergic interneuron",
					"glutamatergic neuron",
					"progenitor cell",
					"Cajal-Retzius cell",
					"forebrain radial glia cell",
					"oligodendrocyte precursor cell",
					"microglial cell",
					"native cell",
					"blood vessel endothelial cell",
					"cerebral cortex endothelial cell")

#split Seurat object into a list of 10 Seurat objects, one for each cell type in the data set
#this step takes several minutes
print("split object")
obj.list <- SplitObject(d, split.by = "cell_type")
print("done")




#create vector for y axis labels for plots
ylabs <- c("local_25th week post-fertilization human stage",
"local_22nd week post-fertilization human stage",
"local_20th week post-fertilization human stage",
"local_19th week post-fertilization human stage",
"local_18th week post-fertilization human stage",
"local_17th week post-fertilization human stage",
"local_16th week post-fertilization human stage",
"local_14th week post-fertilization human stage"
)
ylabs <- rev(ylabs)


#function to convert gene name
get.sym <- function(genes){
return(gconvert(genes,organism="hsapiens",target="ENTREZGENE",filter_na = T, mthreshold = 1)$name)
}




#function to create dot plot for a gene,
#with plot data containing the average and percent expression
#function returns a plot, and plot data is used to create table of data summary
plotGene <- function(h){

p1 <- DotPlot(object = obj.list[[1]], group.by = "development_stage", features = h) +
	ggtitle("progenitor cell") + 
	scale_y_discrete(labels = f) + coord_flip() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p1$data$features.plot <- "progenitor cell"

p2 <- DotPlot(object = obj.list[[2]], group.by = "development_stage", features = h) +
	ggtitle("cerebral cortex GABAergic interneuron") + 
	scale_y_discrete(labels = f) + coord_flip() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p2$data$features.plot <- "cerebral cortex GABAergic interneuron"

p3 <- DotPlot(object = obj.list[[3]], group.by = "development_stage", features = h) +
	ggtitle("native cell") +  scale_x_discrete(labels = "native cell") +
	scale_y_discrete(labels = f) + coord_flip() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p3$data$features.plot <- "native cell"

p4 <- DotPlot(object = obj.list[[4]], group.by = "development_stage", features = h) +
	ggtitle("oligodendrocyte precursor cell") +  
	scale_y_discrete(labels = f) + coord_flip() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p4$data$features.plot <- "oligodendrocyte precursor cell"

p5 <- DotPlot(object = obj.list[[5]], group.by = "development_stage", features = h) +
	ggtitle("forebrain radial glia cell") + 
	scale_y_discrete(labels = f) + coord_flip() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p5$data$features.plot <- "forebrain radial glia cell"

p6 <- DotPlot(object = obj.list[[6]], group.by = "development_stage", features = h) +
	ggtitle("microglial cell") + 
	scale_y_discrete(labels = f) + coord_flip() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p6$data$features.plot <- "microglial cell"

p7 <- DotPlot(object = obj.list[[7]], group.by = "development_stage", features = h) +
	ggtitle("blood vessel endothelial cell") + 
	scale_y_discrete(labels = f) + coord_flip() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p7$data$features.plot <- "blood vessel endothelial cell"

p8 <- DotPlot(object = obj.list[[8]], group.by = "development_stage", features = h) +
	ggtitle("Cajal-Retzius cell") + 
	scale_y_discrete(labels = f) + coord_flip() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p8$data$features.plot <- "Cajal-Retzius cell"

p9 <- DotPlot(object = obj.list[[9]], group.by = "development_stage", features = h) +
	ggtitle("glutamatergic neuron") +
	scale_y_discrete(labels = f) + coord_flip() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p9$data$features.plot <- "glutamatergic neuron"

p10 <- DotPlot(object = obj.list[[10]], group.by = "development_stage", features = h) +
	ggtitle("cerebral cortex endothelial cell") +
	scale_y_discrete(labels = f) + coord_flip() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p10$data$features.plot <- "cerebral cortex endothelial cell"

p <- p1
p$data <- rbind(p1$data, p2$data, p3$data, p4$data, p5$data,
		p6$data, p7$data, p8$data, p9$data, p10$data)
p$data$features.plot <- factor(p$data$features.plot, levels = cells)			

p

}







#get vector of all gene IDs in dataset
all.ensg <- rownames(d)



#for each gene ID in the dataset, get plot data and append to full table of data summary
all.data <- data.frame()
for (x in 1:length(all.ensg)){
	p <- plotGene(all.ensg[[x]])
	p$data$gene <- all.ensg[[x]]
	p$data$sym <- get.sym(all.ensg[[x]])
	all.data <- rbind(all.data, p$data)
}



#write data summary to file
write.csv(quote = F, all.data, file = "all_genes_second_trimester_sc_data.csv")






