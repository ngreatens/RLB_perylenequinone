library(ggplot2)
library(gggenes)
library(dplyr)
library(RColorBrewer)

setwd("/project/fdwsru_fungal/Nick/rlb_light/10_gggenes")

#read into table
genes <- read.csv("genes3.txt", sep = "\t")

########################################################################
# reverse order of genes if on opposite strand to ensure gene arrows are pointed in the opposite direction
# and make sure that the genes that are used for alignemnt are all oriented in the same was 
# through making negative the position values for some of them 


#reverse order of genes if on opposite strand
for (i in 1:nrow(genes)){
  start <- genes[i,2]
  end <- genes[i,3]
  strand <- genes[i,4]
  if (strand == 0){
    genes[i,2] <- end
    genes[i,3] <- start
  }
}
arr <- c()
for (i in 1:nrow(genes)){
  if (genes[i,5] == "g6609.t1"){
      if (genes[i,4] == 1){
        cluster <- genes[i,1]
        arr <- append(arr,cluster)
        }
    }
}
# keep only unique values
arr <- unique(arr)
# make position negative for all clusters in arr. this reverses the genes so they can be aligned later
for (i in 1:length(arr)){
  for (j in 1:nrow(genes)){
    if (arr[i] == genes[j,1]){
      genes[j,2] <- genes[j,2]*(-1)
      genes[j,3] <- genes[j,3]*(-1)
    }
  }
}

# make dummy table to orient genes. per tutorial
dummies <- make_alignment_dummies(
  genes,
  aes(xmin = start, xmax = end, y = molecule, id = prot),
  on = "g6609.t1"
)

######################################################################################################


## get plot of just the genes and save to png.
## was having some problems wiht the viewer in RStudio, so ended up just saving the plot to a file
## every time. 

gene_table<-read.csv(("gene_cluster_functions.csv"), header = TRUE)

# manually set colors. This is the 12 from Set3 plus a brown and white
# unfortunately not colorblind friendly, but not sure how to make it so with so many genes.
# at least can have a table to supplement findings.
# have everything in github and supps, as needed

my_colors=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F", "#cc9963","#FFFFFF")

#get plot
p <- ggplot(genes, aes(xmin = start, xmax = end, y = molecule, fill = prot)) +
  geom_gene_arrow() +
  geom_blank(data = dummies) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  theme_genes() +
  theme(
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12, face = "italic"),
    #legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.title = element_blank()
  )+
  scale_fill_manual(labels = gene_table$Function, values= my_colors)

#save to png
png("example.png", width = 10, height = 4 ,units="in", res = 300)
p
dev.off()

############################################################################################

## now get out gene clusters on tree using ggtree

ggtree(tree, branch.length='none') + 
  geom_tiplab() + xlim_tree(5.5) + 
  geom_facet(mapping = aes(xmin = start, xmax = end, fill = gene),
             data = example_genes, geom = geom_motif, panel = 'Alignment',
             on = 'genE', label = 'gene', align = 'left') +
  scale_fill_brewer(palette = "Set3") + 
  scale_x_continuous(expand=c(0,0)) +
  theme(strip.text=element_blank(),
        panel.spacing=unit(0, 'cm'),
        legend.title = element_blank())

library(ggtree)
library(treeio)
library(glue)
library(dplyr)

tree<- read.iqtree("tree.tree")
df <- read.csv("tree_data.csv", header = TRUE)


#use glue to make new label
df <- dplyr::mutate(df, lab = glue("italic({Genus})~italic({Species})~{nonitalic}~{Assembly_Accession}"))



#get bootstrap values
q <- ggtree(tree)
d <- q$data
#select internal nodes
d <- d[!d$isTip,]


p <- ggtree(tree)  %<+% df +
  geom_tiplab(aes(label = lab), align = TRUE, parse = TRUE, size = 3) + 
  xlim_tree(4)
  #geom_text(aes(label=node))
  #geom_text2(data = d, aes(label = UFboot)) +
p <- ggtree::rotate(p, 16)
p <- ggtree::rotate(p, 18)
# p <- p +geom_text(aes(label=UFboot), size = 2, hjust = -.2, vjust = 1)
p <- p + geom_facet(mapping = aes(xmin = start, xmax = end, fill = prot),
           data = genes, 
           geom = geom_motif, 
           panel = 'Alignment',
           on = "g6609.t1",
           align = 'left') +
  theme(strip.text=element_blank(),
        panel.spacing=unit(0, 'cm'),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=8)) +
  scale_fill_manual(labels = gene_table$Function, values= my_colors) +
  guides(fill=guide_legend(nrow=4 ,byrow=TRUE))

# adjust widths of facets e.g. tree: genes plot is 1:2 
# see: https://yulab-smu.top/treedata-book/chapter12.html

facet_widths(p, widths = c(1, 2))





