library(ggplot2)
library(gggenes)
library(dplyr)
library(RColorBrewer)

setwd("/project/fdwsru_fungal/Nick/rlb_light/13_gggenes_second_go")

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
  if (genes[i,5] == "a6534.t1"){
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
  on = "a6534.t1"
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

my_colors=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#CACACA","#BC80BD","#CCEBC5","#FFED6F", "#cc9963","#787878", "#FFFFFF")

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

library(ggtree)
library(treeio)
library(glue)
library(dplyr)

tree<- read.iqtree("tree.tree")
df <- read.csv("tree_data.csv", header = TRUE)


#use glue to make new label
df <- dplyr::mutate(df, lab1 = glue("italic({Genus})~italic({Species})~{nonitalic}~{Prot_ID}"))
df <- dplyr::mutate(df, lab2 = glue("{Cluster_coords}"))
df <- dplyr::mutate(df, lab3 = glue("{Phaeo_label}"))

#get bootstrap values
q <- ggtree(tree)
d <- q$data
#select internal nodes
d <- d[!d$isTip,]


p <- ggtree(tree)  %<+% df +
  geom_tiplab(aes(label = lab1), parse = TRUE, size = 3, nudge_y = .1) +
  geom_tiplab(aes(label = lab2), parse = TRUE, size = 2.75, nudge_y = -.3) + 
  geom_tiplab(aes(label = lab3), parse = TRUE, size = 3, nudge_y = .1) +
  xlim_tree(4.5)
  #geom_text(aes(label=node))
  #geom_label(data = d, aes(label = UFboot), size = 2.75, hjust = .2, vjust = .2)
p <- ggtree::rotate(p, 20)
p <- ggtree::rotate(p, 22)
# p <- p +geom_text(aes(label=UFboot), size = 2, hjust = -.2, vjust = 1)
p <- p + geom_facet(mapping = aes(xmin = start, xmax = end, fill = prot),
           data = genes, 
           geom = geom_motif, 
           panel = 'Alignment',
           on = "a6534.t1",
           align = 'left') +
  scale_fill_manual(labels = gene_table$Function, values= my_colors) +
  guides(fill=guide_legend(nrow=5 ,byrow=TRUE)) +
  theme(strip.text=element_blank(),
        panel.spacing=unit(0, 'cm'),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=8.75),
        legend.key.size = unit(.4, 'cm'))

# adjust widths of facets e.g. tree: genes plot is 1:2 
# see: https://yulab-smu.top/treedata-book/chapter12.html

p2<- facet_widths(p, widths = c(4, 7))


########################################

df <- read.csv("/project/fdwsru_fungal/Nick/rlb_light/12_expression_chart/perylenequninone_cluster_espression.csv")

blank_df <- data.frame(matrix(nrow=128, ncol = 3.5))
colnames(blank_df) <- c("gene", "treatment", "mapped")
row=1
for (i in 1:nrow((df))){
  for (j in 1:ncol(df)){
    if (endsWith(colnames(df)[j],"light")){
      blank_df[row,1] <- df[i,1]
      blank_df[row,2] <- "light"
      blank_df[row,3] <- df[i,j]
      row <- row + 1
    }
    else if (endsWith(colnames(df)[j], "dark")){
      blank_df[row,1] <- df[i,1]
      blank_df[row,2] <- "dark"
      blank_df[row,3] <- df[i,j]
      row <- row + 1
    }
  }
}

library(ggplot2)
library(dplyr)

mean_df <- blank_df %>% group_by(gene, treatment) %>% summarise(mean = mean(mapped), se = sd(mapped)/sqrt(n()))
mean_df$pvalue=1


light_df <- mean_df[mean_df$treatment == 'light',]
light_df$color <- c("a","b","c","d","e","f","g","h")
dark_df <- mean_df[mean_df$treatment == 'dark',]
dark_df$color <- c("i","j","k","l","m","n","o","p")


light_df$pvalue <- c("*","*","*","*","*","*","","")


library(colorspace) 


# use array of colors from above, part of Set3 from the Rcolorbrewer palettes
# then darken them for the dark conditions

my_colors=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5")
my_colors <- append(my_colors, darken(my_colors, .5))

p3<- ggplot(mean_df, aes(x = gene, y = mean)) + 
  geom_col(data = light_df, aes(fill = color), position = "dodge", just = 1.05,  width = .4, linewidth = .5) +
  geom_errorbar(data = light_df, aes(ymin = mean-se, ymax = mean+se), position = position_nudge(-.225), width = 0.2) +
  geom_col(data = dark_df, aes(fill = color), position = "dodge", just = -.05,  width = .4, linewidth = .5) +
  geom_errorbar(data = dark_df, aes(ymin = mean-se, ymax = mean+se), position = position_nudge(.225), width = .2) +
  geom_label(aes(x = "g6532", y = 14000, label = "Light", angle = 45), size = 3.5, label.size = 0, nudge_x = -.225, ) +
  geom_label(aes(x = "g6532", y = 14000, label = "Dark", angle = 45), size = 3.5, label.size = 0, nudge_x = .225, ) +
  geom_text(data = light_df, aes(label = pvalue, y = mean + se + 100), size = 6) +
  theme_minimal() +
  theme(
    legend.position = "none",
  ) +
  scale_fill_manual(values = my_colors) +
  xlab("Gene") +
  ylab("Normalized count")

######################################################

library(patchwork)

((plot_spacer()) + p2  + plot_spacer() + plot_layout(widths = c(1,48,1))) /((plot_spacer()) + p3  + plot_spacer() + plot_layout(widths = c(1,48,1))) + plot_layout(heights = c(3,1)) + plot_annotation(tag_levels = 'A')
