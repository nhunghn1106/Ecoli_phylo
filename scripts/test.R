# test
nwk <- system.file("extdata", "sample.nwk", package="treeio")

tree <- read.tree(nwk)
ggtree(tree)
circ <- ggtree(tree, layout = "circular")

df <- data.frame(first=c("a", "b", "a", "c", "d", "d", "a", 
                         "b", "e", "e", "f", "c", "f"),
                 second= c("z", "z", "z", "z", "y", "y", 
                           "y", "y", "x", "x", "x", "a", "a"))
rownames(df) <- tree$tip.label

df2 <- as.data.frame(matrix(rnorm(39), ncol=3))
rownames(df2) <- tree$tip.label
colnames(df2) <- LETTERS[1:3]


p1 <- gheatmap(circ, df, offset=.8, width=.2,
               colnames_angle=95, colnames_offset_y = .25) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")


git clong https://github.com/katholt/plotTree.git
setwd("~/Phylogeography/ggtree_tutorial/tree_example_april2015/")

info <- read.csv("info.csv")
tree <- read.tree("tree.nwk")
cols <- c(HCMC='black', Hue='purple2', KH='skyblue2')

p <- ggtree(tree, layout='circular') %<+% info + 
  geom_tippoint(aes(color=location)) + 
  scale_color_manual(values=cols) + 
  geom_tiplab2(aes(label=name), align=T, linetype=NA, 
               size=2, offset=4, hjust=0.5) +
  geom_tiplab2(aes(label=year), align=T, linetype=NA, 
               size=2, offset=8, hjust=0.5)

heatmapData=read.csv("res_genes.csv", row.names=1)
rn <- rownames(heatmapData)
heatmapData <- as.data.frame(sapply(heatmapData, as.character))
rownames(heatmapData) <- rn

heatmap.colours <- c("white","grey","seagreen3","darkgreen",
                     "green","brown","tan","red","orange",
                     "pink","magenta","purple","blue","skyblue3",
                     "blue","skyblue2")
names(heatmap.colours) <- 0:15

gheatmap(p, heatmapData, offset = 10, color=NULL, 
         colnames_position="top", 
         colnames_angle=90, colnames_offset_y = 1, 
         hjust=0, font.size=2) +
  scale_fill_manual(values=heatmap.colours, breaks=0:15)

p2 <- ggtree(tree, layout='circular')

gheatmap(p2, heatmapData, offset = 10, color=NULL, 
         colnames_position="top", 
         colnames_angle=90, colnames_offset_y = 1, 
         hjust=0, font.size=2) +
  scale_fill_manual(values=heatmap.colours, breaks=0:15)

################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtreeExtra")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtreeExtra")

install.packages("ggstar")

library(ggtreeExtra)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(treeio)
library(tidytree)
library(dplyr)
library(ggstar)
library(TDbook)

# load tree_NJIDqgsS and df_NJIDqgsS from TDbook
tr <- tree_NJIDqgsS
metada <- df_NJIDqgsS
metadata <- metada %>%
  select(c("id", "country", "country__colour", 
           "year", "year__colour", "haplotype"))
metadata$haplotype[nchar(metadata$haplotype) == 0] <- NA

countrycolors <- metada %>%
  select(c("country", "country__colour")) %>%
  distinct()

yearcolors <- metada %>%
  select(c("year", "year__colour")) %>%
  distinct()
yearcolors <- yearcolors[order(yearcolors$year, decreasing=TRUE),]

metadata$country <- factor(metadata$country, levels=countrycolors$country) # make the var a factor variable
metadata$year <- factor(metadata$year, levels=yearcolors$year) # make this var a factor variable


p <- ggtree(tr, layout="fan", open.angle=15, size=0.1)

p <- p %<+% metadata

p1 <-p +
  geom_tippoint(
    mapping=aes(colour=country),
    size=1.5,
    stroke=0,
    alpha=0.4
  ) +
  scale_colour_manual(
    name="Country",
    values=countrycolors$country__colour,
    guide=guide_legend(keywidth=0.3,
                       keyheight=0.3,
                       ncol=2,
                       override.aes=list(size=2,alpha=1),
                       order=1)
  ) +
  theme(
    legend.title=element_text(size=5),
    legend.text=element_text(size=4),
    legend.spacing.y = unit(0.02, "cm")
  )

p2 <-p1 +
  geom_fruit(
    geom=geom_star,
    mapping=aes(fill=haplotype),
    starshape=26,
    color=NA,
    size=2,
    starstroke=0,
    offset=0,
  ) +
  scale_fill_manual(
    name="Haplotype",
    values=c("red"),
    guide=guide_legend(
      keywidth=0.3,
      keyheight=0.3,
      order=3
    ),
    na.translate=FALSE
  )

p3 <-p2 +
  new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=year),
    width=0.002,
    offset=0.1
  ) +
  scale_fill_manual(
    name="Year",
    values=yearcolors$year__colour,
    guide=guide_legend(keywidth=0.3, keyheight=0.3, ncol=2, order=2)
  ) +
  theme(
    legend.title=element_text(size=6), 
    legend.text=element_text(size=4.5),
    legend.spacing.y = unit(0.02, "cm")
  )
p3

testcol <- c("#0000FF", "#00008B", "#7CFC00", "#FFD700", "#FF0000", "#FFD700", "#F0F8FF", "#E78AC3", "#FFFFFF")





