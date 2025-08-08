### this script is to choose the outgroups of rooted trees, and their metadata to integrate to the trees
# producing: the assembly list of each outgroup (for downloading sequences), the metadata of 5 seq of each outgroup and the tree metadata with both in- and out- group

# packages

library(data.table)
library(plyr)
library(tidyverse) # this package includes "dplyr"

# now have all the assembly barcodes of selected sequences
# downloaded all the assemblies in fasta format using the python script in Pycharm (modified the codes at the format of output files as fasta) 
# transferred all the downloaded assemblies on the OUCRU server

# run the analysis with tools: IQ tree, prokka (annotation), panaroo (pangenome investigation)

# decide to use ggtree to visualize and annotate the phylogenetics trees
# how to combine metadata (host species, location...) into the trees
# next step is to try to run prokka on the cluster (note to check in which env the tools are installed)
# testing time to run prokka on 191 assemblies with 10 CPU cores on cluster
# upgrade the system to 40CPU
# next step is to write a script to take all the gff files in different directories and put into a directory (done)
# now run prokka to annotate 191 assemblies of ST10_AUS (done, it takes 3hrs to finish annotation)
# move all the gff files to a new directory
# run panaroo on 191 gff files of ST10_AUS -> not generate the core gene alignment, so next step to add options in panaroo command to generate the core gene alignment
# run panaroo on 191 gff files to generate the core_gene_alignment_filter files
# run iqtree on the core_gene_alignment to generate the tree file of ST10_AUS
# visualize the tree file of ST10_AUS using figtree (done)
# now try to generate the tree files of ST11_CAN and ST131_US (done)

# all the phylogenetics trees are created by iqtree 
# unrooted trees, so need to root by an outgroup seq
# since the outgroup (E.fergusonii) in the rooted trees has very long branch length, the trees can be improved by choosing other outgroup sequences of E. coli (as FB suggest)
# same phylogroup but different ST to form the outgroup (ref for the phylogroup and ST: Denamur et al. - 2021 - The population genetics of pathogenic Escherichia)
# ST131: phylo B2, outgroup could be seq from ST73, ST95
# ST11: phylo E, outgroup could be seq from ST32, (or ST280)
# ST10: phylo A, outgroup could be seq from ST93, ST167 (or ST6, ST46)

# search for the outgroup of ST131. Try with ST95

all_seq_E.coli_shigella_2 <- read.csv("data/raw/metadata_E. coli_Shigella.csv", sep = ",")
ST95 <- filter(all_seq_E.coli_shigella_2, ST == "95") # 1861 seq

outgroup_ST95 <- fread("data/outgroup/ST95/ST95_metadata.txt") # download the metadata file from Enterobase first
outgroup_ST95 <- outgroup_ST95[, -40]
outgroup_ST95_merged <- merge(ST95, outgroup_ST95, by = "Uberstrain", all.x = TRUE) # 1861 seq
names(outgroup_ST95_merged)[13] <- "source_type"

table(outgroup_ST95_merged$source_type, useNA = "always") # 7 sources, try to choose 3 seq per source (5 from human)
table(outgroup_ST95_merged$ST_Country, useNA = "always")

# randomly choose 22 seq ST95 from different sources to form the outgroup
# create the assembly list, then download the seq
# annotate the seq using prokka
# exclude the gff file of the outgroup seq of E. fergusonii ESC_LC7308AA_AS from the folder of ST131 and add the 22 gff files of the outgroup
# run the panaroo and iqtree again
# be careful when randomly select the samples as re-run the code can create a different set of samples

ST95_random_swine <- slice_sample(filter(outgroup_ST95_merged, source_type == "Swine"), n = 3)
ST95_random_human <- slice_sample(filter(outgroup_ST95_merged, source_type == "Human"), n = 5)
ST95_random_bovine <- slice_sample(filter(outgroup_ST95_merged, source_type == "Bovine"), n = 3)
ST95_random_Canine <- slice_sample(filter(outgroup_ST95_merged, source_type == "Canine"), n = 3)
ST95_random_poultry <- slice_sample(filter(outgroup_ST95_merged, source_type == "Poultry"), n = 3)
ST95_random_env <- slice_sample(filter(outgroup_ST95_merged, source_type == "Water/River"), n = 3)
ST95_random_avian <- slice_sample(filter(outgroup_ST95_merged, source_type == "Avian"), n = 3)

ST95_outgroup <- bind_rows(ST95_random_swine, ST95_random_human, ST95_random_bovine, ST95_random_Canine, ST95_random_poultry, ST95_random_env, ST95_random_avian)
# write.table(ST95_outgroup, file = "ST95_22outgroup.csv", row.names = FALSE, sep = ",")
assembly_ST95_outgroup <- select(ST95_outgroup, 105)
# write.table(assembly_ST95_outgroup, file = "assemblies.txt", row.names = FALSE, sep = ",", quote = FALSE, col.names = FALSE)

# search for the outgroup of ST11. Try with ST32 (no seq of ST280)
setwd("~/Phylogeography/E. coli_phylogeography/outgroup/ST32/")

ST32 <- filter(all_seq_E.coli_shigella_2, ST == "32") # 739 seq
outgroup_ST32_metadata <- fread("data/outgroup/ST32/ST32_enterobase.txt")
outgroup_ST32_metadata <- outgroup_ST32_metadata[,-40]
outgroup_ST32_merged <- merge(ST32, outgroup_ST32_metadata, by = "Uberstrain", all.x = TRUE) #  739seq
names(outgroup_ST32_merged)[13] <- "source_type"

table(outgroup_ST32_merged$source_type, useNA = "always") # 4 sources, choose 5 seq per source to form the outgroup
table(outgroup_ST32_merged$ST_Country, useNA = "always")

ST32_random_bovine <- slice_sample(filter(outgroup_ST32_merged, source_type == "Bovine"), n = 5)
ST32_random_human <- slice_sample(filter(outgroup_ST32_merged, source_type == "Human"), n = 5)
ST32_random_swine <- slice_sample(filter(outgroup_ST32_merged, source_type == "Swine"), n = 5)
ST32_random_env <- slice_sample(filter(outgroup_ST32_merged, source_type == "Water/River"), n = 5)

ST32_outgroup <- bind_rows(ST32_random_bovine, ST32_random_human, ST32_random_swine, ST32_random_env)
# write.table(ST32_outgroup, file = "ST32_20outgroup.csv", row.names = FALSE, sep = ",")
assembly_ST32_outgroup <- select(ST32_outgroup, 105)
# write.table(assembly_ST32_outgroup, file = "assemblies.txt", row.names = FALSE, sep = ",", quote = FALSE, col.names = FALSE)

# search for the outgroup of ST10. Try with ST93
# ST10: phylo A, outgroup could be seq from ST93, ST167 (or ST6, ST46)

ST93 <- filter(all_seq_E.coli_shigella_2, ST == "93") # 144 seq
ST93_enterobase <- fread("data/outgroup/ST93/ST93_enterobase.txt")
ST93_enterobase <- ST93_enterobase[, -40]

outgroup_ST93_merged <- merge(ST93, ST93_enterobase, by = "Uberstrain", all.x = TRUE) # 144 seq
names(outgroup_ST93_merged)[13] <- "source_type"
table(outgroup_ST93_merged$source_type, useNA = "always") # 7 sources, choose 3seq per source (5 from human)
table(outgroup_ST93_merged$ST_Country, useNA = "always")

ST93_random_env <- slice_sample(filter(outgroup_ST93_merged, source_type == "Water/River"), n = 3)
ST93_random_avian <- slice_sample(filter(outgroup_ST93_merged, source_type == "Avian"), n = 3)
ST93_random_canine <- slice_sample(filter(outgroup_ST93_merged, source_type == "Canine"), n = 3)
ST93_random_bovine <- slice_sample(filter(outgroup_ST93_merged, source_type == "Bovine"), n = 3)
ST93_random_human <- slice_sample(filter(outgroup_ST93_merged, source_type == "Human"), n = 5)
ST93_random_poultry <- slice_sample(filter(outgroup_ST93_merged, source_type == "Poultry"), n = 3)
ST93_random_swine <- slice_sample(filter(outgroup_ST93_merged, source_type == "Swine"), n = 3)

ST93_outgroup <- bind_rows(ST93_random_env, ST93_random_avian, ST93_random_canine, ST93_random_bovine, ST93_random_human, ST93_random_poultry, ST93_random_swine)
# write.table(ST93_outgroup, file = "ST93_outgroup.csv", row.names = FALSE, sep = ",")
assembly_ST93_outgroup <- select(ST93_outgroup, 105)
# write.table(assembly_ST93_outgroup, file = "assemblies.txt", row.names = FALSE, sep = ",", quote = FALSE, col.names = FALSE)

# TREE metadata for ST131SPA and ST131US with 5seq of ST95 as outgroup (just 273 assembly barcodes available for ST131SPA)

metadata_ST95_out <- read.csv("data/outgroup/ST95/ST95_22outgroup.csv")
metadata_ST95_5outgr <- filter(metadata_ST95_out, source_type == "Human")
#write.table(metadata_ST95_5outgr, file = "metadata_ST95_5outgr.csv", row.names = FALSE, sep = ",")
names(metadata_ST95_5outgr)

metadata_ST95_5outgr <- select(metadata_ST95_5outgr, c(1, 14, 29, 45))
metadata_ST131SPA <- select(ST131_SPA_merged, c(106, 13, 28, 44))
metadata_ST131US <- select(ST131_US_merged, c(106, 13, 28, 44))

t <- c("assembly_barcode", "source_type", "species", "ST")
names(metadata_ST95_5outgr) <- t
names(metadata_ST131SPA) <- t
names(metadata_ST131US) <- t

table(metadata_ST131SPA$assembly_barcode, useNA = "always")
metadata_ST131SPA <- filter(metadata_ST131SPA, !is.na(assembly_barcode))
write.table(metadata_ST131SPA, file = "metadata_tree_ST131SPA.csv", row.names = FALSE, sep = ",")
metadata_ST131SPA_5outgr <- bind_rows(metadata_ST131SPA, metadata_ST95_5outgr)
# write.table(metadata_ST131SPA_5outgr, file = "metadata_tree_ST131SPA_5outgr.csv", row.names = FALSE, sep = ",")

metadata_ST131US_5outgr <- bind_rows(metadata_ST131US, metadata_ST95_5outgr)
# write.table(metadata_ST131US_5outgr, file = "metadata_tree_ST131US_5outgr.csv", row.names = FALSE, sep = ",")

# TREE metadata for ST11CAN and ST11NEWZ intergrating with 5 seq of ST32 

metadata_ST32_out <- read.csv("data/outgroup/ST32/ST32_20outgroup.csv")
metadata_ST32_5outgr <- filter(metadata_ST32_out, source_type == "Human")
# write.table(metadata_ST32_5outgr, file = "metadata_ST32_5outgr.csv", row.names = FALSE, sep = ",")

metadata_ST32_5outgr <- select(metadata_ST32_5outgr, c(105, 13, 28, 44))
metadata_ST11NEWZ <- select(ST11_NEWZ_merged, c(106, 13, 28, 44))
metadata_ST11CAN <- select(ST11_CAN_merged, c(106, 13, 28, 44))

names(metadata_ST32_5outgr) <- t
names(metadata_ST11NEWZ) <- t
names(metadata_ST11CAN) <- t

metadata_ST11NEWZ_5outgr <- rbind(metadata_ST11NEWZ, metadata_ST32_5outgr)
# write.table(metadata_ST11NEWZ_5outgr, file = "metadata_tree_ST11NEWZ_5outgr.csv", row.names = FALSE, sep = ",")
metadata_ST11CAN_5outgr <- rbind(metadata_ST11CAN, metadata_ST32_5outgr)
# write.table(metadata_ST11CAN_5outgr, file = "metadata_tree_ST11CAN_5outgr.csv", row.names = FALSE, sep = ",")

# TREE metadata for ST10AUS and ST10CAN intergrating with 5 outgr seq of ST93

metadata_ST93_outgr <- read.csv("data/outgroup/ST93/ST93_outgroup.csv")
metadata_ST93_5outgr <- filter(metadata_ST93_outgr, source_type == "Human")
write.table(metadata_ST93_5outgr, file = "metadata_ST93_5outgr.csv", row.names = FALSE, sep = ",")

names(metadata_ST93_5outgr)

metadata_ST93_5outgr <- select(metadata_ST93_5outgr, c(105, 13, 28, 44))
metadata_ST10CAN <- select(ST10_CAN_merged, c(106, 13, 28, 44))
metadata_ST10AUS <- select(ST10_AUS_merged, c(106, 13, 28, 44))

names(metadata_ST93_5outgr) <- t
names(metadata_ST10CAN) <- t
names(metadata_ST10AUS) <- t

metadata_ST10CAN_5outgr <- rbind(metadata_ST10CAN, metadata_ST93_5outgr)
# write.table(metadata_ST10CAN_5outgr, file = "metadata_tree_ST10CAN_5outgr.csv", row.names = FALSE, sep = ",")
metadata_ST10AUS_5outgr <- rbind(metadata_ST10AUS, metadata_ST93_5outgr)
# write.table(metadata_ST10AUS_5outgr, file = "metadata_tree_ST10AUS_5outgr.csv", row.names = FALSE, sep = ",")




