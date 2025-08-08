### this script is used to select the metadata of sequences of 6 combinations of ST-country with balanced number of host types,
# and then generated the list of assembly barcodes which are used to download the assembled sequences from the enterobase

# packages
library(data.table)
library(plyr)
library(tidyverse) # this package includes "dplyr"

# set working directory
if(grepl(pattern = "blanquart", getwd())) setwd("~/ownCloud/Nhung_Marc_phylogeography/") else setwd("~/Phylogeography/Ecoli_phylo/")

all_seq_E.coli_shigella_2 <- read.csv("data/raw/metadata_E. coli_Shigella.csv", sep = ",")

# the total number of ST of Escherichia is 7197 STs
# 87382 sequences of E. coli and Shigella strains belong to 207 ST-country combinations which has enough data (>100 seq with at least human host type)
# try with one interesting ST-country to see how to download the corresponding sequences
# try ST-131 first (more relevant to MDR)
# some other common STs related to infections: 12, 23, 14, 144, 69, 88, 58...(ref: Burgaya et al. 2023); 11, 131, 73, 10, 95, 21 (Horesh et al. 2021)
# choose ST131_US, ST11_Canada (ST11 associated with the pathotype EHEC)
# ST10_AUS: ST10, a broad-host-range ST

# there are 16 ST131_country in the dataset (means 16 different countries)
ST_131 <- filter(all_seq_E.coli_shigella_2, ST == "131") # there are 16 ST131_country in the dataset
str_extract(tab_e.coli_shigella$ST_Country, "131")
tab_e.coli_shigella <- mutate(tab_e.coli_shigella, is_ST131 = str_extract(tab_e.coli_shigella$ST_Country, "131"))
table(tab_e.coli_shigella$is_ST131)  

# choose a more balanced number of host types
# choose a different country for each ST10, ST11, ST131 (it is easier to choose for ST10 as being in broad hosts)
# choose a more balanced number of host types

tab_e.coli_shigella <- mutate(tab_e.coli_shigella, N_non_HumanHost = apply(tab_e.coli_shigella[,c(3:7)], MARGIN = 1, sum))
tab_e.coli_shigella$N_non_HumanHost
tab_e.coli_shigella_new <- relocate(tab_e.coli_shigella, N_non_HumanHost, .after = N_human_e.coli_shigella)
write.table(tab_e.coli_shigella_new, file = "summary_e.coli_shigella.csv", row.names = FALSE, sep = ",")

# choose ST10_CAN, ST11_NEWZ, ST131_Spain
# ST11_Canada (ST11 associated with the pathotype EHEC, intestinal pathogen): 324 human origin vs 221 bovine vs 75 swine vs 10 water
# ST131_US (extra-intestinal pathogen): 2591 human origin vs 1 bovine vs 36 swine vs 237 poultry vs 112 dog vs 25 water
# ST10_AUS (broad host range, commensal): 114 human origin vs 46 swine vs 6 poultry vs 6 dog vs 19 water

ST131_US <- filter(all_seq_E.coli_shigella_2, ST_Country == "131_United States") # 3002 seq
ST11_CAN <- filter(all_seq_E.coli_shigella_2, ST_Country == "11_Canada") # 630 seq
ST10_AUS <- filter(all_seq_E.coli_shigella_2, ST_Country == "10_Australia") # 191 seq

ST10_CAN <- filter(all_seq_E.coli_shigella_2, ST_Country == "10_Canada") # 167 seq
ST11_NEWZ <- filter(all_seq_E.coli_shigella_2, ST_Country == "11_New Zealand") # 159 seq
ST131_SPA <- filter(all_seq_E.coli_shigella_2, ST_Country == "131_Spain") # 294 seq

# check all are e. coli
table(ST11_CAN$Species)
table(ST131_US$Species)
table(ST10_AUS$Species)
table(ST10_CAN$Species)
table(ST11_NEWZ$Species)
table(ST131_SPA$Species)

# searching for and creating lists of assembly barcodes (ESC...AS) to download sequences
# search the assembly barcodes on enterobase first by searching manually using the criteria like ST, country, then select the experimental data of the assemblies, then export the data in txt file
# the txt files of assembly barcodes will be used to download the seq from enterodatabase (using python script)
# create a metadata for each choosen ST

# ST11_CAN: 630seq

ST11_CAN_enterobase <- fread("data/ST_country/ST11_CAN/ST11_CAN.txt")
ST11_CAN_enterobase <- ST11_CAN_enterobase[, -40]
ST11_CAN_merged <- merge(ST11_CAN, ST11_CAN_enterobase, by = "Uberstrain", all.x = TRUE)
# write.table(ST11_CAN_merged, file = "ST11CAN_merged_database.csv", row.names = FALSE, sep = ",")
names(ST11_CAN_merged)
table(ST11_CAN_merged$Status) # all the assemblies are available
assembly_ST11CAN <- select(ST11_CAN_merged, -c(1:105))
# write.table(assembly_ST11CAN, file = "assemblies.txt", row.names = FALSE, sep = ",", quote = FALSE, col.names = FALSE)

# ST131_US: 3002 seq

ST131_US_enterobase <- fread("data/ST_country/ST131_US/ST131_US.txt")
ST131_US_enterobase <- ST131_US_enterobase[,-40]
ST131_US_merged <- merge(ST131_US, ST131_US_enterobase, by = "Uberstrain", all.x = TRUE)
# write.table(ST131_US_merged, file = "ST131US_merged_database.csv", row.names = FALSE, sep = ",")
assembly_ST131US <- select(ST131_US_merged, -c(1:105))
# write.table(assembly_ST131US, file = "assemblies.txt", row.names = FALSE, sep = ",", quote = FALSE, col.names = FALSE)

# ST10_AUS: 191seq

ST10_AUS_enterobase <- fread("data/ST_country/ST10_AUS/ST10_AUS.txt")
ST10_AUS_enterobase <- ST10_AUS_enterobase[, -40]
ST10_AUS_merged <- merge(ST10_AUS, ST10_AUS_enterobase, by = "Uberstrain", all.x = TRUE)
# write.table(ST10_AUS_merged, file = "ST10AUS_merged_database.csv", row.names = FALSE, sep = ",")
assembly_ST10AUS <- select(ST10_AUS_merged, -c(1:105))
# write.table(assembly_ST10AUS, file = "assemblies.txt", row.names = FALSE, sep = ",", quote = FALSE, col.names = FALSE)

# ST131_SPA: 294 seq (just 273 assembly barcodes available)

ST131_SPA_enterobase <- fread("data/ST_country/ST131_SPA/ST131_SPA.txt") # 293 available to download
ST131_SPA_enterobase <- ST131_SPA_enterobase[, -40]
test3 <- as.data.frame(table(ST131_SPA_enterobase$`Assembly barcode`, useNA = "always"))
ST131_SPA_merged <- merge(ST131_SPA, ST131_SPA_enterobase, by = "Uberstrain", all.x = TRUE)
# write.table(ST131_SPA_merged, file = "ST131SPA_merged_database.csv", row.names = FALSE, sep = ",")
test <- as.data.frame(table(ST131_SPA_merged$`Assembly barcode`, useNA = "always")) # just 273 assembly barcodes available
assembly_ST131SPA <- select(ST131_SPA_merged, -c(1:105))
assembly_ST131SPA <- filter(assembly_ST131SPA, !is.na(barcode))
# write.table(assembly_ST131SPA, file = "assemblies.txt", row.names = FALSE, sep = ",", quote = FALSE, col.names = FALSE)

# ST10_CAN: 167seq 

ST10_CAN_enterobase <- fread("data/ST_country/ST10_CAN/ST10_CAN.txt")
ST10_CAN_enterobase <- ST10_CAN_enterobase[,-40]
ST10_CAN_merged <- merge(ST10_CAN, ST10_CAN_enterobase, by = "Uberstrain", all.x = TRUE)
# write.table(ST10_CAN_merged, file = "ST10CAN_merged_database.csv", row.names = FALSE, sep = ",", quote = FALSE, col.names = FALSE)
table(ST10_CAN_merged$`Assembly barcode`, useNA = "always")
assembly_ST10CAN <- select(ST10_CAN_merged, -c(1:105))
# write.table(assembly_ST10CAN, file = "assemblies.txt", row.names = FALSE, sep = ",", quote = FALSE, col.names = FALSE)

# ST11_NEWZ: 159 seq

ST11_NEWZ_enterobase <- fread("data/ST_country/ST11_NEWZ/ST11_NEWZ.txt")
ST11_NEWZ_enterobase <- ST11_NEWZ_enterobase[,-40]
ST11_NEWZ_merged <- merge(ST11_NEWZ, ST11_NEWZ_enterobase, by = "Uberstrain", all.x = TRUE)
write.table(ST11_NEWZ_merged, file = "ST11NEWZ_merged_database.csv", row.names = FALSE, sep = ",")
table(ST11_NEWZ_merged$`Assembly barcode`, useNA = "always")
assembly_ST11NEWZ <- select(ST11_NEWZ_merged, -c(1:105))
# write.table(assembly_ST11NEWZ, file = "assemblies.txt", row.names = FALSE, sep = ",", quote = FALSE, col.names = FALSE)


