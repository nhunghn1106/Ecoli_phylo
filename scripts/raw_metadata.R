### this script is used to process the raw metadata on enterobase, producing the metadata of E.coli and Shigella sequences,
# which is selected based on the combination of ST-country (each >= 100 seq and at least with human host)

# packages
library(data.table)
library(plyr)
library(tidyverse) # this package includes "dplyr"

rm(list = ls())

# set working directory
if(grepl(pattern = "blanquart", getwd())) setwd("~/ownCloud/Nhung_Marc_phylogeography/") else setwd("~/Phylogeography/Ecoli_phylo/")

# read Enterobase metadata (Enterobase is strain based, including metadata, genomic assemblies (from sequenced reads) and deduced genotyping data (MLST data))
# metadata consists of information act as hyperlinks to the original data sources (isolation year, host, geographical location, unique identifier, strain name)
# experimental data or genomic data (assembly statistics, MLST sequence types), deduced exclusively from the assemblies of sequenced reads
# the negative numbers in MLST scheme represent for missing loci in ST column
# serotype (O antigen and H antigen): prediction, experimental data
# the field "uberstrain" is unique for replicated strains (represent for merging duplicated strains)
# the field "name": unique name given to the strain

# clean raw data
a <- fread("data/raw/all_enterobase_strains.txt")

tab_ST <- table(a$ST)
class(tab_ST)
tab_ST[tab_ST>1000]  # check which ST with number of strains >1000
common_STs <- names(.Last.value)
length(common_STs)
a$is_common_ST <- a$ST %in% common_STs # check if each ST is a common ST >1000 strains
table(a$is_common_ST)

tab_country <- table(a$Country)
tab_country[tab_country>1000]  # check which country with number of strains >1000
common_countries <- names(.Last.value); common_countries <- common_countries[common_countries!=""]
a$is_common_country <- a$Country %in% common_countries
table(a$is_common_country)

# the question: the migration rate (transmission rate) of important clones btw host types (humans vs animals) in different countries
# the concerned variables: ST, Country, Source Type, Source Niche

names(a)
levels(as.factor(a$ST))
levels(as.factor(a$Country))
levels(as.factor(a$`Source Type`))
levels(as.factor(a$`Source Niche`))
levels(as.factor(a$`Source Details`))

length(levels(as.factor(a$ST)))
length(levels(as.factor(a$Country))) # 193 countries
length(levels(as.factor(a$`Source Type`))) # 69 source types
length(levels(as.factor(a$`Source Niche`))) # 14 source niches

table(a$`Source Niche`)      
table(a$`Source Type`)
table(a$`Source Details`)

# clean up sources of strains
# investigate the transmission rates btw human and animals
# so concerned source niches: "human", "companion animal", "environment", "livestock", "poultry"
# concerned source types: "human", "canine", "bovine", "swine", "avian", "poultry"

tab_source <- table(a$`Source Type`, a$`Source Niche`, useNA = "ifany")
which(tab_source > 1000, arr.ind = T)  # show the indices where the sources with more than 1000 strains, have a favor of the combinations of sources with >1000 strains


a$is_bovine <- a$`Source Niche`=="Livestock" & a$`Source Type`=="Bovine" # investigate strains from bovine/livestock
table(a$is_bovine)

a$is_swine <- a$`Source Niche`=="Livestock" & a$`Source Type`=="Swine" # investigate strains from swine/livestock
table(a$is_swine)

a$is_poultry <- (a$`Source Niche`=="Poultry" & a$`Source Type`=="Avian") | (a$`Source Niche`=="Poultry" & a$`Source Type`=="Poultry") # investigate strains from poultry
table(a$is_poultry)

a$is_human <- (a$`Source Niche`=="Human" & a$`Source Type`=="Human") | (a$`Source Niche`=="Human" & a$`Source Type`=="") # investigate strains from humans
table(a$is_human)

a$is_human2 <- (a$`Source Niche`=="Human" & a$`Source Type`=="Human") | (a$`Source Niche`=="Human" & a$`Source Type`=="") | (a$`Source Niche`=="" & a$`Source Type`=="Human")
table(a$is_human2) # should include strains collected in Human as source type but no information about source niches

a$is_dog <- a$`Source Niche`=="Companion Animal" & a$`Source Type` == "Canine" # investigate strains from dogs
table(a$is_dog)

a$is_water <- a$`Source Niche`=="Environment" & a$`Source Type` == "Water/River" # investigate strains from water
table(a$is_water)

a$ST_Country <- paste(a$ST, a$Country, sep = "_") # create the combinations of ST and country

# check how many strains in each concerned host types for each combination of ST_country 

tab_old <- ddply(a[a$Country!="",], .(ST_Country), summarise, N_human = sum(is_human), N_bovine = sum(is_bovine), N_swine = sum(is_swine), N_poultry = sum(is_poultry), N_dog = sum(is_dog), N_water = sum(is_water))
tab <- ddply(a[a$Country!="",], .(ST_Country), summarise, N_human = sum(is_human2), N_bovine = sum(is_bovine), N_swine = sum(is_swine), N_poultry = sum(is_poultry), N_dog = sum(is_dog), N_water = sum(is_water))
count(tab, ST_Country) # unique ST_country

# keep ST_country with >100 sequences, we have a table of ST_country which have more than 100 seq

tab2 <- tab[which(rowSums(tab[, 2:7]) >= 100), ]

# identify the metadata of sequences belonging to ST_country which have enough data

table(tab2$N_human)

tab3 <- filter(tab2, N_human != 0) # exclude which combinations of ST-country with no seq from humans
count(tab3, ST_Country)

# in total we have 207 unique combinations of ST_country, each with >=100 seq (at least with the human host type)
# next is to create a metadata of all seq of these 207 ST-country (tab3 and a)
# take all the seq labelled with these 207 ST_country (88892 sequences)

test <- apply(tab3[, (2:7)], 1, sum)
class(test)
sum(test)  # the total number of seq of eligible ST-country is 88892
test2 <- apply(tab3[, (2:7)], 2, sum)
sum(test2)

table(a$ST_Country)

# firstly need to exclude the strains not from the concerned host types
#  # FB THIS DOES NOT WORK FOR ME   # Nhung: I deleted the code here 
a <- mutate(a, 
            is_concerned_host = ifelse(is_human2 == FALSE & is_bovine == FALSE & is_swine == FALSE & is_poultry == FALSE & is_dog == FALSE & is_water == FALSE, "no", "yes"))
table(a$is_concerned_host)

# create a metadata of all seq of 207 eligible ST_country (88892)

names(tab3)[1] <- "ST_Country"
all_seq <- merge(tab3, a[a$is_concerned_host == "yes", ], by="ST_Country", all.x = TRUE)

table(all_seq$ST_Country)
table(all_seq$Species)  # the number of E. coli is 66029 while Shigella spp is 21353

# create a metadata of only seq of E. coli and Shigella spp. (87382 sequences of 207 ST-country)

all_seq_E.coli_Shigella <- filter(all_seq, Species == "Escherichia coli" | Species == "Shigella boydii" | Species == "Shigella dysenteriae" | Species == "Shigella flexneri" | Species == "Shigella sonnei" | Species == "Shigella sp.")  
length(levels(as.factor(all_seq_E.coli_Shigella$ST_Country)))
table(all_seq_E.coli_Shigella$ST_Country)

# re-correct the number of strains in each host type of each ST_country on the metadata of E.coli and Shigella

tab_e.coli_shigella <- ddply(all_seq_E.coli_Shigella, .(ST_Country), summarise, N_human_e.coli_shigella = sum(is_human2), N_bovine_e.coli_shigella = sum(is_bovine), N_swine_e.coli_shigella = sum(is_swine), N_poultry_e.coli_shigella = sum(is_poultry), N_dog_e.coli_shigella = sum(is_dog), N_water_e.coli_shigella = sum(is_water))

all_seq_E.coli_shigella_2 <- merge(tab_e.coli_shigella, all_seq_E.coli_Shigella, by = "ST_Country", all.x = TRUE)
all_seq_E.coli_shigella_2 <- select(all_seq_E.coli_shigella_2, -(8:13))
all_seq_E.coli_shigella_2 <- relocate(all_seq_E.coli_shigella_2, Uberstrain)

# write.table(all_seq_E.coli_shigella_2, file = "metadata_E. coli_Shigella.csv", row.names = FALSE, sep = ",")





