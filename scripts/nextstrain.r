library(Biostrings)
library(tidyverse)

date_tmp <- "2021-11-30"

fasta_all_ng <- readDNAStringSet("../results/seqs_sub_more.fasta")
fasta_all_ng <- fasta_all_ng[-1]
df_date <- read_tsv("../results/df_date.tsv", col_names=F)

## formating metadata
meta_sample <- read_tsv("~/Downloads/ncov/data/example_metadata.tsv")
# view(meta_sample[1,]) 
meta_sample <- meta_sample[1,]
meta_sample <- meta_sample %>% mutate_all(function(x){return(NA)})
meta_df <- left_join(tibble(strain = names(fasta_all_ng)), meta_sample, "strain")

meta_df$virus  <- rep("ncov", nrow(meta_df))
meta_df$region  <- rep("Asia", nrow(meta_df))
meta_df$location  <- rep("Hong Kong", nrow(meta_df))
meta_df$division  <- rep("Hong Kong", nrow(meta_df))
meta_df$country  <- rep("Hong Kong", nrow(meta_df))
meta_df$date <- sapply(meta_df$strain, function(x){
	if(x=="MN908947_3"){return("2019-12-26")}
	tmp <- df_date$X2[df_date$X1==x]
	if(length(tmp)==0){return("2021-11-10")}else{return(as.character(tmp))}
	# date_tmp = lubridate::ymd(date_tmp)
	# return(as.character(date_tmp))
}, USE.NAMES = F)
meta_df$age <- 10
meta_df$sex <- "M"
meta_df$date_submitted <- date_tmp

## add Wuhan/Hu-1/2019
## add Wuhan/WH01/2019
meta_sample <- read_tsv("~/Downloads/ncov/data/example_metadata.tsv", col_types = cols(.default = "c"))
meta_df <- meta_df %>% mutate_all(as.character)
meta_df <- bind_rows(meta_sample[meta_sample$strain %in% c("Wuhan/Hu-1/2019", "Wuhan/WH01/2019"),], meta_df)

## add info
meta_df$region_exposure <- meta_df$region
meta_df$country_exposure <- meta_df$country
meta_df$division_exposure <- meta_df$division
# meta_df$authors[meta_df$authors == "Epi-I"] <- "Imported"
# meta_df$authors[meta_df$authors == "Epi-L"] <- "Local"
# meta_df$authors[meta_df$authors == "PL"] <- "Local"
# meta_df$authors[meta_df$authors == "Epi-PL"] <- "Local"
# meta_df$authors[meta_df$authors == "I"] <- "Imported"
# meta_df$authors[meta_df$authors == "L"] <- "Local"

meta_df$originating_lab[1:2] <- NA
meta_df$submitting_lab[1:2] <- NA
meta_df$host[1:2] <- NA
meta_df$authors[1:2] <- NA

# meta_df$strain[3:nrow(meta_df)] <- sapply(3:nrow(meta_df), function(i){
# 	df_tmp = meta_df[i,]
# 	paste0("Case_", df_tmp$strain, "/",
# 		"Age_", df_tmp$age, "/",
# 		"Sex_", df_tmp$sex, "/",
# 		"ReportedDate_", df_tmp$date, "/",
# 		"Cluster_", df_tmp$segment, "/",
# 		"Type_", df_tmp$authors, "/",
# 		"Citizen_", df_tmp$originating_lab, "/",
# 		"CitizenDistrict_", df_tmp$submitting_lab, "/")
# })
# meta_df$strain <- gsub(" ", "_", meta_df$strain)
# meta_df$strain <- gsub("&", "", meta_df$strain)

ref_seq <- readDNAStringSet("../../../2020/2020-09-01_COVID_NGS_pipeline/NGS_data_input/reference.fasta")
fasta_all_ng <- c(ref_seq, ref_seq, fasta_all_ng)
meta_df$strain[grepl("EPI_", meta_df$strain)] <- sapply(meta_df$strain[grepl("EPI_", meta_df$strain)], function(x){
	tmp <- strsplit(x, "|", fixed=T)[[1]]
	tmp[grepl("EPI_", tmp)]
})
names(fasta_all_ng) <- meta_df$strain

write_tsv(meta_df, "../results/metadata_nextstrain_hk.tsv")
writeXStringSet(fasta_all_ng, "../results/seq_nextstrain_hk.fasta")

names(fasta_all_ng) %in% meta_df$strain

system("cp -a ../results/metadata_nextstrain_hk.tsv ~/Downloads/ncov/data/metadata.tsv")
system("cp -a ../results/seq_nextstrain_hk.fasta ~/Downloads/ncov/data/sequences.fasta")

# cd ../ncov
# nextstrain build . --cores 8 --use-conda --configfile ./my_profiles/20211108/builds.yaml



