library(tidyverse)
library(Biostrings)

seqs_edit <- readDNAStringSet("../results/../results/seqs_hk_sanger_edit.fasta")
seqs_hk <- seqs_edit[grepl("filled", names(seqs_edit))]
# names(seqs_hk) <- gsub("_filled", "", names(seqs_hk))

# check  genome coverage
cov_check_s <- apply(alphabetFrequency(seqs_edit), 1, function(x){
	sum(x[1:4])
})
cov_check_s <- round(cov_check_s/29903*100, 2)
names(cov_check_s) <- names(seqs_edit)
sort(cov_check_s)
seqs_edit_hq <- seqs_edit[cov_check_s>80]

seq_int <- seqs_hk[grepl("Case_A", names(seqs_hk))][1]

# seqs_compare <- seqs_all[grepl("12098", names(seqs_all))]
seqs_compare <- seqs_edit_hq
num_diff <- sapply(seq_along(seqs_compare), function(i){
	seq_t <- seqs_compare[i]
	check <- strsplit(compareStrings(seq_int, seq_t), "")[[1]]== "?"
	idx <- which(check)
	check1 <- !strsplit(as.character(seq_int), "")[[1]][idx] %in% c("-", "N")
	check2 <- !strsplit(as.character(seq_t), "")[[1]][idx] %in% c("-", "N")
	print(idx[check1 & check2])
	sum(check1 & check2)
})

pos_diff <- sapply(seq_along(seqs_compare), function(i){
	seq_t <- seqs_compare[i]
	check <- strsplit(compareStrings(seq_int, seq_t), "")[[1]]== "?"
	idx <- which(check)
	check1 <- !strsplit(as.character(seq_int), "")[[1]][idx] %in% c("-", "N")
	check2 <- !strsplit(as.character(seq_t), "")[[1]][idx] %in% c("-", "N")
	paste(idx[check1 & check2], collapse=" | ")
})

# subseq(seqs, 6167, 6167)

names(num_diff) <- names(seqs_compare)
head(sort(num_diff))
df_diff <- tibble(sample=names(num_diff), diff_of_nt=num_diff, position = pos_diff) %>% arrange(diff_of_nt)
write_csv(df_diff, "../results/df_diff_to_case_A.csv")
