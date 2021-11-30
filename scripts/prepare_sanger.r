library(Biostrings)

seq_files <- list.files("../data/211126 CPOS/", ".seq$", full.names=T)
file_name <- list.files("../data/211126 CPOS/", ".seq$")
seqs <- lapply(seq_along(seq_files), function(i){
	# x <- seq_files[1]
	x <- seq_files[i]
	tmp <- readLines(x)
	tmp <- paste0(tmp, collapse="")
	tmp <- DNAStringSet(tmp)
	names(tmp) <- file_name[i] 
	if(grepl("_R", x) | grepl("R_", x)){
		tmp <- reverseComplement(tmp)
		names(tmp) <- paste0("RC_", names(tmp))
		}
	return(tmp)	
})
seqs <- do.call(c, seqs)

seq_name <- sapply(file_name, function(x){
	x <- gsub(".seq", "", x, fixed=T)
	tmp <- strsplit(x, "_")[[1]]
	paste0(tmp[grepl("WHP", tmp)], "_", x)
})
names(seqs) <- seq_name
seqs <- seqs[order(names(seqs))]
writeXStringSet(seqs, "../results/sanger_seqs.fasta")

# add others
seqs_hk <- readDNAStringSet("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/latest_genome.fasta")
seqs_int <- seqs_hk[grepl("WHP5043", names(seqs_hk)) | grepl("WHP5047", names(seqs_hk))]
seqs_int <- seqs_int[order(names(seqs_int))]
seqs_add <- c(seqs, seqs_int)
writeXStringSet(seqs_add, "../results/sanger_seqs_add.fasta")


# seq_files <- list.files("../data/211126 CPOS/", ".txt$", full.names=T)
# file_name <- list.files("../data/211126 CPOS/", ".txt$")
# seqs <- lapply(seq_along(seq_files), function(i){
# 	# x <- seq_files[1]
# 	x <- seq_files[i]
# 	tmp <- readDNAStringSet(x)
# 	names(tmp) <- file_name[i] 
# 	if(grepl("_R", x) | grepl("R_", x)){
# 		tmp <- reverseComplement(tmp)
# 		names(tmp) <- paste0("RC_", names(tmp))
# 		}
# 	return(tmp)	
# })
# seqs <- do.call(c, seqs)

# seq_name <- sapply(file_name, function(x){
# 	x <- gsub(".txt", "", x, fixed=T)
# 	tmp <- strsplit(x, "_")[[1]]
# 	paste0(tmp[grepl("WHP", tmp)], "_", x)
# })
# names(seqs) <- seq_name
# seqs <- seqs[order(names(seqs))]
# writeXStringSet(seqs, "../results/sanger_seqs_long.fasta")
