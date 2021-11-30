# compare the mutations between sequences
# and compare the mutations to other variants

library(tidyverse)
library(rvest)
library(ggsci)
library(patchwork)
source("https://raw.githubusercontent.com/Koohoko/Save-ggplot-to-pptx/main/scripts/save_pptx.r")


df_mut_hk <- read_tsv("../data/hk_gisaid_hcov-19_2021_11_27_06.tsv")
df_mut_hk$label <- ifelse(df_mut_hk$`Accession ID`=="EPI_ISL_6716902", "Case A", "Case B") # EPI_ISL_6716902 (12388)
df_mut_hk <- df_mut_hk %>% arrange(label)
df_mut_other <- read_tsv("../data/gisaid_hcov-19_2021_11_26_07.tsv")
df_mut_other$label <- df_mut_other$`Virus name`

df_mut <- bind_rows(df_mut_hk, df_mut_other)
aa_muts <- df_mut$`AA Substitutions`
aa_muts <- gsub("\\(", "", aa_muts)
aa_muts <- gsub("\\)", "", aa_muts)
aa_muts <- strsplit(aa_muts, ",")
df_mut$aa_muts <- aa_muts

tmp <- table(unlist(aa_muts))
common_muts <- names(tmp[tmp>1])

df_mut_sum <- tibble(mut=common_muts)
df_mut_sum$Gene <- sapply(df_mut_sum$mut, function(x){
	strsplit(x, "_")[[1]][1]
})

levels_t <- unique(c(naturalsort::naturalsort(unique(df_mut_sum$Gene))[-(1:3)], "M", "N", "E"))
df_mut_sum$Gene <- factor(df_mut_sum$Gene, levels_t)
df_mut_sum$Mutation <- sapply(df_mut_sum$mut, function(x){
	strsplit(x, "_")[[1]][2]
})
df_mut_sum$Position <- sapply(df_mut_sum$Mutation, function(x){
	as.numeric(gsub("\\D", "", x))
})
df_mut_sum <- df_mut_sum %>% arrange(Gene, Position)

mat_check <- sapply(seq_len(nrow(df_mut_sum)), function(i){
	mut_t <- df_mut_sum$mut[i]
	check <- sapply(seq_len(nrow(df_mut)), function(j){
		mut_t %in% df_mut$aa_muts[[j]]
	})
})
mat_check <- t(mat_check)
mat_check <- ifelse(mat_check, "Y", "")
df_check <- as_tibble(mat_check)
colnames(df_check) <- df_mut$label
df_out <- bind_cols(df_mut_sum, df_check)

## add nt position
df_orf <- read_csv("../data/ORF_SCoV2.csv")
# df_orf$sequence <- gsub("nsp", "NSP", df_orf$sequence)
# df_orf$sequence <- gsub("^S$", "Spike", df_orf$sequence)
df_pos_nt <- lapply(seq_len(nrow(df_out)), function(i){
	# print(i)
	gene_t <- tolower(df_out$Gene[i])
	pos_t_sim <- df_out$Position[i]

	if(gene_t=="nsp12"){
		if(pos_t_sim<=9){gene_t <- "nsp12_1"}else{
			pos_t_sim <- pos_t_sim-9
			gene_t <- "nsp12_2"
			}
	} 
	if(gene_t=="spike"){
		gene_t <- "s"
	} 
	if(sum(tolower(df_orf$sequence) == gene_t)>0){
		ref_data_proteins_t <- df_orf %>% filter(tolower(sequence) == gene_t)
	} else {
		gene_t <- paste0("orf",gene_t)
		ref_data_proteins_t <- df_orf %>% filter(tolower(sequence) == gene_t)
	}
	print(gene_t)
	stopifnot(nrow(ref_data_proteins_t)>0)
	pos_nt_start <- pos_t_sim*3-2 + ref_data_proteins_t$start -1
	pos_nt_stop <- pos_t_sim*3 + ref_data_proteins_t$start -1
	tibble(mut=df_out$mut[i], pos_nt_start=pos_nt_start, pos_nt_stop=pos_nt_stop)
})
df_pos_nt <- bind_rows(df_pos_nt)
df_out <- left_join(df_out, df_pos_nt, "mut")

# in other variants
## from cov-lineages
df_lin_def_snps <- read_csv("../../2021-06-02_lineage_specific_snps/results/df_defining_snps_nt.csv")
unique(df_lin_def_snps$orf)
df_lin_def_snps$orf <- toupper(df_lin_def_snps$orf)
df_lin_def_snps$orf <- gsub("^S", "Spike", df_lin_def_snps$orf)
df_lin_def_snps <- df_lin_def_snps %>% filter(!grepl("+", variant, fixed=T))
df_lin_def_snps <- df_lin_def_snps %>% filter(!grepl("Omicron", variant, fixed=T))
df_covariants <- df_lin_def_snps %>% select(variant:pos_orf)
df_lin_def_snps$variant <- gsub("-like", "", df_lin_def_snps$variant)
df_lin_def_snps$variant <- gsub("B.1.617.1", "Kappa (B.1.617.1)", df_lin_def_snps$variant)
sort(unique(df_lin_def_snps$variant))
## from Covariants
df_covariants <- bind_rows(tibble(variant="Alpha (B.1.1.7)", mut=c("S:H69-","S:V70-","S:Y144-","S:N501Y","S:A570D","S:D614G","S:P681H","S:T716I","S:S982A","S:D1118H")),
tibble(variant="Beta (B.1.351)", mut=c("S:L18F","S:D80A","S:D215G","S:L241-","S:L242-","S:A243-","S:K417N","S:E484K","S:N501Y","S:D614G","S:A701V")),
tibble(variant="Gamma (P.1)", mut=c("S:L18F","S:T20N","S:P26S","S:D138Y","S:R190S","S:K417T","S:E484K","S:N501Y","S:D614G","S:H655Y","S:T1027I","S:V1176F")),
tibble(variant="Delta (B.1.617.2)", mut=c("S:T19R","S:E156-","S:F157-","S:R158G","S:L452R","S:T478K","S:D614G","S:P681R","S:D950N")),
tibble(variant="Kappa (B.1.617.1)", mut=c("S:E154K","S:L452R","S:E484Q","S:D614G","S:P681R","S:Q1071H")),
# tibble(variant="21K\n(B.1.1.529)", mut=c("S:A67V","S:H69-","S:V70-","S:T95I","S:G142-","S:V143-","S:Y144-","S:Y145D","S:N211-","S:G339D","S:S371L","S:S373P","S:S375F","S:K417N","S:N440K","S:G446S","S:S477N","S:T478K","S:E484A","S:Q493R","S:G496S","S:Q498R","S:N501Y","S:Y505H","S:T547K","S:D614G","S:H655Y","S:N679K","S:P681H","S:N764K","S:D796Y","S:N856K","S:Q954H","S:N969K","S:L981F")),
tibble(variant="Eta (B.1.525)", mut=c("S:Q52R","S:A67V","S:H69-","S:V70-","S:Y144-","S:E484K","S:D614G","S:Q677H","S:F888L")),
tibble(variant="Iota (B.1.526)", mut=c("S:L5F", "S:T95I", "S:D253G", "S:E484K", "S:D614G", "S:A701V")),
tibble(variant="Lambda (C.37)", mut=c("S:G75V","S:T76I","S:R246-","S:S247-","S:Y248-","S:L249-","S:T250-","S:P251-","S:G252-","S:D253N","S:L452Q","S:F490S","S:D614G","S:T859N")),
tibble(variant="Mu (B.1.621)", mut=c("S:T95I","S:Y144S","S:Y145N","S:R346K","S:E484K","S:N501Y","S:D614G","S:P681H","S:D950N")),
tibble(variant="B.1.620", mut=c("S:P26S","S:H69-","S:V70-","S:V126A","S:Y144-","S:L241-","S:L242-","S:A243-","S:H245Y","S:S477N","S:E484K","S:D614G","S:P681H","S:T1027I","S:D1118H")))
df_covariants$orf <- "Spike"
df_covariants$pos_orf <- sapply(df_covariants$mut, function(x){
	as.numeric(gsub("\\D", "", x))
})
df_vocvoi_snps <- bind_rows(df_lin_def_snps, df_covariants %>% select(-mut))
df_vocvoi_snps <- df_vocvoi_snps %>% arrange(variant, orf, pos_orf)

df_out_other_Var <- df_out %>% select(mut:Position, pos_nt_start:pos_nt_stop)
df_out_other_Var <- lapply(seq_len(nrow(df_out_other_Var)), function(i){
	gene_i <- df_out_other_Var$Gene[i]
	pos_i <- df_out_other_Var$Position[i]
	tmp1 <- df_vocvoi_snps %>% filter(pos_nt_start==df_out_other_Var$pos_nt_start[i]) %>% .$variant
	tmp2 <- df_vocvoi_snps %>% filter(orf==gene_i, pos_orf==pos_i) %>% .$variant
	tmp <- unique(c(tmp1, tmp2))
	return(tibble(mut=df_out_other_Var$mut[i], found_in=tmp))
	# paste0(unique(tmp), collapse = " | ")
})
df_out_other_Var <- bind_rows(df_out_other_Var)
df_out_other_Var <- left_join(df_out_other_Var, df_out %>% select(mut:Position))
df_out_other_Var$Type="VOC/VOI"

df_plot <- df_out %>% pivot_longer((`Case A`:`hCoV-19/South Africa/NICD-N21607-DX64624/2021`), names_to="found_in") %>% filter(value!="")
df_plot$Type="B.1.1.529"
df_plot <- bind_rows(df_plot, df_out_other_Var)
df_plot <- df_plot %>% filter(!is.na(Gene))
df_plot$gene_facet <- factor(df_plot$Gene, levels=levels(df_plot$Gene), labels=gsub("NSP", "", levels(df_plot$Gene)))

df_plot %>% mutate(value="Y") %>% pivot_wider(id_cols=mut:Position, names_from=found_in, values_from=value) %>% write_csv("../results/mut_table_manual.csv", na="") # https://covariants.org/variants/21G.Lambda manual edit afterwards

df_plot2 <- readxl::read_xlsx("../results/mut_table_manual_edit.xlsx")
df_plot2 <- df_plot2 %>% select(-(`hCoV-19/Botswana/R40B59_BHP_3321001248/2021`:`hCoV-19/Hong Kong/VM21044713/2021` | `Case B`)) %>% pivot_longer((`Case A` | `Lambda (C.37)`:`Epsilon (B.1.427/429)`), names_to="found_in") %>% filter(!is.na(value))
df_plot2$Type <- ifelse(grepl("Case", df_plot2$found_in), "", "VOI")

df_plot2$Gene <- factor(df_plot2$Gene, levels=levels(df_plot$Gene))
levels_t <- df_plot2 %>% arrange(Gene, Position) %>% group_by(Gene, Position) %>% summarise(out=unique(Mutation)) %>% .$out
df_plot2$Mutation <- factor(df_plot2$Mutation, levels_t)

df_plot2$VOC <- grepl("Alpha", df_plot2$found_in) | grepl("Beta", df_plot2$found_in) | grepl("Gamma", df_plot2$found_in) | grepl("Delta", df_plot2$found_in)
df_plot2$Type[df_plot2$VOC] <- "VOC"
df_plot2$VOC <- df_plot2$VOC | grepl("Case", df_plot2$found_in)

df_plot2 %>% filter(found_in=="Case A", Gene=="Spike") %>% .$Mutation %>% unique() %>% length()
df_plot2 %>% filter(found_in=="Case A", Gene!="Spike") %>% .$Mutation %>% unique() %>% length()

df_plot2 %>% filter(found_in!="Case A", Gene=="Spike") %>% .$Mutation %>% unique() %>% length()
15/35
df_plot2 %>% filter(Type=="VOC", Gene=="Spike") %>% .$Mutation %>% unique() %>% length()
11/35

p1 <- ggplot(df_plot2)+
	geom_tile(aes(x=Mutation, y=found_in, fill=VOC), show.legend = FALSE)+
	# facet_wrap(vars(Type), ncol=1, scales="free_y")+
	facet_grid(cols=vars(Gene), rows=vars(Type), scales="free", space="free")+
	theme_classic()+
	scale_fill_manual(values=rev(pal_uchicago()(2)))+
	# theme_minimal()+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position='top', strip.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))+
	ylab("Variants / sequences")+
	xlab("Mutations")+
	ggtitle("B")+
	NULL

p1_rotate <- ggplot(df_plot2)+
	geom_tile(aes(y=Mutation, x=found_in, fill=VOC), show.legend = FALSE)+
	# facet_wrap(vars(Type), ncol=1, scales="free_y")+
	facet_grid(rows=vars(Gene), cols=vars(Type), scales="free", space="free", switch="x")+
	theme_classic()+
	scale_fill_manual(values=rev(pal_uchicago()(2)))+
	scale_x_discrete(position = "top") +
	# theme_minimal()+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position='top', strip.text.y=element_text(angle=0, hjust=0.5, vjust=0.5))+
	xlab("")+
	ylab("Mutations")+
	ggtitle("B")+
	NULL

ggsave("../results/common_muts_all.pdf", width=12, height=12/sqrt(2))
save_pptx("../results/common_muts_all.pptx", width=10*sqrt(2)*0.25, height=10)

load("../results/p2.rdata")
p3 <- p2/p1 + plot_layout(heights = c(0.75, 0.25),guides='collect') &
  theme(legend.position='top')
ggsave("../results/Figure 1.pdf", heigh=13, width = 13/sqrt(2))

p3 <- p2_time+p1_rotate + plot_layout(widths = c(0.75, 0.25))
ggsave("../results/Figure 1_timetree.pdf", height=10, width = 10*sqrt(2))
ggsave("../results/Figure 1_timetree.png", height=10, width = 10*sqrt(2))
save_pptx("../results/Figure 1_timetree.pptx", height=10, width = 10*sqrt(2))

# Parse Covsurver results
rst <- readLines("../data/coronamapBlast.pl.html")
idx <- which.max(sapply(rst, nchar))
rst <- rst[idx]
rst_split <- strsplit(rst, "(){ pop.popOut(", fixed=T)[[1]]
length(rst_split)
rst_split <- rst_split[-1]
head(rst_split)

muts <- sapply(rst_split, function(x){
	tmp <- strsplit(x, "</b></font></td></tr><tr><td align=left colspan=2", fixed=T)[[1]][1]
	tmp <- gsub("'<table><tr><td colspan=2><font size =+1><b>", "", tmp, fixed=T)
}, USE.NAMES = F)
df_covsurver <- tibble(muts=muts)

df_covsurver$Occurrence <- sapply(rst_split, function(x){
	tmp <- strsplit(x, "already occurred ", fixed=T)[[1]][2]
	time_tmp <- strsplit(tmp, "countr", fixed=T)[[1]][1]
	time_tmp <- gsub("\\(.+\\) ", "", time_tmp)
	if(is.na(time_tmp)){return(NA)}
	paste0(time_tmp, "countries")
	# strsplit(rst_split, " times")[[1]][1]
}, USE.NAMES = F)
df_covsurver$Occurrence <- gsub("timein", "time in", df_covsurver$Occurrence)
df_covsurver$Occurrence <- gsub("1 countries", "1 country", df_covsurver$Occurrence)

df_covsurver$Frequency <- sapply(rst_split, function(x){
	tmp <- strsplit(x, "already occurred ", fixed=T)[[1]][2]
	time_tmp <- strsplit(tmp, "countr", fixed=T)[[1]][1]
	time_tmp <- strsplit(time_tmp, " (", fixed=T)[[1]][2]
	time_tmp <- strsplit(time_tmp, " ", fixed=T)[[1]][1]
}, USE.NAMES = F)

df_covsurver$First_collected <- sapply(rst_split, function(x){
	tmp <- strsplit(x, " collected in ", fixed=T)[[1]][2]
	time_tmp <- strsplit(tmp, ",", fixed=T)[[1]][1]
	time_tmp <- strsplit(time_tmp, ".", fixed=T)[[1]][1]
	gsub("\\(.+\\) ", "", time_tmp)
}, USE.NAMES = F)

df_covsurver$color <- sapply(rst_split, function(x){
	tmp <- strsplit(x, "color=", fixed=T)[[1]]
	tmp <- tmp[length(tmp)]
	color_tmp <- strsplit(tmp, ">", fixed=T)[[1]][1]
}, USE.NAMES = F)
unique(df_covsurver$color)
df_covsurver$annotation <- factor(df_covsurver$color, levels=c("black", "blue", "orange", "cyan"), labels=c("No known effects", "Occurring > 100 times",  "Involved in phenotypic effects", "INDELs"))

df_covsurver$mut <- gsub(" ", "_", df_covsurver$muts)

df_out2 <- left_join(df_out %>% select(mut:Position), df_covsurver %>% select(-muts))
# df_out2 <- df_out2 %>% select(mut:Position, Other_varaints:annotation, everything())
writexl::write_xlsx(df_out2 %>% select(-mut, -Position, -color), "../results/mut_check.xlsx")

df_out2 %>% filter(Gene=="Spike" & nchar(Other_varaints)>0)
df_out2 %>% filter(Gene=="Spike" & nchar(Other_varaints)==0) %>% .$annotation %>% table()
df_out2 %>% filter(Gene=="Spike")
16/37*100

other_var_t <- df_out2 %>% filter(Gene=="Spike" & nchar(Other_varaints)>0) %>% .$Other_varaints
tmp <- sort(table(unlist(strsplit(other_var_t, " | ", fixed=T))))
df_vocvoi_num <- tibble(variant=names(tmp), Num_of_common_spike_mutation=tmp)
other_var_t <- df_out2 %>% filter(Gene!="Spike" & nchar(Other_varaints)>0) %>% .$Other_varaints
tmp <- sort(table(unlist(strsplit(other_var_t, " | ", fixed=T))))
df_vocvoi_num <- left_join(df_vocvoi_num, tibble(variant=names(tmp), `Num_of_common_non-spike_mutation`=tmp))
writexl::write_xlsx(df_vocvoi_num, "../results/supplementary table 2.xlsx")

mean_ocu_voc <- df_out2 %>% filter(Gene=="Spike" & nchar(Other_varaints)>0) %>% .$Occurrence
mean_ocu_voc <- sapply(mean_ocu_voc, function(x) {
   as.numeric(strsplit(x, " time")[[1]][1])
})
mean_ocu_nonvoc<- df_out2 %>% filter(Gene=="Spike" & nchar(Other_varaints)==0) %>% .$Occurrence
mean_ocu_nonvoc <- sapply(mean_ocu_nonvoc, function(x) {
   as.numeric(strsplit(x, " time")[[1]][1])
})
mean(mean_ocu_voc)
mean(mean_ocu_nonvoc, na.rm=T)
wilcox.test(mean_ocu_voc, mean_ocu_nonvoc)

# Plot supplementary Figure 1
df_plot <- df_out2 %>% filter(!is.na(First_collected))
df_plot$`In other VOC/VOI` <- nchar(df_plot$Other_varaints)>0
df_plot$gene_facet <- factor(df_plot$Gene, levels=levels(df_plot$Gene), labels=gsub("NSP", "", levels(df_plot$Gene)))
ggplot(df_plot)+
	geom_col(aes(x=Mutation, y=First_collected, fill=`In other VOC/VOI`))+	
	scale_fill_jco()+
	facet_grid(cols=vars(gene_facet), scales="free_x", space="free_x")+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position='top')+
	ylab("First collected (GISAID)")+
	xlab("Mutations of the B.1.1.529 variant")+
	NULL

ggsave("../results/first_collected.pdf", width = 12, height = 6)
