# Figure (Supplementary): Timeline of the two cases
# Figure (Main). Phylogenetic tree of the two cases

library(tidyverse)
library(Biostrings)
library(ggtree)
library(ggsci)
library(scales)
library(lubridate)
library(shadowtext)
library(ggrepel)
library(patchwork)
source("https://raw.githubusercontent.com/Koohoko/Save-ggplot-to-pptx/main/scripts/save_pptx.r")

# Figure (Supplementary): Timeline of the two cases
df_tl <- readxl::read_excel("../data/B.1.1.529_HK.xlsx",2)
df_tl$Name <- gsub("\\", "_", df_tl$Name, fixed=T)
df_tl$Name <- gsub("_n", "\n", df_tl$Name, fixed=T)
df_tl$Date <- lubridate::ymd(df_tl$Time)
df_tl$Type <- NA
df_tl$Type[grepl("^RT-PCR", df_tl$Name)] <- "Tests during quarantine"
df_tl$Type[grepl("^Onset", df_tl$Name)] <- "Symptoms"
df_tl$Type[is.na(df_tl$Type)] <- "Arrival and admission"
df_tl <- df_tl %>% arrange(`Group ID`, Date)
df_tl$direction <- ifelse(df_tl$`Group ID`=="A", 1, -1)
df_tl$position <- ifelse(df_tl$Type=="Tests during quarantine", 0.3, 0.8)
df_tl$position[df_tl$Type=="Symptoms"] <- 0.6
df_tl$position <- df_tl$position*df_tl$direction

text_offset <- 0.1
df_tl <- df_tl %>% group_by(`Group ID`) %>% mutate(date_count=seq(n())) %>% ungroup()
df_tl$text_position <- (df_tl$date_count * text_offset * df_tl$direction) + df_tl$position


df <- df_tl 
timeline_plot<-ggplot(df,aes(x=Date,y=0, col=Type, label=Name))
timeline_plot<-timeline_plot+labs(col="Events")
# timeline_plot<-timeline_plot+scale_color_manual(values=status_colors, labels=status_levels, drop = FALSE)
timeline_plot<-timeline_plot+scale_color_jco()
timeline_plot<-timeline_plot+theme_classic()

# Plot horizontal black line for timeline
timeline_plot<-timeline_plot+geom_hline(yintercept=0, 
                color = "black", size=0.3)

# Plot vertical segment lines for milestones
timeline_plot<-timeline_plot+geom_segment(data=df, aes(y=text_position,yend=0,xend=Date), color='black', size=0.2)

# Plot scatter points at zero and date
timeline_plot<-timeline_plot+geom_point(aes(y=0), size=4, alpha=0.8, position=position_dodge(width = -0.2))

# Don't show axes, appropriately position legend
timeline_plot<-timeline_plot+theme(axis.line.y=element_blank(),
                 axis.text.y=element_blank(),
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.text.x =element_blank(),
                 axis.ticks.x =element_blank(),
                 axis.line.x =element_blank(),
                 legend.position = "bottom"
                )

# Show text for each day
day_buffer <- 2
day_date_range <- seq(min(df$Date) - days(day_buffer), max(df$Date) + days(day_buffer), by='day')
day_format <- format(day_date_range, '%b-%d')
day_df <- data.frame(day_date_range, day_format)

timeline_plot<-timeline_plot+geom_text(data=day_df, aes(x=day_date_range,y=-0.1,label=day_format),size=2.5,vjust=0.5, color='black')

# Show text for each milestone
p1<-timeline_plot+
	geom_shadowtext(aes(y=text_position, label=Name), size=3, bg.color="white", data=. %>% filter(!grepl("^Arr", Type)), show.legend = FALSE)+
	geom_label(aes(y=text_position, label=Name), size=3, bg.color="white", box.padding=0.1, data=. %>% filter(grepl("^Arr", Type)), show.legend = FALSE)+
	guides(color = guide_legend(override.aes = list(alpha=1)))+ 
	NULL

p1 <- p1+
	geom_text(x=ymd("2021-11-08"),y=-0.5,label="Case B", color="Black")+
	geom_text(x=ymd("2021-11-08"),y=0.5,label="Case A", color="Black")+
	# ggtitle("A")+
	NULL
ggsave("../results/timeline.pdf", height = 5, width = 5*sqrt(2))
ggsave("../results/timeline.png", height = 5, width = 5*sqrt(2))

# Figure (Main). Phylogenetic tree of the two cases
## determine helper sequences
df_seq_meta <- read_tsv("../data/metadata.tsv")
df_seq_meta_selected <- df_seq_meta %>% arrange(date) %>% group_by(pango_lineage) %>% mutate(id=seq(n())) %>% filter(id==1) %>% ungroup()
df_seq_meta_selected <- df_seq_meta_selected %>% filter(grepl("^B$", pango_lineage) | grepl("^B.1.1", pango_lineage) | grepl("^B.1.617", pango_lineage)| grepl("^B.1.351", pango_lineage)| grepl("^P.1", pango_lineage)| grepl("^P.2", pango_lineage)| grepl("^P.3", pango_lineage)| grepl("^B.1.621", pango_lineage)| grepl("^C.37", pango_lineage)| grepl("^B.1.526", pango_lineage)| grepl("^B.1.525", pango_lineage)| grepl("^B.1.427", pango_lineage)| grepl("^B.1.429", pango_lineage))
seqs_ref <- readDNAStringSet("../data/aligned.fasta")
seqs_ref <- seqs_ref[names(seqs_ref) %in% df_seq_meta_selected$strain]

## HK  b.1.1.529 sequences
seqs <- readDNAStringSet("../data/seqs_20211125.fasta")
seqs_hk <- seqs[grepl("12388", names(seqs)) | grepl("12404", names(seqs)) ]
names(seqs_hk)[grepl("12388", names(seqs_hk))] <- "Case_A"
names(seqs_hk)[grepl("12404", names(seqs_hk))] <- "Case_B"

writeXStringSet(c(seqs[1], seqs_hk), "../results/seqs_hk_toalign.fasta") ## compare with the Sanger results
file_in <- "../results/sanger_seqs_add.fasta"
file_out <- "../results/seqs_hk_sanger.fasta"
system(paste0("mafft --localpair --maxiterate 1000 --thread 8 --keeplength --addfragments ", file_in, " ../results/ref_seq.fasta > ", file_out)) ## visual inspection
# system(paste0("mafft --localpair --maxiterate 1000 --thread -8 --keeplength --addfragments ", "../results/sanger_seqs_long.fasta", " ../results/ref_seq.fasta > ", "../results/seqs_hk_sanger_long.fasta")) ## visual inspection

seqs_edit <- readDNAStringSet("../results/../results/seqs_hk_sanger_edit.fasta")
seqs_hk <- seqs_edit[grepl("filled", names(seqs_edit))]
names(seqs_hk) <- gsub("_filled", "", names(seqs_hk))

seqs_sub <- c(seqs_ref, seqs_hk)
file_seqs_sub <- "../results/seqs_sub.fasta"
writeXStringSet(seqs_sub, file_seqs_sub)
write_csv(df_seq_meta_selected, "../results/df_seq_sub.csv")

system(paste0("mafft --auto --thread 8 --keeplength --addfragments ", "../data/Omicron_gisaid_hcov-19_2021_12_01_03.fasta", " ../results/ref_seq.fasta > ", "../results/Omicron_global_aligned.fasta")) 
seqs_omic <- readDNAStringSet("../results/Omicron_global_aligned.fasta")
seqs_sub_more <- c(seqs_ref, seqs_hk, seqs_omic)
names(seqs_sub_more)[grepl("Hong_Kong", names(seqs_sub_more))]
seqs_sub_more <- seqs_sub_more[!grepl("Hong_Kong", names(seqs_sub_more))]
file_seqs_sub_more <- "../results/seqs_sub_more.fasta" 
writeXStringSet(seqs_sub_more, file_seqs_sub_more)

df_date <- df_seq_meta_selected %>% select(strain, date)
df_date <- bind_rows(df_date, tibble(strain=c("Case_A", "Case_B"), date=ymd(c("2021-11-13", "2021-11-18"))))
write_tsv(df_date, "../results/df_date.tsv", col_names=F)

system(paste0('~/Documents/iqtree-2.1.3-MacOSX/bin/iqtree2 -T 16 -m GTR+F+R2 --redo -s ', file_seqs_sub, ' -o "Wuhan-Hu-1/2019" --date ../results/df_date.tsv --date-root 2019-12-26 --date-ci 100 --date-options "-l -1"'))
# system(paste0('~/Documents/iqtree-2.1.3-MacOSX/bin/iqtree2 -T 16 -m JC+FQ --redo -s ', file_seqs_sub_more, ' -o "Wuhan-Hu-1/2019"'))

tree <- read.tree("../results/nextstrain__tree.nwk")
p <- ggtree(tree, size=0.3)
df_seq_meta$label <- df_seq_meta$strain
df_seq_meta$label2 <- df_seq_meta$pango_lineage
p$data <- left_join(p$data, df_seq_meta)
p$data$lineage_sim <- ifelse(grepl("^Case_", p$data$label), "B", "Others")
p$data$lineage_sim[p$data$label2 %in% c("B.1.1.7", "B.1.351", "B.1.617.2", "P.1")] <- "B"
p$data$label <-  gsub("Case_", "Case ", p$data$label, fixed=T)
p$data$label2[is.na(p$data$label2) & p$data$isTip] <- p$data$label[is.na(p$data$label2) & p$data$isTip]
p$data$label2[grepl("^Case", p$data$label2)] <- NA

p$data$label2[grepl("^EPI", p$data$label2)] <- sapply(p$data$label2[grepl("^EPI", p$data$label2)], function(x){
	names(seqs_sub_more)[grepl(x, names(seqs_sub_more))]
})
# p$data$text_size <- ifelse(grepl("^Case_", p$data$label), 0.6, 0.2)

font_size=2.2

p2 <- p + 
	geom_tiplab(aes(label=label2), show.legend = FALSE, size=1, color="black")+
	geom_text_repel(aes(label=label), size=5, data=. %>% filter(grepl("^Case", label)), color=pal_uchicago()(1), nudge_x=-5, nudge_y=2, alpha=0.9)+
	geom_tippoint(aes(color=lineage_sim), show.legend = FALSE, alpha=0.8, size=0.3)+
	scale_color_uchicago()+
	theme(legend.position = "bottom")+
	xlim(0, 80)+
	# ggtitle("A")+
	NULL
ggsave("../results/tree.pdf", heigh=7.5*sqrt(2), width = 7.5)
save_pptx("../results/tree.pptx", height=7.5*sqrt(2), width = 7.5)

timetree <- treeio::read.beast("../results/seqs_sub.fasta.timetree.nex")
mrsd <- max(df_date$date, na.rm = T)
p <- ggtree(timetree, mrsd=mrsd, size=0.2) + theme_tree2()
df_seq_meta$label <- df_seq_meta$strain
df_seq_meta$label2 <- df_seq_meta$pango_lineage
p$data <- left_join(p$data, df_seq_meta, "label")
p$data$lineage_sim <- ifelse(grepl("^Case_", p$data$label), "B", "Others")
p$data$lineage_sim[p$data$label2 %in% c("B.1.1.7", "B.1.351", "B.1.617.2", "P.1")] <- "B"
p$data$label <-  gsub("Case_", "Case ", p$data$label, fixed=T)
# p$data$label2[grepl("^Case_", p$data$label)] <- p$data$label[grepl("^Case_", p$data$label)]
# p$data$text_size <- ifelse(grepl("^Case_", p$data$label), 0.6, 0.2)

labels_t <- c("2020-01-01", "2020-07-01", "2021-01-01", "2021-07-01", "2021-12-01")
breaks_t <- Date2decimal(labels_t)

p2_time <- p + 
	geom_tiplab(aes(label=label2), show.legend = FALSE, size=font_size, color="black")+
	# geom_text_repel(aes(label=label2, color=lineage_sim), show.legend = FALSE, size=2)+
	geom_text_repel(aes(label=label), size=5, data=. %>% filter(grepl("^Case", label)), color=pal_uchicago()(1), nudge_x=0.3, alpha=0.9)+
	geom_tippoint(aes(color=lineage_sim), show.legend = FALSE, alpha=0.8, size=0.3)+
	scale_color_uchicago()+
	theme(legend.position = "bottom")+
	scale_x_continuous(limits=c(2019.9, 2022.5), breaks=breaks_t, labels=labels_t)+
	# xlim(0, 0.003)+
	ggtitle("A")+
	NULL

ggsave("../results/timetree.pdf", height=10, width = 6)
save_pptx("../results/timetree.pptx", width=7.5*sqrt(2)*0.75, height=7.5, plot=p2_time)
save(p2, p2_time, file="../results/p2.rdata")

