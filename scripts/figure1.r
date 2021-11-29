# Figure 1A: Timeline of the two cases
# Figure 2B. Phylogenetic tree of the two cases

library(tidyverse)
library(Biostrings)
library(ggtree)
library(ggsci)
library(scales)
library(lubridate)
library(shadowtext)
library(ggrepel)
library(patchwork)

# Figure 1A 
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

p1 <- p1+geom_text(x=ymd("2021-11-08"),y=-0.5,label="Case B", color="Black")+geom_text(x=ymd("2021-11-08"),y=0.5,label="Case A", color="Black")+ggtitle("A")
ggsave("../results/timeline.pdf", height = 5, width = 5*sqrt(2))


# Figure 1B
seqs <- readDNAStringSet("../data/seqs_20211125.fasta")
df_lin <- read_csv("../data/lineage_20211126.csv")
cov_check <- apply(alphabetFrequency(seqs), 1, function(x){
	sum(x[1:4])
})
cov_check <- round(cov_check/29903*100, 2)
df_seqs_ori <- tibble(taxon=names(seqs), coverage=cov_check)
df_seqs_ori <- left_join(df_seqs_ori, df_lin)
df_seqs <- df_seqs_ori %>% filter(grepl("WHP", taxon)) %>% arrange(coverage)
df_seqs <- df_seqs %>% group_by(lineage) %>% mutate(idx=seq(n())) %>% ungroup()
df_seqs_sub_ref <- df_seqs %>% filter(idx==1 & lineage!="None")
df_seqs_sub_ref$label2 <- paste0("HK-case_", df_seqs_sub_ref$lineage)
df_seqs_sub_int <- df_seqs_ori %>% filter(lineage=="B.1.1.529" | grepl("12388", taxon) | grepl("12404", taxon) | grepl("MN908947", taxon))
df_seqs_sub_int$label2 <- df_seqs_sub_int$taxon
df_seqs_sub_int$label2[grepl("12388", df_seqs_sub_int$taxon)] <- "Case A"
df_seqs_sub_int$label2[grepl("12404", df_seqs_sub_int$taxon)] <- "Case B"

df_sub_all <- bind_rows(df_seqs_sub_ref, df_seqs_sub_int)
df_sub_all <- unique(df_sub_all)
df_sub_all$label <- df_sub_all$taxon
df_sub_all$lineage_sim <- "Others"
df_sub_all$lineage_sim[grepl("B.1.1.529", df_sub_all$lineage) | grepl("^Case", df_sub_all$label2)] <- "B.1.1.529"
df_sub_all$label2[grepl("MN908947", df_sub_all$taxon)] <- "Wuhan-Hu-01"
seqs_sub <- seqs[names(seqs) %in% df_sub_all$taxon]
file_seqs_sub <- "../results/seqs_sub.fasta"
writeXStringSet(seqs_sub, file_seqs_sub)

system(paste0('~/softwares/iqtree-2.1.3-MacOSX/bin/iqtree2 -T 16 --redo -s ', file_seqs_sub, ' -o "MN908947_3"')) # using iqtree

tree <- read.tree("../results/seqs_sub.fasta.treefile")
p <- ggtree(tree, size=0.5)
p$data <- left_join(p$data, df_sub_all)

p2 <- p + 
	geom_tiplab(aes(label=label2, color=lineage_sim), show.legend = FALSE, size=1.8)+
	geom_tippoint(aes(color=lineage_sim), show.legend = FALSE, alpha=0.8)+
	scale_color_uchicago()+
	theme(legend.position = "bottom")+
	xlim(0, 0.003)+
	ggtitle("B")

ggsave("../results/tree.pdf", heigh=10, width = 6)

p3 <- p1/p2 + plot_layout(heights = c(0.3, 0.7),guides='collect') &
  theme(legend.position='top')
ggsave("../results/Figure 1.pdf", heigh=13, width = 13/sqrt(2))
