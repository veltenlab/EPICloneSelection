library(ggplot2)
library(ggsci)
dat <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/MB_output/Designer/4135-design-summary_annotated.csv')
types <- gsub('[0-9]', '', dat$Type)
types[grep('high', types)] <- 'DMC'
types[types%in%'IMR'] <- 'IMC'
to_plot <- plyr::count(types)
colnames(to_plot) <- c('Type','Count')
plot_theme <- theme(panel.background = element_blank(),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=10),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_blank(),
                    strip.background = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    #axis.text.x=element_text(angle=45, hjust=1, vjust = 1),
                    axis.text.x=element_blank(),
                    axis.title.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    legend.title=element_blank(),
                    panel.spacing = unit(.1, "lines"))
to_plot$Type <- factor(to_plot$Type, levels=c('Non_cut',
                                              'Always_methylated',
                                              'Always_unmethylated',
                                              'DMC',
                                              'IMC',
                                              'WSH'))
to_plot <- to_plot[order(to_plot$Type), ]
to_plot$prop <- to_plot$Count/sum(to_plot$Count) *100
to_plot$ypos <- cumsum(to_plot$prop)-0.5*to_plot$prop
plot <- ggplot(to_plot, aes(x="", y=Count, fill=Type))+geom_bar(stat="identity", width=1, color="white")+
  coord_polar("y", start=0)+plot_theme+scale_fill_tron()+ylab("")+xlab("")+geom_text(aes(label = Count),
                                                                                     position = position_stack(vjust = 0.5))
ggsave('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/plots/amplicon_type_distribution.png',
       height=150,
       width=150,
       unit='mm')

plot_theme <- theme(panel.background = element_blank(),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=10),
                    axis.text.x = element_text(angle=45, hjust=1, color='black', size=8),
                    strip.background = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    #axis.text.x=element_text(angle=45, hjust=1, vjust = 1),
                    legend.title=element_blank(),
                    panel.spacing = unit(.1, "lines"))

dat <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/MB_output/Designer/4135-design-summary_annotated.csv')
types <- gsub('[0-9]', '', dat$Type)
types[grep('high', types)] <- 'DMC'
types[types%in%'IMR'] <- 'IMC'
dat$Type <- types
tf_names <- colnames(dat)
start_tf <- which(tf_names%in%'CpGCount')+1
end_tf <- which(tf_names%in%"Scl")
tf_names <- tf_names[start_tf:end_tf]
counts_tfs <- unlist(lapply(tf_names, function(x){
  sum(!is.na(dat[, x]))
}))
to_plot <- data.frame(Count=counts_tfs, TF=tf_names)
to_plot$TF <- factor(to_plot$TF, levels=to_plot$TF[order(to_plot$Count, decreasing=T)])
plot <- ggplot(to_plot, aes(x=TF, y=Count))+geom_bar(stat="identity", width=1, color="white")+
  plot_theme
ggsave('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/plots/amplicon_TF_all_distribution.png',
       height=75,
       width=150,
       unit='mm')

to_plot <- c()
for(type in unique(dat$Type)){
  sel_dat <- subset(dat, subset=Type==type)
  counts_tfs <- unlist(lapply(tf_names, function(x){
    sum(!is.na(sel_dat[, x]))
  }))
  to_plot <- rbind(to_plot, data.frame(Type=type, TF=factor(tf_names, levels=tf_names[order(counts_tfs, decreasing=TRUE)]), Count=counts_tfs))
}
plot <- ggplot(to_plot, aes(x=TF, y=Count))+geom_bar(stat="identity", width=1, color="white")+
  plot_theme+facet_wrap(Type~., ncol=3)
ggsave('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/plots/amplicon_TF_stratified.png',
       height=150,
       width=300,
       unit='mm')
