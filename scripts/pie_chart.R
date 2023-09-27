library(ggplot2)
library(dplyr)
panel_info <- read.table('/users/mscherer/cluster/project/Methylome/infos/HSCs/panel_info_dropout_pwm.tsv',
                         sep='\t')
plot_theme <- theme(panel.background = element_blank(),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=20),
                    axis.text=element_blank(),
                    axis.ticks=element_line(color='black', size=.1),
                    strip.background = element_blank(),
                    legend.key=element_rect(color='black', fill=NA),
                    legend.key.size = unit(5, 'mm'),
                    plot.title=element_blank())
to_plot <- panel_info[, 'Type']
to_plot[grepl('high', to_plot)] <- 'DMC'
to_plot <- plyr::count(to_plot)
colnames(to_plot) <- c('Type', 'Count')
to_plot$Type <- factor(to_plot$Type,
                       levels=c('Non_cut',
                                'WSH',
                                'Always_methylated',
                                'DMC',
                                'Always_unmethylated',
                                'IMR'))
to_plot <- to_plot[order(to_plot$Count), ]
to_plot$prop <- to_plot$Count/sum(to_plot$Count) *100
to_plot$ypos <- cumsum(to_plot$prop)-0.5*to_plot$prop
plot <- ggplot(to_plot, aes(x="", y=Count, fill=Type))+geom_bar(stat="identity", width=1, color="white")+
  coord_polar("y", start=0)+plot_theme+ylab("")+xlab("")+geom_text(aes(label = Count),
                                                      position = position_stack(vjust = 0.5),
                                                      size=7.5)+
  scale_fill_manual(values=c('WSH'='#bf8fbf',
                             'IMR'='#bfbf90',
                             'DMC'='#b38456',
                             'Non_cut'='gray90',
                             'Always_methylated'='gray75',
                             'Always_unmethylated'='gray50'))
ggsave('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/plots/amplicon_type_distribution.pdf',
       plot)

tam <- readRDS('/users/mscherer/cluster/project/Methylome/analysis/HSCs/paper/Fig1/old/seurat_lars.RDS')
to_plot <- panel_info[rownames(tam), 'Type']
to_plot[grepl('high', to_plot)] <- 'DMC'
to_plot <- plyr::count(to_plot)
colnames(to_plot) <- c('Type', 'Count')
to_plot$Type <- factor(to_plot$Type,
                       levels=c('Non_cut',
                                'WSH',
                                'Always_methylated',
                                'DMC',
                                'Always_unmethylated',
                                'IMR'))
to_plot <- to_plot[order(to_plot$Count), ]
to_plot$prop <- to_plot$Count/sum(to_plot$Count) *100
to_plot$ypos <- cumsum(to_plot$prop)-0.5*to_plot$prop
plot <- ggplot(to_plot, aes(x="", y=Count, fill=Type))+geom_bar(stat="identity", width=1, color="white")+
  coord_polar("y", start=0)+plot_theme+theme(legend.position='none')+
  scale_fill_manual(values=c('WSH'='#bf8fbf',
                             'IMR'='#bfbf90',
                             'DMC'='#b38456',
                             'Non_cut'='gray90',
                             'Always_methylated'='gray75',
                             'Always_unmethylated'='gray50'))
ggsave('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/plots/amplicon_type_distribution_selected.pdf',
       plot)