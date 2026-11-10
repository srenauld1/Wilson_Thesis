suppressPackageStartupMessages({library(jsonlite);library(ggplot2)})

args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else 'outputs/ves041_blinded_good_walkers_pivots_gt30_20261102'
out <- file.path(root,'talk_figure')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
f <- fromJSON(file.path(root,'fly_summary.json'))
d <- f[f$definition=='primary' & f$subset=='all',]
stopifnot(nrow(d)>0,all(d$cohort=='blinded round 1'),!anyDuplicated(d$fly_key))
control_label <- sprintf('Control (n = %d)',sum(d$group=='control'))
experimental_label <- sprintf('Experimental (n = %d)',sum(d$group=='experimental'))
d$panel <- ifelse(d$group=='control',control_label,experimental_label)
d$panel <- factor(d$panel,levels=c(control_label,experimental_label))
pre <- d;pre$period <- 'Pre-opto';pre$rate <- pre$pre_rate_per_min
on <- d;on$period <- 'Opto-on';on$rate <- on$inhibition_rate_per_min
long <- rbind(pre,on)
long$period <- factor(long$period,levels=c('Pre-opto','Opto-on'))
means <- aggregate(rate~panel+period+group,data=long,FUN=mean)
means$mean_label <- sprintf('%.2f',means$rate)
ymax <- max(8.5,ceiling(max(long$rate)/2)*2+.5)

colors <- c(control='#41464D',experimental='#BA354F')
p <- ggplot(long,aes(period,rate,group=fly_key))+
  geom_line(aes(color=group),alpha=.23,linewidth=.65)+
  geom_point(aes(color=group),alpha=.36,size=2.3)+
  geom_line(data=means,aes(group=group,color=group),linewidth=1.9)+
  geom_point(data=means,aes(group=group,fill=group),shape=21,
             color='white',stroke=1.1,size=5.3)+
  facet_wrap(~panel,nrow=1)+
  scale_color_manual(values=colors)+scale_fill_manual(values=colors)+
  scale_x_discrete(expand=expansion(add=.65))+
  scale_y_continuous(limits=c(0,ymax),breaks=seq(0,ymax,2),expand=expansion(mult=c(0,.01)))+
  labs(title='Reorientations during VES041 inhibition',
       subtitle=NULL,x=NULL,y='Reorientations / min',caption='Bold lines show means')+
  theme_classic(base_size=21,base_family='Arial')+
  theme(
    legend.position='none',
    plot.title=element_text(size=25,face='bold',margin=margin(b=22)),
    plot.title.position='plot',
    strip.background=element_blank(),
    strip.text=element_text(size=22,face='bold',margin=margin(t=4,b=20)),
    panel.spacing=grid::unit(1.3,'cm'),
    axis.title.y=element_text(size=21,margin=margin(r=16)),
    axis.text=element_text(size=19,color='#252525'),
    axis.text.x=element_text(margin=margin(t=12)),
    axis.line=element_line(linewidth=.5,color='#555555'),
    axis.ticks=element_line(linewidth=.5,color='#555555'),
    panel.grid.major.y=element_line(color='#E8EAED',linewidth=.35),
    plot.caption=element_text(size=12,color='#777777',hjust=1,margin=margin(t=20)),
    plot.margin=margin(20,28,15,20),
    plot.background=element_rect(fill='white',color=NA))

ggsave(file.path(out,'ves041_pivot_frequency_talk.png'),p,width=12.8,height=7.2,
       dpi=300,bg='white')
ggsave(file.path(out,'ves041_pivot_frequency_talk.svg'),p,width=12.8,height=7.2,
       bg='white')
cat('Saved talk figure with group means:\n')
print(means[,c('panel','period','rate')],row.names=FALSE)
