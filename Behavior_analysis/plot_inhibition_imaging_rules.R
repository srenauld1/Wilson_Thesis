suppressPackageStartupMessages({library(ggplot2);library(jsonlite)})
root <- 'outputs/ves041_inhibition_imaging_rules_20261104'
d <- fromJSON(file.path(root,'fly_summary.json'))
d <- d[d$definition=='primary' & d$subset %in% c('all','quality_subset'),]
stopifnot(nrow(d)==45)
d$panel <- paste(d$subset,d$group,sep='_')
levels <- c('quality_subset_control','quality_subset_experimental','all_control','all_experimental')
labels <- c('Good walkers\nControl (n = 9)','Good walkers\nExperimental (n = 8)',
            'All usable flies\nControl (n = 17)','All usable flies\nExperimental (n = 11)')
d$panel <- factor(d$panel,levels=levels,labels=labels)
pre <- d;pre$period <- 'Pre-opto';pre$rate <- pre$pre_rate_per_min
on <- d;on$period <- 'Opto-on';on$rate <- on$inhibition_rate_per_min
long <- rbind(pre,on);long$period <- factor(long$period,levels=c('Pre-opto','Opto-on'))
summary <- do.call(rbind,lapply(split(long,list(long$panel,long$period),drop=TRUE),function(q){
 data.frame(panel=q$panel[1],period=q$period[1],group=q$group[1],mean=mean(q$rate),sem=sd(q$rate)/sqrt(nrow(q)))
}))
colors <- c(control='#555555',experimental='#BA354F')
p <- ggplot(long,aes(period,rate,group=fly_key,color=group))+
 geom_line(alpha=.25,linewidth=.55)+geom_point(alpha=.4,size=1.8)+
 geom_errorbar(data=summary,aes(y=mean,ymin=mean-sem,ymax=mean+sem,group=group),width=.13,linewidth=.75)+
 geom_line(data=summary,aes(y=mean,group=group),linewidth=1.5)+
 geom_point(data=summary,aes(y=mean,group=group,fill=group),shape=21,color='white',stroke=.8,size=4)+
 facet_wrap(~panel,ncol=2)+scale_color_manual(values=colors)+scale_fill_manual(values=colors)+
 scale_x_discrete(expand=expansion(add=.55))+
 scale_y_continuous(expand=expansion(mult=c(.03,.09)))+expand_limits(y=0)+
 labs(title='VES041 inhibition: imaging-defined reorientations',
  subtitle='Each line is one fly. Bold lines and error bars show mean +/- SEM.',
  x=NULL,y='Slowdown-associated reorientations / min',
  caption='Matched 5-s pre-opto and opto-on periods. >=15 deg accumulated rotation within 0.5 s after a qualifying slowdown.')+
 theme_classic(base_size=17,base_family='Arial')+
 theme(legend.position='none',strip.background=element_blank(),strip.text=element_text(size=18,face='bold',lineheight=1.15),
  panel.spacing=grid::unit(1,'cm'),plot.title=element_text(size=25,face='bold'),
  plot.subtitle=element_text(size=15,margin=margin(b=18)),axis.text=element_text(color='#333333'),
  panel.grid.major.y=element_line(color='#E9EBED',linewidth=.3),plot.caption=element_text(size=11,color='#666666',margin=margin(t=15)),
  plot.margin=margin(18,20,12,15))
stem <- file.path(root,'imaging_defined_reorientations_paired')
ggsave(paste0(stem,'.svg'),p,width=12.5,height=10.4,device=svglite::svglite,bg='white')
ggsave(paste0(stem,'.png'),p,width=12.5,height=10.4,dpi=180,bg='white')
print(summary)
