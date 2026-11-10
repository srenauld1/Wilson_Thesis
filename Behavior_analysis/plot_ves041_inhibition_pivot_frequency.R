suppressPackageStartupMessages({library(jsonlite);library(ggplot2)})
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else 'outputs/ves041_inhibition_pivot_frequency_20261102'
parameters_path <- file.path(root,'parameters.json')
parameters <- if(file.exists(parameters_path)) fromJSON(parameters_path) else list(primary_angle_deg=45,strict_angle=FALSE)
angle_label <- paste0(ifelse(parameters$strict_angle,'>','>='),parameters$primary_angle_deg)
selection_label <- if(isTRUE(parameters$good_walkers)) 'Previously identified good walkers only.' else 'All readable, genotype-assigned flies are included.'
all <- fromJSON(file.path(root,'fly_summary.json'))
all$group_short <- ifelse(grepl('^experimental',all$group),'Experimental','Control')
all$cohort_label <- paste(all$cohort, paste0(all$duration_s,' s'))
colors <- c(Experimental='#B9364D',Control='#555555')
theme_set(theme_bw(base_size=11))
d <- all[all$definition=='primary' & all$subset=='all' &
  ((all$cohort=='blinded round 1' & all$duration_s==5) |
   (all$cohort=='unblinded RubyACR' & all$duration_s==10)), ]
pre <- d; pre$period <- 'Pre-light';pre$rate <- pre$pre_rate_per_min
on <- d;on$period <- 'Inhibition';on$rate <- on$inhibition_rate_per_min
long <- rbind(pre,on);long$period <- factor(long$period,levels=c('Pre-light','Inhibition'))
p <- ggplot(long,aes(period,rate,group=fly_key,color=group_short)) +
  geom_line(alpha=.45,linewidth=.55)+geom_point(size=2,alpha=.8) +
  facet_grid(group_short~cohort_label) + scale_color_manual(values=colors) +
  labs(title='Pivot-like bout frequency before and during VES041 inhibition',
       subtitle='Each line is one fly; equal-duration pre-light and inhibition windows',
       x=NULL,y='Pivot-like bouts / min',color=NULL,
       caption=paste0('Primary definition: ',angle_label,' degrees directed rotation within 1 s, <=2 mm net displacement, directionality >=0.7.\nOverlapping detections are merged into a single bout. ',selection_label)) +
  theme(legend.position='none',plot.caption=element_text(hjust=0))
ggsave(file.path(root,'paired_pivot_frequency.png'),p,width=10,height=7.5,dpi=200)

p <- ggplot(d,aes(group_short,delta_rate_per_min,color=group_short))+
  geom_hline(yintercept=0,linetype=2,color='grey50')+
  geom_boxplot(outlier.shape=NA,width=.4,alpha=.08)+
  geom_point(position=position_jitter(width=.09,height=0,seed=41),size=2)+
  facet_wrap(~cohort_label)+scale_color_manual(values=colors)+
  labs(title='Within-fly change in pivot-like bout frequency',
       subtitle='Negative values indicate fewer bouts during inhibition',
       x=NULL,y='Inhibition minus pre-light (bouts / min)')+
  theme(legend.position='none')
ggsave(file.path(root,'change_in_pivot_frequency.png'),p,width=9,height=5,dpi=200)

for (label in unique(d$cohort_label)) {
  x <- d[d$cohort_label==label,];x <- x[order(x$group_short,x$delta_rate_per_min),]
  id <- vapply(x$cohort_ids,function(x)paste(x,collapse=','),character(1))
  x$label <- paste0(x$group_short,': ',ifelse(id=='','',paste0('ID ',id,', ')),
                    substr(x$date,6,10),' fly ',x$daily_fly)
  x$label <- factor(x$label,levels=rev(x$label))
  p <- ggplot(x,aes(y=label))+
    geom_segment(aes(x=pre_rate_per_min,xend=inhibition_rate_per_min,yend=label),color='grey60')+
    geom_point(aes(x=pre_rate_per_min),shape=1,size=2.6,color='black')+
    geom_point(aes(x=inhibition_rate_per_min,color=group_short),size=2.4)+
    scale_color_manual(values=colors)+
    labs(title=paste('Individual flies:',label),
         subtitle='Open circles: pre-light; filled circles: inhibition',
         x='Pivot-like bouts / min',y=NULL)+theme(legend.position='none')
  ggsave(file.path(root,paste0(gsub('[^a-z0-9]+','_',tolower(label)),'_individual_flies.png')),
         p,width=10,height=max(5,nrow(x)*.27+1.6),dpi=180)
}

# Sensitivity of the group difference in within-fly changes.
s <- fromJSON(file.path(root,'group_summary.json'))
s <- s[s$subset=='all' & s$cohort %in% c('blinded round 1','unblinded RubyACR') & s$duration_s %in% c(5,10),]
s <- s[!is.na(s$difference_in_changes),]
s$cohort_label <- paste(s$cohort,paste0(s$duration_s,' s'))
p <- ggplot(s,aes(definition,difference_in_changes,group=cohort_label,color=cohort_label))+
  geom_hline(yintercept=0,linetype=2,color='grey50')+geom_line()+geom_point(size=2)+
  labs(title='Sensitivity to the pivot-like bout definition',
       subtitle='Negative values indicate a larger reduction in experimental flies',
       x=NULL,y='Experimental minus control change (bouts / min)',color=NULL)+
  theme(legend.position='bottom')
ggsave(file.path(root,'definition_sensitivity.png'),p,width=10,height=5,dpi=180)

e <- fromJSON(file.path(root,'event_examples.json'))
e <- e[e$cohort %in% c('blinded round 1','unblinded RubyACR') & e$duration_s %in% c(5,10),]
e$group_short <- ifelse(grepl('^experimental',e$group),'Exp','Ctrl')
e$short_cohort <- ifelse(e$cohort=='blinded round 1','Round1','Unblinded')
e$category <- paste(e$short_cohort,e$group_short,e$phase)
chosen <- do.call(rbind,lapply(split(e,e$category),function(z) head(z[order(z$rotation_deg,decreasing=TRUE),],2)))
paths <- do.call(rbind,lapply(seq_len(nrow(chosen)),function(i) {
  z <- chosen[i,];
  data.frame(x=z$x[[1]],y=z$y[[1]],sample=seq_along(z$x[[1]]),
    panel=paste(z$short_cohort,z$group_short,z$phase,substr(z$date,6,10),
                paste0('fly',z$daily_fly),paste0('pulse',z$pulse),paste0(round(z$rotation_deg),' deg')))
}))
p <- ggplot(paths,aes(x,y,color=sample,group=panel))+geom_path(linewidth=.8)+
  facet_wrap(~panel,ncol=4,labeller=label_wrap_gen(width=37))+coord_equal()+scale_color_viridis_c()+
  labs(title='Examples detected as pivot-like bouts',
       subtitle='Largest-rotation examples from each cohort/group/phase; 1 s representative windows',
       x='Relative x (mm)',y='Relative y (mm)')+
  theme(legend.position='none',strip.text=element_text(size=7))
ggsave(file.path(root,'pivot_event_examples.png'),p,width=13,height=10,dpi=160)
