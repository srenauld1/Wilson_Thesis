suppressPackageStartupMessages({library(jsonlite);library(ggplot2)})
root <- 'outputs/ves041_blinded_good_walkers_speed_linked_turns_20261102'
f <- fromJSON(file.path(root,'fly_summary.json'))
d <- f[f$definition=='primary' & f$subset=='all',]
d$group_short <- ifelse(d$group=='experimental','Experimental','Control')
pre <- d;pre$phase <- 'Pre-light';pre$rate <- pre$pre_rate_per_min
on <- d;on$phase <- 'Inhibition';on$rate <- on$inhibition_rate_per_min
long <- rbind(pre,on);long$phase <- factor(long$phase,levels=c('Pre-light','Inhibition'))
colors <- c(Experimental='#BA3C51',Control='#555555')
theme_set(theme_bw(base_size=11))
p <- ggplot(long,aes(phase,rate,group=fly_key,color=group_short))+
  geom_line(alpha=.5)+geom_point(size=2.3)+facet_wrap(~group_short)+
  scale_color_manual(values=colors)+
  labs(title='Speed-linked turns during VES041 inhibition: blinded good walkers',
       subtitle='Each line is one fly; matched 5 s periods; 8 experimental and 9 control flies',
       x=NULL,y='Qualified turns / min',
       caption='Turn >=30 degrees, preceded by >=2 mm/s forward slowdown OR followed by >=2 mm/s forward speed-up.\nSpeed changes compare adjacent 0.5 s windows before/after the turn. No displacement restriction.')+
  theme(legend.position='none',plot.caption=element_text(hjust=0))
ggsave(file.path(root,'paired_speed_linked_turn_frequency.png'),p,width=10,height=5.5,dpi=200)

p <- ggplot(d,aes(group_short,delta_rate_per_min,color=group_short))+
  geom_hline(yintercept=0,linetype=2,color='grey50')+
  geom_boxplot(outlier.shape=NA,width=.4,alpha=.1)+
  geom_point(position=position_jitter(width=.08,height=0,seed=41),size=2.5)+
  scale_color_manual(values=colors)+labs(title='Within-fly change in speed-linked turning',
       x=NULL,y='Inhibition minus pre-light (turns / min)')+theme(legend.position='none')
ggsave(file.path(root,'change_in_speed_linked_turn_frequency.png'),p,width=6,height=5,dpi=200)

s <- fromJSON(file.path(root,'group_summary.json'));s <- s[s$subset=='all',]
labels <- c(primary='Primary: 0.5 s',window025='0.25 s windows',window1='1 s windows',
            yaw10='10 deg/s boundary',yaw30='30 deg/s boundary',turn_vs_flank='During turn vs flanks',
            same_phase_flanks='Flanks within phase',max1s_turn='Turns <=1 s')
s$label <- factor(labels[s$definition],levels=rev(labels))
p <- ggplot(s,aes(difference_in_changes,label))+
  geom_vline(xintercept=0,linetype=2,color='grey50')+geom_point(size=3,color='#215B73')+
  labs(title='Sensitivity of the experimental-control contrast',
       subtitle='Negative values indicate a larger reduction in experimental flies',
       x='Experimental minus control change (turns / min)',y=NULL)
ggsave(file.path(root,'definition_sensitivity.png'),p,width=9,height=5,dpi=180)

e <- fromJSON(file.path(root,'events.json'));e <- e[e$qualifies,]
e$category <- paste(e$group,e$phase)
chosen <- do.call(rbind,lapply(split(e,e$category),function(z) {
  z <- z[order(abs(z$rotation_deg-median(z$rotation_deg))),]
  z[seq_len(min(2,nrow(z))),]
}))
traces <- list();segments <- list();boundaries <- list()
for (i in seq_len(nrow(chosen))) {
  z <- chosen[i,];id <- paste(z$cohort_ids[[1]],collapse=',')
  panel <- paste0(z$group,' ID ',id,', ',z$phase,', pulse ',z$pulse,'\n',
                  round(z$rotation_deg),' deg; slowdown ',round(z$slowdown_mm_s,1),
                  ', speed-up ',round(z$speedup_mm_s,1),' mm/s')
  traces[[i]] <- data.frame(time=z$trace_time_s[[1]],forward=z$trace_forward[[1]],panel=panel)
  segments[[i]] <- data.frame(x=c(-1,-.5,z$duration_s,z$duration_s+.5),
      xend=c(-.5,0,z$duration_s+.5,z$duration_s+1),
      value=c(z$early_pre_mean,z$late_pre_mean,z$early_post_mean,z$late_post_mean),panel=panel)
  boundaries[[i]] <- data.frame(time=c(0,z$duration_s),panel=panel)
}
traces <- do.call(rbind,traces);segments <- do.call(rbind,segments);boundaries <- do.call(rbind,boundaries)
p <- ggplot(traces,aes(time,forward))+
  geom_line(color='#215B73',linewidth=.5)+
  geom_segment(data=segments,aes(x=x,xend=xend,y=value,yend=value),inherit.aes=FALSE,
               color='#D05739',linewidth=1.2)+
  geom_vline(data=boundaries,aes(xintercept=time),color='grey55',linetype=2)+
  facet_wrap(~panel,ncol=2,scales='free_x')+
  labs(title='Detected turns and the forward-velocity criterion',
       subtitle='Blue: forward velocity; orange: adjacent 0.5 s means; dashed: turn start/end',
       x='Time from turn start (s)',y='Forward velocity (mm/s)')+
  theme(strip.text=element_text(size=8))
ggsave(file.path(root,'turn_examples_with_speed_windows.png'),p,width=11,height=11,dpi=170)
