# Scientific figures for VES041 turn direction, using observed signed values.
# Run analyze_turn_direction.py first, then Rscript plot_turn_direction.R.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(jsonlite)})
args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grepl('^--file=',args)][1])
root <- dirname(normalizePath(script));setwd(root)
f <- read.csv('fly_turn_direction.csv')
s <- read.csv('group_turn_direction.csv')
m <- fromJSON('analysis_results.json')
groups <- c('ves041_1_Chrim','ves041_2_Chrim')
colors <- setNames(c('#B85F35','#167D8D'),groups)
f <- f[f$genotype %in% groups,]
s <- s[s$genotype %in% groups,]
stopifnot(nrow(f)==34,all(f$n_pulses==80))
pooled <- m$summary[m$summary$genotype=='VES041 activation pooled',]
theme_set(theme_classic(base_size=15,base_family='Arial')+
 theme(text=element_text(color='#182B3A'),axis.text=element_text(color='#34485A'),
  axis.line=element_line(color='#A7B1BA',linewidth=.35),axis.ticks=element_line(color='#A7B1BA',linewidth=.35),
  axis.title=element_text(size=15),axis.title.x=element_text(margin=margin(t=13)),
  panel.grid.major.x=element_line(color='#E8ECEF',linewidth=.3),
  strip.background=element_blank(),strip.text=element_text(size=18,face='bold',margin=margin(b=16)),
  legend.position='none',panel.spacing=grid::unit(2.0,'lines'),
  plot.margin=margin(15,20,15,15)))
lab <- c('ves041_1_Chrim'='VES041 line 1 + csChrimson  (n = 18)',
         'ves041_2_Chrim'='VES041 line 2 + csChrimson  (n = 16)')
f$genotype <- factor(f$genotype,levels=groups)
s$genotype <- factor(s$genotype,levels=groups)
f <- f[order(f$genotype,f$Fly),]
f$fly_label <- sprintf('Fly %03d',f$Fly)
lev <- c('Group mean',rev(f$fly_label))
f$fly_label <- factor(f$fly_label,levels=unique(lev))
s$fly_label <- factor('Group mean',levels=unique(lev))
group_lines <- data.frame(genotype=factor(groups,levels=groups),
                          y=as.numeric(table(f$genotype))+.5)

build_panel <- function(value,mean_col,low_col,high_col,xlabel,xlim,breaks,ref){
 d <- f;d$value <- d[[value]]
 e <- s;e$mean <- e[[mean_col]];e$lo <- e[[low_col]];e$hi <- e[[high_col]]
 stopifnot(all(is.finite(d$value)),all(d$value>=xlim[1]),all(d$value<=xlim[2]))
 ggplot(d,aes(value,fly_label,color=genotype))+
  geom_vline(xintercept=ref,color='#6F7F8C',linewidth=.55,linetype='dashed')+
  geom_segment(aes(x=ref,xend=value,yend=fly_label),alpha=.35,linewidth=.75)+
  geom_point(size=3.0,alpha=.94)+
  geom_hline(data=group_lines,aes(yintercept=y),inherit.aes=FALSE,color='#D8E0E5',linewidth=.55)+
  geom_errorbar(data=e,aes(x=mean,y=fly_label,xmin=lo,xmax=hi),orientation='y',width=.22,
                inherit.aes=FALSE,
                color='#182B3A',linewidth=1.0)+
  geom_point(data=e,aes(x=mean,y=fly_label),inherit.aes=FALSE,shape=23,
             color='#182B3A',fill='white',size=4.1,stroke=1.05)+
  facet_wrap(~genotype,nrow=1,scales='free_y',labeller=as_labeller(lab))+
  scale_color_manual(values=colors)+
  scale_x_continuous(limits=xlim,breaks=breaks,expand=expansion(mult=c(.04,.04)))+
  scale_y_discrete(limits=function(x)c(rev(sort(x[x!='Group mean'])),'Group mean'),
                   expand=expansion(add=c(.65,.65)))+
  labs(x=xlabel,y=NULL)+theme(axis.text.y=element_text(size=12),axis.text.x=element_text(size=13))
}

signed <- build_panel('mean_signed_light_deg','light_mean_deg','light_ci95_lo','light_ci95_hi',
 'Mean signed heading during light (deg)',c(-55,55),c(-50,-25,0,25,50),0)
signed <- signed+labs(
 title='Turn direction varies between individual VES041 flies',
 subtitle='First 1 s of csChrimson stimulation. Each colored dot is one fly\'s mean across 80 pulses.',
 caption=sprintf('Negative = left; positive = right, following the lab analysis-code convention.\nBlack diamonds: genotype means with 95%% bootstrap intervals across flies.\nPooled 34 flies: mean %.2f deg, 95%% CI %.2f to %.2f deg; two-sided sign-flip p = %.3f.',
  pooled$light_mean_deg,pooled$light_ci95_lo,pooled$light_ci95_hi,m$primary_test$p))+
 theme(plot.title=element_text(size=24,face='bold',margin=margin(b=12)),
             plot.subtitle=element_text(size=15,color='#586776',margin=margin(b=22)),
             plot.caption=element_text(size=12,color='#586776',hjust=0,lineheight=1.3,margin=margin(t=20)),
             plot.margin=margin(20,22,18,18))

f$right_percent <- 100*f$right_fraction
fraction <- build_panel('right_percent','right_mean_percent','right_ci95_lo','right_ci95_hi',
 'Pulses with a rightward net endpoint (%)',c(0,100),seq(0,100,25),50)
fraction <- fraction+labs(
 title='Direction of pulse endpoints for each VES041 fly',
 subtitle='Each colored dot is the percentage of that fly\'s 80 pulses ending with a positive heading shift.',
 caption=sprintf('Vertical reference: 50%% rightward endpoints. No turn-amplitude or walking-speed threshold is applied.\nBlack diamonds: genotype means with 95%% bootstrap intervals across flies.\nPooled 34 flies: %.1f%% rightward, 95%% CI %.1f%% to %.1f%%; %d rightward and %d leftward pulses.',
  pooled$right_mean_percent,pooled$right_ci95_lo,pooled$right_ci95_hi,pooled$right_pulses,pooled$left_pulses))+
 theme(plot.title=element_text(size=23,face='bold',margin=margin(b=12)),
             plot.subtitle=element_text(size=14.5,color='#586776',margin=margin(b=22)),
             plot.caption=element_text(size=12,color='#586776',hjust=0,lineheight=1.3,margin=margin(t=20)),
             plot.margin=margin(20,22,18,18))

save_figure <- function(p,name){
 ggsave(paste0(name,'.png'),p,width=13.2,height=8.9,dpi=220,bg='white',device=ragg::agg_png)
 ggsave(paste0(name,'.svg'),p,width=13.2,height=8.9,bg='white',device=svglite::svglite)
 ggsave(paste0(name,'.pdf'),p,width=13.2,height=8.9,bg='white',device=grDevices::cairo_pdf)
}
save_figure(signed,'VES041_turn_direction_per_fly')
save_figure(fraction,'VES041_rightward_pulse_fraction_per_fly')
writeLines(capture.output(sessionInfo()),'R_session_info.txt')
cat('Saved signed-heading and directional-frequency figures as PNG, SVG and PDF.\n')
