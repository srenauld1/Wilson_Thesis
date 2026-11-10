Sys.setlocale('LC_CTYPE', 'en_US.UTF-8')
library(ggplot2)
library(patchwork)
out <- 'outputs/chrimson_presentation_20261103'
d <- read.csv('/Users/sophiarenauld/Downloads/Wilson Lab/final_ves041_cschrim_dataset.csv')
groups <- list(c(2,7,16,19,21,26,53,58,65,70,72,74,77,79,81,85,119,130),c(1,14,18,24,29,40,44,66,73,75,80,90,100,102,105,108,110),c(3,4,5,8,12,15,25,27,69,99,101,103,111,114,121,124),c(35,52,59,63,89,93,107,113,118,122,123,125,129,132,133,137),c(6,17,20,23,31,34,39,46,50,54,57,62,84,87,88,91,96))
labels <- c('VES041-A\nCsChrimson','VES041-B\nCsChrimson','VES041-A\nGFP','VES041-B\nGFP','Empty\nCsChrimson')
d$Genotype <- NA_character_
for(i in seq_along(groups)) d$Genotype[d$Fly %in% groups[[i]]] <- labels[i]
stopifnot(!anyNA(d$Genotype), nrow(d) == 84, !anyDuplicated(d$Fly))
for (i in seq_along(groups)) stopifnot(all(groups[[i]] %in% d$Fly))
d$Genotype <- factor(d$Genotype,levels=labels)
make_plot <- function(cols,title,ylab,lim,brks){
  a <- data.frame(Fly=d$Fly,Genotype=d$Genotype,Condition='Pre-stim',Value=d[[cols[1]]])
  b <- data.frame(Fly=d$Fly,Genotype=d$Genotype,Condition='Activation',Value=d[[cols[2]]])
  dat <- rbind(a,b); dat$Condition <- factor(dat$Condition,levels=c('Pre-stim','Activation'))
  ann <- data.frame(Genotype=factor(labels[1:2],levels=labels),x=1.5,y=lim[2]*.97,label=if(title=='Rotational speed') c('Change vs GFP: ***','Change vs GFP: ***') else c('Change vs GFP: *','Change vs GFP: ns'))
  ann2 <- ann; ann2$y <- lim[2]*.90
  ann2$label <- if(title=='Rotational speed') c('Change vs Empty: ***','Change vs Empty: *†') else c('Change vs Empty: ns','Change vs Empty: ns')
  ggplot(dat,aes(Condition,Value))+
    stat_summary(aes(fill=Condition),fun=mean,geom='col',width=.56,alpha=.23)+
    geom_line(aes(group=Fly),colour='#7C8793',alpha=.48,linewidth=.45)+
    geom_point(aes(colour=Condition),size=2.1,alpha=.8)+
    stat_summary(fun.data=mean_se,geom='errorbar',width=.19,linewidth=.8,colour='#172431')+
    stat_summary(fun=mean,geom='point',shape=95,size=10,colour='#172431')+
    geom_text(data=rbind(ann,ann2),aes(x=x,y=y,label=label),inherit.aes=FALSE,size=4.3,colour='#172431')+
    facet_wrap(~Genotype,nrow=1)+
    scale_fill_manual(values=c('#4668B2','#DE604A'))+
    scale_colour_manual(values=c('#4668B2','#DE604A'))+
    scale_y_continuous(limits=lim,breaks=brks,expand=expansion(mult=c(0,.04)))+
    labs(title=title,x=NULL,y=ylab,caption=paste0('Paired observations; bars: mean ± SE. Labels compare activation − pre-stim changes; Holm correction across 4 comparisons per outcome.\n* p < 0.05; *** p < 0.001; ns: p ≥ 0.05.',if(title=='Rotational speed') '  † Not significant in the Welch sensitivity test.' else '  Overall genotype × condition interaction: p = 0.0573.'))+
    theme_classic(base_size=19,base_family='Arial')+
    theme(legend.position='none',axis.line.x=element_line(colour='#AAB2B9',linewidth=.4),axis.line.y=element_blank(),axis.ticks=element_blank(),axis.text=element_text(colour='#374151'),axis.text.x=element_text(size=15,margin=margin(t=10)),axis.title.y=element_text(size=19,margin=margin(r=14)),panel.grid.major.y=element_line(colour='#E5E9ED',linewidth=.4),strip.background=element_blank(),strip.text=element_text(size=18,face='bold',margin=margin(b=14)),panel.spacing=grid::unit(1,'lines'),plot.title=element_text(size=28,face='bold',margin=margin(b=22)),plot.caption=element_text(size=12,colour='#64707D',hjust=0,margin=margin(t=18)),plot.margin=margin(20,20,15,15))
}
# A common panel permits brackets spanning the compared genotypes.
rd <- rbind(data.frame(Fly=d$Fly,Genotype=d$Genotype,Condition='Pre-stim',Value=d$Pre.Opto.Rotational.Speed),data.frame(Fly=d$Fly,Genotype=d$Genotype,Condition='Activation',Value=d$Opto.On.Rotational.Speed))
rd$Condition <- factor(rd$Condition,levels=c('Pre-stim','Activation'))
rd$center <- (as.integer(rd$Genotype)-1)*3+1
rd$x <- rd$center+ifelse(rd$Condition=='Pre-stim',-.43,.43)
brackets <- data.frame(x=c(1,4,1,4),xe=c(7,10,13,13),y=c(177,193,209,225),label=c('p < 0.001','p < 0.001','p < 0.001','p = 0.031†'))
r <- ggplot(rd,aes(x,Value))+
 stat_summary(aes(fill=Condition,group=interaction(Genotype,Condition)),fun=mean,geom='col',width=.62,alpha=.23)+
 geom_line(aes(group=Fly),colour='#7C8793',alpha=.48,linewidth=.45)+
 geom_point(aes(colour=Condition),size=2.1,alpha=.8)+
 stat_summary(aes(group=interaction(Genotype,Condition)),fun.data=mean_se,geom='errorbar',width=.22,linewidth=.8,colour='#172431')+
 stat_summary(aes(group=interaction(Genotype,Condition)),fun=mean,geom='point',shape=95,size=10,colour='#172431')+
 geom_segment(data=brackets,aes(x=x,xend=xe,y=y,yend=y),inherit.aes=FALSE,linewidth=.6)+
 geom_segment(data=brackets,aes(x=x,xend=x,y=y-2.5,yend=y),inherit.aes=FALSE,linewidth=.6)+
 geom_segment(data=brackets,aes(x=xe,xend=xe,y=y-2.5,yend=y),inherit.aes=FALSE,linewidth=.6)+
 geom_text(data=brackets,aes(x=(x+xe)/2,y=y+3,label=label),inherit.aes=FALSE,size=6)+
 scale_fill_manual(values=c('#4668B2','#DE604A'))+
 scale_colour_manual(values=c('#4668B2','#DE604A'))+
 scale_x_continuous(breaks=c(1,4,7,10,13),labels=labels,limits=c(-.3,14.3))+
 scale_y_continuous(limits=c(0,237),breaks=seq(0,150,30),expand=expansion(mult=c(0,0)))+
 labs(title='Rotational speed',x=NULL,y='Rotational speed (°/s)',fill=NULL,colour=NULL,caption='Brackets compare activation − pre-stim changes between genotypes; Holm-adjusted planned comparisons from the mixed-design ANOVA.\nBars: mean ± SE; points: individual flies. † Not significant in the Welch sensitivity test.')+
 theme_classic(base_size=19,base_family='Arial')+
 theme(legend.position='top',legend.justification='left',axis.line.y=element_blank(),axis.line.x=element_line(colour='#AAB2B9',linewidth=.4),axis.ticks=element_blank(),axis.text=element_text(colour='#374151'),axis.text.x=element_text(size=16,hjust=0.5,vjust=1,lineheight=1.1,margin=margin(t=12,b=10)),panel.grid.major.y=element_line(colour='#E5E9ED',linewidth=.4),plot.title=element_text(size=28,face='bold',margin=margin(b=10)),plot.caption=element_text(size=12,colour='#64707D',hjust=0,margin=margin(t=18)),plot.margin=margin(20,35,25,25))
fd <- rbind(data.frame(Fly=d$Fly,Genotype=d$Genotype,Condition='Pre-stim',Value=d$Pre.Opto.Forward.Velocity),data.frame(Fly=d$Fly,Genotype=d$Genotype,Condition='Activation',Value=d$Opto.On.Forward.Velocity))
fd$Condition <- factor(fd$Condition,levels=c('Pre-stim','Activation'))
fd$center <- (as.integer(fd$Genotype)-1)*3+1
fd$x <- fd$center+ifelse(fd$Condition=='Pre-stim',-.43,.43)
brackets <- data.frame(x=c(1,4,1,4),xe=c(7,10,13,13),y=c(10.9,11.9,12.9,13.9),label=c('p = 0.026','ns','ns','ns'))
f <- ggplot(fd,aes(x,Value))+
 stat_summary(aes(fill=Condition,group=interaction(Genotype,Condition)),fun=mean,geom='col',width=.62,alpha=.23)+
 geom_line(aes(group=Fly),colour='#7C8793',alpha=.48,linewidth=.45)+
 geom_point(aes(colour=Condition),size=2.1,alpha=.8)+
 stat_summary(aes(group=interaction(Genotype,Condition)),fun.data=mean_se,geom='errorbar',width=.22,linewidth=.8,colour='#172431')+
 stat_summary(aes(group=interaction(Genotype,Condition)),fun=mean,geom='point',shape=95,size=10,colour='#172431')+
 geom_segment(data=brackets,aes(x=x,xend=xe,y=y,yend=y),inherit.aes=FALSE,linewidth=.6)+
 geom_segment(data=brackets,aes(x=x,xend=x,y=y-.15,yend=y),inherit.aes=FALSE,linewidth=.6)+
 geom_segment(data=brackets,aes(x=xe,xend=xe,y=y-.15,yend=y),inherit.aes=FALSE,linewidth=.6)+
 geom_text(data=brackets,aes(x=(x+xe)/2,y=y+.20,label=label),inherit.aes=FALSE,size=6,vjust=0)+
 scale_fill_manual(values=c('#4668B2','#DE604A'))+
 scale_colour_manual(values=c('#4668B2','#DE604A'))+
 scale_x_continuous(breaks=c(1,4,7,10,13),labels=labels,limits=c(-.3,14.3))+
 scale_y_continuous(limits=c(0,14.8),breaks=seq(0,10,2),expand=expansion(mult=c(0,0)))+
 labs(title='Forward velocity',x=NULL,y='Forward velocity (mm/s)',fill=NULL,colour=NULL,caption='Brackets compare activation − pre-stim changes between genotypes; Holm-adjusted planned comparisons from the mixed-design ANOVA.\nBars: mean ± SE; points: individual flies. ns: p ≥ 0.05. Overall genotype × condition interaction: p = 0.0573.')+
 theme_classic(base_size=19,base_family='Arial')+
 theme(legend.position='top',legend.justification='left',axis.line.y=element_blank(),axis.line.x=element_line(colour='#AAB2B9',linewidth=.4),axis.ticks=element_blank(),axis.text=element_text(colour='#374151'),axis.text.x=element_text(size=16,hjust=0.5,vjust=1,lineheight=1.1,margin=margin(t=12,b=10)),panel.grid.major.y=element_line(colour='#E5E9ED',linewidth=.4),plot.title=element_text(size=28,face='bold',margin=margin(b=10)),plot.caption=element_text(size=12,colour='#64707D',hjust=0,margin=margin(t=18)),plot.margin=margin(20,35,25,25))
for(n in c('rotational_speed','forward_velocity')){
 p <- if(n=='rotational_speed') r else f
 ggsave(file.path(out,paste0(n,'.svg')),p,width=14,height=7.875,bg='white')
 ggsave(file.path(out,paste0(n,'.png')),p,width=14,height=7.875,bg='white',dpi=150)
}
ggsave(file.path(out,'both_plots.svg'),r/f,width=16,height=14,bg='white')
