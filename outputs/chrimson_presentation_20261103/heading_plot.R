options(contrasts=c('contr.sum','contr.poly'))
library(car)
library(emmeans)
out <- 'outputs/chrimson_presentation_20261103'
d <- read.csv('/Users/sophiarenauld/Downloads/Wilson Lab/final_ves041_cschrim_dataset.csv')
groups <- list(c(2,7,16,19,21,26,53,58,65,70,72,74,77,79,81,85,119,130),c(1,14,18,24,29,40,44,66,73,75,80,90,100,102,105,108,110),c(3,4,5,8,12,15,25,27,69,99,101,103,111,114,121,124),c(35,52,59,63,89,93,107,113,118,122,123,125,129,132,133,137),c(6,17,20,23,31,34,39,46,50,54,57,62,84,87,88,91,96))
nms <- c('VES041-1 Chrimson','VES041-2 Chrimson','VES041-1 GFP','VES041-2 GFP','Empty Chrimson')
d$Genotype <- NA_character_
for(i in seq_along(groups)) d$Genotype[d$Fly %in% groups[[i]]] <- nms[i]
d$Genotype <- factor(d$Genotype,levels=nms)
stopifnot(nrow(d)==84,!anyDuplicated(d$Fly),!anyNA(d$Genotype))
comparisons <- list('VES041-1 Chrimson minus VES041-1 GFP'=c(1,0,-1,0,0),'VES041-2 Chrimson minus VES041-2 GFP'=c(0,1,0,-1,0),'VES041-1 Chrimson minus Empty Chrimson'=c(1,0,0,0,-1),'VES041-2 Chrimson minus Empty Chrimson'=c(0,1,0,0,-1))

library(ggplot2)
y <- d$Net.Heading.Displacement..deg.
stopifnot(all(is.finite(y)))
m <- aov(y~Genotype,d)
tuk <- as.data.frame(TukeyHSD(m)$Genotype); tuk$comparison<-rownames(tuk)
write.csv(tuk,file.path(out,'heading_tukey_all_comparisons.csv'),row.names=FALSE)
lev<-leveneTest(y~Genotype,d,center=median)
wt<-oneway.test(y~Genotype,d)
ct<-as.data.frame(summary(contrast(emmeans(m,~Genotype),comparisons),infer=c(TRUE,TRUE),adjust='none')); ct$p_holm<-p.adjust(ct$p.value,'holm')
# Match the source figure's Tukey adjustment across all 10 genotype pairs.
get_p<-function(i,j) tuk[paste(nms[j],nms[i],sep='-'),'p adj']
ps<-ct$p_holm
write.csv(ct,file.path(out,'heading_planned_holm_comparisons.csv'),row.names=FALSE)
stars<-function(p) ifelse(p<.001,'***',ifelse(p<.01,'**',ifelse(p<.05,'*','ns')))
ann<-data.frame(x=c(1,2,1,2),xe=c(3,4,5,5),h=c(116,128,140,152),label=ifelse(ps<.001,'p < 0.001',paste0('p = ',formatC(ps,format='f',digits=3))))
sm<-do.call(rbind,lapply(nms,function(g){v<-y[d$Genotype==g];data.frame(Genotype=g,n=length(v),mean=mean(v),SE=sd(v)/sqrt(length(v)))}))
write.csv(sm,file.path(out,'heading_group_summary.csv'),row.names=FALSE)
capture.output(summary(m),cat('\nLevene test\n'),print(lev),cat('\nWelch omnibus\n'),print(wt),cat('\nFour displayed comparisons: Holm adjusted across four planned pairs\n'),print(data.frame(comparison=names(comparisons),p_holm=ps,stars=stars(ps))),file=file.path(out,'heading_statistics.txt'))
d$heading<-y
p<-ggplot(d,aes(Genotype,heading))+
 stat_summary(fun=mean,geom='col',width=.56,fill='#E7B655',alpha=.4)+
 geom_point(position=position_jitter(width=.105,height=0,seed=41),shape=21,fill='#D49B2E',colour='#71521C',stroke=.35,size=2.7,alpha=.85)+
 stat_summary(fun.data=mean_se,geom='errorbar',width=.17,linewidth=.8,colour='#172431')+
 stat_summary(fun=mean,geom='point',shape=95,size=10,colour='#172431')+
 geom_segment(data=ann,aes(x=x,xend=xe,y=h,yend=h),inherit.aes=FALSE,linewidth=.55)+
 geom_segment(data=ann,aes(x=x,xend=x,y=h-2,yend=h),inherit.aes=FALSE,linewidth=.55)+
 geom_segment(data=ann,aes(x=xe,xend=xe,y=h-2,yend=h),inherit.aes=FALSE,linewidth=.55)+
 geom_text(data=ann,aes(x=(x+xe)/2,y=h+2,label=label),inherit.aes=FALSE,size=6)+
 scale_x_discrete(labels=setNames(paste0(c('VES041-1\nChrimson','VES041-2\nChrimson','VES041-1\nGFP','VES041-2\nGFP','Empty\nChrimson'),'\nn = ',sm$n),nms))+
 scale_y_continuous(limits=c(0,160),breaks=seq(0,150,30),expand=expansion(mult=c(0,.02)))+
 labs(title='Net heading displacement',x=NULL,y='Net heading displacement (°)',caption='Points: individual flies; bars: mean ± SE. One-way ANOVA with Holm-adjusted planned comparisons (4 pairs).\nHorizontal jitter separates overlapping observations.')+
 theme_classic(base_size=19,base_family='Arial')+
 theme(axis.line.y=element_blank(),axis.line.x=element_line(colour='#AAB2B9'),axis.ticks=element_blank(),panel.grid.major.y=element_line(colour='#E5E9ED',linewidth=.4),axis.text=element_text(colour='#374151'),axis.text.x=element_text(size=17,margin=margin(t=12)),plot.title=element_text(size=28,face='bold',margin=margin(b=18)),plot.caption=element_text(size=12,hjust=0,colour='#64707D',margin=margin(t=18)),plot.margin=margin(20,25,15,15))
ggsave(file.path(out,'net_heading_displacement.svg'),p,width=14,height=9,bg='white')
ggsave(file.path(out,'net_heading_displacement.png'),p,width=14,height=9,bg='white',dpi=130)
cat(readLines(file.path(out,'heading_statistics.txt')),sep='\n')

w <- do.call(rbind,lapply(seq_along(comparisons),function(i){v<-comparisons[[i]]; tt<-t.test(y[d$Genotype==nms[which(v==1)]],y[d$Genotype==nms[which(v== -1)]]);data.frame(comparison=names(comparisons)[i],p=tt$p.value)}));w$p_holm<-p.adjust(w$p,'holm')
write.csv(w,file.path(out,'heading_welch_sensitivity.csv'),row.names=FALSE)
capture.output(cat('\nWelch planned comparisons; Holm across four\n'),print(w),file=file.path(out,'heading_statistics.txt'),append=TRUE)
print(w)
