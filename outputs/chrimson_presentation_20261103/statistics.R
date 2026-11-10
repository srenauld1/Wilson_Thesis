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
all_anova <- list(); all_contrasts <- list(); all_means <- list(); all_diagnostics <- list()
for(metric in c('Rotational.Speed','Forward.Velocity')){
 pre <- d[[paste0('Pre.Opto.',metric)]]; on <- d[[paste0('Opto.On.',metric)]]
 stopifnot(all(is.finite(pre)),all(is.finite(on)))
 z <- data.frame(Genotype=d$Genotype,delta=on-pre,average=(pre+on)/2)
 md <- lm(delta~Genotype,z); ma <- lm(average~Genotype,z)
 # With two repeated levels, the within-subject interaction is exactly
 # the genotype test on change scores. Type III condition tests the
 # equally weighted mean change across genotypes (sum contrasts).
 ad <- Anova(md,type=3); aa <- anova(ma)
 tab <- data.frame(metric=metric,effect=c('Genotype','Condition','Genotype x Condition'),df1=c(aa$Df[1],ad$Df[1],ad$Df[2]),df2=df.residual(md),F=c(aa$`F value`[1],ad$`F value`[1],ad$`F value`[2]),p=c(aa$`Pr(>F)`[1],ad$`Pr(>F)`[1],ad$`Pr(>F)`[2]))
 tab$partial_eta_squared <- with(tab,F*df1/(F*df1+df2)); all_anova[[metric]]<-tab
 ct <- as.data.frame(summary(contrast(emmeans(md,~Genotype),comparisons),infer=c(TRUE,TRUE),adjust='none'))
 ct$p_holm <- p.adjust(ct$p.value,'holm'); ct$metric<-metric
 all_contrasts[[metric]]<-ct
 all_means[[metric]] <- do.call(rbind,lapply(nms,function(g){ix<-d$Genotype==g; data.frame(metric=metric,Genotype=g,n=sum(ix),pre=mean(pre[ix]),opto=mean(on[ix]),change=mean(on[ix]-pre[ix]),change_SE=sd(on[ix]-pre[ix])/sqrt(sum(ix)))}))
 all_diagnostics[[metric]]<-data.frame(metric=metric,levene_change_p=leveneTest(delta~Genotype,z,center=median)$`Pr(>F)`[1],levene_average_p=leveneTest(average~Genotype,z,center=median)$`Pr(>F)`[1],shapiro_change_residual_p=shapiro.test(residuals(md))$p.value)
}
a<-do.call(rbind,all_anova); cts<-do.call(rbind,all_contrasts); means<-do.call(rbind,all_means); diagnostics<-do.call(rbind,all_diagnostics)
write.csv(a,file.path(out,'mixed_anova_results.csv'),row.names=FALSE)
write.csv(cts,file.path(out,'planned_change_comparisons.csv'),row.names=FALSE)
write.csv(means,file.path(out,'group_means_and_changes.csv'),row.names=FALSE)
write.csv(diagnostics,file.path(out,'assumption_diagnostics.csv'),row.names=FALSE)
sink(file.path(out,'statistics_results.txt'))
cat('Two-way mixed-design ANOVA; 84 independent flies; paired pre/opto.\nGenotype: 5 levels. Condition: 2 levels. Type III condition test, equal genotype weighting.\nTwo repeated levels: change-score genotype F test equals mixed ANOVA interaction.\nPlanned comparisons use common change-score residual variance, df=79.\nHolm adjustment across FOUR planned comparisons WITHIN each outcome.\n95% confidence intervals are pointwise, unadjusted; p-values are Holm-adjusted.\n\n')
print(a,digits=8,row.names=FALSE); print(cts,digits=8,row.names=FALSE); print(means,digits=7,row.names=FALSE); print(diagnostics,digits=7,row.names=FALSE)
sink()
cat(readLines(file.path(out,'statistics_results.txt')),sep='\n')
robust <- list(); omnibus <- list()
for(metric in c('Rotational.Speed','Forward.Velocity')){
 delta<-d[[paste0('Opto.On.',metric)]]-d[[paste0('Pre.Opto.',metric)]]
 wt<-oneway.test(delta~d$Genotype,var.equal=FALSE)
 omnibus[[metric]]<-data.frame(metric=metric,F=unname(wt$statistic),df1=unname(wt$parameter[1]),df2=unname(wt$parameter[2]),p=wt$p.value)
 tab<-do.call(rbind,lapply(seq_along(comparisons),function(i){v<-comparisons[[i]]; x<-delta[d$Genotype==nms[which(v==1)]]; y<-delta[d$Genotype==nms[which(v== -1)]]; tt<-t.test(x,y,var.equal=FALSE); data.frame(metric=metric,contrast=names(comparisons)[i],estimate=mean(x)-mean(y),lower.CL=tt$conf.int[1],upper.CL=tt$conf.int[2],df=unname(tt$parameter),p=tt$p.value)}))
 tab$p_holm<-p.adjust(tab$p,'holm'); robust[[metric]]<-tab
}
w<-do.call(rbind,robust); wo<-do.call(rbind,omnibus)
write.csv(w,file.path(out,'welch_sensitivity_comparisons.csv'),row.names=FALSE)
write.csv(wo,file.path(out,'welch_sensitivity_omnibus.csv'),row.names=FALSE)
capture.output(cat('\nWelch sensitivity tests on changes (unequal variances); Holm within outcome.\n'),print(wo,digits=8,row.names=FALSE),print(w,digits=8,row.names=FALSE),file=file.path(out,'statistics_results.txt'),append=TRUE)
print(wo,digits=8,row.names=FALSE); print(w,digits=8,row.names=FALSE)
