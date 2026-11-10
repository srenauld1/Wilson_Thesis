suppressPackageStartupMessages({library(ggplot2);library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
root <- 'outputs/chrimson_all_flies_mean_sem_20261103'
dir.create(root,recursive=TRUE,showWarnings=FALSE)
destination <- if(length(args)) args[1] else root
stopifnot(dir.exists(destination))
src <- '.build_ves041_activation/reproducibility/source_data/fly_summary_84.csv'
mapping_source <- '.build_ves041_activation/reproducibility/source_data/cschrim_behaviorcrunch.R'
d <- read.csv(src,check.names=FALSE,stringsAsFactors=FALSE)
names(d) <- trimws(names(d))
d$Fly <- suppressWarnings(as.numeric(d$Fly))
d <- d[!is.na(d$Fly),]
stopifnot(nrow(d)==84,!anyDuplicated(d$Fly))
mapping_text <- paste(readLines(mapping_source,warn=FALSE),collapse='\n')
matches <- regmatches(mapping_text,gregexpr('Fly %in% c\\([^)]*\\) ~ "[^"]+"',mapping_text))[[1]]
mapping <- list()
for(m in matches) {
  ids <- as.numeric(strsplit(sub('.*c\\(([^)]*)\\).*','\\1',m),',')[[1]])
  group <- sub('.*~ "([^"]+)"','\\1',m)
  for(id in ids) {stopifnot(is.null(mapping[[as.character(id)]]));mapping[[as.character(id)]] <- group}
}
order <- c('ves041_1_Chrim','ves041_2_Chrim','ves041_1_gfp','ves041_2_gfp','empty_Chrim')
labels <- c('VES041_1 x csChrimson','VES041_2 x csChrimson','VES041_1 x GFP','VES041_2 x GFP','Empty split x csChrimson')
d$genotype <- vapply(d$Fly,function(id)mapping[[as.character(id)]],character(1))
stopifnot(all(d$genotype %in% order))
counts <- table(factor(d$genotype,levels=order))
stopifnot(identical(as.integer(counts),c(18L,17L,16L,16L,17L)))
panel_labels <- setNames(paste0(labels,'\n(n = ',as.integer(counts),')'),order)
columns <- list(
  forward=c('Pre Opto Forward Velocity','Opto On Forward Velocity'),
  rotation=c('Pre Opto Rotational Speed','Opto On Rotational Speed'))
metric_labels <- c(forward='Forward velocity (mm/s)',rotation='Rotational speed (deg/s)')
long <- do.call(rbind,lapply(names(columns),function(metric)do.call(rbind,lapply(1:2,function(j){
  value <- as.numeric(d[[columns[[metric]][j]]]);stopifnot(all(is.finite(value)))
  data.frame(fly=d$Fly,genotype=d$genotype,metric=metric,period=c('Pre-opto','Opto-on')[j],value=value)
}))))
long$period <- factor(long$period,levels=c('Pre-opto','Opto-on'))
long$genotype <- factor(long$genotype,levels=order)
long$metric <- factor(long$metric,levels=c('forward','rotation'))
stopifnot(nrow(long)==336,all(table(long$fly)==4))
summary <- do.call(rbind,lapply(split(long,list(long$genotype,long$metric,long$period),drop=TRUE),function(s){
  data.frame(genotype=s$genotype[1],metric=s$metric[1],period=s$period[1],
             n=nrow(s),mean=mean(s$value),sem=sd(s$value)/sqrt(nrow(s)))
}))
colors <- setNames(c('#B85F35','#167D8D','#8F9AA8','#596A7F','#9D8F72'),order)
p <- ggplot(long,aes(period,value,group=fly,color=genotype))+
  geom_line(alpha=.28,linewidth=.48)+geom_point(alpha=.38,size=1.3)+
  geom_errorbar(data=summary,aes(y=mean,ymin=mean-sem,ymax=mean+sem,group=genotype),
                width=.14,linewidth=.8)+
  geom_line(data=summary,aes(y=mean,group=genotype),linewidth=1.35)+
  geom_point(data=summary,aes(y=mean,group=genotype,fill=genotype),shape=21,
             color='white',stroke=.65,size=3.7)+
  facet_grid(metric~genotype,scales='free_y',switch='y',
             labeller=labeller(metric=metric_labels,genotype=panel_labels))+
  scale_color_manual(values=colors)+scale_fill_manual(values=colors)+
  scale_x_discrete(expand=expansion(add=.42))+
  scale_y_continuous(expand=expansion(mult=c(.06,.1)))+
  labs(title='VES041 activation',
       subtitle='Each line is one fly. Bold lines and error bars show mean +/- SEM.',
       x=NULL,y=NULL,caption='All 84 archived fly summaries. Each fly has equal weight.')+
  theme_classic(base_size=15,base_family='Arial')+
  theme(legend.position='none',strip.background=element_blank(),
        strip.placement='outside',strip.text.x=element_text(size=14,face='bold',lineheight=1.15),
        strip.text.y.left=element_text(angle=90,size=16,margin=margin(r=10)),
        panel.spacing.x=grid::unit(.8,'cm'),panel.spacing.y=grid::unit(.7,'cm'),
        plot.title=element_text(size=25,face='bold'),plot.subtitle=element_text(size=15,margin=margin(b=16)),
        axis.text=element_text(color='#333333',size=13),
        panel.grid.major.y=element_line(color='#E9EBED',linewidth=.3),
        plot.caption=element_text(size=11,color='#666666',margin=margin(t=12)),
        plot.margin=margin(18,20,12,14))
name <- 'csChrimson_all84_forward_rotation_mean_SEM'
dest <- file.path(destination,paste0(name,'.svg'))
stopifnot(!file.exists(dest))
ggsave(dest,p,width=17.8,height=8.5,device=svglite::svglite,bg='white')
ggsave(file.path(root,paste0(name,'.png')),p,width=17.8,height=8.5,dpi=180,bg='white')
write_json(long,file.path(root,'plotted_fly_values.json'),pretty=TRUE,dataframe='rows',digits=NA)
write_json(summary,file.path(root,'group_mean_sem.json'),pretty=TRUE,dataframe='rows',digits=NA)
methods <- c('CSCHRIMSON ACTIVATION: ALL ARCHIVED FLY SUMMARIES',
  paste('Source:',normalizePath(src)),paste('Genotype mapping:',normalizePath(mapping_source)),
  '84 unique flies: VES041_1 x csChrimson 18; VES041_2 x csChrimson 17; VES041_1 x GFP 16; VES041_2 x GFP 16; empty split x csChrimson 17.',
  'All 84 archived rows retained. No new activity filter, outlier exclusion, imputation, or hypothesis tests.',
  'Each line connects the existing pre-opto and opto-on average for one fly. No vertical transformations or jitter.',
  'Group means are arithmetic means of the fly averages. SEM is sample SD across flies divided by sqrt(n), separately for each phase and genotype.',
  'These error bars are SEMs of phase means, not SEMs of paired changes. Fly, not presentation, is the independent unit.',
  'Forward measure retains signed forward velocity as archived; rotational speed is angular-velocity magnitude.',
  'Important source limitation: these are archived fly-level values, not a new all-presentation extraction. The upstream pipeline has a moving-baseline-only summary path and the 84-fly export does not retain per-pulse inclusion metadata. Thus absence of a NEW walking filter does not establish that every stationary presentation contributed to each archived mean.',
  'Do not mix this 84-fly table with the distinct 81-fly pulse-derived dataset. This plot includes IDs 1, 6 and 39 from the archive.',
  paste('SVG:',dest))
writeLines(methods,file.path(root,'methods.txt'))
print(counts);print(summary);cat('Saved',dest,'\n')
