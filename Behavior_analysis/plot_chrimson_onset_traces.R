suppressPackageStartupMessages({library(ggplot2);library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
root <- 'outputs/chrimson_onset_traces_20261103'
out <- if(length(args))args[1] else root
rotation_only <- length(args)>=2 && args[2]=='rotation'
stopifnot(dir.exists(out))
d <- fromJSON(file.path(root,'fly_traces.json'))
s <- fromJSON(file.path(root,'genotype_mean_sem.json'))
r <- fromJSON(file.path(root,'extraction_results.json'))
stopifnot(length(unique(d$fly))==84,nrow(r)==84,all(r$n_pulses_included==80),
          all(is.finite(d$value)),all(is.finite(s$mean)),all(is.finite(s$sem)))
order <- c('ves041_1_Chrim','ves041_2_Chrim','ves041_1_gfp','ves041_2_gfp','empty_Chrim')
labels <- c('VES041_1 x csChrimson','VES041_2 x csChrimson','VES041_1 x GFP','VES041_2 x GFP','Empty split x csChrimson')
counts <- table(factor(r$genotype,levels=order))
stopifnot(identical(as.integer(counts),c(18L,17L,16L,16L,17L)))
panel_labels <- setNames(paste0(labels,'\n(n = ',as.integer(counts),')'),order)
metrics <- c(forward='Forward velocity (mm/s)',rotation='Rotational speed (deg/s)')
if(rotation_only) {d <- d[d$metric=='rotation',];s <- s[s$metric=='rotation',]}
for(name in c('d','s')) {
  z <- get(name);z$genotype <- factor(z$genotype,levels=order)
  z$metric <- factor(z$metric,levels=c('forward','rotation'));assign(name,z)
}
colors <- setNames(c('#B85F35','#167D8D','#8F9AA8','#596A7F','#9D8F72'),order)
if(rotation_only) colors <- setNames(c('#BA354F','#BA354F','#222222','#555555','#888888'),order)
p <- ggplot(d,aes(time,value,group=fly,color=genotype))+
  annotate('rect',xmin=0,xmax=1,ymin=-Inf,ymax=Inf,fill='#B6B6B6',alpha=.16)+
  geom_vline(xintercept=c(0,1),color='#999999',linewidth=.35,linetype='dotted')+
  geom_ribbon(data=s,aes(x=time,y=NULL,ymin=mean-sem,ymax=mean+sem,
                         group=genotype,fill=genotype),color=NA,alpha=.23,inherit.aes=FALSE)+
  geom_line(alpha=.28,linewidth=.4)+
  geom_line(data=s,aes(y=mean,group=genotype),linewidth=1.15)+
  facet_grid(metric~genotype,scales='free_y',switch='y',
             labeller=labeller(metric=metrics,genotype=panel_labels))+
  scale_color_manual(values=colors)+scale_fill_manual(values=colors)+
  scale_x_continuous(limits=c(-1,2),breaks=c(-1,0,1,2),expand=expansion(mult=0))+
  scale_y_continuous(expand=expansion(mult=c(.07,.08)))+
  expand_limits(y=0)+
  labs(title='VES041 activation',
       subtitle='Thin lines: fly averages. Bold lines and bands: mean +/- SEM across flies. Shading: opto-on.',
       x='Time from opto onset (s)',y=NULL,
       caption='84 flies; 80 presentations per fly. All presentations included. Raw velocities averaged in 10-ms bins.')+
  theme_classic(base_size=15,base_family='Arial')+
  theme(legend.position='none',strip.background=element_blank(),strip.placement='outside',
        strip.text.x=element_text(size=14,face='bold',lineheight=1.15),
        strip.text.y.left=element_text(angle=90,size=16,margin=margin(r=10)),
        panel.spacing.x=grid::unit(.7,'cm'),panel.spacing.y=grid::unit(.75,'cm'),
        plot.title=element_text(size=25,face='bold'),
        plot.subtitle=element_text(size=14,margin=margin(b=16)),
        axis.text=element_text(color='#333333',size=13),
        axis.title.x=element_text(size=16,margin=margin(t=12)),
        panel.grid.major.y=element_line(color='#E9EBED',linewidth=.3),
        plot.caption=element_text(size=11,color='#666666',margin=margin(t=12)),
        plot.margin=margin(18,20,12,14))
name <- 'csChrimson_all84_onset_aligned_forward_rotation_mean_SEM'
if(rotation_only) name <- 'csChrimson_all84_onset_aligned_rotation_mean_SEM_red_gray'
dest <- file.path(out,paste0(name,'.svg'));stopifnot(!file.exists(dest))
height <- if(rotation_only)5.3 else 8.5
ggsave(dest,p,width=17.8,height=height,device=svglite::svglite,bg='white')
ggsave(file.path(root,paste0(name,'.png')),p,width=17.8,height=height,dpi=180,bg='white')
cat('Saved',dest,'\n');print(counts)
