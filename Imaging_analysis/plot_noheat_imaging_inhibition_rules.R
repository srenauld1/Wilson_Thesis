suppressPackageStartupMessages({library(ggplot2);library(jsonlite)})
root <- 'outputs/ves041_imaging_inhibition_rules_noheat_20261104'
d <- fromJSON(file.path(root,'fly_traces.json'));s <- fromJSON(file.path(root,'mean_sem.json'))
definitions <- c('imaging','localized')
labels <- c('Imaging slowdown + rotation rule\n259 events; 6 mapped flies','Localized-turn rule from inhibition\n238 events; 6 mapped flies')
metrics <- c('dff','forward','rotation');metric_labels <- c(dff='Original dF/F',forward='Forward velocity (mm/s)',rotation='Rotational speed (deg/s)')
for(name in c('d','s')) {z <- get(name);z$definition <- factor(z$definition,levels=definitions);z$metric <- factor(z$metric,levels=metrics);assign(name,z)}
colors <- c(dff='#BA354F',forward='#222222',rotation='#2767B2')
mask <- data.frame(definition=factor(rep('localized',3),levels=definitions),metric=factor(metrics,levels=metrics),xmin=-.5,xmax=.5,ymin=-Inf,ymax=Inf)
p <- ggplot(d,aes(time,value,group=fly))+
 geom_rect(data=mask,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='#EEEEEE',inherit.aes=FALSE)+
 geom_vline(xintercept=0,linetype='dotted',color='#777777',linewidth=.4)+
 geom_ribbon(data=s,aes(x=time,y=NULL,ymin=mean-sem,ymax=mean+sem,fill=metric,group=metric),inherit.aes=FALSE,alpha=.2)+
 geom_line(color='#999999',alpha=.6,linewidth=.55)+
 geom_line(data=s,aes(y=mean,color=metric,group=metric),linewidth=1.15)+
 facet_grid(metric~definition,scales='free_y',switch='y',labeller=labeller(metric=metric_labels,definition=setNames(labels,definitions)))+
 scale_color_manual(values=colors)+scale_fill_manual(values=colors)+
 scale_x_continuous(limits=c(-2,3),breaks=-2:3,expand=expansion(mult=0))+
 scale_y_continuous(expand=expansion(mult=c(.08,.1)))+
 labs(title='VES041 imaging: no heat',
  subtitle='Matched recordings. Thin lines: fly means. Bold lines and bands: equal-fly mean +/- SEM.',
  x='Time from alignment event (s)',y=NULL,
  caption='Left: strongest deceleration. Right: center of representative 1-s turn window (gray).\nOriginal dF/F; no added baseline subtraction. Six mapped flies include one inferred identity.')+
 theme_classic(base_size=16,base_family='Arial')+
 theme(legend.position='none',strip.background=element_blank(),strip.placement='outside',
  strip.text.x=element_text(size=17,face='bold',lineheight=1.15),strip.text.y.left=element_text(angle=90,size=16,margin=margin(r=10)),
  panel.spacing=grid::unit(.7,'cm'),plot.title=element_text(size=26,face='bold'),plot.subtitle=element_text(size=14,margin=margin(b=15)),
  axis.text=element_text(color='#333333',size=13),axis.title.x=element_text(size=16,margin=margin(t=12)),
  panel.grid.major.y=element_line(color='#E9EBED',linewidth=.3),plot.caption=element_text(size=11,color='#666666',margin=margin(t=12)),
  plot.margin=margin(18,18,12,14))
stem <- file.path(root,'noheat_imaging_original_vs_inhibition_rule')
ggsave(paste0(stem,'.svg'),p,width=12.8,height=10.5,device=svglite::svglite,bg='white')
ggsave(paste0(stem,'.png'),p,width=12.8,height=10.5,dpi=180,bg='white')
