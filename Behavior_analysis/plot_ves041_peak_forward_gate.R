root <- 'outputs/ves041_peak_forward_gate_20261104'
d <- read.csv(file.path(root,'fly_summary.csv'),stringsAsFactors=FALSE)
d <- d[d$definition=='raw_signed',]
plot_breakdown <- function(){
  par(mfrow=c(1,2),family='Arial',mar=c(4.5,4.5,3.2,1.5),oma=c(2,0,3.5,0),cex=1.15)
  mx <- ceiling(max(c(d$pre_rate_per_min,d$light_on_rate_per_min))/5)*5
  for(g in c('experimental','control')){
    z <- d[d$group==g,];n<-nrow(z);yy<-rev(seq_len(n))
    plot(NA,xlim=c(0,mx+1),ylim=c(.5,n+.7),axes=FALSE,xlab='Reorientations / min',ylab='',xaxs='i')
    axis(1,col='#999999',col.axis='#555555',lwd=.6);axis(2,at=yy,labels=z$fly_id,las=1,tick=FALSE,col.axis='#303030')
    mtext('Fly ID',side=2,line=2.5)
    abline(v=axTicks(1),col='#EFEFEF',lwd=.7)
    segments(z$pre_rate_per_min,yy,z$light_on_rate_per_min,yy,col='#BBBBBB',lwd=2)
    points(z$pre_rate_per_min,yy,pch=16,col='#707070',cex=1.2)
    points(z$light_on_rate_per_min,yy,pch=16,col='#963DCE',cex=1.2)
    title(if(g=='experimental')'ACR (8 flies)'else'GFP control (9 flies)',adj=0,cex.main=1.2)
    legend('topright',legend=c('Pre-light','Light on'),pch=16,col=c('#707070','#963DCE'),bty='n',cex=.82)
  }
  mtext('Forward velocity at peak yaw: at most 1 mm/s',outer=TRUE,side=3,line=1.7,font=2,cex=1.25)
  mtext('Same good walkers; signed raw velocity; displacement restriction removed',outer=TRUE,side=3,line=.1,cex=.85)
  mtext('Each dot is one fly; matched pre-light and light-on exposure. >30 deg net turn and directionality >=0.7 retained.',outer=TRUE,side=1,line=.6,cex=.8)
}
svg(file.path(root,'good_walkers_peak_forward_gate.svg'),width=12,height=7,bg='white');plot_breakdown();dev.off()
png(file.path(root,'good_walkers_peak_forward_gate.png'),width=2160,height=1260,res=180,bg='white');plot_breakdown();dev.off()
