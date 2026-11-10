p <- 'outputs/ves041_example_vs_pdf_audit_20261105'
d <- read.csv(file.path(p,'rotation_means_same_axes.csv'))
draw <- function(){
 par(family='Arial',mar=c(4.8,5,3.4,1.2),cex=1.2)
 plot(NA,xlim=c(-2,2),ylim=c(0,180),xlab='Time from strongest deceleration (s)',ylab='Mean rotational speed (deg/s)',bty='l',las=1,xaxs='i',yaxs='i',main='Rotation means: same fly, different event subsets')
 abline(v=0,col='#AAAAAA',lty=3,lwd=1)
 fields<-c('pdf_heat_on','pdf_heat_off','screenshot_control');cols<-c('#9D3C9E','#1454CC','#303030');lty<-c(1,1,2)
 for(i in 1:3){z<-d[d$curve==fields[i],];lines(z$time_s,z$mean_deg_s,col=cols[i],lwd=2.5,lty=lty[i])}
 legend('topleft',legend=c('Screenshot: placed control (15 events)','PDF left: heat off (20 events)','PDF right: heat on (19 events)'),col=c('#303030','#1454CC','#9D3C9E'),lty=c(2,1,1),lwd=2.5,bty='n',cex=.9)
}
svg(file.path(p,'rotation_means_same_axes.svg'),width=10,height=6,bg='white');draw();dev.off()
png(file.path(p,'rotation_means_same_axes.png'),width=1800,height=1080,res=180,bg='white');draw();dev.off()
