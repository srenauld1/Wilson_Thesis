root <- 'outputs/ves041_example_reorientation_timeseries_20261104'
dir.create(file.path(root,'details'),showWarnings=FALSE)
dir.create(file.path(root,'events'),showWarnings=FALSE)
purple <- '#A020F0'; orange <- '#D55E00'; ink <- '#252525'
records <- lapply(c(8,45),function(id) list(id=id,
  label=if(id==8)'VES041_1 x ACR'else'VES041_1 x GFP',
  d=read.csv(file.path(root,paste0('fly',id,'_timeseries.csv'))),
  e=read.csv(file.path(root,paste0('fly',id,'_events.csv'))),
  light=read.csv(file.path(root,paste0('fly',id,'_light_intervals.csv'))),
  scored=read.csv(file.path(root,paste0('fly',id,'_scored_windows.csv')))))
forward_lim <- range(unlist(lapply(records,function(r)r$d$forward_mm_s)))
yaw_max <- ceiling(max(unlist(lapply(records,function(r)abs(r$d$yaw_deg_s))))/100)*100

path_labels <- function(x,y,ids) {
  # Place labels in physical plot coordinates, with leader lines to true positions.
  usr <- par('usr'); pin <- par('pin')
  xx <- (x-usr[1])/diff(usr[1:2])*pin[1]
  yy <- (y-usr[3])/diff(usr[3:4])*pin[2]
  placed <- list()
  for(i in seq_along(ids)) {
    w <- strwidth(as.character(ids[i]),units='inches',cex=.68)+.055
    h <- strheight('0',units='inches',cex=.68)+.045
    candidates <- expand.grid(radius=c(.16,.27,.40,.56,.75,1),angle=seq(0,7)*pi/4)
    best <- NULL; best_cost <- Inf
    for(j in seq_len(nrow(candidates))) {
      cx <- xx[i]+candidates$radius[j]*cos(candidates$angle[j])
      cy <- yy[i]+candidates$radius[j]*sin(candidates$angle[j])
      box <- c(cx-w/2,cx+w/2,cy-h/2,cy+h/2)
      if(box[1]<0 || box[2]>pin[1] || box[3]<0 || box[4]>pin[2])next
      overlaps <- sum(vapply(placed,function(b)box[1]<b[2]&&box[2]>b[1]&&box[3]<b[4]&&box[4]>b[3],logical(1)))
      near_points <- sum(abs(xx-cx)<w/2 & abs(yy-cy)<h/2)
      cost <- overlaps*100+near_points*2+candidates$radius[j]
      if(cost<best_cost){best<-box;best_cost<-cost}
    }
    if(is.null(best))next
    placed[[length(placed)+1]] <- best
    bx <- usr[1]+best[1:2]/pin[1]*diff(usr[1:2]); by <- usr[3]+best[3:4]/pin[2]*diff(usr[3:4])
    segments(x[i],y[i],mean(bx),mean(by),col=adjustcolor(orange,.65),lwd=.6)
    rect(bx[1],by[1],bx[2],by[2],col='white',border=NA)
    text(mean(bx),mean(by),labels=ids[i],cex=.68,col=orange,font=2)
  }
}

draw_path <- function(r,lo=0,hi=1200,whole=TRUE) {
  par(mar=c(4,4,2,1))
  d <- r$d; use <- d$time_s>=lo & d$time_s<=hi
  e <- r$e[r$e$marker_s>=lo & r$e$marker_s<hi,,drop=FALSE]
  if(whole){xl<-c(-600,700);yl<-c(-700,600)}else{
    xl<-range(d$x_mm[use]);yl<-range(d$y_mm[use]);pad<-max(diff(xl),diff(yl))*.12+2
    xl<-xl+c(-pad,pad);yl<-yl+c(-pad,pad)
  }
  plot(d$x_mm[use],d$y_mm[use],type='n',asp=1,xlim=xl,ylim=yl,
       xlab='X position (mm)',ylab='Y position (mm)',bty='l',las=1,
       main=if(whole)'Full flat path'else sprintf('Flat path: %.1f-%.1f s',lo,hi))
  lines(d$x_mm[use],d$y_mm[use],col=ink,lwd=.65)
  ix <- which(use & d$light_on==1)
  xx<-d$x_mm;yy<-d$y_mm;xx[!(use & d$light_on==1)]<-NA;yy[!(use & d$light_on==1)]<-NA
  lines(xx,yy,col=purple,lwd=1.1)
  if(nrow(e))for(i in seq_len(nrow(e))){
    sel<-d$time_s>=e$start_s[i]&d$time_s<=e$end_s[i]
    lines(d$x_mm[sel],d$y_mm[sel],col=orange,lwd=2.6)
  }
  points(d$x_mm[which(use)[1]],d$y_mm[which(use)[1]],pch=21,bg='#32B64A',col='white',cex=.85)
  if(nrow(e)){
    points(e$x_mm,e$y_mm,pch=21,col=orange,bg='white',cex=.65,lwd=1.1)
    path_labels(e$x_mm,e$y_mm,e$event_id)
  }
}

shade_light <- function(r,lo,hi,yl) {
  for(i in seq_len(nrow(r$light)))if(r$light$end_s[i]>lo&&r$light$start_s[i]<hi)
    rect(max(lo,r$light$start_s[i]),yl[1],min(hi,r$light$end_s[i]),yl[2],
         col=adjustcolor(purple,.13),border=NA)
}

draw_raster <- function(r,lo,hi) {
  par(mar=c(.6,5,1.8,1))
  plot(NA,xlim=c(lo,hi),ylim=c(0,1.15),xaxs='i',yaxs='i',axes=FALSE,xlab='',ylab='')
  for(i in seq_len(nrow(r$scored))){
    w<-r$scored[i,];if(w$end_s>lo&&w$start_s<hi)
      rect(max(lo,w$start_s),.08,min(hi,w$end_s),.32,
           col=if(w$phase=='pre')'#D3D3D3'else adjustcolor(purple,.28),border=NA)
  }
  e<-r$e[r$e$marker_s>=lo & r$e$marker_s<hi,,drop=FALSE]
  if(nrow(e))for(i in seq_len(nrow(e))){
    segments(e$start_s[i],.55,e$end_s[i],.55,col=orange,lwd=3)
    points(e$marker_s[i],.55,pch=16,col=orange,cex=.55)
    text(e$marker_s[i],if(i%%2).84 else 1.05,labels=e$event_id[i],cex=.63,col=orange,font=2)
  }
  mtext('Scored / events',side=2,line=1,cex=.8)
  title(sprintf('%.1f-%.1f s | numbered reorientations',lo,hi),cex.main=1,adj=0)
}

draw_trace <- function(r,lo,hi,column,label,yl,bottom=FALSE) {
  par(mar=c(if(bottom)3.7 else .65,5,.3,1))
  d<-r$d;use<-d$time_s>=lo&d$time_s<=hi
  plot(d$time_s[use],d[[column]][use],type='n',xlim=c(lo,hi),ylim=yl,xaxs='i',
       xlab=if(bottom)'Time from recording start (s)'else'',ylab=label,
       bty='l',las=1,xaxt=if(bottom)'s'else'n',cex.lab=.9)
  shade_light(r,lo,hi,par('usr')[3:4])
  abline(h=0,col='#DEDEDE',lwd=.6)
  lines(d$time_s[use],d[[column]][use],col=ink,lwd=.55)
  e<-r$e[r$e$end_s>=lo&r$e$start_s<=hi,,drop=FALSE]
  if(nrow(e))for(i in seq_len(nrow(e))){
    abline(v=e$marker_s[i],col=adjustcolor(orange,.35),lty=3,lwd=.65)
    sel<-d$time_s>=max(lo,e$start_s[i])&d$time_s<=min(hi,e$end_s[i])
    lines(d$time_s[sel],d[[column]][sel],col=orange,lwd=1.5)
  }
}

heading_lim <- function(r,lo,hi){
  limits<-range(r$d$heading_change_deg[r$d$time_s>=lo&r$d$time_s<=hi])
  limits+c(-1,1)*max(10,diff(limits)*.06)
}
header <- function(r,detail=FALSE) {
  mtext(sprintf('Fly %d | %s | %d reorientations (%d pre-light, %d light on)',r$id,r$label,nrow(r$e),
    sum(r$e$phase=='pre'),sum(r$e$phase=='inhibition')),side=3,outer=TRUE,line=2.8,cex=1.25,font=2)
  mtext('Purple: light on   |   Orange: detected 1 s window; matching event IDs on path and traces   |   Gray bar: scored pre-light',
        side=3,outer=TRUE,line=1.1,cex=.82)
  mtext('Existing detections: >30 deg net turn, <=2 mm net displacement, directionality >=0.7. Only 5 s pre-light and light-on windows were scored.',
        side=1,outer=TRUE,line=.6,cex=.8)
}
render <- function(r,mode='overview',lo=0,hi=1200) {
  par(family='Arial',cex=.95,oma=c(3.5,.2,5,.5),mgp=c(2.4,.7,0))
  if(mode=='path'){
    par(oma=c(3,.2,5,.5));draw_path(r);header(r);return(invisible(NULL))
  }
  if(mode=='timeseries')layout(matrix(1:4,ncol=1),heights=c(.65,1,1,1))else
    layout(matrix(c(1,2,1,3,1,4,1,5),ncol=2,byrow=TRUE),widths=c(1,1.65),heights=c(.65,1,1,1))
  if(mode!='timeseries')draw_path(r,lo,hi,whole=(mode=='overview'))
  draw_raster(r,lo,hi)
  draw_trace(r,lo,hi,'forward_mm_s','Forward velocity\n(mm/s)',forward_lim)
  draw_trace(r,lo,hi,'yaw_deg_s','Yaw velocity\n(deg/s)',c(-yaw_max,yaw_max))
  draw_trace(r,lo,hi,'heading_change_deg','Heading change\n(deg)',heading_lim(r,lo,hi),TRUE)
  header(r)
}

for(r in records){
  for(mode in c('overview','timeseries','path')){
    name<-file.path(root,sprintf('VES041_fly%03d_%s_reorientations',r$id,mode))
    w<-if(mode=='path')10 else if(mode=='timeseries')19 else 23
    svg(paste0(name,'.svg'),width=w,height=10,bg='white');render(r,mode);dev.off()
    png(paste0(name,'.png'),width=w*120,height=1200,res=120,bg='white');render(r,mode);dev.off()
  }
  for(lo in seq(0,1000,by=200)){
    name<-file.path(root,'details',sprintf('VES041_fly%03d_%04d_%04ds',r$id,lo,lo+200))
    svg(paste0(name,'.svg'),width=20,height=9,bg='white');render(r,'detail',lo,lo+200);dev.off()
    png(paste0(name,'.png'),width=2000,height=900,res=100,bg='white');render(r,'detail',lo,lo+200);dev.off()
  }
  for(i in seq_len(nrow(r$e))){
    e<-r$e[i,];lo<-max(0,e$start_s-3);hi<-min(1200,e$end_s+3)
    name<-file.path(root,'events',sprintf('VES041_fly%03d_event%03d_%s',r$id,e$event_id,
      if(e$phase=='pre')'pre_light'else'light_on'))
    svg(paste0(name,'.svg'),width=17,height=9,bg='white');render(r,'detail',lo,hi);dev.off()
    # Preview the first event from each phase, without choosing by its shape.
    if(i==which(r$e$phase==e$phase)[1]){
      png(paste0(name,'.png'),width=2040,height=1080,res=120,bg='white');render(r,'detail',lo,hi);dev.off()
    }
  }
  cat('Saved full path, time series, overview and six detail panels for fly',r$id,'\n')
}
