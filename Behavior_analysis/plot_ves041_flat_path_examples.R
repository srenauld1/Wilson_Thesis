suppressPackageStartupMessages(library(jsonlite))
root <- 'outputs/ves041_blinded_good_walkers_pivots_gt30_20261102/flat_path_examples'
records <- fromJSON(file.path(root,'all_flat_paths.json'),simplifyVector=FALSE)
selected <- fromJSON(file.path(root,'selected_examples.json'),simplifyVector=FALSE)
dir.create(file.path(root,'atlas'),showWarnings=FALSE)
values <- function(x) as.numeric(unlist(x))
get_record <- function(key) records[[which(vapply(records,function(r)r$recording_key==key,logical(1)))]]
get_path <- function(r,pulse,phase) r$paths[[which(vapply(r$paths,function(p)p$pulse==pulse && p$phase==phase,logical(1)))]]

crop_path <- function(p,lo,hi) {
  index <- values(p$sample_indices);keep <- index>=lo & index<=hi
  x <- values(p$x)[keep];y <- values(p$y)[keep]
  stopifnot(length(x)>2)
  p$x <- as.list(x-x[1]);p$y <- as.list(y-y[1]);p$sample_indices <- as.list(index[keep])
  p
}
limits <- function(pre,on) {
  xx <- range(c(values(pre$x),values(on$x)));yy <- range(c(values(pre$y),values(on$y)))
  size <- max(diff(xx),diff(yy),1);pad <- .08*size
  list(x=xx+c(-pad,pad),y=yy+c(-pad,pad))
}
draw_path <- function(p,lim,color,title_text='',axes=TRUE,cex_main=1,highlight=TRUE) {
  x <- values(p$x);y <- values(p$y);idx <- values(p$sample_indices)
  plot(x,y,type='n',xlim=lim$x,ylim=lim$y,asp=1,axes=FALSE,xlab='',ylab='',main=title_text,cex.main=cex_main)
  if(axes){axis(1,cex.axis=.85,col.axis='#555555');axis(2,las=1,cex.axis=.85,col.axis='#555555');box(col='#CCCCCC')}
  lines(x,y,col=color,lwd=1.8)
  if(highlight && length(p$event_windows)) for(w in p$event_windows) {
    use <- idx>=w$start & idx<=w$end
    lines(x[use],y[use],col='#7652A2',lwd=2.8)
  }
  points(x[1],y[1],pch=16,cex=.65,col='#202020')
  points(tail(x,1),tail(y,1),pch=2,cex=.8,lwd=1.2,col=color)
}

main_figure <- function(zoom=FALSE,cases=selected) {
  par(mfrow=c(length(cases),2),mar=c(3.1,3.2,4,.8),oma=c(4,2.6,3.6,.5),family='Arial')
  for(i in seq_along(cases)) {
    s <- cases[[i]];r <- get_record(s$recording_key)
    pre <- get_path(r,s$pulse,'pre');on <- get_path(r,s$pulse,'inhibition')
    duration_text <- '5 s';counts <- c(pre$events,on$events)
    if(zoom) {
      start <- min(vapply(pre$event_windows,function(w)w$start,numeric(1)))
      end <- pre$event_windows[[which.min(vapply(pre$event_windows,function(w)w$start,numeric(1)))]]$end
      pre <- crop_path(pre,start,end);on <- crop_path(on,start,end);duration_text <- '1 s'
    }
    lim <- limits(pre,on)
    label <- paste0(ifelse(s$group=='experimental','Experimental','Control'),' ID ',s$fly_id,', pulse ',s$pulse)
    draw_path(pre,lim,'#707070',paste0(label,'\nPre-opto (',duration_text,')',if(!zoom)paste0(': ',counts[1],ifelse(counts[1]==1,' reorientation',' reorientations'))else''))
    draw_path(on,lim,'#BA354F',paste0('Opto-on (',duration_text,')',if(!zoom)paste0(': ',counts[2],ifelse(counts[2]==1,' reorientation',' reorientations'))else''))
  }
  mtext(if(zoom)'Matched flat-path snippets'else'Paired flat paths from individual presentations',outer=TRUE,side=3,line=1.1,cex=1.3,font=2)
  mtext('Relative x (mm)',outer=TRUE,side=1,line=1,cex=1.05)
  mtext('Relative y (mm)',outer=TRUE,side=2,line=.7,cex=1.05)
  mtext(if(zoom)'Same elapsed 1 s in each phase, centered on the first pre-opto reorientation.'else'Each segment is 5 s. Paths are translated only; each paired row uses the same scale.',outer=TRUE,side=1,line=2.3,cex=.76)
  mtext('Purple: detected reorientation window. Dot: start; triangle: end.',outer=TRUE,side=3,line=-.25,cex=.78)
}

for(zoom in c(FALSE,TRUE)) {
  name <- if(zoom)'matched_1s_flat_path_examples'else'paired_5s_flat_path_examples'
  png(file.path(root,paste0(name,'.png')),width=2000,height=2500,res=190,bg='white')
  main_figure(zoom);dev.off()
  svg(file.path(root,paste0(name,'.svg')),width=10.5,height=13.15,bg='white')
  main_figure(zoom);dev.off()
}
two_cases <- selected[vapply(selected,function(s)s$fly_id %in% c('8','30'),logical(1))]
png(file.path(root,'two_baseline_matched_flies.png'),width=2000,height=1450,res=190,bg='white')
main_figure(FALSE,two_cases);dev.off()
svg(file.path(root,'two_baseline_matched_flies.svg'),width=10.5,height=7.63,bg='white')
main_figure(FALSE,two_cases);dev.off()

# Display all 40 paired presentations for every included fly, without trial selection.
index <- c('ALL-PRESENTATION FLAT-PATH ATLAS','Each file contains all 40 presentations from one fly. Gray: pre-opto 5 s; red: opto-on 5 s. Purple marks detected reorientation windows. Each path is translated to start at the origin; no rotation/reflection. Scale is equal within a tile but may differ between tiles. Labels give pulse number and pre/opto reorientation counts.','')
for(r in records) {
  id <- paste(unlist(r$cohort_ids),collapse=',');name <- paste0(r$group,'_ID',id,'_all40.png')
  png(file.path(root,'atlas',name),width=2500,height=1750,res=160,bg='white')
  par(mfrow=c(5,8),mar=c(.4,.4,1.8,.2),oma=c(2,1,3,1),family='Arial')
  for(pulse in 1:40) {
    pre <- get_path(r,pulse,'pre');on <- get_path(r,pulse,'inhibition');lim <- limits(pre,on)
    draw_path(pre,lim,'#8A8A8A',paste0(pulse,': ',pre$events,' / ',on$events),axes=FALSE,cex_main=.8)
    # Add the second measured path in the same coordinate system.
    x <- values(on$x);y <- values(on$y);idx <- values(on$sample_indices)
    lines(x,y,col='#BA354F',lwd=1.3)
    if(length(on$event_windows))for(w in on$event_windows){use <- idx>=w$start & idx<=w$end;lines(x[use],y[use],col='#7652A2',lwd=2.1)}
    points(x[1],y[1],pch=16,cex=.5);points(tail(x,1),tail(y,1),pch=2,cex=.55,col='#BA354F')
  }
  mtext(paste0(ifelse(r$group=='experimental','Experimental','Control'),' ID ',id,': all 40 presentations'),side=3,outer=TRUE,line=1,cex=1.4,font=2)
  mtext('Gray: pre-opto; red: opto-on; purple: reorientation. Tile title: pulse, pre / opto counts. Scale varies between tiles.',side=1,outer=TRUE,line=.4,cex=.85)
  dev.off()
  index <- c(index,paste0(r$group,' ID ',id,': ',file.path(root,'atlas',name)))
}
writeLines(index,file.path(root,'atlas_index.txt'))
cat('Saved paired examples, matched zooms, and all-presentation atlases for ',length(records),' flies.\n',sep='')
