suppressPackageStartupMessages({library(jsonlite); library(ggplot2)})
root <- 'outputs/ves041_inhibition_duration_audit_20261102'
d <- fromJSON(file.path(root, 'behavior_quality_results.json'))
d <- d[d$status == 'ok', ]
d$group_short <- ifelse(grepl('^experimental', d$group), 'Experimental',
                        ifelse(grepl('^control', d$group), 'Control', 'Genotype unassigned'))
d$cohort_id <- vapply(d$cohort_ids, function(x) paste(x, collapse=','), character(1))
d$label <- paste0(ifelse(d$cohort_id == '', '', paste0('ID ', d$cohort_id, ': ')),
                  substr(d$date, 6, 10), ' fly ', d$daily_fly,
                  ' c', d$cell, '/', d$trial)
d$set <- paste(d$cohort, d$duration_s, 's')
colors <- c(Experimental='#C73E4D', Control='#606060', 'Genotype unassigned'='#BF8F32')
theme_set(theme_bw(base_size=10))
for (set in unique(d$set)) {
  x <- d[d$set == set, ]
  x <- x[order(x$baseline_moving_fraction), ]
  x$label <- factor(x$label, levels=x$label)
  p <- ggplot(x, aes(baseline_moving_fraction*100, label, fill=group_short)) +
    geom_col(width=.72) +
    geom_vline(xintercept=c(20, 50), linetype='dashed', colour='grey55', linewidth=.3) +
    geom_text(aes(label=paste0(n_moving_baseline_pulses, '/', n_pulses)),
              x=105, hjust=0, size=2.8) +
    scale_fill_manual(values=colors) +
    scale_x_continuous(limits=c(0, 121), breaks=seq(0,100,20)) +
    labs(title=set, subtitle='Movement during the 1 s immediately before each light onset',
         x='Pre-light samples moving (%)', y=NULL, fill=NULL,
         caption='Moving = total speed >=3 mm/s (existing lab threshold). Right labels: moving-baseline presentations / all presentations.\nA moving baseline has >=50% moving samples and >=95% finite kinematics. Dashed lines are descriptive coverage guides.') +
    theme(legend.position='bottom', plot.caption=element_text(hjust=0), panel.grid.major.y=element_blank())
  slug <- gsub('[^a-z0-9]+','_',tolower(set))
  ggsave(file.path(root,paste0(slug,'_behavior_quality.png')),p,
         width=10,height=max(4,nrow(x)*.24+2),dpi=180)
}

# Time coverage panel for the two main protocols.
for (set in c('blinded round 1 5 s','unblinded RubyACR 10 s')) {
  x <- d[d$set == set, ]; x <- x[order(x$baseline_moving_fraction), ]
  bins <- do.call(rbind, lapply(seq_len(nrow(x)), function(i) {
    b <- x$time_bins[[i]]; b$label <- as.character(x$label[i]); b
  }))
  bins$label <- factor(bins$label, levels=as.character(x$label))
  p <- ggplot(bins,aes(time_s/60,label,fill=moving_fraction*100)) +
    geom_tile() + scale_fill_gradient(low='#F3F3F3',high='#125C80',limits=c(0,100)) +
    labs(title=paste(set,': movement across the recording'),
         subtitle='Each tile summarizes 10 s; includes both light-on and light-off periods',
         x='Recording time (min)',y=NULL,fill='Moving (%)') +
    theme(panel.grid=element_blank(),legend.position='bottom')
  slug <- gsub('[^a-z0-9]+','_',tolower(set))
  ggsave(file.path(root,paste0(slug,'_movement_over_time.png')),p,
         width=10,height=max(5,nrow(x)*.24+2),dpi=180)
}
