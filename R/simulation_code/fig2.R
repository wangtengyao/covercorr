library(putils)

results <- readLines('power.out')
results_spl <- strsplit(results, ' = |, ')
results_df <- do.call(rbind, results_spl)
results_df <- results_df[, c(F,T)]
colnames(results_df) <- results_spl[[1]][c(T,F)]
rownames(results_df) <- NULL
results_df <- as.data.frame(results_df)
print(dim(results_df))
for (j in c(3:17)) results_df[,j] <- as.numeric(results_df[,j])
head(results_df)
results_df[results_df$n==2000, c('p_xi')] <- 1
df_agg <- aggregate(cbind(p_kappa, p_xi, p_dcor, p_hsic, p_kmac, p_USP) ~ type + n + dx + noise, 
                    FUN=function(v){mean(v <= 0.05)}, data=results_df)

palet <- putils::matplotlib_palette()


counter <- 0
types <- c('sinusoidal', 'zigzag', 'circle', 'spiral', 'Lissajous', 'local')
for (type in types){
  counter <- counter + 1
  df_sub <- df_agg[df_agg$type==type, ]
  library(tidyr)
  library(dplyr)
  
  df_long <- df_sub %>%
    pivot_longer(
      cols = starts_with("p_"),           # or cols = 5:10
      names_to = "method",
      names_prefix = "p_",
      values_to = "power"
    )
  
  methods <- unique(df_long$method)
  methods_proper <- c(expression(kappa), expression(xi), 'dCor', 'HSIC', 'KMAc', 'USP')
  ns <- unique(df_long$n)
  pdf(paste0('fig/fig2', letters[counter], '.pdf'), width=6, height=4.5)
  par(tcl=-0.3, mar=c(3.5,3.5,2.5,1), mgp=c(2.1,0.6,0), cex.lab=1.5, cex.axis=1.3, family='serif')
  plot(x=c(0,2), y=c(0,1), type='n', xlab='noise', ylab='power', main=type)
  for (i in seq_along(methods)){
    for (j in seq_along(ns)){
      if (methods[i] == 'xi' && ns[j]==2000) next
      df_small <- df_long[df_long$n==ns[j] & df_long$method==methods[i], ]
      points(df_small$noise, df_small$power, col=palet[i], lty=j, type='b', pch=i, cex=0.8, lwd=2)
    }
  }
  legend(ifelse(counter!=1, 'topright', 'right'), legend=c(methods_proper, 'n=1000','n=2000'), 
         col=c(palet[1:6],rep('black',2)), lty=c(rep(1,6),1:2), pch=c(1:6,NA,NA), pt.cex=0.8, bty='n')
  dev.off()
}


