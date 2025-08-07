results <- readLines('simulation_results/size.out')
results_spl <- strsplit(results, ' = |, ')
results_df <- do.call(rbind, results_spl)
results_df <- results_df[, c(F,T)]
colnames(results_df) <- results_spl[[1]][c(T,F)]
rownames(results_df) <- NULL
results_df <- as.data.frame(results_df)
print(dim(results_df))
for (j in c(1:3)) results_df[,j] <- as.numeric(results_df[,j])
levels <- c(0.01, 0.025, 0.05, 0.1)
for (level in levels){
  results_df[[paste0('percent', level*100)]] <- results_df$pval <= level
}
df_agg <- aggregate(cbind(percent1, percent2.5, percent5, percent10) ~ n+dx, FUN=mean, data=results_df)
se_agg <- df_agg
nreps <- 5000
for (j in 3:6){se_agg[[j]] <- sqrt(df_agg[[j]] * (1 - df_agg[[j]]) / nreps)}
df_agg
se_agg
library(putils)

palet <- putils::matplotlib_palette()


