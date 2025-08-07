library(putils)
source('data_generation.R')

settings <- c('sinusoidal', 'zigzag', 'circle', 'spiral', 'Lissajous', 'local')
noises <- c(0, 0.2, 0.4, 1.0)
counter <- 0
for (setting in settings){
  for (noise in noises){
    counter <- counter + 1
    pdf(paste0('fig/figS1', letters[counter], '.pdf'), width=4, height=4)
    par(tcl=-0.3, mar=c(3.5,3.5,1,1), mgp=c(2.1,0.6,0), cex.lab=1.5, cex.axis=1.3, family='serif')
    bunch(x, y) %=% data_generation(1000, type=setting, noise=noise, dx=1)
    plot(x, y, xlab='X', ylab='Y')
    dev.off()
  }
}
  