library("png")
source("funs.R")
set.seed(123)

# Manhattan

# setup plot
par(mai=rep(0,4)) # no margins

# layout the plots into a matrix
layout(matrix(1:30, ncol=5, byrow=TRUE))

# do the plotting
for(trait in biomarkers) {
  plot(NA,xlim=0:1,ylim=0:1,bty="n",axes=0,xaxs = 'i',yaxs='i')
  img <- readPNG(paste0("data/", trait, "_gwas_man_var.png"))
  rasterImage(img,0,0,1,1)
}

# write to PDF
dev.print(pdf, paste0("data/", "gwas_man_var.pdf"))
dev.off()

# QQ

# setup plot
par(mai=rep(0,4)) # no margins

# layout the plots into a matrix
layout(matrix(1:30, ncol=5, byrow=TRUE))

# do the plotting
for(trait in biomarkers) {
  plot(NA,xlim=0:1,ylim=0:1,bty="n",axes=0,xaxs = 'i',yaxs='i')
  img <- readPNG(paste0("data/", trait, "_gwas_qq_var.png"))
  rasterImage(img,0,0,1,1)
}

# write to PDF
dev.print(pdf, paste0("data/", "gwas_qq_var.pdf"))
dev.off()