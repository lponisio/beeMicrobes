## library(RColorBrewer)
library(viridis)
plotCommDistbyGroup  <- function(dist.mat,
                                 groups,
                                 species.type,
                                  group.name,
                                 f.path="figures/pcoas"){

    f.pcoa <- function(){
        pcoa.comm <- cmdscale(dist.mat)
        dist.m <- as.dist(dist.mat)
        print("betadisper result")
        print(anova(betadisper(dist.m,groups)))

        pcoa.mod <- adonis(dist.mat~groups)
        plot(NA, asp=1,  cex=1.5,
             ylim=range(pcoa.comm[,2]),
             xlim=range(pcoa.comm[,1]),
             xlab='',
             ylab='',
             xaxt='n',
             yaxt='n',
             cex.lab=1.5)

        ## cols <- brewer.pal(length(unique(groups)), "Set3")
        cols <- viridis(length(unique(groups)))
        names(cols) <- unique(groups)
        for(group in unique(groups)){
            ## all points sitting on top of eachother so triple jitter
            points(pcoa.comm[groups == group,],
                   col=cols[group], pch=16, cex=1.5)
            points(pcoa.comm[groups == group,],
                   col="black", pch=1, cex=1.5)
            col.lines <- cols
            col.lines[names(col.lines) != group] <- "#ffffff00"
            ordihull(pcoa.comm, groups, col=col.lines)
        }
        ## ordihull(pcoa.comm, groups, col=cols)

        legend("topright", legend= unique(groups),
               bty="n", col=cols[unique(groups)],
               pch=16, cex=1)
        legend("topright", legend= unique(groups),
               bty="n", col="black", pch=1,
               cex=1)

        mtext('PCoA1', 1, line=2, cex=1.5)
        mtext('PCoA2', 2, line=2, cex=1.5)
        return(pcoa.mod)
    }
    ## function for plotting PcoA axes
    pdf.f(f.pcoa,
          file= file.path(f.path,
                          sprintf("%s_%s_pcoa.pdf",
                                  species.type, group.name)),
          width=7, height=7)
}


