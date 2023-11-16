pdf.f <- function(f, file, ...) {
    cat(sprintf('Writing %s\n', file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
}

## add transparency to named colors
add.alpha <- function(col, alpha=0.2){
    apply(sapply(col, col2rgb)/255, 2,
          function(x)
              rgb(x[1], x[2], x[3],
                  alpha=alpha))
}



## This functions takes site-species-abundance data and creates a
## matrix where the sites are columns and the rows are species.

samp2site.spp <- function(site, spp, abund, FUN=sum) {
  x <- tapply(abund, list(site = site, spp = spp), FUN)
  x[is.na(x)] <- 0
  return(x)
}



phyloHyp <- function(mod, mod.name, mod.gaussian=TRUE){

    ## plot(mod, N = 4, ask = FALSE)
    ## ggsave(sprintf("figures/diagnostics/%s_Diag.pdf", mod.name),
    ##        height=11, width=8.5)

    pdf(sprintf("figures/diagnostics/%s_Diag.pdf", mod.name),
        height=11, width=8.5)
    plot(mod,  N = 4, ask = FALSE)
    dev.off()

    ## if(mod.gaussian){
    ##     hyp <- "sd_GenusSpecies__Intercept^2 / (sd_GenusSpecies__Intercept^2 + sd_GenusSpecies2__Intercept^2 + sigma^2) = 0"
    ##     hyp <- hypothesis(mod, hyp, class = NULL)
    ##     plot(hyp, main="Phylogenetic signal")
    ##     ggsave(sprintf("figures/diagnostics/phyloInt_%s.pdf", mod.name),
    ##            height=4, width=4)
    ## }

  
    if(mod.gaussian){
        hyp <- "sd_GenusSpecies__Intercept^2 / (sd_GenusSpecies__Intercept^2 + sigma^2) = 0"
        hyp <- hypothesis(mod, hyp, class = NULL)
        quartz()
        plot(hyp, main="Phylogenetic signal")
        ggsave(filename=sprintf("figures/diagnostics/phyloInt_%s.pdf", mod.name),
               height=4, width=4)
    }
}

shared_legend <- function(...,
                          ncol = length(list(...)),
                          nrow = 1,
                          position = c("bottom",
                                       "right", "top", "left"),
                          plot = TRUE){
    ## shared legend across plots, modified from lemon package
    plots <- list(...)
    position <- match.arg(position)
    legend <- g_legend(plots[[1]] + theme(legend.position = position,
                                          legend.text =
                                              element_text(size=6),
                                          legend.key.size = unit(3, "point"),
                                          legend.title=element_blank()))
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) {
        if (is.ggplot(x)) {
            x + theme(legend.position = "none")
        }
        else {
            x
        }
    })
    gl <- c(gl, ncol = ncol, nrow = nrow)
    combined <- switch(position, top = arrangeGrob(legend, do.call(arrangeGrob,
        gl), ncol = 1, heights = grid::unit.c(lheight, unit(1, "npc") -
        lheight)), bottom = arrangeGrob(do.call(arrangeGrob,
        gl), legend, ncol = 1, heights = grid::unit.c(unit(1, "npc") -
        lheight, lheight)), left = arrangeGrob(legend, do.call(arrangeGrob,
        gl), ncol = 2, widths = grid::unit.c(lwidth, unit(1, "npc") -
        lwidth)), right = arrangeGrob(do.call(arrangeGrob, gl),
        legend, ncol = 2, widths = grid::unit.c(unit(1, "npc") - lwidth,
            lwidth)))
    if (plot) {
        grid::grid.newpage()
        grid::grid.draw(combined)
    }
    invisible(combined)
}
