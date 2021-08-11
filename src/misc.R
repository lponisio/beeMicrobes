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

