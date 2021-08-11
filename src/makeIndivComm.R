
makeIndivComm <- function(spec, col.sp.names){
    comm.indiv <- spec[, col.sp.names]
    rownames(comm.indiv) <- spec$UniqueID
    ## find those that were not screened to drop them
    not.screened <- apply(comm.indiv, 1, function(x) all(is.na(x)))
    comm.indiv <- comm.indiv[!not.screened, ]

    comm.indiv[is.na(comm.indiv)] <- 0
    comm.indiv  <- bipartite::empty(comm.indiv)
    return(comm.indiv)
}
