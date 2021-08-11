
samp2site.spp <- function(site,spp,abund) {
    x <- tapply(abund, list(site=site,spp=spp), sum)
    x[is.na(x)] <- 0
    return(x)
}

makeCommStruct <- function(spec.dat, type){
    ## prep site by species matrix
    prep.comm <- aggregate(spec.dat[, type],
                           list(site= spec.dat$Site,
                                status= spec.dat$SiteType,
                                sp= spec.dat[, type]), length)

    comm <-  samp2site.spp(site= prep.comm$site,
                           spp= prep.comm$sp, abund=
                                                  prep.comm$x)
    sites <- rownames(comm)
    site.type <- spec.dat$SiteType[match(rownames(comm),
                                         spec.dat$Site)]
    adjsf <- spec.dat$AdjSF[match(rownames(comm),
                                  spec.dat$Site)]

    comm <- bipartite::empty(comm)

    return(list(comm=comm,
                ## sites=sites,
                site.type = site.type,
                adjsf = adjsf))
}


makeStructFromComm <- function(spec, spp){
    spp.pre.comm <- spec[, c("Site", spp)]
    spp.pre.comm <- spp.pre.comm[!apply(spec[, spp], 1,
                                                  function(x) all(is.na(x))),]

    spp.pre.comm <- spp.pre.comm  %>%
        group_by(Site) %>%
        summarise_each(list(mean))

    site.type <- spec$SiteType[match(spp.pre.comm$Site,
                                     spec$Site)]

    adjsf <- spec$AdjSF[match(spp.pre.comm$Site,
                             spec$Site)]

    sites <- spp.pre.comm$Site

    comm <- spp.pre.comm
    comm$Site <- NULL
    comm <- as.matrix(comm)
    rownames(comm) <- spp.pre.comm$Site
    comm[is.na(comm)] <- 0
    comm <- bipartite::empty(comm)

    list(comm=comm,
         ## sites=sites,
         site.type=site.type,
         adjsf = adjsf)
}


getParComm <- function(parasite){
    parasite <- aggregate(list(Parasite=spec[, parasite]),
                          list(GenusSpecies=spec$GenusSpecies,
                               Site=spec$AltSiteName,
                               SiteType=spec$SiteType),
                          function(x) sum(x) / length(x))

    parasite.comm <- samp2site.spp(parasite$Site,
                                   parasite$GenusSpecies,
                                   parasite$Parasite)
    return(parasite.comm)
}
