runMantelSiteType <- function(site.in.type,
                              dis.bee, dist.plant,
                              dist.geo, dist.parasite){

    dist.bee.type <- dist.bee[rownames(dist.bee) %in% site.in.type,
                              colnames(dist.bee) %in% site.in.type]

    dist.plant.type <- dist.plant[rownames(dist.plant) %in% site.in.type,
                                  colnames(dist.plant) %in% site.in.type]

    dist.geo.type <- dist.geo[rownames(dist.geo) %in% site.in.type,
                              colnames(dist.geo) %in% site.in.type]

    dist.parasite.type <- dist.parasite[rownames(dist.parasite) %in% site.in.type,
                                        colnames(dist.parasite) %in% site.in.type]

    out <- MRM(as.dist(dist.parasite.type) ~ as.dist(dist.bee.type) +
                   as.dist(dist.plant.type) + as.dist(dist.geo.type),  nperm=10^4)
    return(out)

}



runMantelBeeSpecies <- function(ind.in.sp,
                                dist.microbe,
                                dist.plant,
                                dist.parasite){
    dist.microbe.sp <- dist.microbe[rownames(dist.microbe) %in% ind.in.sp,
                                    colnames(dist.microbe) %in% ind.in.sp]
    ## print(dim(dist.microbe.sp))
    if(!is.null(dim(dist.microbe.sp))){
        if(dim(dist.microbe.sp)[1] >= 5){
            ## print("here")

            dist.parasite.sp <- dist.parasite[rownames(dist.parasite) %in% ind.in.sp,
                                              colnames(dist.parasite) %in% ind.in.sp]

            dist.plant.sp <- dist.plant[rownames(dist.plant) %in% ind.in.sp,
                                        colnames(dist.plant) %in% ind.in.sp]

            out <- MRM(as.dist(dist.microbe.sp) ~ as.dist(dist.parasite.sp) +
                           as.dist(dist.plant.sp),  nperm=10^4)
        } else{
            out <- NA
        }
    } else{
        out <- NA
    }
    return(out)
}


runMantelBeeSpecies2 <- function(ind.in.sp,
                                dist.plant,
                                dist.microbe,
                                dist.parasite){
  dist.plant.sp <- dist.plant[rownames(dist.plant) %in% ind.in.sp,
                                  colnames(dist.plant) %in% ind.in.sp]
  ## print(dim(dist.microbe.sp))
  if(!is.null(dim(dist.plant.sp))){
    if(dim(dist.plant.sp)[1] >= 5){
      ## print("here")
      
      dist.parasite.sp <- dist.parasite[rownames(dist.parasite) %in% ind.in.sp,
                                        colnames(dist.parasite) %in% ind.in.sp]
      
      dist.microbe.sp <- dist.plant[rownames(dist.microbe) %in% ind.in.sp,
                                  colnames(dist.microbe) %in% ind.in.sp]
      
      out <- MRM(as.dist(dist.plant.sp) ~ as.dist(dist.parasite.sp) +
                   as.dist(dist.microbe.sp),  nperm=10^4)
    } else{
      out <- NA
    }
  } else{
    out <- NA
  }
  return(out)
}

