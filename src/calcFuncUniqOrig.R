library(FD)


## standardize a vector
standardize <- function(x)
  (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)


calcFuncUniqOrigComm <- function(traits, traits.2.keep,
                                 weights,
                                 type = "network",
                                 a,
                                 comm_id = "CommID",
                                 w.abun = TRUE,
                                 scale.out = TRUE,
                                 ...){

  # drop unwanted traits
  these.traits <- traits[, traits.2.keep, drop = FALSE]

  # ---- ALIGN traits and community matrix (this is the key fix) ----
  sp <- intersect(rownames(these.traits), colnames(a))
  if (length(sp) < 2) {
    stop("Need >= 2 shared species between traits rownames and community colnames.")
  }
  these.traits <- these.traits[sp, , drop = FALSE]
  a <- a[, sp, drop = FALSE]
  stopifnot(identical(rownames(these.traits), colnames(a)))
  # ---------------------------------------------------------------

  # convert traits with only 2 categories to binary numeric; others factors; numeric standardized
  for (this.trait in traits.2.keep) {
    out.traits <- unique(these.traits[, this.trait])
    out.traits <- out.traits[!is.na(out.traits)]

    if (!is.numeric(out.traits)) {
      if (length(out.traits) == 2) {
        message("cat == 2 ", this.trait)
        these.traits[, this.trait] <- as.numeric(as.factor(these.traits[, this.trait])) - 1
      } else {
        message("cat > 2 ", this.trait)
        these.traits[, this.trait] <- as.factor(these.traits[, this.trait])
      }
    } else {
      message("numeric ", this.trait)
      these.traits[, this.trait] <- standardize(these.traits[, this.trait])
    }
  }

  site.func.mets <- dbFD(these.traits,
                         a = a,
                         w = weights,
                         corr = "lingoes",
                         print.pco = TRUE,
                         w.abun = w.abun,
                         ...)

  coords <- site.func.mets$x.axes

  coords.comm <- vector(mode = "list", length = nrow(a))
  for (i in seq_len(nrow(a))) {
    # keep only present spp rows (and preserve rownames)
    present <- unlist(a[i, ]) > 0
    spp <- names(which(present))
    coords.comm[[i]] <- coords[spp, , drop = FALSE]
  }

  calcSiteLevelSpMets <- function(coords_sub){
    if (nrow(coords_sub) < 2) {
      out <- data.frame(uniq = NA_real_, originality = NA_real_)
      rownames(out) <- rownames(coords_sub)
      return(out)
    }

    centr <- apply(coords_sub, 2, mean)

    coords2 <- rbind(coords_sub, centr)
    rownames(coords2)[nrow(coords2)] <- "centr"

    dists_centr <- as.matrix(dist(coords2, diag = TRUE, upper = TRUE))
    diag(dists_centr) <- NA

    originality <- dists_centr[rownames(coords_sub), "centr"]

    # nearest neighbour among present species
    nn_mat <- dists_centr[rownames(coords_sub), rownames(coords_sub), drop = FALSE]
    uniq <- apply(nn_mat, 1, min, na.rm = TRUE)

    ## scale within the community unless scale.out = FALSE
    if (scale.out) {
      uniq <- scale(uniq)
      originality <- scale(originality)
    }
    out <- as.data.frame(cbind(uniq, originality))
    colnames(out) <- c("uniq", "originality")
    out
  }

  by.comm.mets <- lapply(coords.comm, calcSiteLevelSpMets)
  names(by.comm.mets) <- rownames(a)

  for (i in seq_along(by.comm.mets)) {
    this.name <- names(by.comm.mets)[i]
    this.comm <- by.comm.mets[[i]]

    this.comm[[comm_id]] <- this.name
    this.comm$GenusSpecies <- rownames(this.comm)
    rownames(this.comm) <- NULL

    by.comm.mets[[i]] <- this.comm
  }

  by.comm.mets <- do.call(rbind, by.comm.mets)
  rownames(by.comm.mets) <- NULL

  list(fd = site.func.mets,
       by.comm.mets = by.comm.mets)
}

make_comm_matrix <- function(df, id_col, sp_col = "GenusSpecies") {
  id_sym <- rlang::sym(id_col)
  sp_sym <- rlang::sym(sp_col)

  df %>%
    filter(!is.na(!!sp_sym), !!sp_sym != "") %>%
    count(!!id_sym, !!sp_sym, name = "N") %>%
    tidyr::pivot_wider(names_from = !!sp_sym, values_from = N, values_fill = 0) %>%
    tibble::column_to_rownames(id_col) %>%
    as.data.frame()
}

add_func_uniq_orig <- function(spec.net, traits,
                               id_col,
                               traits.2.keep, weights,
                               suffix,
                               add_fd = FALSE,
                               w.abun = TRUE,
                               ...) {

  # build comm matrix for this grouping
  comms <- make_comm_matrix(spec.net, id_col = id_col)

  # run FD + uniqueness/originality
  out <- calcFuncUniqOrigComm(traits,
                             traits.2.keep = traits.2.keep,
                             weights = weights,
                             a = comms,
                             comm_id = id_col,
                             w.abun = w.abun,
                             ...)

  # species-by-community uniq/orig
  by_sp <- out$by.comm.mets %>%
    dplyr::rename(
      !!paste0("uniq", suffix) := uniq,
      !!paste0("originality", suffix) := originality
    )

  spec.net <- spec.net %>%
    dplyr::left_join(by_sp, by = c("GenusSpecies", id_col))

  # optional: attach site/community-level FD metrics too
  if (add_fd) {
    fd <- tibble::tibble(
      !!id_col := names(out$fd$FDis),
      !!paste0("BeeFDis", suffix) := unname(out$fd$FDis),
      !!paste0("BeeFEve", suffix) := unname(out$fd$FEve)
    )
    spec.net <- spec.net %>%
      dplyr::left_join(fd, by = id_col)
  }

  spec.net
}


calcFuncUniqOrig <- function(traits, traits.2.keep,
                             weights,
                             type = "network",
                             scale.out = TRUE, ...){
  ## this function calculates trait uniqueness and originality based
  ## coux et al.
  
  these.traits <- traits[, colnames(traits) %in% traits.2.keep]
  
  ## convert traits with only 2 catagories to binary, numeric
  for(this.trait in  traits.2.keep){
    out.traits <- unique(these.traits[, this.trait])
    out.traits <- out.traits[!is.na(out.traits)]
    if(!is.numeric(out.traits)){
      if(length(out.traits) == 2){
        print(paste("cat == 2", this.trait))
        these.traits[, this.trait] <-
          as.numeric(as.factor(these.traits[,this.trait])) -1
      } else {
        print(paste("cat > 2", this.trait))
        these.traits[, this.trait] <-
          as.factor(these.traits[,this.trait])
      }
    } else {
      print(paste("numeric", this.trait))
      these.traits[, this.trait] <-
        standardize(these.traits[,this.trait])
    }
  }
  coords <- dbFD(these.traits, w=weights,
                 corr="lingoes", print.pco=T, ...)$x.axes
  centr <- apply(coords, 2, mean)
  
  ## add centroid coords as last row in dataframe
  coords2 <- data.frame(rbind(coords, centr))
  rownames(coords2)[dim(coords)[1]+1] <- "centr"
  
  ## create a matrix of distances between all species coordinates and
  ## centroid
  dists_centr <- as.matrix(dist(coords2, diag=TRUE, upper=TRUE))
  for (i in 1:dim(dists_centr)[1]) {
    dists_centr[i,i] <- NA
  }
  ## Originality: Distance to centroid of the species present
  originality <- dists_centr[, dim(dists_centr)[1]]
  originality <- originality[-(length(originality))]
  
  ## Uniqueness: nearest neighbour among present species
  ## the the minimum disance for each row
  uniq <- apply(dists_centr[-(dim(dists_centr)[1]),
                            -(dim(dists_centr)[1])],
                1, min, na.rm=TRUE)
  names(uniq) <- names(originality)
  ## scale unless scale.out = FALSE
  if(scale.out){
    uniq <- scale(uniq)
    originality <- scale(originality)
  }
  out <- cbind(uniq, originality)
  colnames(out) <- c("uniq", "originality")
  out <- as.data.frame(out)
  if(is.null(rownames(these.traits)) | type == "all spec"){
    out$GenusSpecies <- traits$GenusSpecies
  } else{
    out$GenusSpecies <- rownames(traits)
  }
  return(out)
}

