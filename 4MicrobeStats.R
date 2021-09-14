## setwd('/Volumes/bombus/Dropbox (University of Oregon)/beeMicrobes')

rm(list=ls())
source("src/initialize.R")
source("src/makeIndivComm.R")
source("src/commPrep.R")

## ***************************************************************
## 16s
##  **************************************************************
microbes <- colnames(spec)[grepl("16s",
                                      colnames(spec))]
comm.microbes.indiv <- makeIndivComm(spec, microbes)

sum.indivs.16s <- apply(comm.microbes.indiv, 2, function(x){
    mean(x > 0)
    })

## number of individuals with each bacteria
sum.indivs.16s <- data.frame(prop=sum.indivs.16s[order(sum.indivs.16s,
                                                  decreasing = TRUE)])

sum.indivs.16s$ASV <- rownames(sum.indivs.16s)
rownames(sum.indivs.16s) <- NULL
sum.indivs.16s$ASV <- factor(sum.indivs.16s$ASV, levels=sum.indivs.16s$ASV)

sum.indivs.16s[1:20,] %>%
  ggplot( aes(x=ASV, y=prop)) +
    geom_bar(stat="identity", fill="darkolivegreen", alpha=.6, width=.4) +
    coord_flip() +
    xlab("") +
    theme_bw() +
    ylim(0,1) + labs(y="Proportion of specimens")


ggsave("figures/top20_16s_indiv_prop.pdf",
       height=4, width=15)


## average proportion of reads
prop.16s <- apply(comm.microbes.indiv, 2, function(x){
    mean(x[x > 0])
    })


prop.16s <- data.frame(prop=prop.16s[order(prop.16s,
                                                  decreasing = TRUE)])

prop.16s$ASV <- rownames(prop.16s)
rownames(prop.16s) <- NULL
prop.16s$ASV <- factor(prop.16s$ASV, levels=prop.16s$ASV)

prop.16s[1:20,] %>%
  ggplot( aes(x=ASV, y=prop)) +
    geom_bar(stat="identity", fill="darkolivegreen", alpha=.6, width=.4) +
    coord_flip() +
    xlab("") +
    theme_bw() +
    ylim(0,1) + labs(y="Mean proportion of reads within individuals")


ggsave("figures/top20_16s_prop_reads.pdf",
       height=4, width=15)

write.csv(sum.indivs.16s, file="saved/16sStats_prop_indiv.csv")
write.csv(prop.16s, file="saved/16sStats_propreads.csv")

##  ****************************************************************
## RBCL
##  ****************************************************************

rbcl <- colnames(spec)[grepl("RBCL", colnames(spec))]
comm.rbcl.indiv <- makeIndivComm(spec, rbcl)

sum.indivs.rbcl <- apply(comm.rbcl.indiv, 2, function(x){
    mean(x > 0)
    })

## number of individuals with each rbcl
sum.indivs.rbcl <- data.frame(prop=sum.indivs.rbcl[order(sum.indivs.rbcl,
                                                  decreasing = TRUE)])

sum.indivs.rbcl$ASV <- rownames(sum.indivs.rbcl)
rownames(sum.indivs.rbcl) <- NULL
sum.indivs.rbcl$ASV <- factor(sum.indivs.rbcl$ASV, levels=sum.indivs.rbcl$ASV)

sum.indivs.rbcl[1:20,] %>%
  ggplot( aes(x=ASV, y=prop)) +
    geom_bar(stat="identity", fill="darkolivegreen", alpha=.6, width=.4) +
    coord_flip() +
    xlab("") +
    theme_bw() +
    ylim(0,1) + labs(y="Proportion of specimens")


ggsave("figures/top20_rbcl_indiv_prop.pdf",
       height=4, width=15)


## average proportion of reads
prop.rbcl <- apply(comm.rbcl.indiv, 2, function(x){
    mean(x[x > 0])
    })


prop.rbcl <- data.frame(prop=prop.rbcl[order(prop.rbcl,
                                                  decreasing = TRUE)])

prop.rbcl$ASV <- rownames(prop.rbcl)
rownames(prop.rbcl) <- NULL
prop.rbcl$ASV <- factor(prop.rbcl$ASV, levels=prop.rbcl$ASV)

prop.rbcl[1:20,] %>%
  ggplot( aes(x=ASV, y=prop)) +
    geom_bar(stat="identity", fill="darkolivegreen", alpha=.6, width=.4) +
    coord_flip() +
    xlab("") +
    theme_bw() +
    ylim(0,1) + labs(y="Mean proportion of reads within individuals")


ggsave("figures/top20_rbcl_prop_reads.pdf",
       height=4, width=15)


write.csv(sum.indivs.rbcl, file="saved/rbclStats_prop_indiv.csv")
write.csv(prop.rbcl, file="saved/rbclStats_propreads.csv")


##  ****************************************************************
## parasites
##  ****************************************************************
## include Apidae as a dummy species to avoid the issue of having
## individuals dropped if they did not have any parasites.

parasite.comm <- spec[, c("UniqueID", "Apidae", parasites)]
parasite.comm <- parasite.comm[parasite.comm$Apidae == 1,]
parasite.comm[, c("Apidae", parasites)] <-
    apply(parasite.comm[, c("Apidae", parasites)], 2, as.numeric)
rownames(parasite.comm) <- parasite.comm$UniqueID
parasite.comm$UniqueID <- NULL
parasite.comm <- as.matrix(parasite.comm)


sum.indivs.parasite <- apply(parasite.comm, 2, function(x){
    mean(x > 0)
    })

## number of individuals with each parasite
sum.indivs.parasite <- data.frame(prop=sum.indivs.parasite[order(sum.indivs.parasite,
                                                  decreasing = TRUE)])

sum.indivs.parasite$parasite <- rownames(sum.indivs.parasite)
rownames(sum.indivs.parasite) <- NULL
sum.indivs.parasite$parasite <- factor(sum.indivs.parasite$parasite,
                                       levels=sum.indivs.parasite$parasite)

sum.indivs.parasite[-1, ] %>%
  ggplot( aes(x=parasite, y=prop)) +
    geom_bar(stat="identity", fill="darkolivegreen", alpha=.6, width=.4) +
    coord_flip() +
    xlab("") +
    theme_bw() +
    ylim(0,1) + labs(y="Proportion of specimens")


ggsave("figures/parasite_indiv_prop.pdf",
       height=3, width=7)

write.csv(sum.indivs.parasite, file="saved/parasiteStats_prop_indiv.csv")
