write.ms.table <- function(mod.output, mod.name){
    sum.mod <- as.data.frame(round(summary(mod.output)$fixed,2))
    samps.mod <- posterior_samples(mod.output,
                                          rownames(sum.mod))[,
                                           paste0("b_",
                                           rownames(sum.mod))]

    sum.mod$Pgt0  <- round(apply(samps.mod, 2, function(x)
        sum(x > 0)/length(x)), 2)

    sum.mod$Plt0  <- round(apply(samps.mod, 2, function(x)
        sum(x < 0)/length(x)),2)


    write.table(sum.mod,
                file=sprintf("saved/tables/%s.txt", mod.name),
                sep="&")

    write.csv(sum.mod,
              file=sprintf("saved/tables/%s.csv", mod.name))
}
