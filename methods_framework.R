library(pbapply)
library(dpcR)

c(1L:5*3, 4L:10*5, 6L:10*10)

get_ints <- function(number_of_exps) {
  adpcr1 <- sim_adpcr(m = 10, n = 765, times = 10000, pos_sums = FALSE, 
                      n_panels = number_of_exps)
  adpcr2 <- sim_adpcr(m = 40, n = 765, times = 10000, pos_sums = FALSE, 
                      n_panels = number_of_exps)
  dat <- binarize(bind_dpcr(adpcr1, adpcr2))
  glmtab <- slot(test_counts(dat, model = "binomial"), "group_coef")
  ratiotab <- slot(test_counts(dat, model = "ratio"), "group_coef")
  dpcrtab <- summary(dat, print = FALSE)[["summary"]]
  list(ints = cbind(dpcrtab[seq(1, nrow(dpcrtab), by = 2), 3:5], 
                    dpcrtab[seq(1, nrow(dpcrtab), by = 2) + 1, 3:5],
                    glmtab[, -1],
                    ratiotab[, -1]),
       real_m = colSums(dat),
       test_obj = test_counts(dat))
}

ints5 <- get_ints(5)


colnames(ints5[["ints"]]) <- rep(c("lambda", "low", "up"), 4)
method_names <- c("Dube", "Bhat", "GLM", "Ratio")
indat <- do.call(rbind, lapply(0L:3, function(i) {
  data.frame(experiment = names(ints5[["real_m"]]),
             name = rep(method_names[i + 1], nrow(ints5[["ints"]])),
             real = ints5[["real_m"]],
             ints5[["ints"]][, 1L:3+3*i])
}))

save(indat, file = "ints_plot.RData")
