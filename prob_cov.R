library(dpcR)
library(binom)
library(reshape2)

coverage <- sapply(c(0.05, 0.1, 1:10/10*2), function(num_mol) {
  number_of_exps <- 1000
  dat <- sim_ddpcr_bkm(num_mol, n_exp = number_of_exps, type = "tnp")
  
  sum <- summary(dat, print = FALSE)[["summary"]]
  
  
  positives <- slot(dat, ".Data")
  total <- slot(dat, "n")
  
  sidak_ci <- dpcR:::fl(binom.confint(as.vector(positives), 
                                      total, 
                                      conf.level = (1 - 0.05)^(1/length(dat)),
                                      "wilson")[, 4L:6])[, 2L:3]
  
  nadj_ci <- dpcR:::fl(binom.confint(as.vector(positives), 
                                     total, 
                                     conf.level = 0.95,
                                     "wilson")[, 4L:6])[, 2L:3]
  
  dube_ci <- sum[sum[["method"]] == "dube", c("lambda.low", "lambda.up")]
  bhat_ci <- sum[sum[["method"]] == "bhat", c("lambda.low", "lambda.up")]
  c(sidak = sum(apply(cbind(rep(num_mol, number_of_exps), sidak_ci), 1, 
                      function(i)
                        abs(i[1]) > abs(i[2]) && i[1] < i[3]))/number_of_exps,
    nadj = sum(apply(cbind(rep(num_mol, number_of_exps), nadj_ci), 1, 
                     function(i)
                       abs(i[1]) > abs(i[2]) && i[1] < i[3]))/number_of_exps,
    dube = sum(apply(cbind(rep(num_mol, number_of_exps), dube_ci), 1, 
                     function(i)
                       abs(i[1]) > abs(i[2]) && i[1] < i[3]))/number_of_exps,
    bhat = sum(apply(cbind(rep(num_mol, number_of_exps), bhat_ci), 1, 
                     function(i)
                       abs(i[1]) > abs(i[2]) && i[1] < i[3]))/number_of_exps)  
})

m_coverage <- melt(coverage)
colnames(m_coverage) <- c("method", "prop", "coverage")
m_coverage[["prop"]] <- unlist(lapply(c(0.05, 0.1, 1:10/10*2), rep, 4))
m_coverage[["prop"]] <- as.factor(format(m_coverage[["prop"]], 1))
levels(m_coverage[["method"]]) <- c("Adjusted", "Non-adjusted", "Dube", "Bhat")

save(m_coverage, file = "coverage.RData")


#simultaneous probability coverage


simultaneous <- sapply(c(0.05, 0.1, 1:10/10*2), function(num_mol) {
  number_of_exps <- 1e4
  dat <- sim_ddpcr_bkm(num_mol, n_exp = number_of_exps, type = "tnp")
  positives <- slot(dat, ".Data")
  total <- slot(dat, "n")
  
  rowSums(sapply(1L:2000, function(dummy) 
    sapply(c(1 - (1 - 0.05)^(1/100), 1:6/20), function(alpha) {
      sidak_ci <- dpcR:::fl(binom.confint(as.vector(sample(positives, 100)), 
                                          rep(20000, 100), 
                                          conf.level = 1 - alpha,
                                          "wilson")[, 4L:6])[, 2L:3]
      
      sum(apply(cbind(rep(num_mol, 100), sidak_ci), 1, 
                function(i)
                  abs(i[1]) > abs(i[2]) && i[1] < i[3]))/100
    }))== 1)/2000
})

simultaneous <- simultaneous[1L:3, ]
rownames(simultaneous) <- c("Adj", "Dube", "Bhat")
m_simultaneous <- melt(simultaneous)
colnames(m_simultaneous) <- c("method", "prop", "coverage")
m_simultaneous[["prop"]] <- unlist(lapply(c(0.05, 0.1, 1:10/10*2), rep, 3))
m_simultaneous[["prop"]] <- as.factor(format(m_simultaneous[["prop"]], 1))
levels(m_simultaneous[["method"]]) <- c("Adjusted", "Dube", "Bhat")

save(m_simultaneous, file = "simultaneous.RData")


m_coverage <- m_coverage[m_coverage[["method"]] != "Non-adjusted", ]
m_coverage[["method"]] <- droplevels(m_coverage[["method"]])
m_coverage[["method"]] <- factor(m_coverage[["method"]],
                                 levels(m_coverage[["method"]])[c(1,3,2)])
m_coverage2 <- cbind(m_coverage, simul = m_simultaneous[, 3])
m_coverage2 <- melt(m_coverage2)
colnames(m_coverage2) <- c("method", "prop", "coverage", "value")

levels(m_coverage2[["coverage"]]) <- c("Average probability coverage", 
                                       "Simultaneous probability coverage")
save(m_coverage2, file = "coverage_poster.RData")
