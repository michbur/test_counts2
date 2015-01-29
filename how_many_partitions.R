library(dpcR)
library(pblapply)

dpcr_size = 20000L

dpcr_data <- function(n, pos)
  c(rep(1, pos), rep(0, n - pos))

sim_res <- t(do.call(cbind, pblapply(1L:8000, function(pos1) 
  sapply(1L:8000, function(pos2) {
    dpcr1 <- create_dpcr(dpcr_data(dpcr_size, pos1), dpcr_size, type = "np", adpcr = TRUE)
    dpcr2 <- create_dpcr(dpcr_data(dpcr_size, pos2), dpcr_size, type = "np", adpcr = TRUE)
    dpcrs <- bind_dpcr(dpcr1, dpcr2)
    
    compb <- test_counts(dpcrs, "binomial")
    compprop <- test_counts(dpcrs, "prop")
    compratio <- test_counts(dpcrs, "ratio")
    
    
    c(dpcr_size, pos1, pos2, sapply(list(compb, compprop, compratio), 
                                    function(single_test)
                                      slot(single_test, "test_res")[, "p_value"]))
  }))))

colnames(sim_res) <- c("total", "pos1", "pos2", "pb", "pp", "pr")
write.csv2(sim_res, file = "sim_res.csv")
