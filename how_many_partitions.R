library(dpcR)
library(pbapply)

dpcr_size = 765

dpcr_data <- function(n, pos)
  c(rep(1, pos), rep(0, n - pos))

sim_res <- t(do.call(cbind, pblapply(1L:600, function(pos1) 
  sapply(1L:600, function(pos2) {
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

sim_res <- read.csv2("sim_res.csv")[, -1]
library(ggplot2)
ggplot(sim_res, aes(x = pos1, y = pos2, fill = pr)) + geom_tile()
