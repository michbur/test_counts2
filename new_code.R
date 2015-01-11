library(dpcR)

adpcr1 <- sim_adpcr(m = 10, n = 765, times = 1e5, pos_sums = FALSE, n_panels = 3)
adpcr2 <- sim_adpcr(m = 50, n = 765, times = 1e5, pos_sums = FALSE, n_panels = 3)
adpcr3 <- sim_adpcr(m = 100, n = 765, times = 1e5, pos_sums = FALSE, n_panels = 3)

adpcrs <- bind_dpcr(adpcr1, adpcr2, adpcr3)

badpcrs <- binarize(adpcrs)




# p_val <- p.adjust(rep(0.05, ncol(badpcrs)), "BH") 
# conf_ints <- dpcR:::fl(binom.confint(colSums(badpcrs), 
#                                      slot(badpcrs, "n"), 
#                                      conf.level = 1 - p_val[1],
#                                      "wilson")[, 4L:6])


#naive approach
positives <- colSums(badpcrs)
total <- slot(badpcrs, "n")
test_ids <- combn(1L:length(total), 2)
all_combns <- apply(test_ids, 2, function(i)
  prop.test(positives[i], total[i]))

p_vals <- p.adjust(sapply(all_combns, function(i)
  i[["p.value"]]), method = "BH")

#res <- t(rbind(test_ids, p_vals > 0.05))

only_signif <- test_ids[, p_vals > 0.05]
groups <- unique(lapply(1L:length(total), function(i)
  sort(unique(as.vector(only_signif[, as.logical(colSums(only_signif == i))])))))


#http://www.stat.ufl.edu/~aa/articles/agresti_bini_bertaccini_ryu_2008.pdf
#http://www.stat.ufl.edu/~aa/cda/R/multcomp/ryu-simultaneous.pdf
#http://cran.r-project.org/web/packages/MCPAN/ - disregard