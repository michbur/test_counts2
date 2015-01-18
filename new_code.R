library(dpcR)

adpcr1 <- sim_adpcr(m = 10, n = 765, times = 1e5, pos_sums = FALSE, n_panels = 3)
adpcr2 <- sim_adpcr(m = 50, n = 765, times = 1e5, pos_sums = FALSE, n_panels = 3)
adpcr3 <- sim_adpcr(m = 100, n = 765, times = 1e5, pos_sums = FALSE, n_panels = 3)

adpcrs <- bind_dpcr(adpcr1, adpcr2, adpcr3)

badpcrs <- binarize(adpcrs)



test_example <- test_counts(badpcrs)

#naive approach
positives <- colSums(badpcrs)
total <- slot(badpcrs, "n")
test_ids <- combn(1L:length(total), 2)
all_combns <- apply(test_ids, 2, function(i)
  prop.test(positives[i], total[i]))

p_vals <- p.adjust(sapply(all_combns, function(i)
  i[["p.value"]]), method = "BH")

statistics <- sapply(all_combns, function(i)
  i[["statistic"]])
#t_res in original test_counts
x_res <- data.frame(X_squared = statistics, p_value = p_vals,
                    row.names = apply(matrix(names(positives)[test_ids], ncol = 2, byrow= TRUE),
                                      1, function(i) paste(i, collapse = " - ")))

#group_coef slot
only_signif <- test_ids[, p_vals > 0.05]
groups <- unique(lapply(1L:length(total), function(i)
  sort(unique(as.vector(only_signif[, as.logical(colSums(only_signif == i))])))))

group_matrix <- sapply(1L:length(total), function(experiment) 
  sapply(groups, function(single_group) experiment %in% single_group))

dimnames(group_matrix) <- list(letters[1L:length(groups)], names(positives))

group_coef<- data.frame(apply(group_matrix, 2, function(i) 
  paste(names(i[which(i)]), collapse = "")), 
  dpcR:::fl(binom.confint(colSums(badpcrs), 
                          slot(badpcrs, "n"), 
                          conf.level = 1 - p.adjust(rep(0.05, ncol(badpcrs)), "BH")[1],
                          "wilson")[, 4L:6]))
colnames(group_coef) <- c("group", "lambda", "lambda.low", "lambda.up")


#http://www.stat.ufl.edu/~aa/articles/agresti_bini_bertaccini_ryu_2008.pdf
#http://www.stat.ufl.edu/~aa/cda/R/multcomp/ryu-simultaneous.pdf
#http://cran.r-project.org/web/packages/MCPAN/ - disregard