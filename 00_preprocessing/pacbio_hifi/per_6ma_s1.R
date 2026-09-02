library(tidyverse)
library(data.table)

adenines <- read_tsv("ref.adenines.tsv",
                     col_names = c("contig", "length", "n_adenines"),
                     show_col_types = FALSE)

sites <- fread(
    cmd    = "zcat sample1.6mA.bed.gz | grep -v '^#'",
    sep    = "\t", header = FALSE,
    select = c(1, 9, 10)
  ) |>
  as_tibble() |>
  set_names(c("contig", "cov", "mod_count"))


N_A      <- sum(adenines$n_adenines)
mean_cov <- mean(sites$cov)
p_fp     <- 1 - (1 - nrow(sites) / N_A)^(1 / mean_cov)

alpha <- 0.05 / N_A                     
sites <- sites |>
  mutate(pval = pbinom(mod_count - 1, cov, p_fp, lower.tail = FALSE))

by_region <- sites |>
  filter(cov >= 10, pval < alpha) |>
  count(contig, name = "n_meth") |>
  right_join(adenines, by = "contig") |>
  mutate(
    n_meth  = replace_na(n_meth, 0),
    pct_6mA = 100 * n_meth / n_adenines
  ) |>
  arrange(desc(pct_6mA))

cat(sprintf("per-read FP rate: %.5f   mean cov: %.1f   alpha: %.2e\n",
            p_fp, mean_cov, alpha))
write_tsv(by_region, "sample1.6mA.pct_by_region.tsv")
print(by_region, n = 40)
