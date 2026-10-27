# 黄瓜杂种优势表型统计：3 个小区（随机区组）
# 输入：每行一个"小区 × 基因型"，共 15 行
#   plot      小区编号（G1/G2/G3）
#   genotype  C / D / S / MD / MS
#   weight_kg 该小区累计采收果实总质量（kg）
#   fruits    该小区累计采收果数
#   plants    该小区株数
# 用法：Rscript heterosis_plot_stats.R input.csv [out_prefix]
# 只依赖 base R；若装了 agricolae，会额外输出 Duncan 字母。

args <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[1] else stop("需要输入 CSV")
outpre <- if (length(args) >= 2) args[2] else "heterosis"

d <- read.csv(infile, stringsAsFactors = FALSE)
need <- c("plot", "genotype", "weight_kg", "fruits", "plants")
stopifnot(all(need %in% names(d)))
d$plot <- factor(d$plot)
d$genotype <- factor(d$genotype, levels = c("C", "D", "S", "MD", "MS"))
stopifnot(!anyNA(d$genotype), nrow(d) == nlevels(d$plot) * 5,
          all(table(d$plot, d$genotype) == 1))   # 完全区组、无缺失

# 小区水平性状
d$yield_g   <- 1000 * d$weight_kg / d$plants   # 单株产量 g
d$fruit_no  <- d$fruits / d$plants             # 单株果数
d$fruit_g   <- 1000 * d$weight_kg / d$fruits   # 单果重 g
traits <- c(yield_g = "Yield_per_plant_g", fruit_no = "Fruits_per_plant", fruit_g = "Fruit_mass_g")

crosses <- list(MD = c("C", "D"), MS = c("C", "S"))
r <- nlevels(d$plot)

# ---- LSD 紧凑字母（等重复时按均值排序的区间法）----
lsd_letters <- function(means, lsd) {
  o <- order(means, decreasing = TRUE); m <- means[o]; k <- length(m)
  ivs <- lapply(seq_len(k), function(i) {
    j <- i; while (j < k && abs(m[i] - m[j + 1]) < lsd) j <- j + 1; c(i, j)
  })
  keep <- !vapply(seq_along(ivs), function(a) any(vapply(seq_along(ivs), function(b)
    b != a && ivs[[b]][1] <= ivs[[a]][1] && ivs[[b]][2] >= ivs[[a]][2] &&
    !(identical(ivs[[b]], ivs[[a]]) && b > a), logical(1))), logical(1))
  ivs <- ivs[keep]
  lab <- character(k)
  for (g in seq_along(ivs)) {
    idx <- ivs[[g]][1]:ivs[[g]][2]
    lab[idx] <- paste0(lab[idx], letters[g])
  }
  setNames(lab, names(m))[names(means)]
}

anova_out <- list(); means_out <- list(); het_out <- list()

for (tr in names(traits)) {
  fit <- aov(as.formula(paste(tr, "~ plot + genotype")), data = d)
  a <- summary(fit)[[1]]
  anova_out[[tr]] <- data.frame(trait = traits[tr], term = trimws(rownames(a)),
                                df = a$Df, SS = a$`Sum Sq`, MS = a$`Mean Sq`,
                                F = a$`F value`, P = a$`Pr(>F)`, row.names = NULL)
  mse <- a["Residuals", "Mean Sq"]; dfe <- a["Residuals", "Df"]

  mu <- tapply(d[[tr]], d$genotype, mean)
  se_g <- sqrt(mse / r)                         # 基于误差均方的标准误
  sd_raw <- tapply(d[[tr]], d$genotype, sd)
  lsd <- qt(0.975, dfe) * sqrt(2 * mse / r)
  let <- lsd_letters(mu, lsd)
  dun <- NA
  if (requireNamespace("agricolae", quietly = TRUE)) {
    dt <- agricolae::duncan.test(fit, "genotype", console = FALSE)$groups
    dun <- dt[names(mu), "groups"]
  }
  means_out[[tr]] <- data.frame(trait = traits[tr], genotype = names(mu),
                                mean = mu, SE_pooled = se_g, SD_raw = sd_raw,
                                LSD_0.05 = let, Duncan_0.05 = dun, row.names = NULL)

  for (h in names(crosses)) {
    p <- crosses[[h]]; f1 <- mu[h]
    mp <- mean(mu[p]); hp_name <- p[which.max(mu[p])]; hp <- mu[hp_name]
    t_mp <- (f1 - mp) / sqrt(1.5 * mse / r)    # Var(F1 - MP) = 1.5σ²/r
    t_hp <- (f1 - hp) / sqrt(2 * mse / r)      # Var(F1 - HP) = 2σ²/r
    het_out[[paste(tr, h)]] <- data.frame(
      trait = traits[tr], hybrid = h, F1 = f1, MP = mp, HP = hp, HP_parent = hp_name,
      MPH_pct = 100 * (f1 / mp - 1), P_MPH = 2 * pt(-abs(t_mp), dfe),
      BPH_pct = 100 * (f1 / hp - 1), P_BPH = 2 * pt(-abs(t_hp), dfe),
      # 显性度 potence ratio: (F1-MP)/(|P1-P2|/2)；>1 超亲，0 加性，±1 完全显性
      # 双亲差异未达 LSD 时显性度无意义，记 NA
      potence = if (abs(diff(mu[p])) >= lsd) (f1 - mp) / (abs(diff(mu[p])) / 2) else NA,
      df_error = dfe, row.names = NULL)
  }
}

A <- do.call(rbind, anova_out); M <- do.call(rbind, means_out); H <- do.call(rbind, het_out)
rownames(A) <- rownames(M) <- rownames(H) <- NULL
write.csv(d, paste0(outpre, "_plot_level.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(A, paste0(outpre, "_anova.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(M, paste0(outpre, "_means_letters.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(H, paste0(outpre, "_heterosis_tests.csv"), row.names = FALSE, fileEncoding = "UTF-8")

op <- options(digits = 4)
cat("\n== 方差分析（plot + genotype）==\n"); print(A)
cat("\n== 均值与多重比较 ==\n"); print(M)
cat("\n== 杂种优势检验（t 检验，误差均方来自 ANOVA）==\n"); print(H)
options(op)
