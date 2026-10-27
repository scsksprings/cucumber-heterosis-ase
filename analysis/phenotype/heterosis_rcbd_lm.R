# 黄瓜杂种优势表型统计：随机区组（允许个别小区缺失）
# 输入 CSV 列：plot, genotype, weight_kg, fruits, plants
#   weight_kg 为 NA 时，该小区的单株产量和单果重按缺失处理（果数保留）
# 用法：Rscript heterosis_rcbd_lm.R input.csv [out_prefix]
# 只依赖 base R

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]; outpre <- if (length(args) >= 2) args[2] else "rcbd"
d <- read.csv(infile, stringsAsFactors = FALSE)
geno <- c("C", "D", "S", "MD", "MS")
d$plot <- factor(d$plot); d$genotype <- factor(d$genotype, levels = geno)
d$yield_g  <- 1000 * d$weight_kg / d$plants
d$fruit_no <- d$fruits / d$plants
d$fruit_g  <- 1000 * d$weight_kg / d$fruits
traits <- c(yield_g = "Yield_per_plant_g", fruit_no = "Fruits_per_plant", fruit_g = "Fruit_mass_g")
crosses <- list(MD = c("C", "D"), MS = c("C", "S"))

holm <- function(p) p.adjust(p, "holm")

# 紧凑字母：基于两两显著性矩阵（插入-吸收法）
cld <- function(means, sig) {
  o <- order(means, decreasing = TRUE); g <- names(means)[o]
  groups <- list(g[1])
  for (x in g[-1]) {
    placed <- FALSE
    for (k in seq_along(groups)) if (all(!sig[x, groups[[k]]])) { groups[[k]] <- c(groups[[k]], x); placed <- TRUE }
    if (!placed) {
      ok <- g[seq_len(which(g == x))]; ok <- ok[!sig[x, ok]]
      groups[[length(groups) + 1]] <- ok
    }
  }
  # 去掉被其他组包含的冗余组
  keep <- !sapply(seq_along(groups), function(i) any(sapply(seq_along(groups), function(j)
    j != i && all(groups[[i]] %in% groups[[j]]) && (length(groups[[j]]) > length(groups[[i]]) || j < i))))
  groups <- groups[keep]
  lab <- setNames(rep("", length(g)), g)
  for (k in seq_along(groups)) lab[groups[[k]]] <- paste0(lab[groups[[k]]], letters[k])
  lab[names(means)]
}

A <- M <- H <- list()
for (tr in names(traits)) {
  dd <- d[!is.na(d[[tr]]), ]
  f <- lm(as.formula(paste(tr, "~ plot + genotype")), data = dd)
  a <- anova(f)   # 基因型在区组之后检验
  A[[tr]] <- data.frame(trait = traits[tr], term = rownames(a), df = a$Df, SS = a$`Sum Sq`,
                        MS = a$`Mean Sq`, F = a$`F value`, P = a$`Pr(>F)`, n_obs = nrow(dd))
  # 最小二乘均值：对区组取平均
  nd <- expand.grid(plot = levels(d$plot), genotype = geno)
  X <- model.matrix(~ plot + genotype, nd)
  L <- t(sapply(geno, function(gg) colMeans(X[nd$genotype == gg, , drop = FALSE])))
  b <- coef(f); V <- vcov(f); dfe <- df.residual(f)
  mu <- drop(L %*% b); se <- sqrt(diag(L %*% V %*% t(L)))
  sig <- matrix(FALSE, 5, 5, dimnames = list(geno, geno))
  for (i in geno) for (j in geno) if (i != j) {
    l <- L[i, ] - L[j, ]; tt <- sum(l * b) / sqrt(drop(t(l) %*% V %*% l))
    sig[i, j] <- 2 * pt(-abs(tt), dfe) < 0.05      # LSD（未校正两两 t 检验）
  }
  M[[tr]] <- data.frame(trait = traits[tr], genotype = geno, lsmean = mu, SE = se,
                        n = as.vector(table(dd$genotype)[geno]), LSD_letter = cld(setNames(mu, geno), sig))
  for (h in names(crosses)) {
    p <- crosses[[h]]; hp <- p[which.max(mu[p])]
    lm_ <- L[h, ] - (L[p[1], ] + L[p[2], ]) / 2; lb <- L[h, ] - L[hp, ]
    tm <- sum(lm_ * b) / sqrt(drop(t(lm_) %*% V %*% lm_)); tb <- sum(lb * b) / sqrt(drop(t(lb) %*% V %*% lb))
    mp <- mean(mu[p]); half <- abs(diff(mu[p])) / 2
    lpar <- L[p[1], ] - L[p[2], ]; tpar <- sum(lpar * b) / sqrt(drop(t(lpar) %*% V %*% lpar))
    H[[paste(tr, h)]] <- data.frame(trait = traits[tr], hybrid = h, F1 = mu[h], MP = mp, HP = mu[hp], HP_parent = hp,
      MPH_pct = 100 * (mu[h] / mp - 1), P_MPH = 2 * pt(-abs(tm), dfe),
      BPH_pct = 100 * (mu[h] / mu[hp] - 1), P_BPH = 2 * pt(-abs(tb), dfe),
      parents_differ_P = 2 * pt(-abs(tpar), dfe),
      potence = if (2 * pt(-abs(tpar), dfe) < 0.05) (mu[h] - mp) / half else NA, df_error = dfe)
  }
}
A <- do.call(rbind, A); M <- do.call(rbind, M); H <- do.call(rbind, H)
# Holm：每个性状内 2 个组合为一个检验家族；另给六项合并的结果供对照
H$P_MPH_Holm_trait <- ave(H$P_MPH, H$trait, FUN = holm)
H$P_BPH_Holm_trait <- ave(H$P_BPH, H$trait, FUN = holm)
H$P_MPH_Holm6 <- holm(H$P_MPH); H$P_BPH_Holm6 <- holm(H$P_BPH)

# 产量 = 果数 × 单果重：按两构成性状中亲值相乘得到的期望产量
mm <- function(tr, g) M$lsmean[M$trait == traits[tr] & M$genotype == g]
Mult <- do.call(rbind, lapply(names(crosses), function(h) {
  p <- crosses[[h]]
  data.frame(hybrid = h, yield_F1 = mm("yield_g", h), yield_MP = mean(sapply(p, mm, tr = "yield_g")),
             yield_expected_from_component_MPs = mean(sapply(p, mm, tr = "fruit_no")) * mean(sapply(p, mm, tr = "fruit_g")))
}))

rownames(A) <- rownames(M) <- rownames(H) <- NULL
for (x in c("A", "M", "H", "Mult")) write.csv(get(x), sprintf("%s_%s.csv", outpre, x), row.names = FALSE, fileEncoding = "UTF-8")
op <- options(digits = 4, width = 200)
print(A); print(M); print(H); print(Mult); options(op)
