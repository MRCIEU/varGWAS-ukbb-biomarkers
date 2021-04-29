model <- function(data, out, chr, pos, oa, ea, rsid) {
    dosage = tryCatch({
      extract_variant_from_bgen(chr, pos, oa, ea)
  }, error = function(error_condition) {
      return(NA)
  })

  # prepare data
  data <- merge(data, dosage, "appieu")
  data <- na.omit(data)
  s <- paste0("chr", chr, "_", pos, "_", oa, "_", ea)
  s <- gsub(" ", "", s, fixed = TRUE)
  data$xsq <- (data %>% pull(!!s))^2
  names(data) <- gsub("-", "_", names(data))
  out <- gsub("-", "_", out)

  # first-stage model
  f <- paste0(out, " ~ ", s, " + sex.31.0.0 + age_at_recruitment.21022.0.0 +", paste0("PC", seq(1, 10), collapse="+"))
  fit1q = tryCatch({
      rq(as.formula(f), tau=0.5, data=data) # quantile
  }, error = function(error_condition) {
      return(NA)
  })
  fit1r <- tidy(lmrob(as.formula(f), data=data)) # SE robust

  # second-stage model
  data$d <- resid(fit1q)^2
  f <- paste0("d ~ ", s, " + xsq")
  fit2 <- lm(as.formula(f), data=data)

  # F-test
  fit0 <- lm(d ~ 1, data=data)
  ftest <- tidy(anova(fit0, fit2))
  fit2t <- tidy(fit2)

  return(data.frame(
      rsid=rsid,
      SNP=s,
      BETA_x=fit2t$estimate[2],
      BETA_xq=fit2t$estimate[3], 
      SE_x=fit2t$std.error[2],
      SE_xq=fit2t$std.error[3],
      Pvar=ftest$p.value[2],
      BETA=fit1r$estimate[2],
      SE=fit1r$std.error[2],
      Pmu=fit1r$p.value[2]
    )
  )
}