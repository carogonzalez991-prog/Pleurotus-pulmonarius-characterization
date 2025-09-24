# placeholder
dirs <- c("figures","outputs","data","scripts")
invisible(lapply(dirs, function(z) if(!dir.exists(z)) dir.create(z)))
set.seed(123)

# ==========================================================
# Pipeline completo — EB (INICIAL vs CONSUMO)
# - Datos simulados (no expone datos reales)
# - Random Forest (tuning CV10, métricas, bootstrap B=500, heatmap)
# - Árbol de Decisión (poda ~8 hojas, CV10 anidada OOF, reglas/varimp, plot inferno invertida)
# - Bland–Altman con NUEVOS DATOS (SOLO modelos INICIALES: RF_ini y DT_ini)
# Salidas:
#   figures/: PNG (heatmaps, árboles, bland-altman)
#   outputs/: CSV/TXT (métricas, importancias, reglas, predicciones, resumen BA)
# ==========================================================

suppressPackageStartupMessages({
  pkgs <- c("ranger","caret","dplyr","tibble","ggplot2","viridisLite",
            "readr","tidyr","rpart","rpart.plot","scales","broom")
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  lapply(pkgs, library, character.only = TRUE)
})

# -----------------------------
# 0) Carpetas (por si no existen)
# -----------------------------
dirs <- c("data","scripts","figures","outputs")
invisible(lapply(dirs, function(z) if (!dir.exists(z)) dir.create(z)))
set.seed(123)

# -----------------------------
# 1) Datos SIMULADOS (misma estructura)
#    Reemplazá estos df si querés usar datos reales.
# -----------------------------
n <- 21
df_ini <- tibble(
  EB = runif(n, 20, 55),
  CA = runif(n, 0.90, 1.80),
  CC = runif(n, 17.0, 45.0),
  L  = runif(n,  3.0, 13.0),
  H  = runif(n,  6.0, 12.5)
)

df_cons <- tibble(
  EB = runif(n, 20, 55),
  CA = runif(n, 0.00, 0.80),
  CC = runif(n,  4.0, 28.5),
  L  = runif(n,  0.0,  4.0),
  H  = runif(n,  0.0,  4.0)
)

# -----------------------------
# 2) Utilitarias
# -----------------------------
R2   <- function(y,yhat) { if (sd(y)==0) NA_real_ else cor(y, yhat)^2 }
RMSE <- function(y,yhat) sqrt(mean((y - yhat)^2))
MAE  <- function(y,yhat) mean(abs(y - yhat))

tg_grid <- expand.grid(
  mtry = 2:4,
  splitrule = "variance",
  min.node.size = c(2,3,4,5)
)

# -----------------------------
# 3) Random Forest (fit + métricas + CV10)
# -----------------------------
fit_rf <- function(df, target = "EB", num.trees = 500){
  set.seed(123)
  idx_tr <- caret::createDataPartition(df[[target]], p = 0.7, list = FALSE)
  tr <- df[idx_tr, , drop = FALSE]
  te <- df[-idx_tr, , drop = FALSE]

  set.seed(123)
  rf_cv <- caret::train(
    stats::as.formula(paste(target, "~ .")),
    data = tr, method = "ranger",
    trControl = caret::trainControl(method="cv", number=10),
    tuneGrid = tg_grid, importance = "permutation",
    num.trees = num.trees, metric = "RMSE"
  )
  best <- rf_cv$bestTune

  set.seed(123)
  rf_fit <- ranger::ranger(
    stats::as.formula(paste(target, "~ .")),
    data = tr,
    mtry = best$mtry,
    splitrule = as.character(best$splitrule),
    min.node.size = best$min.node.size,
    num.trees = num.trees, importance = "permutation", seed=123
  )

  pred_tr <- predict(rf_fit, data=tr)$predictions
  pred_te <- predict(rf_fit, data=te)$predictions

  metrics_split <- dplyr::bind_rows(
    tibble(Set="Train", R2=R2(tr[[target]],pred_tr), RMSE=RMSE(tr[[target]],pred_tr), MAE=MAE(tr[[target]],pred_tr)),
    tibble(Set="Test",  R2=R2(te[[target]],pred_te), RMSE=RMSE(te[[target]],pred_te), MAE=MAE(te[[target]],pred_te))
  )

  set.seed(123)
  rf_cv_all <- caret::train(
    stats::as.formula(paste(target, "~ .")),
    data = df, method="ranger",
    trControl=caret::trainControl(method="cv", number=10),
    tuneGrid=tg_grid, importance="permutation",
    num.trees=num.trees, metric="RMSE"
  )

  best_row <- rf_cv_all$results |>
    dplyr::filter(mtry==rf_cv_all$bestTune$mtry,
                  splitrule==rf_cv_all$bestTune$splitrule,
                  min.node.size==rf_cv_all$bestTune$min.node.size) |>
    dplyr::arrange(RMSE) |>
    dplyr::slice(1)

  metrics_cv10 <- tibble(R2=best_row$Rsquared, RMSE=best_row$RMSE, MAE=best_row$MAE)

  list(model=rf_fit, best=best, metrics_split=metrics_split, metrics_cv10=metrics_cv10)
}

# -----------------------------
# 4) Bootstrap de Importancias (B=500)
# -----------------------------
stab_importance_rf <- function(df, target="EB", B=500, num.trees=500, seed=2024){
  set.seed(seed)
  cv_all <- caret::train(
    stats::as.formula(paste(target, "~ .")),
    data=df, method="ranger",
    trControl=caret::trainControl(method="cv", number=10),
    tuneGrid=tg_grid, importance="permutation",
    num.trees=num.trees, metric="RMSE"
  )
  bt <- cv_all$bestTune

  vars <- setdiff(names(df), target)
  IMP  <- matrix(NA_real_, nrow=B, ncol=length(vars), dimnames=list(NULL, vars))
  RANK <- matrix(NA_integer_, nrow=B, ncol=length(vars), dimnames=list(NULL, vars))

  for (b in 1:B) {
    idx <- sample(seq_len(nrow(df)), replace = TRUE)
    fit_b <- ranger::ranger(
      stats::as.formula(paste(target, "~ .")),
      data=df[idx, ],
      mtry=bt$mtry, splitrule=as.character(bt$splitrule),
      min.node.size=bt$min.node.size,
      num.trees=num.trees, importance="permutation",
      seed=seed+b
    )
    impb <- ranger::importance(fit_b)
    IMP[b, names(impb)] <- impb
    RANK[b, ] <- rank(-IMP[b, ], ties.method="min")
  }

  as_tibble(IMP) |>
    tidyr::pivot_longer(everything(), names_to="Variable", values_to="Imp") |>
    dplyr::group_by(Variable) |>
    dplyr::summarise(Mean=mean(Imp,na.rm=TRUE), SD=sd(Imp,na.rm=TRUE), .groups="drop") |>
    dplyr::left_join(
      as_tibble(RANK) |>
        tidyr::pivot_longer(everything(), names_to="Variable", values_to="rk") |>
        dplyr::group_by(Variable) |>
        dplyr::summarise(Top1_pct=mean(rk==1,na.rm=TRUE)*100,
                         Top2_pct=mean(rk<=2,na.rm=TRUE)*100, .groups="drop"),
      by="Variable"
    ) |>
    dplyr::mutate(CI95_low=Mean-1.96*SD, CI95_high=Mean+1.96*SD) |>
    dplyr::arrange(desc(Top1_pct), desc(Mean))
}

# -----------------------------
# 5) Heatmap 2D (Top-2 por estabilidad)
# -----------------------------
heatmap_from_stability <- function(df, target="EB", top2_vars, titulo, ngrid=180){
  stopifnot(length(top2_vars)==2)
  set.seed(123)
  rf_all <- caret::train(
    stats::as.formula(paste(target, "~ .")),
    data=df, method="ranger",
    trControl=caret::trainControl(method="cv", number=10),
    tuneGrid=tg_grid, importance="permutation",
    num.trees=500, metric="RMSE"
  )
  best <- rf_all$bestTune

  fit <- ranger::ranger(
    stats::as.formula(paste(target, "~ .")),
    data=df,
    mtry=best$mtry, splitrule=as.character(best$splitrule),
    min.node.size=best$min.node.size,
    num.trees=500, seed=999
  )

  v1 <- top2_vars[1]; v2 <- top2_vars[2]
  rest <- setdiff(names(df)[names(df)!=target], c(v1, v2))

  grid <- expand.grid(
    x1 = seq(min(df[[v1]]), max(df[[v1]]), length.out=ngrid),
    x2 = seq(min(df[[v2]]), max(df[[v2]]), length.out=ngrid)
  )
  names(grid) <- c(v1, v2)

  add0 <- as.data.frame(matrix(0, nrow=nrow(grid), ncol=length(rest)))
  names(add0) <- rest
  newdat <- cbind(grid, add0)
  names(newdat) <- c(v1, v2, rest)

  newdat$Ypred <- predict(fit, data=newdat)$predictions

  p <- ggplot(newdat, aes(x=.data[[v1]], y=.data[[v2]], fill=Ypred)) +
    geom_tile() +
    stat_contour(aes(z=Ypred), color="white", linewidth=0.35, bins=12, show.legend=FALSE) +
    scale_fill_viridis_c(option="inferno", direction=-1, name=paste0(target,"\npred.")) +
    labs(title=titulo,
         subtitle=paste0("Top-2 (B=500): ", v1, " & ", v2, " | Otras=0"),
         x=v1, y=v2) +
    theme_minimal(base_size=12) +
    theme(plot.title=element_text(face="bold"))
  p
}

# -----------------------------
# 6) Árbol de Decisión (poda ~8 hojas, CV10 OOF)
# -----------------------------
fit_pruned_tree <- function(df, target="EB", leaves_target=8L){
  set.seed(123)
  idx_te <- caret::createDataPartition(df[[target]], p=0.30, list=FALSE)
  te <- df[idx_te, , drop = FALSE]
  tr <- df[-idx_te, , drop = FALSE]

  ctrl_big <- rpart.control(cp=1e-4, minsplit=3, minbucket=2, maxdepth=30, xval=0)
  fit_big  <- rpart::rpart(stats::as.formula(paste(target, "~ .")), data=tr, method="anova", control=ctrl_big)

  tab <- fit_big$cptable
  if (is.null(tab) || nrow(tab)==0) {
    fit_8 <- fit_big
  } else {
    i <- which(tab[, "nsplit"] == (leaves_target - 1L))
    if (!length(i)) i <- which.min(abs(tab[, "nsplit"] - (leaves_target - 1L)))
    fit_8 <- prune(fit_big, cp=tab[i, "CP"])
  }

  pred_tr <- predict(fit_8, newdata=tr)
  pred_te <- predict(fit_8, newdata=te)
  metrics_split <- dplyr::bind_rows(
    tibble(Set="Train", R2=R2(tr[[target]],pred_tr), RMSE=RMSE(tr[[target]],pred_tr), MAE=MAE(tr[[target]],pred_tr)),
    tibble(Set="Test",  R2=R2(te[[target]],pred_te), RMSE=RMSE(te[[target]],pred_te), MAE=MAE(te[[target]],pred_te))
  )

  # CV10 anidada con OOF
  set.seed(123)
  folds <- caret::createFolds(df[[target]], k=10, returnTrain=FALSE)
  oof <- rep(NA_real_, nrow(df)); fold_rows <- list()
  for (k in seq_along(folds)) {
    idx_val <- folds[[k]]
    dtr <- df[-idx_val, , drop = FALSE]
    dva <- df[ idx_val, , drop = FALSE]
    fb <- rpart::rpart(stats::as.formula(paste(target, "~ .")), data=dtr, method="anova", control=ctrl_big)
    cpt <- fb$cptable
    if (is.null(cpt) || nrow(cpt)==0) {
      f8 <- fb
    } else {
      ii <- which(cpt[, "nsplit"] == (leaves_target - 1L))
      if (!length(ii)) ii <- which.min(abs(cpt[, "nsplit"] - (leaves_target - 1L)))
      f8 <- prune(fb, cp=cpt[ii, "CP"])
    }
    p_tr <- predict(f8, newdata=dtr)
    p_va <- predict(f8, newdata=dva)
    oof[idx_val] <- p_va
    fold_rows[[k]] <- tibble(
      fold=k,
      leaves=sum(f8$frame$var=="<leaf>"),
      R2_train=R2(dtr[[target]],p_tr),
      RMSE_train=RMSE(dtr[[target]],p_tr),
      MAE_train=MAE(dtr[[target]],p_tr),
      R2_val=R2(dva[[target]],p_va),
      RMSE_val=RMSE(dva[[target]],p_va),
      MAE_val=MAE(dva[[target]],p_va)
    )
  }
  cv10_folds   <- dplyr::bind_rows(fold_rows)
  cv10_summary <- cv10_folds %>% dplyr::summarise(
    folds = dplyr::n(),
    leaves_mean = mean(leaves), leaves_sd = sd(leaves),
    R2_val_mean = mean(R2_val, na.rm=TRUE),
    RMSE_val_mean = mean(RMSE_val, na.rm=TRUE),
    MAE_val_mean = mean(MAE_val, na.rm=TRUE)
  )
  oof_tbl <- tibble(Obs=df[[target]], Pred_OOF=oof)
  oof_metrics <- tibble(R2_OOF=R2(oof_tbl$Obs,oof_tbl$Pred_OOF),
                        RMSE_OOF=RMSE(oof_tbl$Obs,oof_tbl$Pred_OOF),
                        MAE_OOF=MAE(oof_tbl$Obs,oof_tbl$Pred_OOF))
  list(model=fit_8, metrics_split=metrics_split, cv10_folds=cv10_folds, cv10_summary=cv10_summary,
       oof_tbl=oof_tbl, oof_metrics=oof_metrics)
}

# -----------------------------
# 7) Plot Árbol con inferno invertida + contraste
# -----------------------------
plot_tree_inferno <- function(fit, file_png, title="Árbol de Decisión", dpi=200){
  pal <- viridisLite::inferno(100, direction=-1)
  yvals <- fit$frame$yval
  idx <- round(scales::rescale(yvals, to=c(1,100)))
  box_cols <- pal[idx]

  rgbm <- t(col2rgb(box_cols))/255
  lum  <- 0.2126*rgbm[,1] + 0.7152*rgbm[,2] + 0.0722*rgbm[,3]
  text_cols <- ifelse(lum < 0.45, "white", "black")

  png(file_png, width=1800, height=1200, res=dpi)
  rpart.plot::rpart.plot(
    fit, type=5, extra=101, under=TRUE,
    fallen.leaves=TRUE, branch.lty=1,
    shadow.col=rgb(0,0,0,0.15),
    box.col=box_cols, border.col="grey15",
    col=text_cols, roundint=FALSE, tweak=1.05,
    main=title
  )
  dev.off()
}

save_rules_varimp <- function(fit, prefix){
  rules_txt <- capture.output(rpart.plot::rpart.rules(fit, roundint=FALSE, cover=TRUE))
  writeLines(rules_txt, paste0("outputs/", prefix, "_rules.txt"))
  imp <- tibble::enframe(fit$variable.importance, name="Variable", value="Importance") %>% arrange(desc(Importance))
  readr::write_csv(imp, paste0("outputs/", prefix, "_varimp.csv"))
}

# -----------------------------
# 8) Correr RF y DT (INICIAL y CONSUMO)
# -----------------------------
rf_ini <- fit_rf(df_ini, target="EB")
rf_con <- fit_rf(df_cons, target="EB")
readr::write_csv(rf_ini$metrics_split, "outputs/RF_EB_INICIAL_metrics_split.csv")
readr::write_csv(rf_ini$metrics_cv10,  "outputs/RF_EB_INICIAL_metrics_cv10.csv")
readr::write_csv(rf_con$metrics_split, "outputs/RF_EB_CONSUMO_metrics_split.csv")
readr::write_csv(rf_con$metrics_cv10,  "outputs/RF_EB_CONSUMO_metrics_cv10.csv")

stab_ini <- stab_importance_rf(df_ini, target="EB", B=500)
stab_con <- stab_importance_rf(df_cons, target="EB", B=500)
readr::write_csv(stab_ini, "outputs/RF_EB_INICIAL_importance_bootstrap.csv")
readr::write_csv(stab_con, "outputs/RF_EB_CONSUMO_importance_bootstrap.csv")

top2_ini <- stab_ini$Variable[1:2]
top2_con <- stab_con$Variable[1:2]
p_hm_ini <- heatmap_from_stability(df_ini, target="EB", top2_vars=top2_ini, titulo="RF — EB (INICIAL)")
p_hm_con <- heatmap_from_stability(df_cons, target="EB", top2_vars=top2_con, titulo="RF — EB (CONSUMO)")
ggsave("figures/RF_EB_INICIAL_heatmap_top2.png", p_hm_ini, width=7.5, height=6, dpi=600)
ggsave("figures/RF_EB_CONSUMO_heatmap_top2.png", p_hm_con, width=7.5, height=6, dpi=600)

dt_ini <- fit_pruned_tree(df_ini, target="EB", leaves_target=8L)
dt_con <- fit_pruned_tree(df_cons, target="EB", leaves_target=8L)
readr::write_csv(dt_ini$metrics_split, "outputs/DT_EB_INICIAL_metrics_split.csv")
readr::write_csv(dt_ini$cv10_folds,    "outputs/DT_EB_INICIAL_cv10_folds.csv")
readr::write_csv(dt_ini$cv10_summary,  "outputs/DT_EB_INICIAL_cv10_summary.csv")
readr::write_csv(dt_ini$oof_tbl,       "outputs/DT_EB_INICIAL_oof_predictions.csv")
readr::write_csv(dt_ini$oof_metrics,   "outputs/DT_EB_INICIAL_oof_metrics.csv")

readr::write_csv(dt_con$metrics_split, "outputs/DT_EB_CONSUMO_metrics_split.csv")
readr::write_csv(dt_con$cv10_folds,    "outputs/DT_EB_CONSUMO_cv10_folds.csv")
readr::write_csv(dt_con$cv10_summary,  "outputs/DT_EB_CONSUMO_cv10_summary.csv")
readr::write_csv(dt_con$oof_tbl,       "outputs/DT_EB_CONSUMO_oof_predictions.csv")
readr::write_csv(dt_con$oof_metrics,   "outputs/DT_EB_CONSUMO_oof_metrics.csv")

plot_tree_inferno(dt_ini$model, "figures/DT_EB_INICIAL_8leaves.png", "Árbol de Decisión — EB (INICIAL, 8 hojas)")
plot_tree_inferno(dt_con$model, "figures/DT_EB_CONSUMO_8leaves.png", "Árbol de Decisión — EB (CONSUMO, 8 hojas)")
save_rules_varimp(dt_ini$model, "DT_EB_INICIAL")
save_rules_varimp(dt_con$model, "DT_EB_CONSUMO")

# -----------------------------
# 9) Bland–Altman con NUEVOS DATOS (SOLO modelos INICIALES)
#    Reemplazá este bloque "new_tbl" con tus datos observados reales.
# -----------------------------
new_tbl <- tibble(
  CA = runif(15, 0.90, 1.80),
  CC = runif(15, 17.0, 45.0),
  L  = runif(15,  3.0, 13.0),
  H  = runif(15,  6.0, 12.5),
  EB_obs = runif(15, 20, 55)
)

# Predicciones SOLO con modelos INICIALES
pred_rf_ini <- predict(rf_ini$model, data=new_tbl)$predictions
pred_dt_ini <- as.numeric(predict(dt_ini$model, new_tbl))

pred_tbl <- tibble(
  CA=new_tbl$CA, CC=new_tbl$CC, L=new_tbl$L, H=new_tbl$H,
  EB_obs=new_tbl$EB_obs,
  RF_ini=pred_rf_ini,
  DT_ini=pred_dt_ini
)
readr::write_csv(pred_tbl, "outputs/Predicciones_EB_nuevos_datos_INICIAL.csv")

# Función Bland–Altman
bland_altman <- function(obs, pred, model_name, file_png){
  diff <- pred - obs
  avg  <- (pred + obs) / 2
  bias <- mean(diff, na.rm=TRUE)
  sd_d <- sd(diff, na.rm=TRUE)
  loa_low  <- bias - 1.96*sd_d
  loa_high <- bias + 1.96*sd_d
  R2v   <- R2(obs, pred)
  RMSEv <- RMSE(obs, pred)
  MAEv  <- MAE(obs, pred)

  dfp <- tibble(avg=avg, diff=diff)
  p <- ggplot(dfp, aes(x=avg, y=diff)) +
    geom_point(alpha=0.8) +
    geom_hline(yintercept=bias, color="black", linetype="solid") +
    geom_hline(yintercept=loa_low,  color="red", linetype="dashed") +
    geom_hline(yintercept=loa_high, color="red", linetype="dashed") +
    labs(title=paste0("Bland–Altman — ", model_name),
         subtitle=sprintf("Bias=%.2f | LoA=[%.2f, %.2f] | R2=%.3f | RMSE=%.2f | MAE=%.2f",
                          bias, loa_low, loa_high, R2v, RMSEv, MAEv),
         x="Media (pred + obs)/2", y="Diferencia (pred − obs)") +
    theme_minimal(base_size=12)
  ggsave(file_png, p, width=7.5, height=5.5, dpi=600)

  tibble(Model=model_name, N=sum(!is.na(diff)), Bias=bias,
         LoA_low=loa_low, LoA_high=loa_high, R2=R2v, RMSE=RMSEv, MAE=MAEv)
}

ba_rf_ini <- bland_altman(new_tbl$EB_obs, pred_tbl$RF_ini, "RF — Inicial", "figures/BA_RF_inicial.png")
ba_dt_ini <- bland_altman(new_tbl$EB_obs, pred_tbl$DT_ini, "DT — Inicial", "figures/BA_DT_inicial.png")

ba_summary <- dplyr::bind_rows(ba_rf_ini, ba_dt_ini)
readr::write_csv(ba_summary, "outputs/BlandAltman_summary_nuevos_datos_INICIAL.csv")

# -----------------------------
# 10) LOG final
# -----------------------------
cat("\n== Archivos generados ==\n",
    "- outputs/RF_EB_INICIAL_metrics_split.csv | RF_EB_INICIAL_metrics_cv10.csv\n",
    "- outputs/RF_EB_CONSUMO_metrics_split.csv | RF_EB_CONSUMO_metrics_cv10.csv\n",
    "- outputs/RF_EB_INICIAL_importance_bootstrap.csv | RF_EB_CONSUMO_importance_bootstrap.csv\n",
    "- figures/RF_EB_INICIAL_heatmap_top2.png | RF_EB_CONSUMO_heatmap_top2.png\n",
    "- outputs/DT_EB_INICIAL_* (metrics_split/cv10_folds/cv10_summary/oof_*)\n",
    "- outputs/DT_EB_CONSUMO_* (metrics_split/cv10_folds/cv10_summary/oof_*)\n",
    "- figures/DT_EB_INICIAL_8leaves.png | DT_EB_CONSUMO_8leaves.png\n",
    "- outputs/Predicciones_EB_nuevos_datos_INICIAL.csv\n",
    "- figures/BA_RF_inicial.png | BA_DT_inicial.png\n",
    "- outputs/BlandAltman_summary_nuevos_datos_INICIAL.csv\n")
