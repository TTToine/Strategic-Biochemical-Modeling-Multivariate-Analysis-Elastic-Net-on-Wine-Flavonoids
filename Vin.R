library(gRbase)      
library(ggplot2)
library(gridExtra)
library(knitr)
library(broom)       
library(car)         
library(ggcorrplot)
library(qgraph)
library(dplyr)       
library(tidyr)       
library(reshape2)    
library(glmnet)      

# ------------------- RESET GRAFICO -------------------
graphics.off()

data(wine)
colnames(wine) <- c("Type", "Alcohol", "Malic", "Ash", "Alcalinity",
                    "Magnesium", "Phenols", "Flavanoids", "NonFlavanoids",
                    "Proanthocyanins", "Color", "Hue", "Dilution", "Proline")
wine$Type <- as.factor(wine$Type)

str(wine)

desc_stats <- wine %>%
  pivot_longer(cols = -Type, names_to = "Variable", values_to = "Value") %>%
  group_by(Variable, Type) %>%
  summarise(
    Mean = mean(Value),
    Median = median(Value),
    SD = sd(Value),
    IQR = IQR(Value),
    .groups = "drop"
  )

desc_wide <- desc_stats %>%
  pivot_wider(
    names_from = Type,
    values_from = c(Mean, Median, SD, IQR),
    names_glue = "{.value}_{Type}"
  )

kable(head(desc_wide, 6), digits = 3,
      caption = "Statistiche marginali a confronto tra i tre vitigni")
 
mq0 <- lm(Flavanoids ~ ., data = wine)

aic_0 <- step(mq0, direction = "both", k = 2, trace = 0)
bic_0 <- step(mq0, direction = "both", k = log(nrow(wine)), trace = 0)

mod_full_quad <- lm(Flavanoids ~ Type + Malic + I(Malic^2) + Color + I(Color^2) + 
                      Proline + I(Proline^2), data = wine)
mod_full_all  <- lm(Flavanoids ~ (.)^2, data = wine)

mod_quad_aic <- step(mod_full_quad, direction = "both", k = 2, trace = 0)
mod_quad_bic <- step(mod_full_quad, direction = "both", k = log(nrow(wine)), trace = 0)
mod_full_aic <- step(mod_full_all, direction = "both", k = 2, trace = 0)
mod_full_bic <- step(mod_full_all, direction = "both", k = log(nrow(wine)), trace = 0)

modelli <- list(
  "Completo"               = mq0,
  "AIC stepwise base"      = aic_0,
  "BIC stepwise base"      = bic_0,
  "Quadratici AIC"         = mod_quad_aic,
  "Quadratici BIC"         = mod_quad_bic,
  "Satura interazioni AIC" = mod_full_aic,
  "Satura interazioni BIC" = mod_full_bic
)

risultati <- data.frame(
  Modello   = names(modelli),
  AIC       = sapply(modelli, AIC),
  BIC       = sapply(modelli, BIC),
  adj_R2    = sapply(modelli, function(m) summary(m)$adj.r.squared),
  n_coef    = sapply(modelli, function(m) length(coef(m)))
)
risultati <- risultati[order(risultati$BIC), ]
rownames(risultati) <- NULL

kable(risultati, digits = 3, caption = "Confronto modelli (ordinati per BIC)")

x_vars <- model.matrix(Flavanoids ~ . - 1, data = wine)
y_var  <- wine$Flavanoids

set.seed(42)
cv_en <- cv.glmnet(x_vars, y_var, alpha = 0.5, nfolds = 10)

plot(cv_en, main = "Elastic Net: Devianza vs Log(Lambda)\n")

coef_en <- coef(cv_en, s = "lambda.min")
print(coef_en)

dir.create("grafici_wine", showWarnings = FALSE)

# 3.1 Scatter matrix + loess
numeric_vars <- setdiff(names(wine)[sapply(wine, is.numeric)], "Flavanoids")
plot_list <- list()
for (v in numeric_vars) {
  p <- ggplot(wine, aes_string(x = v, y = "Flavanoids", color = "Type")) +
    geom_point(alpha = 0.7, size = 1) +
    geom_smooth(method = "loess", se = FALSE, span = 0.9) +
    labs(title = v, x = NULL, y = "Flavanoids") +
    theme_minimal() +
    theme(legend.position = "none")
  plot_list[[v]] <- p
}
scatter_grid <- grid.arrange(grobs = plot_list, ncol = 3,
                             top = "Relazione tra Flavanoids e variabili chimiche (per cultivar)")
ggsave("grafici_wine/01_scatter_matrix.png", plot = scatter_grid,
       width = 14, height = 12, dpi = 300, bg = "white")

# 3.2 Focus su Flavanoids (box + density affiancati - come nel report)
p_box_flav <- ggplot(wine, aes(x = Type, y = Flavanoids, fill = Type)) +
  geom_boxplot() +
  labs(title = "Boxplot Flavanoids per cultivar", x = "Cultivar", y = "Flavanoids") +
  theme_minimal() +
  theme(legend.position = "none")

p_dens_flav <- ggplot(wine, aes(x = Flavanoids, fill = Type)) +
  geom_density(alpha = 0.6) +
  labs(title = "Densità di Flavanoids per cultivar", x = "Flavanoids", y = "Densità") +
  theme_minimal() +
  theme(legend.position = "right")

grid.arrange(p_box_flav, p_dens_flav, ncol = 2)
ggsave("grafici_wine/02_focus_flavanoids_grid.png", 
       plot = grid.arrange(p_box_flav, p_dens_flav, ncol = 2),
       width = 10, height = 5, dpi = 300, bg = "white")

# 3.3 Boxplot faceted per TUTTE le variabili
df_long <- wine %>%
  pivot_longer(cols = -Type, names_to = "variable", values_to = "value")

p_box_facet <- ggplot(df_long, aes(x = Type, y = value, fill = Type)) +
  geom_boxplot(outlier.size = 0.6) +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal() +
  labs(title = "Boxplot di tutte le variabili per Cultivar", x = "", y = "Valore") +
  theme(legend.position = "bottom")
print(p_box_facet)
ggsave("grafici_wine/13_boxplot_faceted_all_vars.png", plot = p_box_facet,
       width = 14, height = 10, dpi = 300, bg = "white")

# 3.4 Density plot faceted per TUTTE le variabili
p_density_facet <- ggplot(df_long, aes(x = value, fill = Type, color = Type)) +
  geom_density(alpha = 0.25) +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal() +
  labs(title = "Density plot di tutte le variabili per Cultivar", x = "Valore", y = "Densità") +
  theme(legend.position = "bottom")
print(p_density_facet)
ggsave("grafici_wine/14_density_faceted_all_vars.png", plot = p_density_facet,
       width = 14, height = 10, dpi = 300, bg = "white")

# 3.5 Barplot BIC
bic_bar <- ggplot(risultati, aes(x = reorder(Modello, BIC), y = BIC, fill = BIC)) +
  geom_col() + coord_flip() + theme_minimal() +
  labs(title = "Confronto BIC dei modelli", x = NULL, y = "BIC") +
  scale_fill_gradient(low = "steelblue", high = "tomato")
print(bic_bar)
ggsave("grafici_wine/04_bic_comparison.png", plot = bic_bar,
       width = 9, height = 7, dpi = 300, bg = "white")

# 3.6 Heatmap differenza correlazioni (Type 2 vs Type 1)
cor_type1 <- cor(wine[wine$Type == "v1", -1])
cor_type2 <- cor(wine[wine$Type == "v2", -1])
cor_diff  <- cor_type2 - cor_type1
cor_diff_m <- reshape2::melt(cor_diff)

p_diff_heat <- ggplot(cor_diff_m, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  labs(title = "Differenza correlazioni (Type 2 - Type 1)", x = "", y = "")
print(p_diff_heat)
ggsave("grafici_wine/15_correlation_difference_Type2_vs_Type1.png", plot = p_diff_heat,
       width = 10, height = 9, dpi = 300, bg = "white")

# 3.7 Scatter + retta per le 3 variabili più correlate
cor_with_flav <- cor(wine[, sapply(wine, is.numeric)])["Flavanoids", ]
top_vars <- names(sort(abs(cor_with_flav), decreasing = TRUE))[2:4]

for (v in top_vars) {
  p_lm <- ggplot(wine, aes_string(x = v, y = "Flavanoids", color = "Type")) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_minimal() +
    labs(title = paste(v, "vs Flavanoids (con retta di regressione)"),
         x = v, y = "Flavanoids")
  print(p_lm)
  ggsave(paste0("grafici_wine/16_scatter_lm_", v, ".png"), plot = p_lm,
         width = 8, height = 6, dpi = 300, bg = "white")
}

final_model <- bic_0

# VIF
vif_table <- car::vif(final_model)
vif_df <- as.data.frame(vif_table)
vif_df$Variable <- rownames(vif_df)
rownames(vif_df) <- NULL
names(vif_df)[1] <- "VIF"
if ("GVIF" %in% colnames(vif_df)) {
  vif_df$VIF <- vif_df$GVIF ^ (1 / (2 * vif_df$Df))
  vif_df <- vif_df[, c("Variable", "VIF")]
}

vif_bar <- ggplot(vif_df, aes(x = reorder(Variable, VIF), y = VIF)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red") +
  labs(title = "VIF del modello BIC selezionato", x = NULL, y = "VIF")
print(vif_bar)
ggsave("grafici_wine/05_vif_barplot.png", plot = vif_bar,
       width = 9, height = 7, dpi = 300, bg = "white")

# Coefficient plot
coef_plot <- tidy(final_model, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  ggplot(aes(x = estimate, y = reorder(term, estimate))) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Coefficienti del modello BIC con IC 95%", x = "Stima", y = NULL) +
  theme_minimal()
print(coef_plot)
ggsave("grafici_wine/06_coefficients.png", plot = coef_plot,
       width = 9, height = 7, dpi = 300, bg = "white")

# Fitted vs Observed
fitted_obs <- data.frame(Observed = wine$Flavanoids, Fitted = fitted(final_model))
fitted_obs_plot <- ggplot(fitted_obs, aes(x = Observed, y = Fitted)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Fitted vs Observed", x = "Flavanoids osservati", y = "Flavanoids previsti") +
  theme_minimal()
print(fitted_obs_plot)
ggsave("grafici_wine/07_fitted_vs_observed.png", plot = fitted_obs_plot,
       width = 8, height = 6, dpi = 300, bg = "white")

# Matrici correlazione e parziali
numeric_data <- wine[, sapply(wine, is.numeric)]
corr_matrix <- cor(numeric_data)
corr_plot <- ggcorrplot(corr_matrix, method = "square", type = "lower", lab = TRUE,
                        title = "Matrice di Correlazione", ggtheme = theme_minimal())
print(corr_plot)
ggsave("grafici_wine/08_correlation_matrix.png", plot = corr_plot,
       width = 10, height = 9, dpi = 300, bg = "white")

R_inv <- solve(corr_matrix)
pcor <- -cov2cor(R_inv)
diag(pcor) <- 0
pcor_thresh <- pcor
pcor_thresh[abs(pcor_thresh) < 0.15] <- 0

pcor_plot <- ggcorrplot(pcor, method = "square", type = "lower", lab = TRUE,
                        title = "Matrice Correlazioni Parziali", ggtheme = theme_minimal())
print(pcor_plot)
ggsave("grafici_wine/09_partial_correlation.png", plot = pcor_plot,
       width = 10, height = 9, dpi = 300, bg = "white")

# Grafo correlazioni parziali
png("grafici_wine/10_partial_correlation_graph.png", width = 1200, height = 1000, res = 150)
qgraph(pcor_thresh, layout = "spring", labels = colnames(pcor_thresh),
       edge.labels = TRUE, edge.color = "blue",
       title = "Grafo correlazioni parziali (|r| > 0.15)")
dev.off()

# Diagnostica residui
pdf("grafici_wine/11_diagnostic_plots.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))
plot(final_model, which = 1:4, col = "darkblue", pch = 16, cex = 0.8)
dev.off()

# Partial residual plots
pdf("grafici_wine/12_partial_residual_plots.pdf", width = 10, height = 8)
crPlots(final_model, terms = ~ Color + Phenols + Ash, 
        main = "Partial Residual Plots")
dev.off()

# Intervalli di confidenza
conf_intervals <- confint(final_model, level = 0.95)
kable(conf_intervals, digits = 3, 
      caption = "Intervalli di confidenza 95% - Modello BIC finale")