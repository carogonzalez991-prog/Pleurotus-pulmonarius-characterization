# scripts/main.R — pipeline principal
req <- c("ranger","caret","dplyr","tibble","ggplot2","viridisLite",
         "readr","tidyr","rpart","rpart.plot","scales","broom")
to_install <- setdiff(req, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

# asegurar carpetas de salida
dirs <- c("data","scripts","figures","outputs")
invisible(lapply(dirs, function(z) if (!dir.exists(z)) dir.create(z)))

# correr tu pipeline (ajustá el nombre si usás otro)
source("scripts/01_pipeline_EB.R")

cat("\n✅ Listo. Ver figuras/ y outputs/\n")
