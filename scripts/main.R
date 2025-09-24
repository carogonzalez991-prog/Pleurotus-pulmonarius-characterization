# main.R — corre todo el pipeline de punta a punta

# 0) Paquetes (instala faltantes)
req <- c("ranger","caret","dplyr","tibble","ggplot2","viridisLite",
         "readr","tidyr","rpart","rpart.plot","scales")
to_install <- setdiff(req, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

# 1) Asegurar carpetas
dirs <- c("data","scripts","figures","outputs")
invisible(lapply(dirs, function(z) if(!dir.exists(z)) dir.create(z)))

# 2) Correr el pipeline (guarda PNG en figures/ y CSV/TXT en outputs/)
source("scripts/01_pipeline_EB.R")

cat("\n✅ Listo. Mirá resultados en:\n- figures/ (PNG)\n- outputs/ (CSV/TXT)\n")
