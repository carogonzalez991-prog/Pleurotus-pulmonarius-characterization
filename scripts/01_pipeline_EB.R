# placeholder
dirs <- c("figures","outputs","data","scripts")
invisible(lapply(dirs, function(z) if(!dir.exists(z)) dir.create(z)))
set.seed(123)

