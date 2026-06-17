suppressPackageStartupMessages({
  library(MSdev)
  library(dplyr)
})
Sys.setenv(TRACE_PKG_ROOT = normalizePath("."))
devtools::load_all(".", quiet = TRUE)

object <- MSdev_load("C:/Users/91879/OneDrive/Documents/YLF_Lab/Project/2025.10.10.PAVE/PAVE-Matlab/example/raw_neg-20251010T060553Z-1-001/MSdev_2026_06_03.Rdata")
if (is.null(TRACE:::.trace_get_temp(object, "Negative", "cn.seed.assign"))) {
  object <- TRACE_network_assignment(object, i.pol = 0)
}

cache_dir <- "cache"
demos <- TRACE_cache_all_network_split_demos(
  object,
  i.pol = 0,
  cache_dir = cache_dir,
  min_group_size = 4L
)

cat("\nExported", length(demos), "demos:\n")
idx <- readRDS(file.path(cache_dir, "demo_index.rds"))$index
print(idx)
message("DONE_CACHE")
