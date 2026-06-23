suppressPackageStartupMessages({
  library(MSdev)
  library(dplyr)
})
Sys.setenv(TRACE_PKG_ROOT = normalizePath("."))
devtools::load_all(".", quiet = TRUE)

object <- MSdev_load("C:/Users/91879/OneDrive/Documents/YLF_Lab/Project/2025.10.10.PAVE/PAVE-Matlab/example/raw_neg-20251010T060553Z-1-001/MSdev_2026_06_03.Rdata")
cpdb <- "d:/data/2025.12.26.PAVE2/trace.cp.db.xlsx"
message("Re-running TRACE_network_assignment with current blocking strategy...")
object <- TRACE_network_assignment(object, i.pol = 0)
message("Re-running TRACE_annotate...")
object <- TRACE_annotate(object, i.pol = 0, cpdb = cpdb)

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
