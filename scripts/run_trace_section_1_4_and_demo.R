suppressPackageStartupMessages({
  library(MSdev)
  library(dplyr)
})
Sys.setenv(TRACE_PKG_ROOT = normalizePath("."))
devtools::load_all(".", quiet = TRUE)

object_path <- "C:/Users/91879/OneDrive/Documents/YLF_Lab/Project/2025.10.10.PAVE/PAVE-Matlab/example/raw_neg-20251010T060553Z-1-001/MSdev_2026_06_03.Rdata"
export_path <- "c:/Users/91879/OneDrive/Documents/YLF_Lab/Project/2025.10.10.PAVE/PAVE-Matlab/example/TRACE_export.xlsx"
cpdb <- "d:/data/2025.12.26.PAVE2/trace.cp.db.xlsx"
cache_dir <- "cache"

message("Loading MSdev object...")
object <- MSdev_load(object_path)

message("TRACE_network_assignment (i.pol = 0)...")
object <- TRACE_network_assignment(object, i.pol = 0)

message("TRACE_annotate (i.pol = 0)...")
object <- TRACE_annotate(object, i.pol = 0, cpdb = cpdb)

message("TRACE_export...")
TRACE_export(object, file = export_path)

message("MSdev_save...")
MSdev_save(object)

message("Regenerating network split demos...")
demos <- TRACE_cache_all_network_split_demos(
  object,
  i.pol = 0,
  cache_dir = cache_dir,
  min_group_size = 4L
)

cat("\nExported", length(demos), "demos:\n")
idx <- readRDS(file.path(cache_dir, "demo_index.rds"))$index
print(idx)
message("DONE")
