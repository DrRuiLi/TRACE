open_PAVE <- function(){

  pave.file <- MSdev:::get_dir_expand_from_onedrive("Documents/YLF_Lab/Project/2025.10.10.PAVE/code/PAVE_data_Analysis.R")
  rstudioapi::documentOpen(pave.file)
  setRStudioDir()
  return(pave.file)
}
