RenameGenesSeurat <- function(obj , newnames ) {
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA
  names(newnames) <- rownames(obj)

  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]] <- newnames
    if (length(RNA@scale.data)) rownames(RNA@scale.data) <- newnames[rownames(RNA@scale.data)]
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}


