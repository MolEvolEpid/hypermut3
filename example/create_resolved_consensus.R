library(ape)

ex <- as.matrix(read.FASTA("example.fasta"))

# depending on the sequences, you might have to add other conditionals for other multistate characters
resolved_consensus <- sapply(1:ncol(ex), function(x){
  char <- sort(unname(unique(as.character(ex[,x]))))
  char <- paste0(char[char != '-'], collapse = '')
  if(nchar(char) > 1){
    if(char == 'ag') char = 'r'
    if(char == 'ct') char = 'y'
    if(char == 'at') char = 'w'
    if(char == 'ac') char = 'm'
    if(char == 'gt') char = 'k'
  } 
  return(char)
})

# we don't need to save it here because the resolved consensus is already included
# write.FASTA(rbind(ex, resolved_consensus = as.DNAbin(resolved_consensus)), "example.fasta")