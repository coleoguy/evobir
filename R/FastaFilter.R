# Simple function takes fasta and returns
# fitlered fasta with only the longest sequence 
# for each species.

FastaFilter <- function(f.in, folder = F, prefix="pruned"){
  # f.in if we are processing just a single file this should be
  # a text string holding the name or name and path to the file
  
  # folder if T then all files in the folder will be 
  # processed by the function
  
  # prefix as a default we save new copies with the prefix
  # pruned but we can put whatever we want here
  
  # first we have a function for a single file
  single.file <- function(f.in, prefix){
    # this reads in our file
    seqs <- read.dna(f.in, as.character = T, format="fasta", as.matrix = F)
    # gets a lits of the unique taxa
    taxa.present <- unique(names(seqs))
    # this will hold our final set of sequences
    final.list <- list()
    # here we find the longest sequence for each species and store it
    for(i in 1:length(taxa.present)){
      hits <- which(names(seqs) == taxa.present[[i]])
      possible <- seqs[hits]
      z <- which.max(unlist(lapply(possible, length)))
      final.list[i] <- possible[z]
      names(final.list)[i] <- names(possible)[z]
    }
    # save the result as a fasta file
    write.dna(x=final.list, file=paste(prefix, ".", f.in, sep = ""), 
              format = "fasta")
  }
  # here we call our function for a single file
  if(folder == F){
    single.file(f.in, prefix)
  }
  # here we get a list of all files and call our
  # function itterativly to process each file
  if(folder == T){
    files <- list.files()
    for(i in 1:length(files)){
      single.file(f.in = files[i], prefix = prefix)
    }
  }
}