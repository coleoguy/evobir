# re-order alignment base on physical mapping 
# Sean Chien 
# count - until hit the first nucleotide
# re-order them by number of '-' before the first nucleotide
# reference sequences go to front. vector number or character?
# in order to run this function, seqinr is required. 

# check if references are in the file 
# print references 
# library(seqinr)
######################################################################################
# Intro
# re-ordering sequences based on physical mapping in order to visualize easier.
# file for original file's path and name. newfile for new file's path and name. 
# ref is optional. It is a vector and can contain either numeric or character. 
# ref default is off. If it is on, the reference sequence(s) will be in the beginning.
######################################################################################

ReOrderAlignment <- function(file, newfile, ref=NULL){
  dat <- read.fasta(file)
  ### no reference
  if (is.null(ref)){
    result = c()
    print("No reference")
    print("Re-ordering...")
    for (j in 1:length(dat)){
      dash = 0
      for (i in 1:length(dat[[j]])){
        if (dat[[j]][i] == "-"){
          dash <- dash + 1
        }
        else if (!dat[[j]][i] == "-") break
      }
      result[j]<- dash
    } 
    newlist<-dat[order(result)]
    write.fasta(newlist, names = names(newlist), file.out = newfile)
    print("Done!")
  }
  
  ### with reference(s)
  # find where is the references and take it away?
  # check if the references are in the file first.
  if (!is.null(ref)){
    for (i in 1:length(ref)){
      if ((names(dat[ref][i])) %in% names(dat)){
        print(paste("The reference you picked: ", names(dat[ref][i])))
        if (i == length(ref)){
          print("Re-ordering...")
          # taking off the reference(s) for the original data and run re-order
          norefname<-names(dat)[-c(match(names(dat[ref]),names(dat)))]
          noreferdat<-dat[norefname]
          # run 
          result = c()
          for (j in 1:length(noreferdat)){
            dash = 0
            for (i in 1:length(noreferdat[[j]])){
              if (noreferdat[[j]][i] == "-"){
                dash <- dash + 1
              }
              else if (!noreferdat[[j]][i] == "-") break
            }
            result[j]<- dash
          } 
          newlist<-noreferdat[order(result)]
          newlist<-c(dat[c(match(names(dat[ref]),names(dat)))],newlist)
          write.fasta(newlist, names = names(newlist), file.out = newfile)
          print("Done!")
        }
      }
      # if the references are not in the file
      else if (!(names(dat[ref][i])) %in% names(dat)){
        print(paste("Oops! the reference,'",ref[i],"',is not in your file. Please check again."))
        break
        # stop the whole process
      }
    }
  }
}


