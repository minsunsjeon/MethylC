############################################################
# CPGCluster #
############################################################

# algorithm to be run on input csv file to combine islands by distance between position via the meanshift algorithm
# relies on the C++ methylation package (mean_shift0324.cpp)
# only works on linux or mac machines (for the parallel package to work)

# INPUT
#
# file - file (path of csv file) \
# input data generated from MethylDackel
# data format: [Chromosome, Start Position, End Position, Sample 1_Methylated, Sample 1_Unmethylated, ...]
# expects the columns, "Chromosome," "Start Position", "End Position"
#
# setpsize - size of window for the meanshift algorithm (int)
# default value of 50
#

# OUTPUT
#
# file - meanshift_result.csv (file in directory)
# [Chromosome, Start Position, End Position, Sets, Sample 1_Methylated, Sample 1_Unmethylated, ...]
# Sets is the number of sets combined into each row
#

############################################################
# FUNCTIONS #
############################################################

######## sum the columns ########
computeplus <- function(x){
  x <- na.omit(x)
  if (length(x)==0) {return(NA)}
  else {return(sum(x))}
}

######## get the start of the column ########
getstart <- function(x){
  return(x[1])
}

######## get the end of the column ########
getend <- function(x){
  return(x[length(x)])
}


####### apply meanshift #######
# [ returns list of lists of indices that make up islands for each island ]
# x ----------- dataframe (in this case the global var data)
#               format: [Chromosome, Start Position, End Position, Sample 1_Methylated, Sample 1_Unmethylated, ..., ...]
# stepsize ---- integer (window size that defines range single sets
#               can differ from the island center to be in the island)
applyMeanshift <- function(x, stepsize=50){

  tmppos <- x$`Start Position`

  # Ensure the C++ function is available
  outpos = mean_shift(tmppos, stepsize)
}



########## formatting ##########
# [ returns final data with combined islands and unchanged single sets ]
# x - list of dataframe for each chromosome, and list of island indices
#     columns: [Chromosome1, [[20, 21]]]
formatIslands <- function(x) {
  chrom <- x[[1]]
  islands <- x[[2]]

  islands_unlisted <- unlist(islands)
  if (length(islands_unlisted) != 0) {
    finalobj <- chrom[-islands_unlisted, ] # all rows but those in islands
    finalobj <- finalobj %>%
      mutate(Sets=1, .after=`End Position`)

    summed_rows <- lapply(islands, function(x) c(chrom[getstart(x),"Chromosome"], chrom[getstart(x), "Start Position"], chrom[getend(x), "End Position"], length(x), dplyr::summarise_all(chrom[getstart(x):getend(x),4:ncol(chrom)],computeplus))) # list of summed rows
    summed_rows <- lapply(summed_rows, function(x) setNames(x, colnames(finalobj))) # change it so the column names match
    summed_rows <- rbindlist(summed_rows, use.names=TRUE, fill=TRUE, idcol=NULL) # put the summed_rows into DF

    finalobj <- rbind(finalobj, summed_rows) # add to final DF
    finalobj <- finalobj[order(finalobj$`Start Position`),] # put them in the correct order

    finalobj <- data.table(finalobj)
  }
  else {
    finalobj <- chrom
  }

  return(finalobj)
}


############################################################
# USING FUNCTIONS ON DATA #
############################################################

#' Run Meanshift Clustering
#'
#' @param file Path to the CSV file containing data.
#' @param stepsize Integer window size for clustering (default: 50).
#' @return A data frame with meanshift clustering results.
#' @export
run_meanshift <- function(file, stepsize = 50) {
  # Load C++ function
  sourceCpp(system.file("src", "mean_shift0324.cpp", package = "CpGCluster"))

  message("\nLoad data")
  data <- data.table::fread(file)

  # Get number of cores for parallel processing
  numCores <- parallel::detectCores()

  ############
  # APPLYING
  ############

  message("\nSeparate the data by chromosome")
  ori.list <- list()
  for(i in c(1:22,"X","Y")){
    ori.list[[i]] <- data[which(data$Chromosome == paste(i)),]
  }

  message("\nGet islands")
  # List of indices in islands
  island.list <- parallel::mclapply(ori.list, applyMeanshift, stepsize, mc.cores = numCores)

  # List of dataframe for each chromosome and indices in the island
  format.list <- list()
  for(i in c(1:22,"X","Y")){
    format.list[[i]] <- list()
    format.list[[i]][[1]] <- ori.list[[i]]
    format.list[[i]][[2]] <- island.list[[i]]
  }

  message("\nStart formatting")
  out.list <- parallel::mclapply(format.list, formatIslands, mc.cores = numCores)

  # Create an empty data frame to store results
  outframe <- data.frame(matrix(0, nrow = 0, ncol = 20,
                                dimnames = list(NULL, paste0(colnames(ori.list[[1]])))))

  class(outframe$Chromosome) <- "character"
  result <- data.table::setDT(outframe)

  result <- rbindlist(out.list, use.names=TRUE, fill=TRUE, idcol=NULL)

  return(result)
}
