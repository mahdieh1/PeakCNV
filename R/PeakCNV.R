#' Runs all or selected steps in PeakCNV analysis.
#'
#' @param run.to Numeric vector defining which steps to run, default = 1, 2, 3.
#' @param working.dir Working directory, default = NULL will use current working directory.

PeakCNV = function(run.to =c(1,2,3),working.dir = NULL)
{

  # Set working directory
  if (is.null(working.dir)) {

    working.dir = getwd()

  }

  ## check format of working directory ##
  if (substr(working.dir, nchar(working.dir), nchar(working.dir)) != "/") {

    working.dir = paste(working.dir, "/", sep = "")

  }
  setwd(working.dir)
  # Load all dependencies ##
  if("BiocManager" %in% rownames(installed.packages()) == FALSE) {

    print("install [BiocManager]")
    install.packages("BiocManager")
    print("install [HelloRanges]")
    BiocManager::install("HelloRanges")
  }

  if("data.table" %in% rownames(installed.packages()) == FALSE) {

    print("install [data.table]")
    install.packages("data.table")

  }

  if("dbscan" %in% rownames(installed.packages()) == FALSE) {

    print("install [dbscan]")
    install.packages("dbscan")

  }
if("magrittr" %in% rownames(installed.packages()) == FALSE) {

    print("install [magrittr]")
    install.packages("magrittr")

  }
if("dplyr" %in% rownames(installed.packages()) == FALSE) {

    print("install [dplyr]")
    install.packages("dplyr")

  }
  ## Load all required packages ##
  library(HelloRanges)
  library(data.table)
  #library(NbClust)
  library(dbscan)
  library(dplyr) 
  library(magrittr)
  
  if (1 %in% run.to) {
    print("step1: Building CNVRs ...")
    BuildingCNVRs()
  }
  if (2 %in% run.to) {
    print("step2: Genrating Input Matrix for clustering ...")
    InputMatrix()
    print("step2: Clustering ...")
    Clustering()
  }
  if (3 %in% run.to) {
    print("step3: Selection Process ...")
    Selection()
  }
}


