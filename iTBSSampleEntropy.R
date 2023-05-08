#calculate Sample Entropy of COP data in AP/ML directions
#created 3/4/2022 by Jack Manning
#Packages: astsa, seewave, pracma, zoo, rlist
#C:\Users\jackpmanning\OneDrive - Texas A&M University\Documents\Projects\Manuscript JM_ES\Manning DT Data
#Purpose: The code sets the working directory, loops through subfolders, reads in data from text files, 
#creates new data frames, downsamples the data, calculates center of pressure (COP) positions, and calls the 
#function COPAnalysis.

#1. housekeeping print execution time, whats currently loaded, and then load the packages

Sys.time()
ls()
rm(list=ls())
ls()
starttime <- Sys.time()
#sink("C:\\Users\\jackpmanning\\OneDrive - Texas A&M University\\Documents\\Projects\\Manuscript JM_ES\\nohardcodeloop.RData", split=TRUE)
#sink("\\nohardcodeloop.RData", split=TRUE)

#load in appropriate packeges
library(astsa)
library(seewave)
library(pracma)
library(zoo)
library(rlist)
library(writexl)
#library(sazedR)

#_______________________________________________________________________________
#Create Functions

#create functions for downsample and COPAnalysis
#downsample source code
downsample <- function(data, window_size)
{
  return(ts(as.ts(rollapply(zoo(data), width=window_size, by=window_size, FUN=mean)), frequency=1))
}

# This function is for manipulating text objects
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#_______________________________________________________________________________

#Desktop in 332
#setwd("C:\\Users\\jackpmanning\\Documents\\MSThesisData")
#setwd("005")	#change this each time to the correct folder

#Work Laptop
#1. set the working directory and the output files
setwd("C:\\Users\\jackpmanning\\OneDrive - Texas A&M University\\Documents\\Projects\\AshwiniiTBS\\Subjects")
upone <- file.path("C:\\Users\\jackpmanning\\OneDrive - Texas A&M University\\Documents\\Projects\\AshwiniiTBS")
outputs <- file.path(paste(upone, "\\", "Outputs", sep=""))
xlout <- file.path(paste(upone, "\\", "xlout", sep=""))


#2. get a list of all the subfolders
foldernames <- as.list(list.dirs(getwd(), full.names = FALSE, recursive = FALSE))

# for(i in 1:length(foldernames)){
#   #7. Make new folders for each data output
#   new_folder_name <- paste(foldernames[i], "pdf", sep="")
#   # set the path of the new subfolder
#   new_folder_path <- file.path(outputs, "\\", new_folder_name)
#   new_folder_path <- gsub("/", "\\", new_folder_path)
#   # create the new subfolder
#   dir.create(new_folder_path, recursive = TRUE)
# }

senoutputDF <- data.frame(matrix(ncol = 3, nrow = 0))

p <- 0
#1.cloop through a folder that contains folders which contain data in .txt files
for(i in 1:length(foldernames)){
  
  #1a. create a list with the names of all data files in current working directory
  filename <- as.list(list.files((paste(getwd(), "/", foldernames[i], sep="")), pattern="\\.txt$", full.names=TRUE))
  for(j in 1:length(filename)){
    
    #clean up the file name string
    current.file <- gsub("/", "", substrRight(filename[j], 13))
    current.file <- gsub(".txt", "", current.file)
    print(current.file)
    
    #Load current data file which is filename[[j]]
    data1 <- read.delim(filename[[j]], header = F, sep = ",")
    
    #rename the data frame columns
    names(data1)[3] <- "forcez"
    names(data1)[4] <- "momentx"
    names(data1)[5] <- "momenty"
    
    #4. create vectors of 6000 0's
    forcez <- c(rep(0,6000))
    momentx <- c(rep(0,6000))
    momenty <- c(rep(0,6000))
    
    #create data frame with zero vectors
    dfzeros <- data.frame(forcez, momentx, momenty)
    
    #downsample data then replace dfzeros vectors with downsampled data
    dfzeros$forcez <- downsample(as.ts(ifelse(data1[["forcez"]] == 0, data1[["forcez"]][1], data1[["forcez"]])), window_size = 2)
    dfzeros$momentx <- downsample(as.ts(data1[["momentx"]]), window_size = 2)
    dfzeros$momenty <- downsample(as.ts(data1[["momenty"]]), window_size = 2)
    
    #create x and y position vectors
    dfzeros$xpos <- as.ts(dfzeros[["momentx"]]/dfzeros[["forcez"]])
    dfzeros$ypos <- as.ts(dfzeros[["momenty"]]/dfzeros[["forcez"]])
    
    #Sample Entropy Analysis  
    SEnX.2 <- sample_entropy(dfzeros$xpos, edim=2, r = 0.2*sd(dfzeros$xpos), tau=1)
    SEnY.2 <- sample_entropy(dfzeros$ypos, edim=2, r = 0.2*sd(dfzeros$ypos), tau=1)
    
    #create data frame with sen values
    senoutput <- data.frame(SEnX = SEnX.2, SenY = SEnY.2)
    
    # Overwrite copoutput with current file's cop analysis
    senoutput <- cbind(as.data.frame(current.file), senoutput)
    
    # Add cop analysis from current file to data frame containing all data
    senoutputDF <- rbind(senoutputDF, senoutput)
    
  }
  
  # Write data frame to excel file
  write_xlsx(senoutputDF, paste(xlout, "\\", foldernames[[i]], "NonlinearAnalysis.xlsx", sep=""))
  
  # re-initialize output data frame
  senoutputDF <- data.frame(matrix(ncol = 3, nrow = 0))
  
  #count loop iterations
  p <- p+1
  print(p)
  
}


#print runtime
end_time <- Sys.time()

end_time - starttime



