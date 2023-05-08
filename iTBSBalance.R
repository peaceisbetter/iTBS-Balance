#calculate linear analysis of COP data in AP/ML directions
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


#create functions for downsample and COPAnalysis
#downsample source code
downsample <- function(data, window_size)
{
  return(ts(as.ts(rollapply(zoo(data), width=window_size, by=window_size, FUN=mean)), frequency=1))
}

#COP Analysis code
COPAnalysis <- function(PathName, FileName, sampfreq, COPx, COPy) {
  library(writexl)
  xlout <- file.path("C:\\Users\\jackpmanning\\OneDrive - Texas A&M University\\Documents\\Projects\\Manuscript JM_ES\\xlout")
  
  # calculate number of data points
  n <- length(COPx)
  
  # create time vector
  t <- (0:(n-1)) * (1/sampfreq)
  
  # get current date and time for filename
  datetime <- format(Sys.time(), "%Y_%m_%d_%H-%M")
  
  ## Begin analysis of the segments
  
  # Normalizing x and y data.
  COPxnorm <- COPx - mean(COPx)
  COPynorm <- COPy - mean(COPy)
  
  # Finds radius of COP path.
  COPd <- sqrt(COPxnorm^2 + COPynorm^2)
  
  # Plot COP data, x and y versus the time.
  #send graphical output to pdf
 
  
  ## Linear analyses
  
  # Root Mean Square (x,y,d)
  rms_x <- sqrt((1/length(COPxnorm)) * sum(COPxnorm^2))
  rms_y <- sqrt((1/length(COPynorm)) * sum(COPynorm^2))
  rms_d <- sqrt((1/length(COPd)) * sum(COPd^2))
  
  # Range of mediolateral and anteroposterior motion
  rangeml <- max(COPxnorm) - min(COPxnorm)
  rangeap <- max(COPynorm) - min(COPynorm)
  
  # Sway Path
  swaypathd <- sum(abs(diff(COPd)))
  swaypathx <- sum(abs(diff(COPxnorm)))
  swaypathy <- sum(abs(diff(COPynorm)))
  
  # This version of the calculation provides the total tangential distance.
  i <- 1:(n-1)
  swaypathtangential <- sum(sqrt((COPxnorm[i+1] - COPxnorm[i])^2 + 
                                   (COPynorm[i+1] - COPynorm[i])^2))
  
  # Area of 95% Confidence Circle
  # (Notes) Finds the area of a circle that includes 95% of radii(COPd).
  # This is a one sided test so the z score has a value of 1.645. It
  # assumes a normal distribution of radii, as does Prieto. However the
  # distribution of radii is not normally distributed. For now use 1.645
  # from the normal pdf. Chi-square pdf needed?
  
  confcir95 <- mean(COPd) + 1.645 * sd(COPd)
  Acir <- pi * confcir95^2
  
  # Area of 95% Confidence Ellipse
  # (Notes) This is from Prieto (1996) and from Sokal and Rohlf
  # (1995) Biometry p589, also cited by Prieto.
  
  sML <- sqrt((1/length(COPxnorm)) * sum(COPxnorm^2))  # Prieto
  sAP <- sqrt((1/length(COPynorm)) * sum(COPynorm^2))  # Prieto
  sAPML <- (1/length(COPxnorm)) * sum(COPxnorm * COPynorm)  # Prieto
  
  prieto <- list()
  
  prieto$f <- 3
  
  # (Notes) Prieto has missed squaring the sum of the first two terms in his
  # equation(16) page 959. This equation is actually from Sokal and Rohlf
  # (1995).
  
  prieto$d <- sqrt((sAP^2 + sML^2)^2 - 4 * (sAP^2 * sML^2 - sAPML^2))
  
  prieto$ellipRA <- sqrt(prieto$f * (sAP^2 + sML^2 + prieto$d))  # Prieto eq 14, (1995)
  prieto$ellipRB <- sqrt(prieto$f * (sAP^2 + sML^2 - prieto$d))  # Prieto eq 15
  prieto$ellipA <- 2 * pi * prieto$f * sqrt(sAP^2 * sML^2 - sAPML^2)  # Prieto eq 18
  
  prieto$lambda <- (sAP^2 + sML^2 + prieto$d) / 2  # from Sokal and Rohlf
  prieto$slope <- prieto$ellipRA / prieto$ellipRB
  # prieto$slope <- sAPML / (prieto$lambda - sML^2) # from Sokal and Rohlf
  prieto$angle <- atan(prieto$slope)

 
  # # Save figure
  # # saveas(powfig,fullfile(PathName,[FileName '_Power_Spectrum']),'jpg')
  # 
  # # Median Frequency
  # # (Notes) Prieto discards the first two points in the spectra
  # # and only uses data to 5 Hz. We have frequencies up at 7, so
  # # cutoff here is at 10 Hz - to avoid 60 Hz noise
  # # when doing power spectral densities. See left
  # 
  # freqindex <- 1
  # while (f[freqindex] < 10 && freqindex != length(f)){
  #   freqindex <- freqindex + 1
  # }
  # analysisspectrum <- Spectrum[3:freqindex]
  # analysisfrequency <- f[3:freqindex]
  # 
  # CumSumPower <- cumsum(analysisspectrum)
  # FindMedian <- which(CumSumPower > 0.5 * sum(analysisspectrum))
  # MedianIndex <- min(FindMedian)
  # 
  # medianfreq <- (0.5 / length(Spectrum)) * sampfreq * MedianIndex

  # #Frequency dispersion
  # Mu0 <- (1 / length(analysisspectrum)) * sum(analysisspectrum)
  # Mu1 <- (1 / length(analysisspectrum)) * sum(analysisspectrum * analysisfrequency)
  # Mu2 <- (1 / length(analysisspectrum)) * sum(analysisspectrum * analysisfrequency^2)
  # 
  # freqdisp <- sqrt(1 - Mu1^2 / (Mu0 * Mu2))

  # #Add new variables to output  
  # copoutput <- data.frame(rms_x, rms_y, rms_d, 
  #                         rangeml, rangeap, 
  #                         swaypathd, swaypathx, 
  #                         swaypathy, swaypathtangential)
  
  copoutput <- data.frame(
    'rms_x' = rms_x,
    'rms_y' = rms_y,
    'rms_d' = rms_d,
    'rangeap' = rangeap,
    'rangeml' = rangeml,
    'swaypathd' = swaypathd,
    'swaypathx' = swaypathx,
    'swaypathy' = swaypathy,
    'swaypathtangential' = swaypathtangential,
    'Acir' = Acir,
    'prieto.ellipA' = prieto$ellipA
    # 'medianfreq' = medianfreq,
    # 'freqdisp' = freqdisp
    )
  
  assign("copoutput", copoutput, envir = .GlobalEnv)
  
  # write_xlsx(copoutput, paste(xlout, current.file, ".xlsx", sep=""))
  
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#________________________________________________________________________________________________


#_______________________________________________________________________________________________________________

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

copoutputDF <- data.frame(matrix(ncol = 12, nrow = 0))

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
    
    #COPAnalysis <- function(PathName, FileName, sampfreq, COPx, COPy)  
    COPAnalysis(paste(outputs, "\\",foldernames[[i]], sep=""), current.file, sampfreq = 50, COPx =  dfzeros$xpos, dfzeros$ypos)
    
    # Overwrite copoutput with current file's cop analysis
    copoutput <- cbind(as.data.frame(current.file), copoutput)
    
    # Add cop analysis from current file to data frame containing all data
    copoutputDF <- rbind(copoutputDF, copoutput)
    
  }
  
  # Write output data frame to excel file
  write_xlsx(copoutputDF, paste(xlout, "\\", foldernames[[i]], "LinearAnalysis.xlsx", sep=""))
  
  # re-initialize output data frame
  copoutputDF <- data.frame(matrix(ncol = 12, nrow = 0))
  
  #count loop iterations
  p <- p+1
  print(p)
  
}


#print runtime
end_time <- Sys.time()

end_time - start_time





