####################################################################################################################
#                                                                                                                  #
#                     Jasper's Chip-seq pipeline                                                                   #
#                                                                                                                  #
####################################################################################################################

cat("\014")  # Clears the console.

library("devtools") 
library("fastqcr")
library("rJava")
library("tools")
library("R.oo") 
library("R.methodsS3") 
library("methods")
library("R.utils") 
library("reticulate")
library("stringr")

###  Constants  ####################################################################################################

directory <- "/home/jasper/Documents/ChIP/test/"  # Set directory to retreive files from.
Reference <- list.files(pattern = "D39.fasta") # Set reference genome.
genome_size <- "2.03e+6" # Set the size of your reference genome in basepairs (for MACS2)
fastqc_path <- "/usr/bin/fastqc" # Set you path to the program fastqc (already set to standard)
trimmomatic_location <- "~/Downloads/Trimmomatic-0.38/trimmomatic-0.38.jar" # Set the path to your trimmomatic jar file.

# Please set your sequence files as <type>-<samplename>-<replicate#>.fastq.gz, where type is input or IP 
# (e.g. IP-D39-R1.fastq.gz or input-ParB-R2.fastq.gz)
# If you don't want to change the names, MACS2 will not function (all other functions will).

##  Functions  #####################################################################################################

# Function to move files from one folder (oldSubDir) to new one (newSubDir) based on extention of files
Move.Files <- function(oldSubDir, extension, newSubDir){
  FilesList <- list.files(path = oldSubDir, pattern = extension)
  for (file in FilesList){ # Move all files to a new folder
    file.rename(from = paste(oldSubDir, "/", file, sep = ""),  to = paste(newSubDir, "/", file, sep = ""))
  }
}

# Function to unzip input files from gunzip format to fastq format
Unzip.Files <- function(){
  setwd(directory) # Set location of files just to be sure.
  subDir <- "1-ZippedFiles" # Create a subdir for input files.
  dir.create(file.path(subDir), showWarnings = FALSE) # Give a warning when the dir already excists, but we dont care about that.
  subDir2 <- "2-FastqFiles" # Also create a subdir for the output files
  dir.create(file.path(subDir2), showWarnings = FALSE) 
  
  files <- list.files(path = subDir, pattern = "\\.gz$", ignore.case=TRUE) # List al gunzipped files
  for (file in files){
    print(paste("unzipping", file, "..."))
    gunzip(paste(subDir, "/", file, sep = "")) # Gunzip file in first dir
  }
  Move.Files(oldSubDir = subDir, extension = "fastq", newSubDir = subDir2) # Use move file function to move new fastq files to the second folder
  print("Done unzipping")
}

# Function that uses Trimmomatic to trim the fastq files to specified qualitiy settings and creates equal forward en reverse paired-end files.
Trim.Files <- function(){
  setwd(directory) # Set location of files just to be sure.
  subDir <- "2-FastqFiles"
  FastqSequences <- list.files(path = subDir, pattern = ".fastq", ignore.case=TRUE)
  setwd(paste(directory, subDir, sep = "/")) # Go into the fastq folder to find the files.
  Listlength <- length(FastqSequences)
  for (i in seq(from=1, to=Listlength, by=2)) {  # Use by=2 to skip one file every loop. This way the forward and reverse files can be trimmed together to delete unpaired sequences.
    # Set up strings for Trimmomatic and concatenate to one string.
    # Settings in First and Last, feel free to edit.
    TrimStringFirst <- paste("java -jar ", trimmomatic_location, " PE", sep="")
    TrimStringLast <- "ILLUMINACLIP:TruSeq3-PE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100"
    # Variable File names, input is FastqSecuence <file_name.fastq, output will be Trimmed-<file_name>.fastq
    #TrimmedSequence <- paste("Trimmed", String, sep = "_")
    R1_P <-  paste(sub('\\.fastq$', '', FastqSequences[i]), "_paired.fastq", sep = "") # removes .fastq
    R1_UP <-  paste(sub('\\.fastq$', '', FastqSequences[i]), "_unpaired.fastq", sep = "") # removes .fastq
    R2_P <-  paste(sub('\\.fastq$', '', FastqSequences[i+1]), "_pairedasd.fastq", sep = "") # removes .fastq
    R2_UP <-  paste(sub('\\.fastq$', '', FastqSequences[i+1]), "_unpaired.fastq", sep = "") # removes .fastq
    # Concatenate all strings together as input command.
    TrimmomaticString <- paste(TrimStringFirst, FastqSequences[i], FastqSequences[i+1], R1_P, R1_UP, R2_P, R2_UP, TrimStringLast, collapse = " ")
    # Input the concatenated Trimmomatic string to system
    system(TrimmomaticString)
  }
  file.remove(list.files(pattern = "unpaired")) # unpaired files are not of use anymore and removed to save space.
  setwd(directory)
  subDir2 <- "3-TrimmedFiles" # Also create a subdir for the output files
  dir.create(file.path(subDir2), showWarnings = FALSE) # Give a warning when the dir already excists, but we dont care about that.
  newExtension <- "paired"
  Move.Files(oldSubDir = subDir, extension = newExtension, newSubDir = subDir2)
  oldNames <- list.files(path = subDir2, pattern = "paired")
  newNames <- paste(sep = "", sub('\\_paired.fastq$', '', oldNames), ".fastq")
  
  for (i in 1:length(oldNames)){
    file.rename(paste(subDir2, "/", oldNames[i], sep = ""), paste(subDir2, "/", newNames[i], sep = ""))
  }
}

# This function performs a qc report of the fastq files after trimming.
Make.QC.Report <- function(){
  subDir <- "3-TrimmedFiles" # Set to subdir where trimmed files are located.
  fastqc(fq.dir = subDir, qc.dir = NULL, threads = 4, fastqc.path = fastqc_path) # Perform the Fastqc in the right subdirectory and set to 4 cores.
  }

# This function uses bowtie2 (better for sequences above the 50 bp's) to align the fastq (forward and reverse) files to the reference genome (set in the constants tab)
Make.Alignment <- function(){
  # Fist, build index out of the reference genome, named Reference. Output is called: "reference.index".
  system(paste("bowtie2-build", Reference, "reference.index"))
  subDir <- "3-TrimmedFiles/"  
  FastqSequences <- list.files(path = subDir, pattern = "fastq")
  subDir2 <- "4-Alignments" # Also create a subdir for the output files
  dir.create(file.path(subDir2), showWarnings = FALSE) # Give a warning when the dir already excists, but we dont care about that.
  Listlenght <- length(FastqSequences)
  for (i in seq(from=1, to=Listlenght, by=2)) {
    SamAligned <- paste(sub('\\fastq$', '', FastqSequences[i]), ".sam", sep = "")
    # Now call bowtie2, -x for index, -U for fastq file and -S for output sam file
    # system(paste("bowtie2", "-x", "D39.index", "-U", TrimmedSequence, "-S", SamAligned)) ---- Single end bowtie (-U means unpaired)
    system(paste("bowtie2 -x reference.index -1", paste(subDir, "/", FastqSequences[i], sep = ""), "-2", paste(subDir, "/", FastqSequences[i+1], sep = "")
                 , "-S", paste(subDir2, "/", SamAligned, sep = ""), "-p 8", collapse = " "))
  }      
}

# This function uses samtools to compress the sam files to bam files, saving space and enabling the files to be used for MACS2 (requires bam files for PE runs)
Compress.Files <- function(){
  setwd(paste(directory, "4-Alignments", sep="/"))
  SamFiles <- list.files(pattern = "sam")
  for (file in SamFiles){
    system(paste("samtools view -Sb", file, " > ", paste(sub('\\.sam$', '', file), ".bam", sep = "")))  # Use 'samtools view -Sb' as command for sam > bam conversion
  }
}

# This function makes a file list of all aligned files for later use in MACS2
Make.File.List <- function(){
  setwd(paste(directory, "4-Alignments", sep="/"))
  BamSequences <- list.files(pattern = "\\.bam$")
  FileBase <- list() # Remove the extentions and save basename to FileBase
  x = 1
  for (sequence in BamSequences){ # Reomove the extentions
    FileBase[[x]] <- sub('\\.bam$', '', sequence) # removes .fastq
    x = x + 1
  }

  SampleNames <- list() # Create a list to save all filenames
  x = 1
  for (name in FileBase){ # Recognize the name of a sample, between the 2 dashes.
    SampleNames[x] <- str_match(name, "-(.*?)-") # Saves as "-Name-"
    SampleNames[x] <- substr(SampleNames[x], 2, nchar(SampleNames[x])-1) # Removes dashes to save as "Name" 
    x = x + 1
  }
  SampleNames <- unique(SampleNames) # Removes duplicates from list
  
  FileStrings <- list() # Creates a list to save all names in a nice structure for later use. list -> FileName -> Type -> strings...
  for (name in SampleNames){
    inputNames <- paste("input", name, sep = "-") # Put togther the inputs: input-Name
    IPNames <- paste("IP", name, sep = "-") # And the IP's
    
    input <- intersect(list.files(pattern = inputNames), list.files(pattern = "bam")) # Now filter those strings containing the input-Name and .sam extention and save in input vector
    IP <- intersect(list.files(pattern = IPNames), list.files(pattern = "bam"))
    
    FileStrings[[name]] <- list(input, IP) # Save the two vectors to the list under the sample name
    names(FileStrings[[name]]) <- c("input", "IP") # Give the vectors the right name
  }
  return(FileStrings) # Return for MACS2
}

# This function performs peak calling using MACS2 and the file list of all IP and input samples. It compares the inputs to IP's.
Call.Peaks <- function(File_list){
  for (Sample in File_list){
    # Get the file lists in order. SampleNamesList((inputSampleNameList-R#),(IPSampleNameList-R#))
    # MACS2: -t = sample, -c = mock control, add --broad for broad peaks, -g = genome size (bp's)
    firststring <- "macs2 callpeak -t"
    laststring <- paste("-f BAMPE --broad -g ", genome_size, sep="")
    
    input <- vector() # Initializing the vector, for the appending function later on. Clears it every loop.
    for (i in 1:length((Sample$input))){
      input <- c(input, Sample$input[[i]]) # appending all values into the vector for the MACS2
    }
    IP <- vector()
    for (i in 1:length((Sample$IP))){
      IP <- c(IP, Sample$IP[[i]])
    }
    
    macs2_string <- paste(c(firststring, IP, "-c", input, "-n", Sample$input[[1]], laststring), collapse=" ")
    system(paste(macs2_string))
  }
}

###  Program  flow ###############################################################################################

Unzip.Files()
Trim.Files()
Make.QC.Report()
Make.Alignment()
Compress.Files()
File_list <- Make.File.List()
Call.Peaks(File_list)
