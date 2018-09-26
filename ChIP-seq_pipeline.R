################################################################################
#                                                                              #
#                       Jasper's Chip-seq pipeline                             #
#                                                                              #
#   This pipeline processes paired-end chip-seq data using several programs.   #
#     First check the constants section before using the pipeline please.      #
#                                                                              #   
################################################################################

cat("\014")  # Clears the console.

library("fastqcr")
library("R.utils") 
library("stringr") 

###  Global Constants  #########################################################

directory <- "/home/jasper/Documents/ChIP/"  # Set directory to retreive files from.
Reference <- list.files(pattern = "D39.fasta") # Set reference genome.
genome_size <- "2.03e+6" # Set the size of your reference genome in basepairs (for MACS2)
fastqc_path <- "/usr/bin/fastqc" # Set you path to the program fastqc (already set to standard)
trimmomatic_location <- "~/Downloads/Trimmomatic-0.38/trimmomatic-0.38.jar" # Set the path to your trimmomatic jar file.

# Please set your sequence files as <type>-<samplename>-<replicate#>.fastq.gz, where type is input or IP 
# (e.g. IP-D39-R1.fastq.gz or input-ParB-R2.fastq.gz)
# If you don't want to change the names, MACS2 will not function (all other functions will).

# Subdirectories that will be used:
GunzippedFiles  <- "1-GunzippedFiles"
UnzippedFiles   <- "2-FastqFiles"
TrimmedFiles    <- "3-TrimmedFastqFiles"
QcReports       <- "4-QcReports"
SamFiles        <- "5-SamFiles"
BamFiles        <- "6-BamFiles"
PeakFiles       <- "7-Peaks"

##  Functions  #################################################################

# Function used to create all subdirectories (if not excistant yet).
Create.dirs <- function() {
  dir.create(file.path(GunzippedFiles), showWarnings = FALSE) # Warnings will be shown if subdirectory already excists, but we don't care about that.
  dir.create(file.path(UnzippedFiles), showWarnings = FALSE) # So, warnings set to FALSE.
  dir.create(file.path(TrimmedFiles), showWarnings = FALSE)
  dir.create(file.path(QcReports), showWarnings = FALSE)
  dir.create(file.path(SamFiles), showWarnings = FALSE)
  dir.create(file.path(BamFiles), showWarnings = FALSE)
  dir.create(file.path(PeakFiles), showWarnings = FALSE)
}

# Function to move files from one folder (oldSubDir) to new one (newSubDir) based on extention of files
Move.Files <- function(oldSubDir, extension, newSubDir) {
  FilesList <- list.files(path = oldSubDir, pattern = extension)
  for (file in FilesList) { # Move all files to a new folder
    file.rename(from = paste(oldSubDir, "/", file, sep = ""),  
                to = paste(newSubDir, "/", file, sep = ""))
  }
}

# Function to copy files from one folder (oldSubDir) to new one (newSubDir) based on extention of files without overwriting excisting files in new directory.
Copy.Files <- function(oldSubDir, extension, newSubDir) {
  FilesList <- list.files(path = oldSubDir, pattern = extension)
  for (file in FilesList) { # Copy all files to a new folder
    file.copy(from = paste(oldSubDir, "/", file, sep = ""),  
              to = paste(newSubDir, "/", file, sep = ""), overwrite = FALSE)
  }
}

# Function to remove files based on pattern
Remove.Files <- function(extension) {
  FilesList <- list.files(pattern = extension)
  for (file in FilesList) { # Copy all files to a new folder
    file.remove(file)
  }
}

# Function to unzip input files from gunzip format to fastq format
Unzip.Files <- function() {
  files <- list.files(path = GunzippedFiles, pattern = "\\.gz$", 
                      ignore.case = TRUE) # List al gunzipped files
  for (file in files) {
    print(paste("unzipping", file, "..."))
    gunzip(paste(GunzippedFiles, "/", file, sep = "")) # Gunzip file in first dir
  }
  Move.Files(oldSubDir = GunzippedFiles, extension = "fastq", 
             newSubDir = UnzippedFiles) # Use move file function to move new fastq files to the second folder
  print("Done unzipping")
}

# Function that uses Trimmomatic to trim the fastq files to specified qualitiy settings and creates equal forward en reverse paired-end files.
Trim.Files <- function() {
  FastqSequences <- list.files(path = UnzippedFiles, pattern = ".fastq", 
                               ignore.case = TRUE)
  setwd(paste(directory, UnzippedFiles, sep = "/")) # Go into the fastq folder to find the files.
  Listlength <- length(FastqSequences)
  for (i in seq(from = 1, to = Listlength, by = 2)) {  # Use by=2 to skip one file every loop. This way the forward and reverse files can be trimmed together to delete unpaired sequences.
    # Set up strings for Trimmomatic and concatenate to one string.
    # Settings in First and Last, feel free to edit.
    TrimStringFirst <- paste("java -jar ", trimmomatic_location, " PE", 
                             sep = "")
    TrimStringLast <- "ILLUMINACLIP:TruSeq3-PE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100" # Adjust to required settings.
    # Variable File names, input is FastqSecuence <file_name.fastq, output will be Trimmed-<file_name>.fastq
    #TrimmedSequence <- paste("Trimmed", String, sep = "_")
    R1_P <-  paste(sub('\\.fastq$', '', FastqSequences[i]), 
                   "_paired.fastq", sep = "") # removes .fastq
    R1_UP <-  paste(sub('\\.fastq$', '', FastqSequences[i]), 
                    "_unpaired.fastq", sep = "") # removes .fastq
    R2_P <-  paste(sub('\\.fastq$', '', FastqSequences[i+1]), 
                   "_paired.fastq", sep = "") # removes .fastq
    R2_UP <-  paste(sub('\\.fastq$', '', FastqSequences[i+1]), 
                    "_unpaired.fastq", sep = "") # removes .fastq
    # Concatenate all strings together as input command.
    TrimmomaticString <- paste(TrimStringFirst, FastqSequences[i], 
                               FastqSequences[i + 1], R1_P, R1_UP, R2_P, R2_UP, 
                               TrimStringLast, collapse = " ")
    # Input the concatenated Trimmomatic string to system
    system(TrimmomaticString)
  }
  file.remove(list.files(pattern = "unpaired")) # unpaired files are not of use anymore and removed to save space.
  setwd(directory) # Revert back to working directory.
  newExtension <- "paired"
  Move.Files(oldSubDir = UnzippedFiles, extension = newExtension, 
             newSubDir = TrimmedFiles)
  oldNames <- list.files(path = TrimmedFiles, pattern = "paired") # Remove paired from file extention.
  newNames <- paste(sep = "", sub('\\_paired.fastq$', '', oldNames), ".fastq")
  
  for (i in 1:length(oldNames)) {
    file.rename(paste(TrimmedFiles, "/", oldNames[i], sep = ""), 
                paste(TrimmedFiles, "/", newNames[i], sep = ""))
  }
}

# This function performs a qc report of the fastq files after trimming.
Make.QC.Report <- function(){
  fastqc(fq.dir = TrimmedFiles, qc.dir = NULL, threads = 4, 
         fastqc.path = fastqc_path) # Perform the Fastqc in the right subdirectory and set to 4 cores.
  }

# This function uses bowtie2 (better for sequences above the 50 bp's) to align the fastq (forward and reverse) files to the reference genome (set in the constants tab)
Make.Alignment <- function() {
  # Fist, build index out of the reference genome, named Reference. Output is called: "reference.index".
  system(paste("bowtie2-build", Reference, "reference.index"))
  FastqSequenceList <- list.files(path = TrimmedFiles, pattern = "fastq")
  Listlenght <- length(FastqSequenceList)
  for (i in seq(from = 1, to = Listlenght, by = 2)) {
    SamAligned <- paste(sub('.fastq', '', FastqSequenceList[i]), "sam", 
                        sep = ".")
    # Now call bowtie2, -x for index, -U for fastq file and -S for output sam file
    # system(paste("bowtie2", "-x", "D39.index", "-U", TrimmedSequence, "-S", SamAligned)) ---- Single end bowtie (-U means unpaired)
    system(paste("bowtie2 -x reference.index -1", paste(TrimmedFiles, "/", 
                                                        FastqSequenceList[i], 
                                                        sep = ""), "-2", 
                 paste(TrimmedFiles, "/", FastqSequenceList[i + 1], sep = "")
                 , "-S", paste(SamFiles, "/", SamAligned, sep = ""), 
                 "-p 8", collapse = " "))
  }      
}

# This function uses samtools to compress the sam files to bam files, saving space and enabling the files to be used for MACS2 (requires bam files for PE runs)
Compress.Files <- function() {
  Copy.Files(SamFiles, "\\.sam$", BamFiles) # First copy sam files to bam folder
  setwd(paste(directory, BamFiles, sep = "/")) # Go into BamFiles dir
  SamFileList <- list.files(pattern = "sam") # List all sam files to convert
  for (file in SamFileList) { # Convert all files in loop
    system(paste("samtools view -Sb", file, " > ", paste(sub('\\.sam$', '', 
                                                    file), ".bam", sep = "")))  # Use 'samtools view -Sb' as command for sam > bam conversion
  }
  Remove.Files("sam") # Now remove the sam files from the folder
  setwd(directory) # Revert directory back to standard
}

# This function makes a file list of all aligned files for later use in MACS2
Make.File.List <- function() {
  setwd(paste(directory, BamFiles, sep = "/"))
  BamSequences <- list.files(pattern = "\\.bam$")
  FileBase <- list() # Remove the extentions and save basename to FileBase
  x = 1
  for (sequence in BamSequences) { # Reomove the extentions
    FileBase[[x]] <- sub('\\.bam$', '', sequence) # removes .fastq
    x = x + 1
  }

  SampleNames <- list() # Create a list to save all filenames
  x = 1
  for (name in FileBase) { # Recognize the name of a sample, between the 2 dashes.
    SampleNames[x] <- str_match(name, "-(.*?)-") # Saves as "-Name-"
    SampleNames[x] <- substr(SampleNames[x], 2, nchar(SampleNames[x]) - 1) # Removes dashes to save as "Name" 
    x = x + 1
  }
  SampleNames <- unique(SampleNames) # Removes duplicates from list
  
  FileStrings <- list() # Creates a list to save all names in a nice structure for later use. list -> FileName -> Type -> strings...
  for (name in SampleNames) {
    inputNames <- paste("input", name, sep = "-") # Put togther the inputs: input-Name
    IPNames <- paste("IP", name, sep = "-") # And the IP's
    
    input <- intersect(list.files(pattern = inputNames), 
                       list.files(pattern = "bam")) # Now filter those strings containing the input-Name and .sam extention and save in input vector
    IP <- intersect(list.files(pattern = IPNames), list.files(pattern = "bam"))
    
    FileStrings[[name]] <- list(input, IP) # Save the two vectors to the list under the sample name
    names(FileStrings[[name]]) <- c("input", "IP") # Give the vectors the right name
  }
  setwd(directory)
  return(FileStrings) # Return for MACS2
}

# This function performs peak calling using MACS2 and the file list of all IP and input samples. It compares the inputs to IP's.
Call.Peaks <- function(File_list) {
  setwd(paste(directory, BamFiles, sep = "/"))
  for (Sample in File_list) {
    # Get the file lists in order. SampleNamesList((inputSampleNameList-R#),(IPSampleNameList-R#))
    # MACS2: -t = sample, -c = mock control, add --broad for broad peaks, -g = genome size (bp's)
    firststring <- "macs2 callpeak -t"
    laststring <- paste("-f BAMPE --broad -g ", genome_size, sep = "")
    
    input <- vector() # Initializing the vector, for the appending function later on. Clears it every loop.
    for (i in 1:length((Sample$input))) {
      input <- c(input, Sample$input[[i]]) # appending all values into the vector for the MACS2
    }
    IP <- vector()
    for (i in 1:length((Sample$IP))) {
      IP <- c(IP, Sample$IP[[i]])
    }
    
    macs2_string <- paste(c(firststring, IP, "-c", input, "-n", 
                            Sample$input[[1]], laststring), collapse = " ")
    system(paste(macs2_string))
  }
  setwd(directory)
  Move.Files(BamFiles, "peaks", PeakFiles)
}

###  Program  flow  ############################################################

setwd(directory) # Set your working directory of the files.
Create.dirs() # Create new subdirectories to save temporary files in.
Move.Files(getwd(), "\\.gz$", GunzippedFiles) # Move gunzipped files to first directory.
Unzip.Files() # Commence rest of program flow.
Trim.Files()
Make.QC.Report()
Make.Alignment()
Compress.Files()
File_list <- Make.File.List()
Call.Peaks(File_list)
