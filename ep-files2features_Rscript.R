#!/usr/bin/Rscript
library("data.table")
library("pbapply")
library("optparse")

rm(list=ls(all=TRUE))
invisible(gc())
options(warn=1)

#################### DEFINITION OF REQUIRED FUNCTIONS###################################################

# Function to read textfile content into character vector, removing punctuation and lowercasing tokens
ReadTextfile <- function(filename) {
  current.text <- scan(filename,
                       what = "char",
                       strip.white = TRUE,
                       quote = "",
                       quiet = TRUE)
  
  current.text <- current.text[intersect(grep("\\w", current.text), grep("[^”]", current.text))] # instead of matching \\w (= word character) only, match \\w  AND not ” (this quotation mark frequently appears in corpus) by intersecting regex patterns "\\w" and "[^”]" 
  current.text <- tolower(current.text)
}


GetWordlistLongFormat <- function(texts) {
  # To Do: Documentation
  cat("\nCreating Matrix of Word Counts In Long Format ...")
  # Calculate absolute word type frequencies for each text and store them in list
  freqs.abs.l <- lapply(texts, table)
  cat("\n\tfreqs.abs.l created!")
  # Create list of data frames in long format with absolute word counts for each text
  counts.dfs.l <- mapply(data.frame,
                         ID = seq_along(texts),
                         ABSFRQ = freqs.abs.l,
                         SIMPLIFY = FALSE,
                         MoreArgs = list(stringsAsFactors = FALSE))
  cat("\n\tcounts.dfs.l created!")
  rm(freqs.abs.l)
  
  # Convert list of data frames with absolute word counts to one large data frame of type data.table
  counts.dt <- rbindlist(counts.dfs.l)
  colnames(counts.dt) <- c("ID", "WORD", "FREQ_ABS") # Renaming columns
  cat("\n\tcounts.dt created!")
  rm(counts.dfs.l)
  # Relative Frequencies
  counts.dt$FREQ_REL <- counts.dt$FREQ_ABS / ave(counts.dt$FREQ_ABS, counts.dt$ID, FUN = sum)
  cat("\n  --> DONE!\n\n")
  return(counts.dt)
}



GetMostFrequentWords <- function(wordlist_long, method, parameter) {
  # To Do: Documentation
  cat("Getting Most Frequent Words Across Whole Corpus ...")
  cat("   Method: ", method, " - Parameter: ", parameter)
  mfws <- "" # Initializing character vector to be returned
  if (method == "wordlist") {
    # if external file list that exists at specified path is supplied
    wordlist.v <- scan(parameter,
                       what = "char",
                       strip.white = TRUE,
                       quote = "",
                       quiet = TRUE,
                       sep="\n")
    mfws <- wordlist.v # External list read in
  } else if (method == "relative") {
    nr.texts.c <- length(unique(wordlist_long$ID))
    cat("\n\tGenerating list of mean relative frequencies over whole corpus ...\n")
    freqs.rel.mean.df <- aggregate(FREQ_REL ~ WORD,
                                   data = wordlist_long[ , -3],
                                   function(x) { sum(x) / nr.texts.c })

    cat("\tFrom list of mean relative frequencies over whole corpus, get word types...\n")
    mfws <- as.character(freqs.rel.mean.df[freqs.rel.mean.df$FREQ_REL >= parameter, "WORD"])
    if (length(mfws) < 10) {
      #if relative criteria yield too short a list --> fallback to absolute list
      freqs_abs.df <- aggregate(FREQ_ABS ~ WORD, data = wordlist_long, sum)
      freqs_abs.df <- freqs_abs.df[order(freqs_abs.df$FREQ_ABS, decreasing = TRUE), ]
      mfws <- as.character(freqs_abs.df$WORD[1 : 100])
    }
  } else {
    freqs_abs.df <- aggregate(FREQ_ABS ~ WORD, data = wordlist_long, sum)
    freqs_abs.df <- freqs_abs.df[order(freqs_abs.df$FREQ_ABS, decreasing = TRUE), ]
    mfws <- as.character(freqs_abs.df$WORD[1 : parameter])
  }
  rm(wordlist_long)
  cat("  --> DONE identifying", length(mfws), "Most Frequent Words!\n\n")
  return(mfws)
}



CalculateSTTR <- function(text, segment.length) {
  # Function to calculate Standardised Type-Token-Ratio (STTR)
  # To Do: Documentation
  # Further details on lexical diversity measures, see:
  # http://corpora.ids-mannheim.de/libac/doc/libac-addOn-LexikalVielfalt.pdf
  
  segments.txt <- split(text, ceiling(seq_along(text) / segment.length))
  segments.vocabulary <- lapply(segments.txt, unique)
  segments.ttr <- lapply(segments.txt, function(x) {
    length(unique(x)) / length(x)
    })
  ttr.mean <- mean(unlist(segments.ttr))
  return(ttr.mean)
}


CalculateMTLD <- function(txt, threshold) {
  # Function to calculate Measure of Textual Lexical Diversity (cf. McCarthy, P.M. & Jarvis, S. (2010): MTLD, vocd-D, and HD-D: A validation study of sophisticated approaches to lexical diversity assessment. In: Behavior Research Methods 42/2, pp. 381-392. doi:10.3758/BRM.42.2.381 )
  # To Do: Documentation
  mtld.score <- 0.0 # Initializing variable to be returned
  RunMTLD <- function(txt, threshold) {
    # Inner function called within CalculateMTLD(txt, threshold)
    # To Do: Documentation
    factors <- 0
    start.token.position <- 1 # Start calculations from 1st token in text
    for (i in 1:length(txt)) {
      # Calculation from start token to current token
      ttr <- length(unique(txt[start.token.position : i])) / length(txt[start.token.position : i])
      if (ttr < threshold) {
        factors <- factors +1 # increase factor by 1 if ttr below threshold
        ttr <- 1
        start.token.position <- i + 1 # setting new start token to token after current one
      }
    }
    remainder = 1 - ttr
    excess.range = 1 - threshold
    partial.factor = remainder / excess.range
    factors <- factors + partial.factor
    return(factors)
  } # End of inner function definition
  
  mtld.forward <- RunMTLD(txt, threshold) # MTLD calculated twice: from start to end and vice versa
  mtld.reverse <- RunMTLD(rev(txt), threshold)
  val1 <- length(txt) / mtld.forward
  val2 <- length(txt) / mtld.reverse
  mtld.score <- mean(c(val1, val2))
  if (mtld.forward == 0 | mtld.reverse == 0) {
    mtld.score <- 0
  }
  return(mtld.score)
}


ttr_4chars <- function(txt, threshold) {
  # This function is currently not required in conversion of ep-files to feature representation
  # To Do: Documentation
  ttrs <- c(1:length(txt))
  start.token.position <- 1
  for (i in 1:length(txt)) {
    ttr <- length(unique(txt[start.token.position:i]))/length(txt[start.token.position:i])
    cat(txt[i], ":", ttr, "\n", sep=" ")
    ttrs[i] <- ttr
  }
  return(ttrs)
}


GetDocumentTermMatrixAbsolute <- function(wordlist.long, mostFrequentWords) {
  # Function for creation of Document-Term-Matrix (raw counts) from list of text vectors and set of most frequent words
  # To Do: Documentation
  cat("Creating Document-Term-Matrix of Absolute Frequencies ...")
  counts.mfwsonly.long.dt <- wordlist.long[wordlist.long$WORD %in% mfws.c,]
  # Reshape data frame of counts from long to wide format
  counts.mfwsonly.wide.dt <- dcast(counts.mfwsonly.long.dt,
                                   ID ~ WORD,
                                   value.var = "FREQ_ABS",
                                   drop=TRUE)
  rm(counts.mfwsonly.long.dt)
  # Identify texts, which have no word counts at all and therefore do not appear in data frame after conversion to wide format
  missingtexts.v <- setdiff(c(1:length(textlengths.tokens.v)), counts.mfwsonly.wide.dt$ID)
  if (length(missingtexts.v) > 0) {
    missingtexts.m <- matrix(rep(NA, length(missingtexts.v) * (length(colnames(counts.mfwsonly.wide.dt)))),
                             nrow = length(missingtexts.v),
                             dimnames = list(missingtexts.v,colnames(counts.mfwsonly.wide.dt)))

  missingtexts.m[ , "ID"] <- missingtexts.v
  counts.mfwsonly.wide.dt <- rbind(counts.mfwsonly.wide.dt, missingtexts.m)
  counts.mfwsonly.wide.dt <- counts.mfwsonly.wide.dt[order(counts.mfwsonly.wide.dt[,"ID"]),]
  }
  colnames(counts.mfwsonly.wide.dt) <- toupper(colnames(counts.mfwsonly.wide.dt))
  colnames(counts.mfwsonly.wide.dt)[-1] <- 
    paste("FREQ-ABS_",
          colnames(counts.mfwsonly.wide.dt)[-1],
          sep = "")
  cat("  --> DONE!\n\n")
  return(counts.mfwsonly.wide.dt)
}


GetDocumentTermMatrixTFIDF <- function(dtm_rel) {
  # To Do: Documentation
  cat("Creating  TF-IDF Weighted Document-Term-Matrix ...  ")
  idfs <- log(nrow(dtm_rel) / colSums(dtm_rel > 0, na.rm = TRUE)) # IDF(t) = log(Total number of documents / Number of documents with term t in it)
  for (i in seq_along(dtm_rel)) {
    set(dtm_rel,
        i = which(is.na(dtm_rel[[i]])),
        j = i, value = 0)
  }
  tfidf.m <- as.matrix(dtm_rel) %*% diag(idfs)
  colnames(tfidf.m) <- colnames(dtm_rel)
  rownames(tfidf.m) <- c(1:nrow(tfidf.m))
  tfidf.m <- cbind(tfidf.m, ID = rownames(tfidf.m))
  tfidf.m <- tfidf.m[,c(ncol(tfidf.m), 1 : (ncol(tfidf.m) - 1))] # move last column (ID) to first
  colnames(tfidf.m) <- gsub("FREQ-REL_",
                            "TFIDF_",
                            colnames(tfidf.m))
  cat("  --> DONE!\n\n")
  return(tfidf.m)
}


CalculateMeanWordRank <- function(wordlist.long) {
  # Caluclates mean word rank using the whole corpus as a reference
  # To Do: Documentation
  # To Do: Option to provide frequency list from external reference corpus
  counts.dt <- wordlist.long
  mwr <- 0 # Initialize variable to be returned
  cat("Calculating Mean Word Rank for Each Text in Corpus ... ")
  # Calculate word type frequencies (absolute counts) over whole corpus
  freqs.abs.t <- sort(xtabs(FREQ_ABS ~ WORD, data = wordlist.long.dt),
                       decreasing = TRUE)
  # Rank words by absolute frequency without gaps after ties ("dense rank")
  ranks.v <- frankv(-freqs.abs.t,
                    ties.method = "dense")
  names(ranks.v) <- names(freqs.abs.t)
  # Since counts.dt shows types and not tokens, rank of each type needs to be multiplied by its frequency
  counts.dt$RANKSUM <- ranks.v[as.character(counts.dt$WORD)] * counts.dt$FREQ_ABS
  # Finally, mean word rank for each text can be computed by dividing a text's agregate ranksum by total number of tokens in the text 
  ranks_aggregate <- aggregate(RANKSUM ~ ID, data = counts.dt, sum)
  tokens_total <- aggregate(FREQ_ABS ~ ID, data = counts.dt, sum)
  mwr <- ranks_aggregate / tokens_total
  mwr <- mwr$RANKSUM
  rm(counts.dt)
  cat("  --> DONE!\n\n")
  return(mwr)
}


ConfigureMFWExtraction <- function(method, parameter, lang) {
  # Parameter supplied to selected MFW extraction method, corresponding to
  #  1) threshold for the mean relative frequency of words across whole corpus if selected method is 'relative', or
  #  2) n most frequent words in sorted frequency list across whole corpus if selected method is 'absolute', or
  #  3) path to an external wordlist in plain text (one word per line) if
  #     selected method is 'wordlist'.
  
  if (method == "relative") { # if method = relative, set parameter to default value if parameter contains letters or is >= 1
    if (grepl("[:alpha:]", parameter)) {
      parameter <- 0.001
    } else if (as.numeric(parameter >= 1)) {
      parameter <- 0.001
    }
    parameter <- as.numeric(parameter)
  } else if (method == "absolute") { # if method = absolute, set parameter to default if parameter contains letters or if < 1
    if (grepl("[:alpha:]", parameter)) {
      parameter <- 100
    } else if (as.numeric(parameter) < 1) {
      parameter <- 100
    }
    parameter <- as.numeric(parameter)
  } else { # if method = wordlist, set method to relative with default threshold if no .txt files found in wordlist directory
    wordlist.path <- trimws(paste0(parameter, "/wordlist_", tolower(lang), ".txt"))
    if (!file.exists(wordlist.path)) {
      method <- "relative"
      parameter <- 0.001
    } else {
      parameter <- wordlist.path
    }
  }
  return(list(method, parameter))
}


########################################################### END OF FUNCTION DEFINITIONS

########################################################### PARSE CLI INPUT

lang.choices <- c("BG", "CS", "DA", "DE", "EL", "EN", "ES", "ET", "FI",
                      "FR", "HU", "IT", "LT", "LV", "NL", "PL", "PT", "RO",
                      "SK", "SL", "SV")
mfw.choices <- c("relative", "absolute", "wordlist")

cli.options <- list(
  make_option(c("-i", "--inputDir"),
              type = "character", default = NULL, 
              help = "Path to input directory containing tokenised EuroParl files",
              metavar = "path"),
  make_option(c("-o", "--outputDir"),
              type = "character",
              default = "analyses/", 
              help = "Path to output directory [default %default]",
              metavar = "path"),
  make_option(c("-l", "--languages"), type = "character", 
              help = paste0("One or more languages, choose from {",
                           paste(lang.choices, collapse="|"),
                           "}, comma-separated without blanks"),
              metavar = "choice"),
  make_option(c("-m", "--mfwMethod"),
              type = "character",
              default = "relative",
              help = paste0("Method to Calculate MFWs, choose from [",
                            paste(mfw.choices, collapse = "|"),
                            "]"),
              metavar = "choice"),
  make_option(c("-p", "--mfwParameter"),
              type = "character",
              default = "0.001",
              help = "Parameter supplied to chosen method for MFW calulation: threshold of mean relative frequency if method is 'relative' [default 0.001], number of words if method is 'absolute' [default 100], or path to folder with wordlists (wordlist_xy.txt) if method is 'wordlist'",
              metavar = "character"),
  make_option(c("-t", "--tfidf"),
              action = "store_true",
              default = "FALSE",
              help = "Use TF-IDF weighting instead of relative frequency counts in document-term-matrix of Most Frequent Words [default %default].",
              metavar = "boolean"))


# $ Rscript analyses/ep-files2features.R -i output_test/ -o testfiles_mini/ -m relative -p 7 -l DE,PL,FR

opt_parser = OptionParser(usage = "%prog [options]",
                          option_list = cli.options,
                          prog = "ep-files2features.R",
                          description = "\nScript for converting tokenised EuroParl files to feature representations in matrix format.",
                          epilogue = "Usage example: ep-files2features.R -i inputs/ -o outputs/\n")

cli.args = parse_args(opt_parser)

# Terminate program if agruments -i or -l are not speficied
if (is.null(cli.args$inputDir) | is.null(cli.args$languages)) {
  print_help(opt_parser)
  stop("Mandatory argument -inputFolder and/or -languages is missing! See help page (--help)")
}

# Parse path to input folder and terminate program if folder does not exist
indir <- cli.args$inputDir
if (!dir.exists(indir)) {
  print_help(opt_parser)
  stop("Specified input folder does not exist! See help page (--help)")
}

# Make sure output folder exists
outdir <- cli.args$outputDir
if (!dir.exists(outdir)) {
  dir.create(outdir)
}

# Parse argument -l and terminate program if unsupported language specified
langs <- unlist(strsplit(cli.args$languages, ","))
if (any(!(langs %in% lang.choices))) { # if any of specified languages not in lang.choices
  print_help(opt_parser)
  stop("One or more unsupported languages selected! See help page (--help)")
}

# Parse argument -m that specifies method for calculation of Most Frequent Words
# and terminate program if unsupported method specified
mfw.method <- cli.args$mfwMethod
if (!mfw.method %in% mfw.choices) { # if any of specified languages not in lang.choices
  print_help(opt_parser)
  stop("Unsupported MFW method selected! See help page (--help)")
}

mfw.parameter <- cli.args$mfwParameter

tfidf.weighting <- cli.args$tfidf

cat("\n~~~~ PARSING CLI ARGUMENTS ~~~~\n\n")
cat("Indir\t\t", indir, "\n")
cat("Outdir\t\t", outdir, "\n")
cat("Languages\t", langs, "\n")
cat("MFWMethod\t", mfw.method, "with parameter: ", mfw.parameter, "\n")
cat("TFIDF\t\t", tfidf.weighting, "\n")


#langs <- c("PL") #DE# languages to be converted to feature representations # 2DO: Read in language codes interactively
#lang <- "PL" # make this in loop
cat("\nRetrieving list of EuroParl files in specified input folder ...  ")
files.alllanguages.v <- dir(indir, recursive = T,  full.names = T) #"../output_full/comparable"
cat("  --> DONE,", length(files.alllanguages.v), "files in input folder\n\n")

for (lang in langs) {

  # Create vector of all filenames in input folder
  # From vector of filenames files.alllanguages.v, select chosen language as specified by language code only 
  files.v <- files.alllanguages.v[grep(paste0("/", lang, "/"),
                                       files.alllanguages.v)]

  
  cat("###########################   ", lang, "   ###########################\n\n")
  cat(length(files.v), "files in input folder for", lang, "\n\n")
  
  if (length(files.v) < 1) {
    next
  }
  
  # Apply ReadTextfile(fn) to each file to return a list, where each text is
  # the lowercased sequence of tokens (punctuation removed)
  cat("Reading ", length(files.v), "text files into memory ...\n") ## Read text files into list of text vectors for each text
  texts.l <- pblapply(files.v, function(filename) ReadTextfile(filename)) # if no progress bar is required, use lapply instead of pblapply
  # Removing empty texts from texts.l
  cleanup.v <- sapply(texts.l, function(x) { length(x) > 0 })
  files.v <- files.v[cleanup.v]
  texts.l <- texts.l[cleanup.v]
  
  cat("   --> DONE,", length(texts.l), "read into memory\n\n")
  # Count number of tokens and types in each text
  textlengths.tokens.v <- sapply(texts.l,
                                 length,
                                 simplify = TRUE)
  textlengths.types.v <- sapply(texts.l,
                                function(text.x) { length(unique(text.x)) })
  
  cat("\nCalculating Standardised Type-Token-Ratio for Each Text  ...  ")
  sttr.v <- unlist(lapply(texts.l,
                          function(txt) { CalculateSTTR(txt, 100) }))
  cat("  --> DONE!\n\n\n")
  
  cat("Calculating Measure of Textual Lexical Diversity for Each Text ... \n")
  mtld.l <- pblapply(texts.l, function(txt) { CalculateMTLD(txt, 0.72)}) # ORIG unlist(pblapply(texts.l, function(txt) { CalculateMTLD(txt, 0.72)}))
  cat("   --> DONE!\n\n")
  
  
  # Convert list of texts to long format data table of
  # 3 columns (Text-ID, Word Type, Count)
  wordlist.long.dt <- GetWordlistLongFormat(texts.l)
  rm(texts.l)
  
  # Retrieve most frequent words in corpus
  mfw.config.method <- ConfigureMFWExtraction(mfw.method, mfw.parameter, lang)[[1]]
  mfw.config.parameter <- ConfigureMFWExtraction(mfw.method, mfw.parameter, lang)[[2]]
  mfws.c <- GetMostFrequentWords(wordlist.long.dt, mfw.config.method, mfw.config.parameter) # could be optimised by passing long format word list instead of list of texts as argument

  # Create Document-Term-Matrix with absolute frequency counts for each of the most frequent words
  dtm.raw.dt <- GetDocumentTermMatrixAbsolute(wordlist.long.dt, mfws.c)

  # From Document-Term-Matrix of raw counts, create DTM of relative frequency counts
  cat("Calculating Document-Term-Matrix of Relative Frequencies ... ")
  dtm.rel.dt <- dtm.raw.dt[,-1] / textlengths.tokens.v
  colnames(dtm.rel.dt) <- gsub("FREQ-ABS_", "FREQ-REL_", colnames(dtm.rel.dt))
  cat("  --> DONE!\n\n")
  
  
  # From Document-Term-Matrix of relative counts, creat TF-IDF weighted DTM
  if (tfidf.weighting == TRUE) {
    dtm.tfidf.m <- GetDocumentTermMatrixTFIDF(dtm.rel.dt)
  }
  
  
  
  # Calculate mean word ranks
  mean_word_ranks.v <- CalculateMeanWordRank(wordlist.long.dt)
  
  
  # VERSCHIEBEN?
  # For column 'Corpus Type', create vector of length of files.v with constant value 'comp'
  corpustype.v <- c(rep("comp", length(textlengths.tokens.v)))
  
  
  # VERSCHIEBEN?
  # For column 'Corpus Section', create vector indicating whether filename is from (non-)translated subsection
  corpussection.v <- grepl("non-translated", files.v) # TRUE at each position where filelist matches 'non-translated'
  length(corpussection.v)
  corpussection.v <- gsub("TRUE", "O", corpussection.v)
  corpussection.v <- gsub("FALSE", "T", corpussection.v)
  corpussection.v <- corpussection.v[cleanup.v] # TEST THIS
  
  # Get source language SL from filename
  # First, extract 2 caps before optional - optional cap optional cap followed by /digit (digit marks beginning of base file name)
  # Second, use substr() to strip extacted string after - and before /digit
  sl.c <- regmatches(files.v, regexpr("([A-Z]{2}-?[A-Z]?[A-Z]?/[0-9])", files.v))
  sl.c <- substr(regmatches(files.v, regexpr("([A-Z]{2}-?[A-Z]?[A-Z]?/[0-9])", files.v)),1, 2)
  sl.c <- sl.c[cleanup.v]
  
  # Create vector indicating target language (PL) whose length corresponds to number of files in PL comparable subcorpus
  tl.c <- rep("PL", length(textlengths.tokens.v))
  
  ### Get year from filenames
  year.c <- regmatches(files.v, regexpr("(/[0-9]{2})", files.v)) ##regmatches(files.v, regexpr("([0-9]{2}-?[A-Z]?[A-Z]?/[0-9])", files.v))
  year.c <- substr(year.c, 2,3)
  year.c <- sapply(year.c, function(year) {
    if (grepl("^9", year)) { paste("19", year, sep="") }
    else { paste("20", year, sep="") }
  })
  year.c <- year.c[cleanup.v]
                                                                                                                                                                                                                                                                                          
  # Create Data Frame
  cat("Writing Data Frame ...  ")
  features.comp.df <- data.frame(TYPE = corpustype.v,
                                 SECTION = corpussection.v,
                                 SL = sl.c,
                                 TL = tl.c,
                                 YEAR = year.c,
                                 FILENAME = files.v,
                                 LEN_TOK = textlengths.tokens.v,
                                 LEN_TYPES = textlengths.types.v,
                                 TTR = textlengths.types.v / textlengths.tokens.v, # features.comp.df$TTR <- apply(features.comp.df, 1, function(row.x) (strtoi(row.x[8])/strtoi(row.x[7])))
                                 STTR = sttr.v,
                                 MTLD = unlist(mtld.l),
                                 MEANWORDRANK = mean_word_ranks.v)
                                # missing: STTR & DTM
  if (tfidf.weighting == FALSE) {
    features.comp.df <- cbind(features.comp.df, dtm.rel.dt)
  } else {
    features.comp.df <- cbind(features.comp.df, dtm.tfidf.m)
    }         
                            
  cat("  --> DONE!")
  write.csv(features.comp.df,
            file = paste0(outdir, "/ep_comparable_", tolower(lang), ".csv"))
  
  #write.csv(features.comp.df,
  #          file = paste0("ep_comparable_"pl.csv"))
  #write.csv(wordvectors.df, file="wordvectors.csv")
  #write.table(features.comp.df, file="ep_comparable_pl.csv", sep="\t")
  save(features.comp.df,
       file = paste0(outdir, "/ep_comparable_", tolower(lang), ".RData"))
  
  if(length(langs) > 1) {
    cat("Cleaning up environment prior to starting next language ...  ")
    rm(wordlist.long.dt)
    rm(dtm.raw.dt)
    rm(dtm.rel.dt)
    if (tfidf.weighting == TRUE) {
      rm(dtm.tfidf.m)
    }
    rm(mtld.l)
    rm(features.comp.df)
    cat("  --> DONE!\n\n")
    gc()

  }
  cat("\n\n\tEXTRACTION FOR LANGUAGE", lang, "COMPLETED!\n\n")
}



