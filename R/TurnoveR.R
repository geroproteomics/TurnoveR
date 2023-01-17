#Test Change
#test change 2
#Written by Cameron Wehrfritz
#and Natan Basisty, PhD
#Schilling Lab, Buck Institute for Research on Aging
#Novato, California, USA
#March, 2020
#updated: January 28, 2021


# PROTEIN TURNOVER ANALYSIS
# STEP 1
# PROCESS DATA FROM SKYLINE
# CORRECT FOR THE NATURALLY OCCURING HEAVY ISOTOPES OF COMMON ELEMENTS
# CALCULATE PRECURSOR POOL AND AVERAGE TURNOVER SCORE
# FILTER DATA (OPTIONAL)

# OUTPUT: DATA TABLES AND PDFs OF PLOTS

######################
#### Begin Program ###
######################



#------------------------------------------------------------------------------------
# START CODE FOR RUNNING IN RSTUDIO (comment out if running from TurnoveR)
#------------------------------------------------------------------------------------


#
# #------------------------------------------------------------------------------------
# #set working directory
# setwd("/Volumes/GibsonLab/users/Cameron/2020_0814_Skyline_Turnover_Tool/Turnover_R_scripts") # VPN mac
# setwd("//bigrock/GibsonLab/users/Cameron/2020_0814_Skyline_Turnover_Tool/Turnover_R_scripts") # VPN windows
#setwd("C:/Users/olingerba/Desktop/TurnoveR Rework/test4/")
#setwd("C:/Users/olingerba/OneDrive - National Institutes of Health/Desktop/TurnoveR Rework/test5")
#setwd("C:/Users/olingerba/OneDrive - National Institutes of Health/Desktop/TurnoveR Rework/package2_part4")
setwd("C:/Users/olingerba/OneDrive - National Institutes of Health/Desktop/TurnoveR Rework/Nick")
# #------------------------------------------------------------------------------------
#

#libraries
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)
library(purrr)
library(reshape2)
library(gridExtra)
library(forcats)
library(pracma)
library(seqinr)
library(hablar)
library(regexPipes)
library(rlang)
library(tidyverse)

#test run of TuronoveR
TurnoveR(df.input, AA = "leucine", total_HE = 3, HE = "H", num_HE = 3, second_element = FALSE, HE2 = NULL, num_HE2 = 0,
 AA_B = NA, total_HE_B = 0, HE_B = NULL, num_HE_B = 0, second_element_B = FALSE,  HE2_B = NULL, num_HE2_B = 0,
 NIAC = TRUE, PP_Enrich = TRUE, diet.enrichment = 0.999999, Detection.Qvalue.threshold = 1, ATS.threshold=0.7,
 reference.cohort = "OCon", IDP.threshold = 0)

###
###
#Main Function-------------------------------------------------------------------
TurnoveR <- function(df.input, AA, total_HE, HE, num_HE, second_element = TRUE, HE2 = NULL, num_HE2 = 0,
                     AA_B = NA, total_HE_B = 0, HE_B = NULL, num_HE_B = 0, second_element_B = TRUE,  HE2_B = NULL, num_HE2_B = 0,
                     NIAC = TRUE, PP_Enrich = TRUE, diet.enrichment = 0.999999, Detection.Qvalue.threshold = 1, ATS.threshold=0.7,
                     reference.cohort = "OCon", IDP.threshold = 0){

  #turn off isotopic envelope deconvolution if the mass difference between the natural and heavy AA is 4+
  if(total_HE >= 4){
    NIAC = FALSE
  }
  #remove_AA_sym = FALSE, remove_num_HE = NULL,

  #check for second heavy labeled amino acid; set cycle count variable to two if second label is present
  counter_AA <- 1
  if(!is.na(AA_B)){
    counter_AA <- 2
  }

  #Iterate first through AA, then again for AA_B
  for (val2 in 1:counter_AA){

    #for first iteration, keep vals the same, assign AA_B to be removed
    if (val2 == 1){
      if (!is.na(AA_B)){
        remove_AA_sym = AA_B
        remove_num_HE = total_HE_B
      }
    }

    #for second iteration, change AA_B values to regular values, calculate turnover using these instead and remove AA
    if (val2 == 2){
      if (!is.na(AA_B)){
        #remove AA from dataset
        remove_AA_sym = AA
        remove_num_HE = total_HE
        #substitute all values with second heavy AA
        AA = AA_B
        total_HE = total_HE_B
        HE = HE_B
        num_HE = num_HE_B
        second_element = second_element_B
        HE2 = HE2_B
        num_HE2 = num_HE2_B
      }
    }



    # #------------------------------------------------------------------------------------
    # # Set Default Values
    #
    # ## these may be user defined in the future
    min.abundance <- 10**(-5) # minimum abundance
    resolution <- 0.1 # resolution for distinguishing peaks
    p.tolerance <- 0.05 # tolerance for combining masses in observed data
    Detection.Qvalue.threshold <- 1 # value for filtering on Detection.Qvalue where 0 is the most stringent and 1 is least stringent
    #
    # # diet enrichment
    diet.enrichment <- 99.9999 # percent enrichment of heavy AA in diet - Update to user specified value
    # # if user specifies diet.enrichment of 100%, change to 99.9999% (since 100% will not work in FBC step)
    # diet.enrichment <- ifelse(diet.enrichment==100, 99.9999, diet.enrichment)
    diet.enrichment <- diet.enrichment/100 # transform percent diet enrichment from 0-100 % to 0.0 - 1.0
    # #------------------------------------------------------------------------------------


    #------------------------------------------------------------------------------------
    # END CODE FOR RUNNING IN RSTUDIO
    #------------------------------------------------------------------------------------




    #------------------------------------------------------------------------------------
    # Preliminary Filters and Cleaning:

    df <- df.input %>%
      filter(Is.Decoy == "FALSE") %>% # filter out Decoys (which are used for training algorithm in Skyline)
      filter(!Protein=="Biognosys|iRT-Kit_WR_fusion") %>% # filter out Biognosys rows, since these are spiked in to the sample for Quality Control
      filter(!Fragment.Ion=="precursor [M-1]") %>% # filter out [M-1] precursor observations (since these cause an issue with building Matrix A in the FBC step)
      mutate_at(vars(Timepoint), list(~as.numeric(.))) %>% # convert Timepoint variable to numeric
      mutate_at(vars(Detection.Q.Value), list(~as.numeric(.))) %>% # convert Detection.Q.Value variable to numeric
      mutate_at(vars(Detection.Q.Value), list(~ifelse(is.na(.), 0.00123, .))) %>% # FOR TESTING PURPOSES ONLY -- if Qvalue is missing replace NA with value=0.00123 -- TO BE REMOVED IN OFFICIAL TOOL -- FOR TESTING PURPOSES ONLY
      filter(!is.na(Detection.Q.Value)) # filter out observations which do not have numerical Detection.Q.Value
    #------------------------------------------------------------------------------------

    #------------------------------------------------------------------------------------
    #bao change: remove all AAs containing arginine, so that testing can be done with just heavy leucine
    #string2 = "hey"
    #string1 = "heythere"
    #grepl(string1,string2, fixed=TRUE)


    #conduct following only if second heavy labeled amino acid is present and needs to be removed for first roung
    if (!is.na(AA_B)){
      ###
      ###
      #for loop to cycle through each row, find rows with heavy arginine, then remove these from df

      #convert AA to be removed to AA symbol
      if (paste(remove_AA_sym) == "leucine"){
        remove_AA_sym = "L"
      }
      if (paste(remove_AA_sym) == "valine"){
        remove_AA_sym = "V"
      }
      if (paste(remove_AA_sym) == "lysine"){
        remove_AA_sym = "K"
      }
      if (paste(remove_AA_sym) == "arginine"){
        remove_AA_sym = "R"
      }

      #define expression to locate heavy AA to be removed
      remove_expression <- paste(paste(remove_AA_sym), "[+",paste(remove_num_HE),"]", sep="")
      remove_expression

      row_list = c()
      for (i in 1:nrow(df)){
        if (grepl(df[i,"Modified.Sequence"],remove_expression,fixed=TRUE)){
          row_list = c(row_list,i)
        }
        #bool_list = c(bool_list,!grepl(df[i,14],"R[+10]",fixed=TRUE))
      }
      #use row_list to remove all rows containing R[+10] from df
      df = df[-c(row_list),]
    }#end check for second heavy-labeled amino acid

    ###select first 4 proteins just while testing pipeline, remove this after testing
    protein_4 <- unique(df$Protein.Gene)[1:4]
    #subset df for only first 4 proteins to speed up analysis
    row_list2 = c()
    for (i in 1:nrow(df)){
      if (df[i,"Protein.Gene"] %in% protein_4){
        row_list2 = c(row_list2,i)
      }
      #bool_list = c(bool_list,!grepl(df[i,14],"R[+10]",fixed=TRUE))
    }
    #use row_list to keep only rows with first 4 proteins
    df = df[c(row_list2),]
    #remove line 167-179 after testing
    ###
    ###
    #------------------------------------------------------------------------------------
    #end filtering out arginine

    # # need updated test data set with numerical Qvalues before using this chunk
    #------------------------------------------------------------------------------------
    # Detection Qvalue Filter
    # Detection.Qvalue.threshold can be between [0,1) where smaller values are more stringent
    # The default should be 1, corresponding to no filter, thereby retaining all of the data
    df <- df %>%
      filter(Detection.Q.Value < Detection.Qvalue.threshold)
    #------------------------------------------------------------------------------------

    #-----------------------------------------------------
    #BAO change:
    #1)function to select amino acid string
    #2)if statements to select amino acid symbol
    #Intake user defined heavy amino acid
    #take input for amino acid
    #take input for mass increase of heavy AA (example: 3 in indicates normal heavy leucine with +3 deuterium, 6 indicates L-arg:HCL which has 6 heavy carbons)
    #need modifications for this in order to take user input from skyline


    #create amino acid symbol based on amino acid
    if (paste(AA) == "leucine"){
      AA_symbol = "L"
    }
    if (paste(AA) == "valine"){
      AA_symbol = "V"
    }
    if (paste(AA) == "lysine"){
      AA_symbol = "K"
    }
    if (paste(AA) == "arginine"){
      AA_symbol = "R"
    }


    #define heavy isotopes from user-selected elements
    if (HE == "H"){
      HE_iso = "D"
    }
    if (HE == "C"){
      HE_iso = "X"
    }
    if (HE == "O"){
      HE_iso = "Y"
    }
    if (HE == "N"){
      HE_iso = "Z"
    }


    #repeat HE_iso selection from user-defined heavy element (only if second element is selected)
    if (second_element == TRUE){
      if (HE2 == "H"){
        HE_iso2 = "D"
      }
      if (HE2 == "C"){
        HE_iso2 = "X"
      }
      if (HE2 == "O"){
        HE_iso2 = "Y"
      }
      if (HE2 == "N"){
        HE_iso2 = "Z"
      }
      #end name second heavy isotope
    }



    #------------------------------------------------------------------------------------
    # Set Default Values for Mass Distributions
    # Elements: Hydrogen, Carbon, Oxygen, Nitrogen, Sulfur

    # Create Tibble with Elements, masses and abundances
    tbl.elements <- tibble("Element"=c("Hydrogen", "Hydrogen", "Carbon", "Carbon", "Oxygen", "Oxygen", "Oxygen", "Nitrogen", "Nitrogen", "Sulfur", "Sulfur", "Sulfur", "Sulfur"),
                           "Symbol"=c("H", "H", "C", "C", "O", "O", "O", "N", "N", "S", "S", "S", "S"),
                           "Mass"= c(1.0078246, 2.0141021, 12.000000, 13.0033554, 15.9949141, 16.9991322, 17.9991616, 14.0030732, 15.0001088, 31.972070, 32.971456, 33.967866, 35.96708),
                           "Abundance"= c(0.999855,  0.000145, 0.98916,  0.01084, 0.997576009706, 0.000378998479, 0.002044991815,  0.99633, 0.00366, 0.95021, 0.00745, 0.04221, 0.00013)
    )


    #bao change: add row to tbl.elements depending on heavy element, also initialize variable to hold molecular weight of the element for later (line 883)
    if (HE == "H"){
      # Add Deuterium with user specified amount for diet enrichment
      tbl.elements <- add_row(tbl.elements, Element="Deuterium", Symbol="D", Mass=1.0078246, Abundance=1-diet.enrichment) %>% # hydrogen
        add_row(Element="Deuterium", Symbol="D", Mass=2.0141021, Abundance=diet.enrichment) # deuterium
      HE_weight = 1.0078246
    }
    if (HE == "C"){
      # Add Heavy Carbon with user specified amount for diet enrichment
      tbl.elements <- add_row(tbl.elements, Element="hCarbon", Symbol="X", Mass=12.000000, Abundance=1-diet.enrichment) %>% # carbon
        add_row(Element="hCarbon", Symbol="X", Mass=13.0033554, Abundance=diet.enrichment) # heavy carbon
      HE_weight = 12.000000
    }
    if (HE == "O"){
      # Add Heavy Oxygen with user specified amount for diet enrichment
      tbl.elements <- add_row(tbl.elements, Element="hOxygen", Symbol="Y", Mass=15.9949141, Abundance=1-diet.enrichment) %>% # oxygen
        add_row(Element="hOxygen", Symbol="Y", Mass=16.9991322, Abundance=diet.enrichment) # heavy oxygen
      HE_weight = 15.9949141
    }
    if (HE == "N"){
      # Add Heavy Oxygen with user specified amount for diet enrichment
      tbl.elements <- add_row(tbl.elements, Element="hNitrogen", Symbol="Z", Mass=14.0030732, Abundance=1-diet.enrichment) %>% # oxygen
        add_row(Element="hNitrogen", Symbol="Z", Mass=15.0001088, Abundance=diet.enrichment) # heavy nitrogen
      HE_weight = 14.0030732
    }

    #repeat addition of row for second heavy element (if present)
    if (second_element == TRUE){

      #bao change: add row to tbl.elements depending on heavy element, also initialize variable to hold molecular weight of the element for later (line 883)
      if (HE2 == "H"){
        # Add Deuterium with user specified amount for diet enrichment
        tbl.elements <- add_row(tbl.elements, Element="Deuterium", Symbol="D", Mass=1.0078246, Abundance=1-diet.enrichment) %>% # hydrogen
          add_row(Element="Deuterium", Symbol="D", Mass=2.0141021, Abundance=diet.enrichment) # deuterium
        HE_weight = 1.0078246
      }
      if (HE2 == "C"){
        # Add Heavy Carbon with user specified amount for diet enrichment
        tbl.elements <- add_row(tbl.elements, Element="hCarbon", Symbol="X", Mass=12.000000, Abundance=1-diet.enrichment) %>% # carbon
          add_row(Element="hCarbon", Symbol="X", Mass=13.0033554, Abundance=diet.enrichment) # heavy carbon
        HE_weight = 12.000000
      }
      if (HE2 == "O"){
        # Add Heavy Oxygen with user specified amount for diet enrichment
        tbl.elements <- add_row(tbl.elements, Element="hOxygen", Symbol="Y", Mass=15.9949141, Abundance=1-diet.enrichment) %>% # oxygen
          add_row(Element="hOxygen", Symbol="Y", Mass=16.9991322, Abundance=diet.enrichment) # heavy oxygen
        HE_weight = 15.9949141
      }
      if (HE2 == "N"){
        # Add Heavy Oxygen with user specified amount for diet enrichment
        tbl.elements <- add_row(tbl.elements, Element="hNitrogen", Symbol="Z", Mass=14.0030732, Abundance=1-diet.enrichment) %>% # oxygen
          add_row(Element="hNitrogen", Symbol="Z", Mass=15.0001088, Abundance=diet.enrichment) # heavy nitrogen
        HE_weight = 14.0030732
      }

    }

    # create nested list-column from this tibble for each element
    ntbl.elements <- tbl.elements %>%
      group_by(Element, Symbol) %>%
      nest()
    #------------------------------------------------------------------------------------


    # FUNCTIONS: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    #------------------------------------------------------------------------------------
    # Small and Simple functions:

    # function to normalize vector; sum to one; as in probabilities
    sum.to.one.fun <- function(x){x/sum(x)}

    # function to normalize vector; as in unit vector
    normalize.fun <- function(x){x / sqrt(sum(x^2))}

    # function to count number of heavy-labeled leucines
    leucine.tally.fun <- function(x){
      tally <- length(unlist(str_split(x, "L\\[\\+3\\]", simplify=FALSE)))-1
    }

    # function to count number of leucines
    leucine.count.fun <- function(x){
      count <- length(unlist(str_split(x, "L", simplify=FALSE)))-1
    }

    # function to strip away the leucine modification but leave all other modifications
    # when L[+3] occurs replace it with L
    peptide.modification.fun <- function(x){
      new.sequence <- paste0(unlist(str_split(x, "L\\[\\+3\\]", simplify=FALSE)), collapse="L")
    }

    #BAO change:
    #3)function to count number of valine, lysine, arginine
    #4)function to count the number of heavy AA
    #4)function to create string splitting pattern based on chosen amino acid
    #5)function to strip away AA modifications


    #create string splitting pattern as string
    #can update this later to change the number based on added weight to heavy AA
    split_pattern <- paste(paste(AA_symbol),"\\[\\+",paste(total_HE),"\\]", sep="")

    # function to count number of heavy-labeled amino acid
    AA.heavy.tally.fun <- function(x){
      tally <- length(unlist(str_split(x, paste(split_pattern), simplify=FALSE)))-1
    }

    # function to count number of amino acid
    AA.count.fun <- function(x){
      count <- length(unlist(str_split(x, paste(AA_symbol), simplify=FALSE)))-1
    }

    # function to strip away the AA modification but leave all other modifications
    # when L[+3] occurs replace it with L
    aa.modification.fun <- function(x){
      new.sequence <- paste0(unlist(str_split(x, paste(split_pattern), simplify=FALSE)), collapse=paste(AA_symbol))
    }


    # function to count number of charges (plus signs "+")
    charge.count.fun <- function(x){
      charge.count <- length(unlist(str_split(x, "\\+", simplify=FALSE)))-1
    }
    #------------------------------------------------------------------------------------



    # Larger more complex functions:

    #------------------------------------------------------------------------------------
    # Mass Distribution - Add Function
    # Function adds two mass distributions together by adding masses together and multiplies abundances.

    add.fun <- function(self, addition){
      map <- data.frame(Mass=as.numeric(), Abundance=as.numeric()) # initiate data frame to write 'mass' and 'abundance' to
      for(mass1 in seq_along(self[["Mass"]])){
        for(mass2 in seq_along(addition[["Mass"]])){
          mass <- self$Mass[mass1] + addition$Mass[mass2] # new mass
          if(mass %in% map$Mass){   # update frequency if it already exists in map for this mass
            map$Abundance[map$Mass %in% mass] = map$Abundance[map$Mass %in% mass] + self$Abundance[mass1]*addition$Abundance[mass2] # update frequency for when it previously exists
          } else{ # calculate frequency and add new row if this mass does not already exist in map
            frequency <- self$Abundance[mass1]*addition$Abundance[mass2] # product of abundances
            map <- add_row(map, Mass=mass, Abundance=frequency) # append mass and frequency to the map tibble
          }
        } #end for
      } #end for
      return(map)
    } #end function
    #------------------------------------------------------------------------------------


    #------------------------------------------------------------------------------------
    # Mass Distribution - Multiply
    # Function returns the result of adding this Mass Distribution to itself the specified number of times.
    # i.e. get the mass distribution of C27 from C.

    multiply.fun <- function(self, factor){
      if(factor==0){
        return(result) # Return mass distribution; this should be the terminating condition
      }

      if(factor==1){
        result <- self
        return(result)
      }

      result <- add.fun(self, self) # calls the addition function to build result onto itself (eg C2 becomes C4, C4 becomes C8, etc)
      # note: 'self' becomes 'result' and for large factors 'result' continues to build onto 'result'
      if(factor>=4){
        result <- multiply.fun(result, floor(factor/2)) # recursive; reduces factor by two and calls function again
      }

      if(factor%%2 !=0){ # if factor is odd, we still have to add one more self to result
        result <- add.fun(self, result)
        return(result)
      }

      return(result)
    } #end function
    #------------------------------------------------------------------------------------


    #------------------------------------------------------------------------------------
    # Modify Molecular Formula Function
    # Function to Create Data Frame with modified molecular formulas.
    # Inputs are Molecular Formula (as a string) and Number of Heavy AA
    # Produces data frames with modified molecular formula and
    # the corresponding number of heavy AA
    # N number of natural isotopes (H, C, O, or N) are removed based on user defined number and the same number of isotoptes are added to the molecular formula for
    # each additional heavy AA
    modify.mol.formula.fun <- function(molecular.formula, num.heavy.AA, heavy_element,isotope,num_iso){
      # define elements vector to search against
      element.vector <- c("H","C","N","O","S",paste(isotope)) # all elements (including heavy isotope)
      element.vector.noH <- c("H","C","N","O","S") # vector of elements, next line removed heavy element from this list
      element.vector.noH <- element.vector.noH[! element.vector.noH %in% heavy_element]    # remove Heavy element or heavy isotope

      if(num.heavy.AA==0){ # make no changes if there are no heavy AA
        # make data frame with molecular formula and the number of heavy AA
        df.mod.mol.formula <- data.frame(Molecular.formula=as.character(molecular.formula), Number.Heavy.AA=as.numeric(num.heavy.AA))
        return(df.mod.mol.formula)
      } else { # modify molecular formula by reducing number of heavy element (ie H) and increasing number of heavy element isotope (ie D) in equal measure
        # initialize data frame to store modified molecular formulas
        df.mod.mol.formula <- dplyr::as_tibble(matrix(NA, ncol = 2, nrow = num.heavy.AA+1, # initializing a tibble with 2 named columns
                                                      dimnames=list(NULL, c("Molecular.Formula", "Number.Heavy.AA")))) # names of columns

        #convert was superceded by mutate_at in Dplyr
        #df.mod.mol.formula <- df.mod.mol.formula %>%
        #  convert(chr(Molecular.Formula), int(Number.Heavy.AA)) # convert column types: Molecular.Formula to character, Number.Heavy.AA to integer

        #convert formula column to character
        df.mod.mol.formula <- df.mod.mol.formula %>%
          mutate_at(c('Molecular.Formula'), as.character)

        #convert num AA column to integer
        df.mod.mol.formula <- df.mod.mol.formula %>%
          mutate_at(c('Number.Heavy.AA'), as.integer)

        df.mod.mol.formula[1,] <- list(molecular.formula, 0) # write to first row of tibble, which is 0 heavy AA and thus not molecular.formula is not modified

        mol.formula.pieces <- unlist(strsplit(molecular.formula, "")) # molecular formula pieces
        HE.index <- which(mol.formula.pieces==paste(heavy_element)) # get index of heavy element in the molecular formula
        # define remaining molecular formula pieces after hydrogen
        remaining.mol.formula.pieces <- mol.formula.pieces[(HE.index+1):length(mol.formula.pieces)] # keep only HE and remaining pieces after hydrogen in the molecular formula
        # get HE factor from molecular formula
        HE.factor <- ifelse(any(remaining.mol.formula.pieces %in% element.vector.noH),  # are there any remaining elements after the HE?
                            as.numeric(unlist(strsplit(paste(replace(remaining.mol.formula.pieces, which(remaining.mol.formula.pieces %in% element.vector), ","), collapse=""), ","))[1]), # if YES there are remaining elements after hydrogen
                            as.numeric(paste(remaining.mol.formula.pieces, sep="", collapse=""))) # if NO there are not any remaining elements after hydrogen
        # get HE factor indices from molecular formula; used to modify the molecular formula
        HE.factor.indeces <- seq(from=HE.index+1, to=HE.index+nchar(HE.factor)) # indeces where HE factor is stored in molecular formula
        # define first part of the molecular formula up to and including HE
        first.part.mol.formula <- paste(mol.formula.pieces[1:HE.index], sep="", collapse="")
        # define last part of the molecular formula pieces after HE
        last.part.mol.formula <- paste(remaining.mol.formula.pieces[-seq(from=1, to=nchar(HE.factor))], sep="", collapse="")

        # Add 3 HE isotope (ie Deuterium) and Remove 3 HE (ie Hydrogen) in each iteration according to the number of heavy AA
        #note: num_iso is user defined increase in molecular weight (ie +3 D or +6 heavy C, etc)
        for(i in 1:num.heavy.AA){
          num.HE_iso <- num_iso*i # increase in increments of 3/6/etc for each iteration (ie each heavy leucine has three deuteriums)
          num.HE <- HE.factor-num_iso*i # decrease number of hydrogen by 3 for each iteration, deuterium is replacing hydrogen
          if(num.HE>0){ # positive number of HE, which should be the case
            mod.molecular.formula.temp <- paste(first.part.mol.formula, num.HE, last.part.mol.formula, sep="", collapse="") # replace HE factor with the newly modified number of hydrogen
            mod.molecular.formula <- paste(mod.molecular.formula.temp, paste(isotope), num.HE_iso, sep="", collapse="") # add HE_iso (ie deteurium etc) to end of molecular formula
            df.mod.mol.formula[1+i,] <- list(as.character(mod.molecular.formula), as.numeric(i)) # store results
          } else{ # zero or negative number of HE (ie hydrogens), which should hopefully never be the case
            df.mod.mol.formula[1+i,] <- list(NA, i) # store NA for molecular formula and "i" for the number of heavy AA
          }
        } #end for

        return(df.mod.mol.formula)
      }
    } #end function
    #------------------------------------------------------------------------------------


    # Modify Molecular Formula Function (only when second heavy element is specific (ie, the heavy-labeled amino acid contains both heavy C and N, for example))
    # Function to Create Data Frame with modified molecular formulas.
    # Inputs are Molecular Formula (as a string) and Number of Heavy AA
    # Produces data frames with modified molecular formula and
    # the corresponding number of heavy AA
    # N number of natural isotopes (H, C, O, or N) are removed based on user defined number and the same number of isotoptes are added to the molecular formula for
    # each additional heavy AA
    modify.mol.formula.fun2 <- function(molecular.formula, num.heavy.AA, heavy_element,isotope,num_iso, heavy_element2, isotope2,num_iso2){
      # define elements vector to search against
      element.vector <- c("H","C","N","O","S",paste(isotope)) # all elements (including heavy isotope)
      element.vector.noH <- c("H","C","N","O","S") # vector of elements, next line removed heavy element from this list
      element.vector.noH <- element.vector.noH[! element.vector.noH %in% heavy_element]    # remove Heavy element or heavy isotope

      #repeat with second heavy element
      # define elements vector to search against
      element.vector2 <- c("H","C","N","O","S",paste(isotope2)) # all elements (including heavy isotope)
      element.vector.noH2 <- c("H","C","N","O","S") # vector of elements, next line removed heavy element from this list
      element.vector.noH2 <- element.vector.noH2[! element.vector.noH2 %in% heavy_element2]    # remove Heavy element or heavy isotope

      if(num.heavy.AA==0){ # make no changes if there are no heavy AA
        # make data frame with molecular formula and the number of heavy AA
        df.mod.mol.formula <- data.frame(Molecular.formula=as.character(molecular.formula), Number.Heavy.AA=as.numeric(num.heavy.AA))
        return(df.mod.mol.formula)
      } else { # modify molecular formula by reducing number of heavy element (ie H) and increasing number of heavy element isotope (ie D) in equal measure
        # initialize data frame to store modified molecular formulas
        df.mod.mol.formula <- dplyr::as_tibble(matrix(NA, ncol = 2, nrow = num.heavy.AA+1, # initializing a tibble with 2 named columns
                                                      dimnames=list(NULL, c("Molecular.Formula", "Number.Heavy.AA")))) # names of columns

        #convert is superceded by mutate_at in dplyr
        #df.mod.mol.formula <- df.mod.mol.formula %>%
        #  convert(chr(Molecular.Formula), int(Number.Heavy.AA)) # convert column types: Molecular.Formula to character, Number.Heavy.AA to integer

        #convert formula column to character
        df.mod.mol.formula <- df.mod.mol.formula %>%
          mutate_at(c('Molecular.Formula'), as.character)

        #convert num AA column to integer
        df.mod.mol.formula <- df.mod.mol.formula %>%
          mutate_at(c('Number.Heavy.AA'), as.integer)

        df.mod.mol.formula[1,] <- list(molecular.formula, 0) # write to first row of tibble, which is 0 heavy AA and thus not molecular.formula is not modified

        mol.formula.pieces <- unlist(strsplit(molecular.formula, "")) # molecular formula pieces

        HE.index <- which(mol.formula.pieces==paste(heavy_element)) # get index of heavy element in the molecular formula
        # define remaining molecular formula pieces after hydrogen
        remaining.mol.formula.pieces <- mol.formula.pieces[(HE.index+1):length(mol.formula.pieces)] # keep only HE and remaining pieces after hydrogen in the molecular formula
        # get HE factor from molecular formula
        HE.factor <- ifelse(any(remaining.mol.formula.pieces %in% element.vector.noH),  # are there any remaining elements after the HE?
                            as.numeric(unlist(strsplit(paste(replace(remaining.mol.formula.pieces, which(remaining.mol.formula.pieces %in% element.vector), ","), collapse=""), ","))[1]), # if YES there are remaining elements after hydrogen
                            as.numeric(paste(remaining.mol.formula.pieces, sep="", collapse=""))) # if NO there are not any remaining elements after hydrogen
        # get HE factor indices from molecular formula; used to modify the molecular formula
        HE.factor.indeces <- seq(from=HE.index+1, to=HE.index+nchar(HE.factor)) # indeces where HE factor is stored in molecular formula
        # define first part of the molecular formula up to and including HE
        first.part.mol.formula <- paste(mol.formula.pieces[1:HE.index], sep="", collapse="")
        # define last part of the molecular formula pieces after HE
        last.part.mol.formula <- paste(remaining.mol.formula.pieces[-seq(from=1, to=nchar(HE.factor))], sep="", collapse="")


        # Add 3 HE isotope (ie Deuterium) and Remove 3 HE (ie Hydrogen) in each iteration according to the number of heavy AA
        #note: num_iso is user defined increase in molecular weight (ie +3 D or +6 heavy C, etc)
        for(i in 1:num.heavy.AA){
          num.HE_iso <- num_iso*i # increase in increments of 3/6/etc for each iteration (ie each heavy leucine has three deuteriums)
          num.HE <- HE.factor-num_iso*i # decrease number of hydrogen by 3 for each iteration, deuterium is replacing hydrogen

          if(num.HE>0){ # positive number of HE, which should be the case
            mod.molecular.formula.temp <- paste(first.part.mol.formula, num.HE, last.part.mol.formula, sep="", collapse="") # replace HE factor with the newly modified number of hydrogen
            mod.molecular.formula <- paste(mod.molecular.formula.temp, paste(isotope), num.HE_iso, sep="", collapse="") # add HE_iso (ie deteurium etc) to end of molecular formula

            ##
            #repeat total process with second heavy element
            mol.formula.pieces <- unlist(strsplit(mod.molecular.formula, "")) # molecular formula pieces
            #repeat indexing for second heavy element
            HE.index2 <- which(mol.formula.pieces==paste(heavy_element2)) # get index of heavy element in the molecular formula
            # define remaining molecular formula pieces after hydrogen
            remaining.mol.formula.pieces2 <- mol.formula.pieces[(HE.index2+1):length(mol.formula.pieces)] # keep only HE and remaining pieces after hydrogen in the molecular formula
            # get HE factor from molecular formula
            HE.factor2 <- ifelse(any(remaining.mol.formula.pieces2 %in% element.vector.noH2),  # are there any remaining elements after the HE?
                                 as.numeric(unlist(strsplit(paste(replace(remaining.mol.formula.pieces2, which(remaining.mol.formula.pieces2 %in% element.vector2), ","), collapse=""), ","))[1]), # if YES there are remaining elements after hydrogen
                                 as.numeric(paste(remaining.mol.formula.pieces2, sep="", collapse=""))) # if NO there are not any remaining elements after hydrogen
            # get HE factor indices from molecular formula; used to modify the molecular formula
            HE.factor.indeces2 <- seq(from=HE.index2+1, to=HE.index2+nchar(HE.factor2)) # indeces where HE factor is stored in molecular formula
            # define first part of the molecular formula up to and including HE
            first.part.mol.formula2 <- paste(mol.formula.pieces[1:HE.index2], sep="", collapse="")
            # define last part of the molecular formula pieces after HE
            last.part.mol.formula2 <- paste(remaining.mol.formula.pieces2[-seq(from=1, to=nchar(HE.factor2))], sep="", collapse="")

            #second HE isotope counts
            num.HE_iso2 <- num_iso2*i # increase in increments of 3/6/etc for each iteration (ie each heavy leucine has three deuteriums) #second HE
            num.HE2 <- HE.factor2-num_iso2*i # decrease number of hydrogen by 3 for each iteration, deuterium is replacing hydrogen #second HE

            mod.molecular.formula.temp <- paste(first.part.mol.formula2, num.HE2, last.part.mol.formula2, sep="", collapse="") # replace HE factor with the newly modified number of hydrogen
            mod.molecular.formula <- paste(mod.molecular.formula.temp, paste(isotope2), num.HE_iso2, sep="", collapse="") # add HE_iso (ie deteurium etc) to end of molecular formula

            ##
            #End Repeat with second HE

            df.mod.mol.formula[1+i,] <- list(as.character(mod.molecular.formula), as.numeric(i)) # store results
          } else{ # zero or negative number of HE (ie hydrogens), which should hopefully never be the case
            df.mod.mol.formula[1+i,] <- list(NA, i) # store NA for molecular formula and "i" for the number of heavy AA
          }
        } #end for


        return(df.mod.mol.formula)
      }
    } #end function

    #------------------------------------------------------------------------------------
    # Molecular Formula Nest Function
    # Function which takes a data frame with Modified Molecular Formulas with accompanying number of heavy leucines as input
    # and returns a nested data frame containing the element-factor information for each modified molecular formula

    df.mol.formula.fun <- function(df.mod.mol.formula){
      # define element and numeric vector to search against
      element.vector <- c("H", "C", "N", "O", "S", paste(HE_iso)) # vector of elements which we are correcting distribution for: Hydrogen, Carbon, Nitrogen, Oxygen, Sulfur, HE_iso (is deuterium, etc)
      #add second isotope if user-defined
      if(second_element != FALSE){
        element.vector <- c(c(element.vector), paste(HE_iso2)) # vector of elements which we are correcting distribution for: Hydrogen, Carbon, Nitrogen, Oxygen, Sulfur, HE_iso (is deuterium, etc)
      }

      numeric.vector <- c(0:9) # vector of numbers 0 to 9

      # first create empty data frame to write results to
      df.lc.mod.mol.formula <- data.frame(Molecular.Formula=as.character(), Number.Heavy.AA=as.numeric(),
                                          Element=as.character(), Factor=as.numeric())

      for(i in seq_along(df.mod.mol.formula[["Molecular.Formula"]])){   # iterate through all molecular formulas
        # mol.formula <- df.mod.mol.formula[i, "Molecular.Formula"] %>%   # take the molecular formula and convert to character string, such as "C62H108N18O25"
        #   as.character(.)
        mol.formula <- df.mod.mol.formula$Molecular.Formula[[i]]
        mol.formula.pieces <- unlist(strsplit(mol.formula, "")) # split and unlist the molecular formula into individual single element characters

        # ELEMENTS
        element.logical <- is.element(unlist(strsplit(mol.formula, "")), element.vector) # check which elements of molecular formula are in our element list
        elements.vector <- mol.formula.pieces[element.logical] # elements in molecular.formula in the order they appear

        # FACTORS
        factors <- paste(replace(mol.formula.pieces, element.logical, ","), collapse="") # collapse factors seperated by commas
        factors.vector <- unlist(strsplit(factors, ","))[-1] # unlist factors and remove the first element which is blank, now we have a vector of factors in the order they appear

        # Data Frame of elements with their factors
        df.mol.formula <- data.frame(Element=as.character(elements.vector), Factor=as.numeric(factors.vector)) # make data frame of the molecular formula

        # create a temporary data frame to store molecular formula information
        temp <- df.mol.formula %>%
          cbind("Molecular.Formula"=mol.formula, "Number.Heavy.AA"=as.numeric(df.mod.mol.formula$Number.Heavy.AA[[i]]))

        # add temp data frame to df.lc.mod.mol.formula which exists outside of for loop which will be nested after running this loop
        df.lc.mod.mol.formula <- add_row(df.lc.mod.mol.formula, "Molecular.Formula"=temp$Molecular.Formula, "Number.Heavy.AA"=temp$Number.Heavy.AA,"Element"=temp$Element, "Factor"=temp$Factor)
      } #end for

      # nest for each modified molecular formula
      df.lc.mod.mol.formula <- df.lc.mod.mol.formula %>%
        group_by(Molecular.Formula, Number.Heavy.AA) %>%
        nest()

      # rename list column 'data' to 'df.mol.formula' since the molecular formula is stored there as a data frame
      names(df.lc.mod.mol.formula)[names(df.lc.mod.mol.formula)=="data"] <- "df.mol.formula"

      return(df.lc.mod.mol.formula)
    } #end function
    #------------------------------------------------------------------------------------


    #------------------------------------------------------------------------------------
    # Molecular Formula Distribution Calculator Function
    # Function which calculates Mass Distribution using Data Frame of elements with their factors.
    # Input is a data frame with molecular formula stored in an a nested data frame with columns "Element" and "Factor".
    # Output is a data frame with list column "data" which stores each element-factor mass distribution.
    df.mol.formula.distribution.fun <- function(df.lc.mod.mol.formula){

      # first create empty data frame to write results to
      df.lc.mod.mol.formula.dist <- data.frame(Molecular.Formula=as.character(), Number.Heavy.AA=as.numeric(), Mass=as.numeric(), Abundance=as.numeric())

      for(j in seq_along(df.lc.mod.mol.formula$Molecular.Formula)){
        df.mol.formula <- df.lc.mod.mol.formula$df.mol.formula[[j]] # define molecular formula for given modified molecular formula

        # Calculate distribution for the molecular formula and write out results to df.lc.mol.formula
        for(i in seq_along(df.mol.formula$Element)){
          element <- df.mol.formula$Element[[i]] # declare element
          self <- ntbl.elements[[ which(element==ntbl.elements[["Symbol"]]), "data"]] %>% # declare naturally occuring isotopic element distribution, which we call 'self'
            as.data.frame() # self needs to be data frame for proper references in add.fun
          factor <- df.mol.formula$Factor[[i]] # declare factor
          distribution <- multiply.fun(self, factor) # calculate mass distribution for this element-factor combination

          # Filter for Minimum Abundance
          distribution <- filter(distribution, distribution$Abundance>min.abundance)

          # temporary data frame with mass distribution - used to add onto external data frame to store results
          df.temp <- distribution %>%
            cbind("Molecular.Formula"=df.lc.mod.mol.formula$Molecular.Formula[[j]], "Number.Heavy.AA"=as.numeric(df.lc.mod.mol.formula$Number.Heavy.AA[[j]]),
                  "Element"=element, "Factor"=factor)

          # row bind onto external data frame
          df.lc.mod.mol.formula.dist <- rbind(df.lc.mod.mol.formula.dist, df.temp)
        } #end for

      } #end for

      # nest by molecular formula and number heavy leucines
      df.lc.mod.mol.formula.dist <-  df.lc.mod.mol.formula.dist %>%
        group_by(Molecular.Formula, Number.Heavy.AA) %>%
        nest()

      # create empty data frame to write results to
      test.out <- data.frame(Molecular.Formula=NA, Number.Heavy.AA=NA, Element=NA, Factor=NA, data=NA)

      # for each molecular formula nest by element and factor the already nested data frame again
      for(k in seq_along(df.lc.mod.mol.formula.dist$Molecular.Formula)){
        df.distribution <- df.lc.mod.mol.formula.dist$data[[k]] %>%
          group_by(Element, Factor) %>%
          nest()


        # add molecular formula and number of heavy leucines back onto the distribution, so we can nest after for loop
        test <- cbind("Molecular.Formula"=rep(as.character(df.lc.mod.mol.formula.dist$Molecular.Formula[[k]]), length=nrow(df.distribution)),
                      "Number.Heavy.AA"=rep(df.lc.mod.mol.formula.dist$Number.Heavy.AA[[k]], length=nrow(df.distribution)),
                      df.distribution)

        # row bind to test.out which exists outside of for loop
        # bind_rows is a workaround since we are binding a nested data frame, this will throw warnings but not an error
        test.out <- test.out %>%
          bind_rows(test)
      } #end for

      # remove first row of NA's which were there as placeholder
      test.out <- na.omit(test.out)

      # nest
      test.out <- test.out %>%
        group_by(Molecular.Formula, Number.Heavy.AA) %>%
        nest()

      # Return molecular formula data frame with list column which stores the mass distribution of each element-factor combination
      return(test.out)
    } #end function
    #------------------------------------------------------------------------------------


    #------------------------------------------------------------------------------------
    # Adder Function
    # Function to Calculate Full Molecular Formula Mass Distribution from individual mass distributions stored in list-column of data frame.
    # Input is the list column named "data" which contains the distributions which need to be added together, which we do by recursively using itself and the add.function.
    # Output is the full molecular formula mass distribution.

    adder.fun <- function(x){
      y <- x$data # y is the mass distributions which need to be added together
      if (nrow(x)==1){ return ( (y[[1]]))} # terminating condition, return first element of list-column when only one row remains
      else {return( add.fun( (y[[1]]), adder.fun(x[-1,])))} # recursively add first distribution with all of the rest of the distributions
    } #end function
    #------------------------------------------------------------------------------------


    #------------------------------------------------------------------------------------
    # Help Adder Function
    # Helper function which prepares the nested data for input into adder function
    # and calls on the adder function to calculate the full mass distribution for any given molecular formula.
    # Input is a data frame with molecular formulas, which contains a nested data frame with the molecular formula stored as element-factor columns with the individual distributions in another nested data frame
    help.adder.fun <- function(x){
      # create empty data frame to write results to
      df.final <- data.frame(Molecular.Formula=as.character(), Number.Heavy.AA=as.numeric(), Mass=as.numeric(), Abundance=as.numeric())

      # calculate the full mass distribution for each modified molecular formula
      for(i in seq_along(x$Molecular.Formula)){
        z <- x$data[[i]] # access the df.lc.mod.mol.formula structure which is the first nested data frame
        df.distribution <- adder.fun(z) # calculate full mass distribution for this given molecular formula using the adder.function
        # add molecular formula and number of heavy leucines onto the distribution so we can seperate by group later
        df.final <- df.distribution %>%
          cbind("Molecular.Formula"=rep(x$Molecular.Formula[i], length=nrow(df.distribution)), "Number.Heavy.AA"=rep(x$Number.Heavy.AA[i], length=nrow(df.distribution))) %>%
          rbind(df.final)
      } #end for

      # nest the distributions by molecular formula
      df.lc.final <- df.final %>%
        group_by(Molecular.Formula, Number.Heavy.AA) %>%
        nest()

      # rename nested data frame to 'df.distribution'
      names(df.lc.final)[names(df.lc.final)=="data"] <- "df.distribution"
      return(df.lc.final)
    } #end function
    #------------------------------------------------------------------------------------


    #------------------------------------------------------------------------------------
    # Resolution Function
    # Function to combine masses within the mass resolution limit by weighted average
    # Inputs are x and y
    # x is the mass distribution - the initial input of the function
    # y is the reduced mass distribution - the result of the function
    # Output is the reduced mass distribution data frame, using resolution
    # 'resolution' is a globally defined value; we can make it a user defined value if desired

    weighted.average.fun <- function(x, y){
      if(nrow(x)==0){ # terminating condition: if input distribution is empty then we are done, return y
        return(y) # reduced mass distribution, computed when finished
      } else{
        cur.mass <- x$Mass[1] # grab first mass
        df.cur <- x %>% # peel off all masses within resolution of current mass
          filter(abs(x$Mass-cur.mass)<resolution)
        x <- x %>% # reduce df.distribution, retain all masses which are not within resolution of current mass
          filter(!abs(x$Mass-cur.mass)<resolution)

        # Store Results
        if(nrow(df.cur)==1){ # if df.cur only has one mass write this to results df
          y <- rbind(y, df.cur)
        } else{ # calculate weighted average of masses and the sum of the abundances
          weighted.avg.mass <- sum(df.cur$Mass*df.cur$Abundance)/sum(df.cur$Abundance)
          new.abundance <- sum(df.cur$Abundance)
          # write to results df
          y <- add_row(y, "Mass"=weighted.avg.mass, "Abundance"=new.abundance)
        }
        return(weighted.average.fun(x, y)) # recursive call
      }
    } #end function
    #------------------------------------------------------------------------------------


    #------------------------------------------------------------------------------------
    # Help Resolution Function
    # Helper function which prepares the nested data for input into adder function
    # and calls on the adder function to calculate the full mass distribution for any given molecular formula.
    # Input is a data frame with molecular formulas (df.lc.final), which contains nested data frames of mass distributions for each modified molecular formula
    help.weighted.avg.fun <- function(x){

      # initialize empty data frame to write results to
      df.simplified.out <- data.frame(Molecular.Formula=as.character(), Number.Heavy.AA=as.numeric(), Mass=as.numeric(), Abundance=as.numeric())

      # initialize empty data frame to write reduced mass distribution to. This will be passed to weighted.average.function
      df.distribution.reduced <- data.frame(Mass=as.numeric(), Abundance=as.numeric())

      # simplify the mass distribution by resolution for each modified molecular formula
      for(i in seq_along(x$Molecular.Formula)){
        df.simplified <- weighted.average.fun(x=x$df.distribution[[i]], y=df.distribution.reduced) %>% # simplify the mass distribution using resolution
          arrange(Mass) # and arrange by Mass
        # add molecular formula and number of heavy AA onto the distribution so we can seperate groups and nest mass distributions later
        df.simplified.out <- df.simplified %>%
          cbind("Molecular.Formula"=rep(x$Molecular.Formula[i], length=nrow(df.simplified)), "Number.Heavy.AA"=rep(x$Number.Heavy.AA[i], length=nrow(df.simplified))) %>%
          rbind(df.simplified.out)
      } #end for

      # Filter for Minimum Abundance
      df.simplified.out <- filter(df.simplified.out, df.simplified.out$Abundance>min.abundance)

      # nest the distributions by molecular formula
      df.lc.simplified.out <- df.simplified.out %>%
        group_by(Molecular.Formula, Number.Heavy.AA) %>%
        nest()

      # rename nested data frame to 'df.distribution.reduced'
      names(df.lc.simplified.out)[names(df.lc.simplified.out)=="data"] <- "df.distribution.reduced"
      return(df.lc.simplified.out)
    } #end function
    #------------------------------------------------------------------------------------


    #------------------------------------------------------------------------------------
    # Build Matrix P Function
    # combining observations within mass tolerance
    # Function to combine masses within a tolerance region
    # Inputs are x and y
    # x is the data frame containing observations - the initial input of the function
    # y is the reduced data frame - the result of the function

    # p.tolerance <- 0.05 # tolerance for combining masses in observed data together,
    # this is many times the resolution yet small enough to distinguish between
    # M+0, M+1, M+2, etc... peaks when not overlapping due to Heavy Leucine contribution (+3)

    p.matrix.build.fun <- function(x, y){
      if(nrow(x)==0){ # terminating condition: if input distribution is empty then we are done, return y
        return(y)
      } else{
        # this chunk could go into a helper function, since it only needs to be run once on the onset
        x <- x %>%
          select(Mass.recalc, Area) # keep only Mass and Area columns
        x$Area <- as.numeric(x$Area)

        cur.mass <- x$Mass.recalc[1] # grab first mass
        df.cur <- x %>% # peel off all masses within tolerance of current mass
          filter(abs(x$Mass.recalc-cur.mass)<p.tolerance)

        x <- x %>% # reduce df.distribution, retain all masses which are not within tolerance of current mass
          filter(!abs(x$Mass.recalc-cur.mass)<p.tolerance)

        # Store Results
        if(nrow(df.cur)==1){ # if df.cur only has one mass write this to results df
          y <- rbind(y, df.cur)
        } else{ # calculate mean of masses and mean of areas
          avg.mass <- mean(as.numeric(df.cur$Mass.recalc))
          new.area <- mean(as.numeric(df.cur$Area))
          # write to results df
          y <- add_row(y, "Mass.recalc"=avg.mass, "Area"=new.area)
        }
        return(p.matrix.build.fun(x, y)) # recursive call
      }
    } #end function
    #------------------------------------------------------------------------------------

    # end FUNCTIONS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



    #------------------------------------------------------------------------------------
    # Add Columns using defined functions from above

    # Add number of heavy AA column
    df <- df %>%
      mutate("Number.Heavy.AA"= map_dbl(Modified.Sequence, AA.heavy.tally.fun)) # add number of heavy AA column

    # Add a new modified peptide sequence column, which strips away all Leucine modifications
    df <- df %>%
      mutate("Modified.Peptide.Seq"= map_chr(Modified.Sequence, aa.modification.fun))

    # Move most important columns to the left of the data frame, and keep all other columns on the right
    df <- df %>%
      dplyr::select(Protein.Accession, Protein.Gene, Peptide, Modified.Peptide.Seq, Replicate.Name, Condition, Timepoint, Number.Heavy.AA, Area, everything()) # NO COHORT - removed cohort from below on 7/24/2020 before first test run ~4:15pm
    #------------------------------------------------------------------------------------


    #------------------------------------------------------------------------------------
    # Find Best Combination - FBC
    # Find best linear combination of non-heavy and heavy leucine peptides
    # from the observed data that match the theoretical mass distribution with correcting for naturally occuring heavy isotopes.

    # all proteins
    proteins <- unique(df$Protein.Gene) # all proteins, once the script works the first time we can try running through all proteins

    # # test case - just one protein for testing
    # # if we want to test a shorter version then just run one protein
    # proteins <- unique(df$Protein.Gene)[1] # just one protein

    # create df.test data frame to run FBC on
    df.test <- df %>%
      filter(df$Protein.Gene %in% proteins)

    # Initialize Solutions Data Frame for storing the Find Best Combination solutions
    # long format with Number.Heavy.Leucine and individual Isotope.Dot.Product (IDP) values
    df.solutions <- data.frame(matrix(NA,
                                      nrow=length(unique(df.test$Replicate.Name))*length(unique(df.test$Peptide))*length(unique(df.test$Product.Charge)*max(df.test$Number.Heavy.AA)),
                                      ncol=13))

    # name columns
    names(df.solutions)[1:13] <- c("Protein.Gene", "Protein.Accession", "Peptide", "Modified.Peptide", "Replicate.Name", "Condition",
                                   "Timepoint", "Product.Charge", "Number.Heavy.AA", "Detection.Q.Value", "Total.Area.MS1",
                                   "Isotope.Dot.Product", "FBC.Solution")

    # set counter before loop
    counter <- 1

    # loop through PROTEINS
    for(k in seq_along(proteins)){
      # subset for protein
      df.protein <- filter(df, df$Protein.Gene==proteins[k])

      # vector of unique peptides
      peps <- unique(df.protein$Peptide)

      # loop through PEPTIDES
      for(n in seq_along(peps)){
        # subset for peptide
        df.peptide <- filter(df.protein, df.protein$Peptide==peps[n])

        # vector of unique modified.peptide.sequences # these are variable modifications that can exist
        mod.peps <- unique(df.peptide$Modified.Peptide.Seq)

        # loop through MODIFIED.PEPTIDE.SEQUENCES
        for(p in seq_along(mod.peps)){
          # subset for modified peptide sequence
          df.mod.peptide <- filter(df.peptide, df.peptide$Modified.Peptide.Seq==mod.peps[p])

          # vector of unique replicates
          reps <- unique(df.mod.peptide$Replicate.Name)

          # add 1 to the end of Molecule.Formula column if the last element is not a number
          df.mod.peptide <- df.mod.peptide %>%
            mutate_at( vars(Molecule.Formula), list(~ifelse(unlist(str_split(unique(df.mod.peptide$Molecule.Formula), ""))[nchar(unique(df.mod.peptide$Molecule.Formula))] %in% c("H","C","N","O","S"),
                                                            paste(c(unlist(str_split(unique(df.mod.peptide$Molecule.Formula), "")), "1"), collapse=""), unique(df.mod.peptide$Molecule.Formula))))

          # Find Best Combination if the data have both: heavy AA peptides and zero heavy AA peptides -- otherwise cannot find best combination
          # df.mod.peptide must have the full sequence of number.heavy.AA from 0 to the maximum value -- this will be checked for again at the df.charge level, but if it fails here then we can exit out of find.best.combination
          if(all(0:max(unique(df.mod.peptide$Number.Heavy.AA)) %in% unique(df.mod.peptide$Number.Heavy.AA), # must have the full sequence of number.heavy.AA from 0 to the maximum value
                 max(df.mod.peptide$Number.Heavy.AA)>0 # maximum Number.Heavy.AA must be greater than Zero
          )){
            # Matrix A:
            # Generate distributions for Matrix A based on Molecular Formula and max Number Heavy AA in df.mod.peptide
            #1
            #check for second element
            if (second_element == TRUE){
              df.mod.mol.formula <- modify.mol.formula.fun2(unique(df.mod.peptide$Molecule.Formula), max(df.mod.peptide$Number.Heavy.AA),HE,HE_iso,num_HE,HE2,HE_iso2,num_HE2)
            }
            if (second_element == FALSE){
              df.mod.mol.formula <- modify.mol.formula.fun(unique(df.mod.peptide$Molecule.Formula), max(df.mod.peptide$Number.Heavy.AA),HE,HE_iso,num_HE)
            }
            # 2
            df.lc.mod.mol.formula <- df.mol.formula.fun(df.mod.mol.formula)
            # 3
            df.lc.mod.mol.formula.dist <- df.mol.formula.distribution.fun(df.lc.mod.mol.formula)
            # 4
            df.lc.final <- help.adder.fun(df.lc.mod.mol.formula.dist)
            # 5
            df.lc.simplified.out <- help.weighted.avg.fun(df.lc.final) ## Theoretical data for Matrix A
            # Unnest df.lc.simplified.out
            df.A <- unnest(df.lc.simplified.out, cols = c(df.distribution.reduced))

            # loop through each REPLICATE.NAME
            for(m in seq_along(reps)){
              # subset for replicate
              df.rep <- filter(df.mod.peptide, df.mod.peptide$Replicate.Name==reps[m])

              # vector of product.charges
              charges <- unique(df.rep$Product.Charge)

              # loop through each PRODUCT.CHARGE
              for(l in seq_along(charges)){
                # subset for product.charge
                df.charge <- df.rep %>%
                  filter(df.rep$Product.Charge==charges[l])

                # check to see if full sequence of number.heavy.AA exists in df.charge ... it must in order to carry on with Find.Best.Combination
                if(all(0:max(unique(df.charge$Number.Heavy.AA)) %in% unique(df.charge$Number.Heavy.AA))){

                  # Keep only unique entries (drop duplicate isomers)
                  # back transform from Product.Mz to Precursor.Neutral.Mass, so that masses in observed matrix P are comparable to masses in theoretical matrix A
                  # arrange ascending number of heavy AA
                  df.charge <- df.charge %>%
                    filter(! duplicated(df.charge$Product.Mz)) %>% # keep observations with unique Product.Mz - remove duplicates, which come from isomers
                    mutate(Mass.recalc=Product.Mz*Precursor.Charge-unique(Precursor.Charge)*HE_weight) %>% # multiply product.mz by charge and subtract a charge number of the radiolabeled element (ie H)
                    mutate("Number.Heavy.AA"= map_dbl(Modified.Sequence, AA.heavy.tally.fun)) %>% # add number of heavy AA column
                    arrange(Number.Heavy.AA) # arrange ascending number of heavy AA
                  #bao change to incorporate heavy amino acids
                  #mutate(c("Number.Heavy.",paste(AA))= map_dbl(Modified.Sequence, AA.heavy.tally.fun)) %>% # add number of heavy AA column
                  #arrange(c("Number.Heavy.",paste(AA))) # arrange ascending number of heavy AA


                  # P Matrix
                  # Observed Areas - this is the mass spec data
                  x <- df.charge %>% select(Mass.recalc, Area) # input is the recalculated masses
                  y <- data.frame(Mass.recalc=as.numeric(), Area=as.numeric()) # empty data frame for results
                  P <- p.matrix.build.fun(x,y)
                  # Build Matrix P: column matrix of normalized observed Abundances (aka Areas, Intensities)
                  P.mtx <- P %>%
                    select(Area) %>%
                    apply(2, sum.to.one.fun)
                  # define length of P matrix for number of masses which matrix A should also have
                  p.length <- nrow(P.mtx)


                  # For each mass in P find the closest match in A and if it is within tolerance (fixed value) then keep that mass in A
                  # ensure a one to one correspondence for P to A
                  df.A.keep <- data.frame(Molecular.Formula=c(), Number.Heavy.AA=c(), Mass=c(), Abundance=c())
                  for(i in seq_along(df.charge$Mass.recalc)){
                    mass <- df.charge$Mass.recalc[i] # get mass
                    no.AA <- df.charge$Number.Heavy.AA[i] # get number of heavy AA
                    df.A.temp <- subset(df.A, Number.Heavy.AA==no.AA) # subset for corresponding number of heavy AA
                    mass.delta <- abs(mass-df.A.temp$Mass) # mass differences
                    best.A <- ungroup(filter(df.A.temp, min(mass.delta)==mass.delta)) # grab row with closest mass, ungroup for row binding
                    df.A.keep <- rbind(df.A.keep, best.A) # save to external data frame
                  } #end for

                  if(NIAC == FALSE){
                    df.A.keep$Abundance <- rep(1, nrow(df.A.keep))
                  }

                  # using df.A.keep, build matrix A
                  df.A.keep.zero <- filter(df.A.keep, Number.Heavy.AA==0)
                  df.A.keep.nonzero <- filter(df.A.keep, !Number.Heavy.AA==0)
                  df.reference <- df.A.keep.zero # define reference data frame to look back at, initially the one with 0 Heavy AA

                  # Calculate starting rows for 1,2,3,... heavy AA columns
                  rows <- c() # vector to store starting rows for 1,2,3,etc heavy AA columns
                  for(i in seq_along(unique(df.A.keep.nonzero$Number.Heavy.AA))){
                    df.cur <- subset(df.A.keep.nonzero, Number.Heavy.AA==i) %>% # define current data frame
                      arrange(Mass) # arrange from low to high mass
                    #BAO change, adding conditional, testing to see if masses overlap
                    #if masses dont overlap (ie the heavy AA is much heavier than the naturally occurring AA, thus no isotope envelope overlap)
                    #then assign row number 1 greater than the row of the heaviest mass of the previous iteration for num.heavy.AA
                    if (min(abs(df.cur$Mass[1]-df.reference$Mass)) <0.1){
                      row.match <- which(min(abs(df.cur$Mass[1]-df.reference$Mass))==abs(df.cur$Mass[1]-df.reference$Mass)) # take first mass in df.cur and find the row for the closest-matching mass in df.reference

                    }
                    else{
                      row.match <- length(df.reference$Mass) + 1
                    }
                    rows <- append(rows, row.match) # append matching row number to rows vector
                    if(i>1){
                      rows[i] <- rows[i] + rows[i-1]-1 # adjust current element of rows vector by adding the previous element of rows vector less one
                    }
                    df.reference <- df.cur # update df.reference to df.current prior to next iteration

                    if(i==max(unique(df.A.keep.nonzero$Number.Heavy.AA))){ # calculate the number of rows needed to build Matrix A
                      num.rows <- rows[i] + nrow(df.cur) - 1 # current element of rows vector plus the size of the current data less by one
                    }
                  } #end for
                  rows <- append(rows, 1, after=0) # append the value 1 to the beginning of rows vector, corresponding Zero Heavy AA column being positioned starting at the first row in Matrix A

                  # Build template for Matrix A now that we know starting row numbers
                  A.template <- data.frame(matrix(NA, nrow = 2, ncol = max(unique(df.A$Number.Heavy.AA))+1))
                  #A.column.names <- c("Zero", "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight") # column names from which to take and use in naming A Template
                  A.column.names <- c("0", "1", "2", "3", "4", "5", "6", "7", "8") # column names from which to take and use in naming A Template, numeric, for easy subsetting

                  names(A.template) <- A.column.names[1:ncol(A.template)] # name columns with the required names based on number of columns needed
                  row.names(A.template) <- c("Above", "Below") # rows for how much to pad above and below with zeros
                  for(i in seq_along(unique(df.A.keep$Number.Heavy.AA))){
                    num.heavy.AA <- i-1 # number of heavy AA is i less one
                    df.cur <- subset(df.A.keep, Number.Heavy.AA==num.heavy.AA) %>% # define current data frame for populating template
                      arrange(Mass) # arrange from low to high mass
                    top.pad <- rows[i]-1 # number of zeros to pad above in the template
                    bottom.pad <- p.length - top.pad - length(df.cur$Mass) # number of zeros to pad below in the template; top.pad + nrow(df.cur) + bottom.pad = constant (constant=length of matrix P)
                    A.template[,names(A.template)==num.heavy.AA] <- c(top.pad, bottom.pad) # write pad lengths to template
                  } #end for

                  # Building Matrix A:
                  # Now that we have the template for Matrix A, let's build it!
                  A <- data.frame(matrix(NA, nrow=p.length, ncol=max(unique(df.A$Number.Heavy.AA))+1))
                  names(A) <- A.column.names[1:ncol(A)]
                  for(i in seq_along(unique(df.A.keep$Number.Heavy.AA))){
                    num.heavy.AA <- i-1 # number of heavy AA is i less one
                    df.cur <- subset(df.A.keep, Number.Heavy.AA==num.heavy.AA) %>% # define current data frame for populating matrix
                      arrange(Mass) # arrange from low to high mass
                    top.pad <- A.template[,names(A.template)==num.heavy.AA][1] # grab top pad length
                    bottom.pad <- A.template[,names(A.template)==num.heavy.AA][2] # grab bottom pad length
                    A[,names(A)==num.heavy.AA] <- c(rep(0, length=top.pad), df.cur$Abundance, rep(0, length=bottom.pad))
                  } #end for

                  A <- as.matrix(A) # convert to matrix

                  # Solutions of Find Best Combination (FBC):
                  # Matrix Math from equation 5 of Brauman 1966 'Least Squares Analysis and Simplification of Multi-Isotope Mass Spectra'
                  A.tA <- t(A) %*% A # A.transpose times A
                  #bao change, add tolerance to solve function to allow inversion of matrix, turned off but would be: A.tA.inv <- solve(A.tA, tol = 1e-100)
                  A.tA.inv <- solve(A.tA) # inverse of (A.transpose times A)
                  A.math <- A.tA.inv %*% t(A) # inverse of (A.transpose times A) times A.transpose (this is all of the matrix math left of matrix P in equation 5 of Brauman 1966)
                  solutions <- A.math %*% P.mtx # P.mtx is the column matrix of observed abundances
                  solutions <- solutions %>%
                    apply(2, sum.to.one.fun) # normalize solutions


                  # now write these solutions out in a long format
                  rows.write.out <- counter:(counter+length(unique(df.charge$Number.Heavy.AA))-1)
                  counter <- max(rows.write.out) + 1 # increment counter
                  print(counter)

                  # Write Out to df.solutions:
                  df.solutions[rows.write.out, "FBC.Solution"] <- solutions
                  # Replicate Name
                  df.solutions[rows.write.out, "Replicate.Name"] <- unique(df.charge$Replicate.Name)
                  # protein.accession
                  df.solutions[rows.write.out, "Protein.Accession"] <- unique(df.charge$Protein.Accession)
                  # protein
                  df.solutions[rows.write.out, "Protein.Gene"] <- unique(df.charge$Protein.Gene)
                  # peptide
                  df.solutions[rows.write.out, "Peptide"] <- unique(df.charge$Peptide)
                  # modified.peptide
                  df.solutions[rows.write.out, "Modified.Peptide"] <- unique(df.charge$Modified.Peptide.Seq)
                  # Condition
                  df.solutions[rows.write.out, "Condition"] <- unique(df.charge$Condition)
                  # Timepoint
                  df.solutions[rows.write.out, "Timepoint"] <- as.numeric(unique(df.charge$Timepoint))
                  # Product.Charge
                  df.solutions[rows.write.out, "Product.Charge"] <- as.numeric(unique(df.charge$Product.Charge))
                  # Number Heavy AA
                  df.solutions[rows.write.out, "Number.Heavy.AA"] <- as.numeric(unique(df.charge$Number.Heavy.AA))
                  # Detection Qvalue
                  df.solutions[rows.write.out, "Detection.Q.Value"] <- as.numeric(unique(df.charge$Detection.Q.Value)) ### is this right? example data set has #N/A entries...check with good values
                  # Total Area MS1
                  df.solutions[rows.write.out, "Total.Area.MS1"] <- as.numeric(unique(df.charge$Total.Area.MS1))
                  # Isotope.Dot.Product - individual IDP values; one unique per peptide by Number of Heavy AA
                  # first create a vector of unique values based on Number Heavy AA ... then write them out
                  IDP <- c() # initialize IDP vector
                  nhl <- unique(df.charge$Number.Heavy.AA)
                  for(i in seq_along(unique(df.charge$Number.Heavy.AA))){
                    IDP.temp <- df.charge %>%
                      filter(Number.Heavy.AA==nhl[i]) %>%
                      pull(Isotope.Dot.Product) %>%
                      unique() %>%
                      as.numeric() %>%
                      replace_na(.,0) %>%
                      median() %>% # take the median here ensures IDP.temp is only one single quantity, which makes a one-to-one correspondence with number heavy AA
                      as.numeric()
                    IDP <- append(IDP, IDP.temp)
                  }
                  df.solutions[rows.write.out, "Isotope.Dot.Product"] <- IDP
                } else{ ## full sequence of number.heavy.AA is not contained in the data -- exited at the CHARGE level

                  # write out in a long format
                  rows.write.out <- counter:(counter+length(unique(df.charge$Number.Heavy.AA))-1)
                  counter <- max(rows.write.out) + 1 # increment counter
                  print(counter)

                  # Write Out to df.solutions:
                  #df.solutions[rows.write.out, "FBC.Solution"] <- solutions
                  # Replicate Name
                  df.solutions[rows.write.out, "Replicate.Name"] <- unique(df.charge$Replicate.Name)
                  # protein.accession
                  df.solutions[rows.write.out, "Protein.Accession"] <- unique(df.charge$Protein.Accession)
                  # protein
                  df.solutions[rows.write.out, "Protein.Gene"] <- unique(df.charge$Protein.Gene)
                  # peptide
                  df.solutions[rows.write.out, "Peptide"] <- unique(df.charge$Peptide)
                  # modified.peptide
                  df.solutions[rows.write.out, "Modified.Peptide"] <- unique(df.charge$Modified.Peptide.Seq)
                  # Condition
                  df.solutions[rows.write.out, "Condition"] <- unique(df.charge$Condition)
                  # Timepoint
                  df.solutions[rows.write.out, "Timepoint"] <- as.numeric(unique(df.charge$Timepoint))
                  # Product.Charge
                  df.solutions[rows.write.out, "Product.Charge"] <- as.numeric(unique(df.charge$Product.Charge))
                  # Number Heavy AA
                  df.solutions[rows.write.out, "Number.Heavy.AA"] <- as.numeric(unique(df.charge$Number.Heavy.AA))
                  # Detection Qvalue
                  df.solutions[rows.write.out, "Detection.Q.Value"] <- as.numeric(unique(df.charge$Detection.Q.Value)) ### is this right? example data set has #N/A entries...check with good values
                  # Total Area MS1
                  df.solutions[rows.write.out, "Total.Area.MS1"] <- as.numeric(unique(df.charge$Total.Area.MS1))
                  # Isotope.Dot.Product - individual IDP values; one unique per peptide by Number of Heavy AA
                  # first create a vector of unique values based on Number Heavy AA ... then write them out
                  IDP <- c() # initialize IDP vector
                  nhl <- unique(df.charge$Number.Heavy.AA)
                  for(i in seq_along(unique(df.charge$Number.Heavy.AA))){
                    IDP.temp <- df.charge %>%
                      filter(Number.Heavy.AA==nhl[i]) %>%
                      pull(Isotope.Dot.Product) %>%
                      unique() %>%
                      as.numeric() %>%
                      replace_na(.,0) %>%
                      median() %>% # take the median here ensures IDP.temp is only one single quantity, which makes a one-to-one correspondence with number heavy AA
                      as.numeric()
                    IDP <- append(IDP, IDP.temp)
                  }
                  df.solutions[rows.write.out, "Isotope.Dot.Product"] <- IDP
                } #end else


              } #end for - product charge level
            } #end for - replicate level
          } else { #end if, else -- no solutions to be found, exited at the PEPTIDE level

            # write out in a long format
            rows.write.out <- counter:(counter+length(unique(df.mod.peptide$Number.Heavy.AA))-1) # replace df.charge with df.mod.peptide since for this iteration df.mod.peptide is the furthest along (df.charge does not exist)
            counter <- max(rows.write.out) + 1 # increment counter
            print(counter)

            # Write Out to df.solutions:
            #df.solutions[rows.write.out, "FBC.Solution"] <- solutions
            # Replicate Name
            df.solutions[rows.write.out, "Replicate.Name"] <- unique(df.mod.peptide$Replicate.Name)[1] # take just the first element in case there are multiple values
            # protein.accession
            df.solutions[rows.write.out, "Protein.Accession"] <- unique(df.mod.peptide$Protein.Accession)
            # protein
            df.solutions[rows.write.out, "Protein.Gene"] <- unique(df.mod.peptide$Protein.Gene)
            # peptide
            df.solutions[rows.write.out, "Peptide"] <- unique(df.mod.peptide$Peptide)
            # modified.peptide
            df.solutions[rows.write.out, "Modified.Peptide"] <- unique(df.mod.peptide$Modified.Peptide.Seq)
            # Condition
            df.solutions[rows.write.out, "Condition"] <- unique(df.mod.peptide$Condition)[1] # take just the first element in case there are multiple values
            # Timepoint
            df.solutions[rows.write.out, "Timepoint"] <- as.numeric(unique(df.mod.peptide$Timepoint))[1] # take just the first element in case there are multiple values
            # Product.Charge
            #bao change: use min charge, ask nate about this:
            #df.solutions[rows.write.out, "Product.Charge"] <- as.numeric(unique(df.mod.peptide$Product.Charge))
            df.solutions[rows.write.out, "Product.Charge"] <- min(as.numeric(unique(df.mod.peptide$Product.Charge)))
            # Number Heavy AA
            df.solutions[rows.write.out, "Number.Heavy.AA"] <- as.numeric(unique(df.mod.peptide$Number.Heavy.AA))
            # Detection Qvalue
            df.solutions[rows.write.out, "Detection.Q.Value"] <- as.numeric(unique(df.mod.peptide$Detection.Q.Value))[1] # take just the first element in case there are multiple values
            # Total Area MS1
            df.solutions[rows.write.out, "Total.Area.MS1"] <- as.numeric(unique(df.mod.peptide$Total.Area.MS1))[1] # take just the first element in case there are multiple values
            # Isotope.Dot.Product - individual IDP values; one unique per peptide by Number of Heavy AA
            # first create a vector of unique IDP values based on Number Heavy AA ... then write them out
            IDP <- c() # initialize IDP vector
            nhl <- unique(df.mod.peptide$Number.Heavy.AA)
            for(i in seq_along(unique(df.mod.peptide$Number.Heavy.AA))){
              IDP.temp <- df.mod.peptide %>%
                filter(Number.Heavy.AA==nhl[i]) %>%
                pull(Isotope.Dot.Product) %>%
                unique() %>%
                as.numeric() %>%
                replace_na(.,0) %>%
                median() %>% # take the median here ensures IDP.temp is only one single quantity, which makes a one-to-one correspondence with number heavy AA
                as.numeric()
              IDP <- append(IDP, IDP.temp)
            }
            df.solutions[rows.write.out, "Isotope.Dot.Product"] <- IDP
          } # end else -- no solution

        } #end for - modified.peptide.sequence level
      } #end for - peptide level
    } #end for - protein level

    # trim any extra rows in the data frame past the counter, since there may be extranneous rows at the time of the initialization of df.solutions
    # this will leave a row of NA's at the end of df.solutions, since we've increased the counter at the end of each iteration of the loop.
    # we can trim this last row off also by doing na.omit() or we can instead increase the counter at the beginning of each iteration of the loop.
    df.solutions <- df.solutions[1:(counter-1),]

    # write out
    write.csv(df.solutions, file="df_solutions_date.csv", row.names = FALSE)
    #------------------------------------------------------------------------------------


    #------------------------------------------------------------------------------------
    # Clean Up df.solutions

    # remove rows which have 'NA' for FBC solution
    df.solutions <- df.solutions %>%
      filter(!is.na(FBC.Solution))

    # clean up df.solutions - for successive R scripts
    df.solutions <- df.solutions %>%
      dplyr::filter(Number.Heavy.AA<5) %>% # keep rows with at most 4 heavy AA
      dplyr::rename(Cohort=Condition) %>% # rename Condition to Cohort, this matches the nomenclature of the successive R scripts
      dplyr::rename(Modified.Peptide.Seq=Modified.Peptide) %>% # adjust this name for future steps
      dplyr::mutate(Condition=paste0(Cohort, "_D", Timepoint)) %>% # create new Condition column, by merging cohort with timepoint and including "D" for day
      dplyr::mutate(Total.Replicate.Name=paste0(Condition, "_", Replicate.Name)) %>% # create new Total.Replicate.Name column
      dplyr::mutate(Area.Number.Heavy.AA=paste0("Area", Number.Heavy.AA)) %>% # create new column for casting to wide format later
      dplyr::mutate(FBC.Solution=ifelse(FBC.Solution<0, 0, ifelse(FBC.Solution>1, 1, FBC.Solution))) # FBC solution must be [0,1]: negative values are changed to 0, greater than 1 are changed to 1

    # store max number of heavy AA in df.solutions
    max.num.heavy.AA <- max(df.solutions$Number.Heavy.AA)

    # column names of heavy labeled areas
    # retain up to Area4 at most (any heavier labelled areas are unnecessary)
    # if max.num.heavy.AA is greater than 4 then retain Area0 to Area4
    # else retain up to max.num.heavy.AA
    if(max.num.heavy.AA>4){
      heavy.label.areas <- paste0("Area", 1:4)
    } else{heavy.label.areas <- paste0("Area", 1:max.num.heavy.AA)}
    #------------------------------------------------------------------------------------


    #------------------------------------------------------------------------------------
    # Isotope Dot Product Filter
    # IDP.threshold can be between [0,1) where 1 is most stringent.
    # The default should be 0, corresponding to no filter, thereby retaining all of the data.
    IDP.threshold <- 0.0 # value for filtering by Isotope Dot Product
    df.solutions.filtered <- df.solutions %>%
      filter(Isotope.Dot.Product > IDP.threshold)
    #------------------------------------------------------------------------------------


    #------------------------------------------------------------------------------------
    # Rearrange variables in df.solutions.filtered
    # and write out this filtered version
    df.solutions.filtered <- df.solutions.filtered %>%
      select(Protein.Gene, Protein.Accession, Peptide, Modified.Peptide.Seq, Total.Replicate.Name, Condition, Cohort,
             Replicate.Name, Timepoint, Product.Charge, Number.Heavy.AA, Area.Number.Heavy.AA,
             everything()) # all other variables

    # write out
    write.csv(df.solutions.filtered, file=paste("df_solutions_filtered_",Sys.Date(),"_",paste(AA),".csv",sep=""), row.names = FALSE)
    #------------------------------------------------------------------------------------



    #------------------------------------------------------------------------------------
    # Create Area data frame with FBC solutions cast in wide format for each peptide
    df.areas.charge <- dcast(data = df.solutions.filtered, formula = Protein.Accession + Protein.Gene + Replicate.Name + Total.Replicate.Name + Cohort + Condition + Timepoint + Modified.Peptide.Seq + Product.Charge ~ Area.Number.Heavy.AA , value.var=c("FBC.Solution")) %>% # wide format cast
      mutate("Number.AA"= map_dbl(Modified.Peptide.Seq, AA.count.fun)) %>% # add Number of AA column
      filter(!is.na(Area0)) %>% # keep only rows with a non-NA Area0
      mutate(Percent.Label = rowSums(select(., all_of(heavy.label.areas)), na.rm=TRUE)/rowSums(select(., contains("Area")), na.rm=TRUE)) # create new column `Percent.Label`=Label Area/Total Area

    # split off the 1 AA data; we won't be able to calculate individual `precursor pool` with this data but we can
    # back calculate the `% new synthesized` once we have calculated the `median precursor pools` for each Condition with the other data containing multiple AA
    df.areas.one.l <- df.areas.charge %>%
      filter(Number.AA==1 & !is.na(Area0) & !is.na(Area1)) %>% # keep all peptides with 1 AA where both Area0 and Area1 are not NA
      select(-c(heavy.label.areas[-1])) # drop heavy AA columns above Area1

    # for the rest of the data, keep peptides with 2, 3 or 4 AA
    # this data will be used to calculate individual `precursor pool` for each peptide
    # and then we will calculate median precursor pool by condition
    df.areas.charge <- df.areas.charge %>%
      filter(Number.AA>1 & Number.AA<5) # only keep data with 2-4 heavy AA


    #if the max number of AA is 1, then turn off precursor pool enrichment
    if(max.num.heavy.AA < 2){
      PP_Enrich <- FALSE
    }


    #This version of the pipeline runs if either the maximum number of heavy AA < 2, or if Precursor Pool Enrichment is Off
    if(PP_Enrich == FALSE){
      # Calculate Percent Newly Synthesized for each Peptide
      # using peptides with 2,3 or 4 AA.
      # Average.Turnover.Score: 0 best, 1 worst.


      #run this if max.num.AA > 1
      if(max.num.heavy.AA > 1){

        # initialize new data frame to contain % newly synthesized values without precursor pool correction
        df.precursor.pool <- df.areas.charge %>%
          cbind("Precursor.Pool"=NA)

        # Loop
        perc.new.synth <- rep(NA, length(df.precursor.pool$Condition)) # initialize vector for storing % New Synthesized
        #avg.turn.score <- rep(NA, length(df.precursor.pool$Condition)) # initialize vector for storing Average Turnover Score
        for(i in 1:length(df.precursor.pool$Condition)){
          df.loop <- df.precursor.pool[i, ] # subset data for this iteration
          num.AA <- df.loop$Number.AA # number of AA for control-flow
          print(i, num.AA)


          if(num.AA==2){ # two AA
            # Percent New Synthesied:
            # calculate new component of the Observed Unlabeled-Area: Obs.Area0.new = Th0*(Obs1+Obs2)/(Th1+Th2) .... where Obs is observed and Th is theoretical
            #O0new <- df.binom2[median.pp, 1]*(df.loop[, "Area1"] + df.loop[, "Area2"])/(df.binom2[median.pp, 2] + df.binom2[median.pp, 3])
            # calculate percent new synthesized: % new synthesized = (Obs1 + Obs2)/(Obs0 + Obs1 + Obs2)
            perc.new.synth[i] <- (df.loop[, "Area1"] + df.loop[, "Area2"])/(df.loop[, "Area0"] + df.loop[, "Area1"] + df.loop[, "Area2"])

          } else if(num.AA==3){ # three AA
            # Percent New Synthesied:
            #O0new <- df.binom3[median.pp, 1]*(df.loop[, "Area1"] + df.loop[, "Area2"] + df.loop[, "Area3"])/(df.binom3[median.pp, 2] + df.binom3[median.pp, 3] + df.binom3[median.pp, 4])
            perc.new.synth[i] <- (df.loop[, "Area1"] + df.loop[, "Area2"] + df.loop[, "Area3"])/(df.loop[, "Area0"] + df.loop[, "Area1"] + df.loop[, "Area2"] + df.loop[, "Area3"])

          } else if(num.AA==4){ # else four AA
            # Percent New Synthesied:
            #O0new <- df.binom4[median.pp, 1]*(df.loop[, "Area1"] + df.loop[, "Area2"] + df.loop[, "Area3"] + df.loop[, "Area4"])/(df.binom4[median.pp, 2] + df.binom4[median.pp, 3] + df.binom4[median.pp, 4] + df.binom4[median.pp, 5])
            perc.new.synth[i] <- (df.loop[, "Area1"] + df.loop[, "Area2"] + df.loop[, "Area3"] + df.loop[, "Area4"])/(df.loop[, "Area0"] + df.loop[, "Area1"] + df.loop[, "Area2"] + df.loop[, "Area3"] + df.loop[, "Area4"])

          } else{ # ELSE - this case should not happen, since there should only be 2, 3 or 4 AA in the df.precursor.pool data
            print("This should not happen !!!")
            perc.new.synth[i] <- NA
            avg.turn.score[i] <- NA
          }
        } # end for


        # add new columns onto df.precursor.pool
        df.precursor.pool <- df.precursor.pool %>%
          cbind("Perc.New.Synth"=perc.new.synth)
        #cbind("Avg.Turnover.Score"=1-avg.turn.score) # taking the complement of avg.turn.score; now 1 is the best score, 0 is the worst


        # write out
        write.csv(df.precursor.pool, paste("Step1_Data_Output_Skyline_multiAA_peps_",Sys.Date(),"_",paste(AA),".csv", sep=""), row.names = FALSE)
        #end loop checking if max.num.heavy.AA > 1
      }



      # Calculate Percent Newly Synthesized for Peptides with 1 AA
      # we are doing this now that the median precursor pool for each Condition has been established
      # since 1 AA peptides cannot be used to calculate individual precursor pool

      perc.new.synth.one.l <- rep(NA, length(df.areas.one.l$Condition)) # initialize vector for storing % New Synthesized
      for(i in 1:length(df.areas.one.l$Condition)){
        print(i)

        df.loop <- df.areas.one.l[i, ] # subset data for this iteration

        # Percent New Synthesized:
        num.AA <- df.loop$Number.AA
        if(num.AA==1){ # one AA
          # calculate percent new synthesized: % new synthesized = (Obs1 + Obs2)/(Obs0 + Obs1 + Obs2)
          perc.new.synth.one.l[i] <- (df.loop[, "Area1"])/(df.loop[, "Area0"] + df.loop[, "Area1"])
          print((df.loop[, "Area1"])/(df.loop[, "Area0"] + df.loop[, "Area1"]))
        }
      }

      # add column onto data frame
      df.areas.one.l <- df.areas.one.l %>%
        cbind("Perc.New.Synth"=perc.new.synth.one.l)

      ####arrange depreciated in Dplyr

      #df.areas.one.l <- df.areas.one.l %>%
      # arrange(Timepoint)

      #arrange by timepoint using order
      df.areas.one.l <- df.areas.one.l[order(df.areas.one.l$Timepoint),]

      # write out
      #write.csv(df.areas.one.l, "Step1_Data_Output_Skyline_singleAA_peps_date.csv", row.names = FALSE)
      #write out with heavy label identified
      write.csv(df.precursor.pool, paste("Step1_Data_Output_Skyline_singleAA_peps_",Sys.Date(),"_",paste(AA),".csv", sep=""), row.names = FALSE)
      #end Percent Newly Synthesized without PP-Enrichment
    }




    #find percent newly synthesized with PP-Enrichment
    if(PP_Enrich == TRUE){
      # Calculate Percent Newly Synthesized and Average Turnover Score for each Peptide
      # using peptides with 2,3 or 4 AA.
      # Average.Turnover.Score: 0 best, 1 worst.

      ###
      ###
      ###Begin calculation of Precursor Pool 1/16/23

      # Theoretical Distrubtion which we'll be matching our observed data against in order to calculate Precursor Pool
      # Make the binomial probablity distribution in R [for peptides with 1-4 AA, which is most of the data]
      # n is the number of trials (aka number of AA in a peptide which will either be Unlabeled or Labeled)

      # Make Binomial Distributions

      # n=1, corresponds to peptide with 1 AA
      df.binom1 <- as.data.frame( matrix(data=NA, nrow=100, ncol=2))
      colnames(df.binom1) <- c("1L_L", "1L_1H")
      for(p in 1:100){
        df.binom1[p, ] <- dbinom(0:1, 1, p/100) # 1 trials, probability of success= p/100 = p%
      }

      # n=2, corresponds to peptide with 2 heavy AA
      df.binom2 <- as.data.frame( matrix(data=NA, nrow=100, ncol=3))
      colnames(df.binom2) <- c("2L_L", "2L_1H", "2L_2H")
      for(p in 1:100){
        df.binom2[p, ] <- dbinom(0:2, 2, p/100) # 2 trials, probability of success= p/100 = p%
      }

      # n=3, corresponds to peptide with 3 AA
      df.binom3 <- as.data.frame( matrix(data=NA, nrow=100, ncol=4))
      colnames(df.binom3) <- c("3L_L", "3L_1H", "3L_2H", "3L_3H")
      for(p in 1:100){
        df.binom3[p, ] <- dbinom(0:3, 3, p/100) # 3 trials, probability of success= p/100 = p%
      }

      # n=4, corresponds to peptide with 3 AA
      df.binom4 <- as.data.frame( matrix(data=NA, nrow=100, ncol=5))
      colnames(df.binom4) <- c("4L_L", "4L_1H", "4L_2H", "4L_3H", "4L_4H")
      for(p in 1:100){
        df.binom4[p, ] <- dbinom(0:4, 4, p/100) # 3 trials, probability of success= p/100 = p%
      }


      # Normalize the full distributions (both the unlabeled and labeled parts)
      df.binom1.norm <- as.matrix(apply(df.binom1, 2, as.numeric)) %>%
        apply(1, sum.to.one.fun) %>%
        t()

      df.binom2.norm <- as.matrix(apply(df.binom2, 2, as.numeric)) %>%
        apply(1, sum.to.one.fun) %>%
        t()

      df.binom3.norm <- as.matrix(apply(df.binom3, 2, as.numeric)) %>%
        apply(1, sum.to.one.fun) %>%
        t()

      df.binom4.norm <- as.matrix(apply(df.binom4, 2, as.numeric)) %>%
        apply(1, sum.to.one.fun) %>%
        t()


      # Normalize the heavy-components of the distributions (only the heavy-labeled parts)
      df.binom2.heavynorm <- df.binom2 %>%
        select(-c(1)) %>% # drop unlabeled column
        apply(1, sum.to.one.fun) %>%
        t()

      df.binom3.heavynorm <- df.binom3 %>%
        select(-c(1)) %>% # drop unlabeled column
        apply(1, sum.to.one.fun) %>%
        t()

      df.binom4.heavynorm <- df.binom4 %>%
        select(-c(1)) %>% # drop unlabeled column
        apply(1, sum.to.one.fun) %>%
        t()

      # For each Replicate.Name calculate the Precursor Pool across all proteins
      # For a given Replicate.Name, proteins are unique by two factors: modified.peptide.sequence and the number of charges ("+" signs)

      # initialize vectors for looping
      reps <- unique(df.areas.charge$Replicate.Name)

      # initialize new data frame for writing out precursor pool
      df.precursor.pool <- df.areas.charge %>%
        cbind("Precursor.Pool"=NA)

      for(i in 1:length(unique(df.areas.charge$Replicate.Name))){
        df.iloop <- filter(df.areas.charge, Replicate.Name==reps[i]) # subset data by replicate
        peps <- unique(df.iloop$Modified.Peptide.Seq)
        for(j in 1:length(peps)){
          df.jloop <- filter(df.iloop, Modified.Peptide.Seq==peps[j]) # subset data from i-loop by peptide
          charges <- sort(unique(df.jloop$Product.Charge)) # vector of charges for looping through in k-loop
          for(k in 1:length(charges)){
            # row index, writing precursor pool result out to the correct row in data frame df.areas.charge
            row.index <- which(df.areas.charge$Replicate.Name==reps[i] & df.areas.charge$Modified.Peptide.Seq==peps[j] & df.areas.charge$Product.Charge==charges[k])
            print(c(i,j,k,row.index))

            df.loop <- filter(df.jloop, Product.Charge==charges[k]) # subset data from j-loop based on charge state

            num.AA <- df.loop$Number.AA # store number of AA - for use in if statements

            # Declare P, the heavy-labelled data which will match against theoretical distribution using distmat (Euclidean distance function)
            P <- df.loop %>%
              select(heavy.label.areas) %>% # keep only heavy label areas
              discard(is.na(.)) %>% # drop all NA elements
              sum.to.one.fun() %>% # normalize data
              as.numeric(.) %>% # force numeric
              as.vector(.) # force vector

            if( (sum(!is.na(P)))!=num.AA){ # do not proceed with turnover calculation if any components of the Data P are NA
              indx <- NA
              print("Not Enough Data to proceed with Precusor Pool matching")
            } else if(num.AA==2){ # for 2 AA, calculate precursor pool, using 2-AA distribution
              Q <- df.binom2.heavynorm # Q is the theoretical distribution to match against
              distances <- distmat(Q, P)  # Euclidean distance between Q and P, returns a matrix, the columns are the distances of each slice of Q to P
              indx <- which(distances==min(distances))
            } else if(num.AA==3){ # for 3 AA, calculate precursor pool, using 3-AA distribution
              Q <- df.binom3.heavynorm # Q is the theoretical distribution to match against
              distances <- distmat(Q, P)  # Euclidean distance between Q and P, returns a matrix , the columns are the distances of each slice of Q to P
              indx <- which(distances==min(distances))
            } else if(num.AA==4) { # for 4 AA, calculate precursor pool, using 4-AA distribution
              Q <- df.binom4.heavynorm # Q is the theoretical distribution to match against
              distances <- distmat(Q, P) # Euclidean distance between Q and P, returns a matrix , the columns are the distances of each slice of Q to P
              indx <- which(distances==min(distances))
            } else {
              print("Else???") # should never be any 'else' cases
            }
            df.precursor.pool[row.index, "Precursor.Pool"] <- as.numeric(indx) # write out the Precursor Pool to data frame
          } # end for k-loop
        } # end for j-loop
      } # end for i-loop

      # add up all Areas (Area0 to max Area)
      df.precursor.pool <- df.precursor.pool %>%
        mutate(Total.Area = rowSums(select(., contains("Area")), na.rm=TRUE))

      # write out Precursor Pool Data frame
      write.csv(df.precursor.pool, "Precursor_Pool_date.csv", row.names = FALSE)

      # Calculate Precursor Pool % for each Condition.
      # individual precursor pool has already been calculated and stored in the data frame df.precursor.pool
      # here we will group by Replicate.Name and calculate the median precursor pool for each group
      # then we will group by Condition and calculate the final median precursor pool for each condition

      df.pp.medians <- df.precursor.pool %>%
        dplyr::group_by(Total.Replicate.Name) %>% # first group by Total Replicate Name and calculate median precursor pool
        mutate(Total.Replicate.Median.PP=median(Precursor.Pool, na.rm = TRUE)) %>%
        dplyr::group_by(Condition) %>% # group by Condition and calculate median precursor pool
        mutate(Precursor.Pool=median(Total.Replicate.Median.PP)) %>%
        ungroup() %>%
        arrange(Timepoint) %>% # arrange in ascending order by timepoint
        select(Condition, Precursor.Pool) %>% # keep only these variables
        unique() # keep only the unique conditions and their median precursor pool values

      # write out Precursor Pool Data frame
      write.csv(df.pp.medians, "Precursor_Pool_PPmedians_date.csv", row.names = FALSE)

      ###
      ###
      ###End calculation of Precursor Pool 1/16/23

      # Loop
      perc.new.synth <- rep(NA, length(df.precursor.pool$Condition)) # initialize vector for storing % New Synthesized
      avg.turn.score <- rep(NA, length(df.precursor.pool$Condition)) # initialize vector for storing Average Turnover Score
      for(i in 1:length(df.precursor.pool$Condition)){
        df.loop <- df.precursor.pool[i, ] # subset data for this iteration
        num.AA <- df.loop$Number.AA # number of AA for control-flow
        print(i, num.AA)


        # get median Precursor Pool based on Condition
        median.pp <- df.pp.medians %>%
          dplyr::filter(as.character(df.loop$Condition)==df.pp.medians[,"Condition"]) %>% # filter based on Condition
          dplyr::pull(Precursor.Pool) %>% # take only the value
          round(digits=0) # round to nearest whole number since we will use median.pp to grab from a certain row which should be an integer


        if(num.AA==2){ # two AA
          # Percent New Synthesied:
          # calculate new component of the Observed Unlabeled-Area: Obs.Area0.new = Th0*(Obs1+Obs2)/(Th1+Th2) .... where Obs is observed and Th is theoretical
          O0new <- df.binom2[median.pp, 1]*(df.loop[, "Area1"] + df.loop[, "Area2"])/(df.binom2[median.pp, 2] + df.binom2[median.pp, 3])
          # calculate percent new synthesized: % new synthesized = (Obs0new + Obs1 + Obs2)/(Obs0 + Obs1 + Obs2)
          perc.new.synth[i] <- (O0new + df.loop[, "Area1"] + df.loop[, "Area2"])/(df.loop[, "Area0"] + df.loop[, "Area1"] + df.loop[, "Area2"])

          # Average Turnover Score:
          # heavy labelled components of observed data
          P <- select(df.loop, c("Area1", "Area2")) %>%
            sum.to.one.fun() %>% # normalize data
            as.numeric(.) %>%
            as.vector(.)
          # theoretical distribution -  heavy labelled normalized only
          Q <- df.binom2.heavynorm[median.pp, c(1,2)]
          avg.turn.score[i] <- distmat(Q, P)

        } else if(num.AA==3){ # three AA
          # Percent New Synthesied:
          O0new <- df.binom3[median.pp, 1]*(df.loop[, "Area1"] + df.loop[, "Area2"] + df.loop[, "Area3"])/(df.binom3[median.pp, 2] + df.binom3[median.pp, 3] + df.binom3[median.pp, 4])
          perc.new.synth[i] <- (O0new + df.loop[, "Area1"] + df.loop[, "Area2"] + df.loop[, "Area3"])/(df.loop[, "Area0"] + df.loop[, "Area1"] + df.loop[, "Area2"] + df.loop[, "Area3"])

          # Average Turnover Score:
          # heavy labelled components of observed data
          P <- select(df.loop, c("Area1", "Area2", "Area3")) %>%
            sum.to.one.fun() %>% # normalize data
            as.numeric(.) %>%
            as.vector(.)
          # theoretical distribution -  heavy labelled normalized only
          Q <- df.binom3.heavynorm[median.pp, c(1,2,3)]
          avg.turn.score[i] <- distmat(Q, P)

        } else if(num.AA==4){ # else four AA
          # Percent New Synthesied:
          O0new <- df.binom4[median.pp, 1]*(df.loop[, "Area1"] + df.loop[, "Area2"] + df.loop[, "Area3"] + df.loop[, "Area4"])/(df.binom4[median.pp, 2] + df.binom4[median.pp, 3] + df.binom4[median.pp, 4] + df.binom4[median.pp, 5])
          perc.new.synth[i] <- (O0new + df.loop[, "Area1"] + df.loop[, "Area2"] + df.loop[, "Area3"] + df.loop[, "Area4"])/(df.loop[, "Area0"] + df.loop[, "Area1"] + df.loop[, "Area2"] + df.loop[, "Area3"] + df.loop[, "Area4"])

          # Average Turnover Score:
          # heavy labelled components of observed data
          P <- select(df.loop, c("Area1", "Area2", "Area3", "Area4")) %>%
            sum.to.one.fun() %>% # normalize data
            as.numeric(.) %>%
            as.vector(.)
          # theoretical distribution -  heavy labelled normalized only
          Q <- df.binom4.heavynorm[median.pp, c(1,2,3,4)]
          avg.turn.score[i] <- distmat(Q, P)

        } else{ # ELSE - this case should not happen, since there should only be 2, 3 or 4 AA in the df.precursor.pool data
          print("This should not happen !!!")
          perc.new.synth[i] <- NA
          avg.turn.score[i] <- NA
        }
      } # end for

      # add new columns onto df.precursor.pool
      df.precursor.pool <- df.precursor.pool %>%
        cbind("Perc.New.Synth"=perc.new.synth) %>%
        cbind("Avg.Turnover.Score"=1-avg.turn.score) # taking the complement of avg.turn.score; now 1 is the best score, 0 is the worst

      # write out
      write.csv(df.precursor.pool, "Step1_Data_Output_Skyline_multiAA_peps_date.csv", row.names = FALSE)



      # Calculate Percent Newly Synthesized for Peptides with 1 AA
      # we are doing this now that the median precursor pool for each Condition has been established
      # since 1 AA peptides cannot be used to calculate individual precursor pool

      perc.new.synth.one.l <- rep(NA, length(df.areas.one.l$Condition)) # initialize vector for storing % New Synthesized
      for(i in 1:length(df.areas.one.l$Condition)){
        print(i)

        df.loop <- df.areas.one.l[i, ] # subset data for this iteration

        # get median Precursor Pool based on Condition
        median.pp <- df.pp.medians %>%
          dplyr::filter(as.character(df.loop$Condition)==df.pp.medians[,"Condition"]) %>% # filter based on Condition
          dplyr::pull(Precursor.Pool) %>% # take only the value
          round(digits=0) # round to nearest whole number since we will use median.pp to grab from a certain row

        # Percent New Synthesied:
        # calculate new component of the Observed Unlabeled-Area: Obs.Area0.new = Th0*(Obs1)/(Th1) .... where Obs is observed and Th is theoretical
        O0new <- df.binom1[median.pp, 1]*df.loop[, "Area1"]/df.binom1[median.pp, 2]
        # calculate percent new synthesized: % new synthesized = (Obs0new + Obs1 + Obs2)/(Obs0 + Obs1 + Obs2)
        perc.new.synth.one.l[i] <- (O0new + df.loop[, "Area1"])/(df.loop[, "Area0"] + df.loop[, "Area1"])
      }

      # add column onto data frame
      df.areas.one.l <- df.areas.one.l %>%
        cbind("Perc.New.Synth"=perc.new.synth.one.l)


      ####
      df.areas.one.l <- df.areas.one.l %>%
        arrange(Timepoint)

      # write out
      write.csv(df.areas.one.l, "Step1_Data_Output_Skyline_singleAA_peps_date.csv", row.names = FALSE)


    } #End Calculation of Newly Synthesized Peptides using precursor pool enrichment




    # first relevel factors for Condition variable, for the data frames which we will use to plot data:

    if(max.num.heavy.AA > 1){
      # unique cohorts
      cohorts <- unique(df.precursor.pool$Cohort)

      #pull unique and ordered timepoints
      order_df <- df.precursor.pool[order(df.precursor.pool$Timepoint),]
      timepoints <- unique(order_df$Timepoint)

      #this block no longer works, duplicate name error so replaced with block above 1/16/23
      # unique timepoints - in ascending order
      #timepoints <- df.precursor.pool %>%
      # arrange(Timepoint) %>%
      #pull(Timepoint) %>%
      #unique()
    }

    if(max.num.heavy.AA < 2){
      # unique cohorts
      cohorts <- unique(df.areas.one.l$Cohort)

      # unique timepoints - in ascending order
      timepoints <- df.areas.one.l %>%
        arrange(Timepoint) %>%
        pull(Timepoint) %>%
        unique()
    }

    conditions.relevel <- c() # initialize vector
    # loop through timepoints, creating conditions.relevel vector with all cohorts along the way
    for(i in timepoints){
      conditions.relevel <- append(conditions.relevel, paste0(cohorts, sep="_D", i))
    }

    # relevel Condition variable in these three data frames:
    if (max.num.heavy.AA > 1){
      df.pp.medians <- df.pp.medians %>%
        mutate(Condition = fct_relevel(Condition, conditions.relevel))

      df.precursor.pool <- df.precursor.pool %>%
        mutate(Condition = fct_relevel(Condition, conditions.relevel))
    }
    df.areas.one.l <- df.areas.one.l %>%
      mutate(Condition = fct_relevel(Condition, conditions.relevel))
    ###



    # Plots:

    ### Peptides with one AA:
    # Box Plots of Percent Newly Synthesized by Condition
    # reorder Conditions using {Forcats} function fct_relevel
    boxplot.oneAA.percentnewsynth <- df.areas.one.l %>%
      ggplot(aes(y=Perc.New.Synth, fill=Condition)) + # optional: use linetype=group to use different linetypes
      geom_boxplot(aes(y=Perc.New.Synth, col=Condition, alpha=0.1)) +
      labs(title=paste("Percent Newly Sythnesized Distribution by Condition - Peptides with 1 ",paste(AA), sep=""), x="Condition", y="Percent Newly Synthesized") +
      theme_bw()

    # save plot
    ggsave(paste("Boxplot_Percent-Newly-Synthesized_single-AA-peptides_",paste(AA),".pdf", sep=""),
           plot = boxplot.oneAA.percentnewsynth,
           width = 7, height = 5,
           units = "in", # inches
           dpi = 300)



    #### Peptides with multiple AA:
    if (max.num.heavy.AA > 1 & PP_Enrich == TRUE){
      # Density curves - Precursor Pool
      # all Conditions in the same plot
      density.groups <- df.precursor.pool %>%
        ggplot(aes(x=Precursor.Pool, fill=Condition)) +
        geom_histogram(aes(y=..density..), alpha=0.2, col="black", position='identity') +
        geom_density(alpha=0.2) +
        geom_vline(data=df.pp.medians, aes(xintercept=Precursor.Pool, col=Condition), linetype="dashed", show.legend=FALSE) +
        labs(title="Precursor Pool Distrubtion by Condition", x="Precursor Pool", y="Density") +
        theme_bw()

      # Density curves - Precursor Pool
      # Facet by Condition
      density.precursor.pool <- df.precursor.pool %>%
        ggplot(aes(x=Precursor.Pool, fill=Condition)) +
        geom_histogram(aes(y=..density..), alpha=0.2, col="black", position='identity') +
        geom_density(alpha=0.2) +
        geom_vline(data=df.pp.medians, aes(xintercept=Precursor.Pool, col=Condition), linetype="dashed", show.legend=FALSE) +
        facet_wrap(~ Condition, ncol = 2, scales = "fixed") +
        labs(title="Precursor Pool Distrubtion by Condition", x="Precursor Pool", y="Density") +
        theme_bw()

      # save plot
      ggsave(paste("Density_Precursor-Pool_",paste(AA),".pdf", sep=""),
             plot = density.precursor.pool,
             width = 7, height = 5,
             units = "in", # inches
             dpi = 300)


      #change Na values in percent new synth column of df.precursor.pool to 0
      df.precursor.pool[["Perc.New.Synth"]][is.na(df.precursor.pool[["Perc.New.Synth"]])] <-0
    }

    if (max.num.heavy.AA > 1) {
      # Density Curves - Percent Newly Synthesized
      # Facet by Condition
      density.percent.new.synthesized <- df.precursor.pool %>%
        ggplot(aes(x=Perc.New.Synth, fill=Condition)) +
        geom_histogram(aes(y=..density..), alpha=0.2, col="black", position='identity') +
        geom_density(alpha=0.2) +
        #geom_vline(data=df.pp.medians, aes(xintercept=Precursor.Pool, col=Condition), linetype="dashed", show.legend=FALSE) +
        facet_wrap(~ Condition, ncol = 2, scales = "fixed") +
        labs(title="Percent Newly Synthesized Distrubtion by Condition", x="Percent Newly Synthesized", y="Density") +
        theme_bw()

      # save plot
      ggsave(paste("Density_Percent-Newly-Synthesized_",paste(AA),".pdf", sep=""),
             plot = density.percent.new.synthesized,
             width = 7, height = 5,
             units = "in", # inches
             dpi = 300)
    }

    if (max.num.heavy.AA > 1 & PP_Enrich == TRUE){
      # Box Plots  Precursor Pool by Condition
      boxplot.precursor.pool <- df.precursor.pool %>%
        ggplot(aes(y=Precursor.Pool, fill=Condition)) + # optional: use linetype=group to use different linetypes
        geom_boxplot(aes(y=Precursor.Pool, col=Condition, alpha=0.1)) +
        labs(title="Precursor Pool Distribution by Condition", x="Condition", y="Precursor Pool") +
        theme_bw()

      # save plot
      ggsave(paste("Boxplot_Precursor-Pool_",paste(AA),".pdf", sep=""),
             plot = boxplot.precursor.pool,
             width = 7, height = 5,
             units = "in", # inches
             dpi = 300)

      # Box Plots - Average Turnover Score by Condition
      boxplot.avg.turnover.score <- df.precursor.pool %>%
        ggplot(aes(y=Avg.Turnover.Score, fill=Condition)) + # optional: use linetype=group to use different linetypes
        geom_boxplot(aes(y=Avg.Turnover.Score, col=Condition, alpha=0.1)) +
        labs(title="Average Turnover Score Distribution by Condition", x="Condition", y="Average Turnover Score") +
        theme_bw()

      # save plot
      ggsave(paste("Boxplot_Average-Turnover-Score_",paste(AA),".pdf", sep=""),
             plot = boxplot.avg.turnover.score,
             width = 7, height = 5,
             units = "in", # inches
             dpi = 300)
    }

    if (max.num.heavy.AA > 1) {
      # Box Plots - Percent Newly Synthesized by Condition
      boxplot.percent.newly.synthesized <- df.precursor.pool %>%
        ggplot(aes(y=Perc.New.Synth, fill=Condition)) + # optional: use linetype=group to use different linetypes
        geom_boxplot(aes(y=Perc.New.Synth, col=Condition, alpha=0.1)) +
        labs(title="Percent Newly Sythnesized Distribution by Condition", x="Condition", y="Percent Newly Synthesized") +
        theme_bw()

      # save plot
      ggsave(paste("Boxplot_percent-newly-synthesized_",paste(AA),".pdf", sep=""),
             plot = boxplot.percent.newly.synthesized,
             width = 7, height = 5,
             units = "in", # inches
             dpi = 300)
    }

    if(max.num.heavy.AA > 1 & PP_Enrich == TRUE){
      # Scatterplot - Percent New Synthesized vs. Precursor Pool
      # scatterplot.groups <- df.precursor.pool %>%
      #   ggplot(aes(x=Precursor.Pool, y=Perc.New.Synth, fill=Condition)) + # optional: use linetype=group to use different linetypes
      #   geom_point(aes(x=Precursor.Pool, y=Perc.New.Synth, col=Condition, alpha=0.1)) +
      #   geom_vline(data=df.pp.medians, aes(xintercept=Precursor.Pool, col=Condition), linetype="dashed", show.legend=FALSE) +
      #   labs(title="Percent New Synthesized vs. Precursor Pool", x="Precursor Pool", y="Percent New Synthesized") +
      #   theme_bw()

      # Scatterplot - Percent New Synthesized vs. Precursor Pool
      # Facet by Condition
      scatterplot.percent.new.synth.vs.precursor.pool <- df.precursor.pool %>%
        ggplot(aes(x=Precursor.Pool, y=Perc.New.Synth, fill=Condition)) + # optional: use linetype=group to use different linetypes
        geom_point(aes(x=Precursor.Pool, y=Perc.New.Synth, col=Condition, alpha=0.1)) +
        geom_vline(data=df.pp.medians, aes(xintercept=Precursor.Pool, col=Condition), linetype="dashed", show.legend=FALSE) +
        facet_wrap(~ Condition, ncol = 2, scales = "fixed") +
        labs(title="Percent New Synthesized vs. Precursor Pool", x="Precursor Pool", y="Percent New Synthesized") +
        theme_bw()

      # save plot
      ggsave(paste("Scatterplot_Percent-Newly-Synthesized_vs_Precursor-Pool_",paste(AA),".pdf", sep=""),
             plot = scatterplot.percent.new.synth.vs.precursor.pool,
             width = 7, height = 5,
             units = "in", # inches
             dpi = 300)

      # Percent New Synthesized vs. Average Turnover Score
      # avg.turnover.plot <- df.precursor.pool %>%
      #   ggplot(aes(x=Avg.Turnover.Score, y=Perc.New.Synth, fill=Condition)) + # optional: use linetype=group to use different linetypes
      #   geom_point(aes(x=Avg.Turnover.Score, y=Perc.New.Synth, col=Condition, alpha=0.1)) +
      #   labs(title="Percent New Synthesized vs. Average Turnover Score", x="Average Turnover Score", y="Percent New Synthesized") +
      #   theme_bw()

      # Percent New Synthesized vs. Average Turnover Score
      # Facet by Condition
      scatterplot.percent.new.synth.vs.avg.turnover.score <- df.precursor.pool %>%
        ggplot( aes(x=Avg.Turnover.Score, y=Perc.New.Synth, fill=Condition)) + # optional: use linetype=group to use different linetypes
        geom_point(aes(x=Avg.Turnover.Score, y=Perc.New.Synth, col=Condition, alpha=0.1)) +
        geom_hline(aes(yintercept=1)) +
        facet_wrap(~ Condition, ncol = 2, scales = "fixed") + # can do fixed or free scales
        labs(title="Percent New Synthesized vs. Average Turnover Score", x="Average Turnover Score", y="Percent New Synthesized") +
        theme_bw()

      # save plot
      ggsave(paste("Scatterplot_Percent-Newly-Synthesized_vs_Avg-Turnover-Score_",paste(AA),".pdf", sep=""),
             plot = scatterplot.percent.new.synth.vs.avg.turnover.score,
             width = 7, height = 5,
             units = "in", # inches
             dpi = 300)



      # Average Turnover Score Filter

      # first see the histogram of Average Turnover Score
      hist(df.precursor.pool$Avg.Turnover.Score, breaks=100, main="Average Turnover Score", xlab="Average Turnover Score")

      # average turnover score filter
      # between [0,1) where 1 is most stringent
      # the default should be 0
      ATS.threshold <- 0.70 # average turnover score value, 70% is a typically a good starting place

      df.pp.ats.filtered <- df.precursor.pool %>%
        filter(Avg.Turnover.Score>ATS.threshold)

      density.percnew <- df.pp.ats.filtered %>%
        # mutate(Condition = fct_relevel(Condition, "OCon_D3", "OCR_D3", "OCon_D7", "OCR_D7", "OCon_D12", "OCR_D12", "OCon_D17", "OCR_D17")) %>%
        mutate(Condition = fct_relevel(Condition, conditions.relevel)) %>%
        ggplot(aes(x=Perc.New.Synth, fill=Condition)) +
        geom_histogram(aes(y=..density..), alpha=0.2, col="black", position='identity') +
        geom_density(alpha=0.2) +
        facet_wrap(~ Condition, ncol = 2, scales = "fixed") +
        labs(title=paste("Percent Newly Synthesized by Condition; Average Turnover Score >", ATS.threshold) , x="Percent Newly Synthesized", y="Density") +
        theme_bw()
    }

    #end for loop
  }
}







