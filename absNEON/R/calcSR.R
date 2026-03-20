################################################################################
#' @title calcSR

#' @author
#' Robert Hensley \email{hensley@battelleecology.org} \cr

#' @description This function calculates the spectral slope ratio (SR), 275-295 nm
#' divided by 350-400 nm from National Ecological Observatory Network (NEON) water
#' chemistry data.

#' @param absorbanceData User input of the table of NEON absorbance data [dataframe]
#' @param concentrationData User input of the table of NEON concentration data [dataframe]
#' @param correctFe User input of whether correction for overlapping absorbption of Fe(III)
#' should also be included. See README for more information. Defaults to FALSE. [boolean]

#' @import plyr
#' @import tidyverse

#' @return This function returns a table of SR values

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords NEON, DOC, SUVA, UV Absorbance

#' @examples
#' #Using an example file
#' #outputData <- calcSR(
#' absorbanceData=swc_externalLabAbsorbanceScan,
#' concentrationData=swc_externalLabDataByAnalyte,
#' correctFe=FALSE)

#' @export

# changelog and author contributions / copyrights
################################################################################
calcSR<-function(
    absorbanceData,
    concentrationData,
    correctFe=FALSE){
  
  if(nrow(absorbanceData) < 1){
    print("Error: No absorbance data")
    return(NULL)
  }
  
  #' Averages replicate absorbance scans
  absorbanceData$sampleID.wavelength<-paste(absorbanceData$sampleID, absorbanceData$wavelength, sep=".")
  absorbanceData<-plyr::ddply(absorbanceData,c("sampleID.wavelength"),summarise,sampleID=unique(sampleID),domainID=unique(domainID),siteID=unique(siteID),
                              collectDate=unique(collectDate),wavelength=unique(wavelength),absorbance=mean(decadicAbsorbance)) 
  absorbanceData$sampleID.wavelength<-NULL
  
  #' No concentration required if Fe correction not performed
  combinedData<-absorbanceData
  
  #' Adds Fe concentrations if input condition is true
  if(correctFe == TRUE){
    Fe<-concentrationData[(concentrationData$analyte=="Fe"),]
    if(nrow(Fe) < 1){
      print("Error: No concentration data")
      return(NULL)
    }
    Fe<-Fe[,c("sampleID","analyteConcentration")]
    colnames(Fe)<-c("sampleID","Fe")
    combinedData<-merge(absorbanceData, Fe,by.x="sampleID",by.y="sampleID")
    combinedData$absorbanceFe<-combinedData$Fe*((-0.00000044*(combinedData$wavelength^2))-(0.00007755*combinedData$wavelength)+0.11337671)
    combinedData$absorbanceCorrected<-combinedData$absorbance-combinedData$absorbanceFe
  }
  
  #' Performs calculation without any corrections
  #' Calculates slope from 275-295 nm
  slopes275 <- combinedData %>% 
    filter(wavelength >= 275 & wavelength <= 295) %>% 
    group_by(sampleID) %>% 
    nest() %>% 
    mutate(S = map(data, ~lm(log(absorbance) ~ wavelength, data = .x))) %>% 
    mutate(slope = map(S, ~tidy(.x))) %>% 
    unnest(slope) %>% 
    select(sampleID, term, estimate) %>% 
    pivot_wider(names_from = term, values_from = estimate) %>% 
    dplyr::rename("intercept275" = `(Intercept)`,"slope275" = wavelength) %>% 
    ungroup() %>% 
    mutate(slope275 = slope275*-1)
  #' Calculates slope from 350-400 nm
  slopes350 <- combinedData %>% 
    filter(wavelength >= 350 & wavelength <= 400) %>% 
    group_by(sampleID) %>% 
    nest() %>% 
    mutate(S = map(data, ~lm(log(absorbance) ~ wavelength, data = .x))) %>% 
    mutate(slope = map(S, ~tidy(.x))) %>% 
    unnest(slope) %>% 
    select(sampleID, term, estimate) %>% 
    pivot_wider(names_from = term,values_from = estimate) %>% 
    dplyr::rename("intercept350" = `(Intercept)`,"slope350" = wavelength) %>% 
    ungroup() %>% 
    mutate(slope350 = slope350*-1)
  #' Calculates spectral slope ratio
  SlopeRatios <- merge(x = slopes275, y = slopes350,by = "sampleID") %>%
    mutate(SR = slope275/slope350)
  
  #' Performs calculation with Fe correction
  if(correctFe == TRUE){
    #' Calculates slope from 275-295 nm using corrected absorbance
    slopes275Corrected <- combinedData %>% 
      filter(wavelength >= 275 & wavelength <= 295) %>% 
      group_by(sampleID) %>% 
      nest() %>% 
      mutate(S = map(data, ~lm(log(absorbanceCorrected) ~ wavelength, data = .x))) %>% 
      mutate(slope = map(S, ~tidy(.x))) %>% 
      unnest(slope) %>% 
      select(sampleID, term, estimate) %>% 
      pivot_wider(names_from = term, values_from = estimate) %>% 
      dplyr::rename("intercept275" = `(Intercept)`,"slope275" = wavelength) %>% 
      ungroup() %>% 
      mutate(slope275 = slope275*-1)
    #' Calculates slope from 350-400 nm using corrected absorbance
    slopes350Corrected <- combinedData %>% 
      filter(wavelength >= 350 & wavelength <= 400) %>% 
      group_by(sampleID) %>% 
      nest() %>% 
      mutate(S = map(data, ~lm(log(absorbanceCorrected) ~ wavelength, data = .x))) %>% 
      mutate(slope = map(S, ~tidy(.x))) %>% 
      unnest(slope) %>% 
      select(sampleID, term, estimate) %>% 
      pivot_wider(names_from = term,values_from = estimate) %>% 
      dplyr::rename("intercept350" = `(Intercept)`,"slope350" = wavelength) %>% 
      ungroup() %>% 
      mutate(slope350 = slope350*-1)
    #' Calculates spectral slope ratio
    SlopeRatiosCorrected <- merge(x = slopes275Corrected, y = slopes350Corrected,by = "sampleID") %>%
      mutate(SRCorrected = slope275/slope350)
  }
  
  #' Formats output table
  outputTable<-unique(absorbanceData[,c("domainID","siteID","sampleID","collectDate")])
  outputTable<-merge(outputTable,SlopeRatios,by.x="sampleID",by.y="sampleID",all.x=T,all.y=F)
  outputTable<-outputTable[,c("domainID","siteID","sampleID","collectDate","SR")]

  if(correctFe == TRUE){
    outputTable<-merge(outputTable,SlopeRatiosCorrected,by.x="sampleID",by.y="sampleID",all.x=T,all.y=F)
    outputTable<-outputTable[,c("domainID","siteID","sampleID","collectDate","SR","SRCorrected")]
  }
  
  return(outputTable)
}
