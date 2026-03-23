################################################################################
#' @title calcSuva254

#' @author
#' Robert Hensley \email{hensley@battelleecology.org} \cr

#' @description This function calculates the specific ultra-violet absorbance
#' (SUVA) at 254 nm from National Ecological Observatory Network (NEON) water
#' chemistry data.

#' @param absorbanceData User input of the table of NEON absorbance data [data.frame]
#' @param concentrationData User input of the table of NEON concentration data [data.frame]
#' @param correctFe User input of whether correction for overlapping absorbption of Fe(III)
#' should also be included. See README for more information. Defaults to FALSE. [logical]
 
#' @import tidyverse

#' @return This function returns a table of SUVA254 values

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords NEON, DOC, SUVA, UV Absorbance

#' @examples
#' #Using an example file
#' #outputData <- calcSuva254(
#' #absorbanceData=swc_externalLabAbsorbanceScan,
#' #concentrationData=swc_externalLabDataByAnalyte,
#' #correctFe=FALSE)

#' @export

# changelog and author contributions / copyrights
################################################################################
calcSuva254<-function(
    absorbanceData,
    concentrationData,
    correctFe=FALSE){
  
  # Keeps only wavelengths from full-spectrum scan required for specific function
  absorbanceData<-absorbanceData[(absorbanceData$wavelength=="254"),]
  if(nrow(absorbanceData) < 1){
    print("Error: No absorbance data")
    return(NULL)
  }
  
  # Averages replicate absorbance scans
  absorbanceAveraged<-absorbanceData
  absorbanceAveraged$sampleID.wavelength<-paste(absorbanceAveraged$sampleID, absorbanceAveraged$wavelength, sep=".")
  absorbanceAveraged<-absorbanceAveraged |> dplyr::group_by(sampleID.wavelength) |> dplyr::summarize(
                sampleID=unique(sampleID),domainID=unique(domainID),siteID=unique(siteID),
                collectDate=unique(collectDate),wavelength=unique(wavelength),absorbance=mean(decadicAbsorbance)) 

  absorbanceAveraged$sampleID.wavelength<-NULL
  
  #' Formats and adds older, discrete absorbance values from concentration table
  absorbanceDiscrete<-concentrationData[(concentrationData$analyte=="UV Absorbance (254 nm)"),]
  absorbanceDiscrete$wavelength<-254
  absorbanceDiscrete<-absorbanceDiscrete[,c("sampleID","domainID","siteID","collectDate","wavelength","analyteConcentration")]
  colnames(absorbanceDiscrete)<-c("sampleID","domainID","siteID","collectDate","wavelength","absorbance")
  absorbanceDiscrete <- absorbanceDiscrete[!(absorbanceDiscrete$sampleID %in% absorbanceAveraged$sampleID), ]
  absorbanceAveraged<-rbind(absorbanceDiscrete,absorbanceAveraged)
  
  
  # Adds DOC concentration data 
  DOC<-concentrationData[(concentrationData$analyte=="DOC"),]
  if(nrow(DOC) < 1){
    print("Error: No concentration data")
    return(NULL)
  }
  DOC<-DOC[,c("sampleID","analyteConcentration")]
  colnames(DOC)<-c("sampleID","DOC")
  combinedData<-merge(absorbanceAveraged, DOC,by.x="sampleID",by.y="sampleID")
  # Adds Fe concentrations if input condition is true
  if(correctFe == TRUE){
    Fe<-concentrationData[(concentrationData$analyte=="Fe"),]
    Fe<-Fe[,c("sampleID","analyteConcentration")]
    colnames(Fe)<-c("sampleID","Fe")
    combinedData<-merge(combinedData, Fe,by.x="sampleID",by.y="sampleID")
  }
  
  # Calculates specific ultra-violet absorbance (units - L/mg-m)
  combinedData$suva254<-combinedData$absorbance/combinedData$DOC*100 
  # Corrects for Fe absorbance and calculates suva if input condition is true
  if(correctFe == TRUE){
    combinedData$absorbanceFe<-combinedData$Fe*((-0.00000044*(combinedData$wavelength^2))-(0.00007755*combinedData$wavelength)+0.11337671)
    combinedData$absorbanceCorrected<-combinedData$absorbance-combinedData$absorbanceFe
    combinedData$suva254Corrected<-combinedData$absorbanceCorrected/combinedData$DOC*100
  }
  
  # Formats output table
  outputTable<-combinedData[,c("domainID","siteID","sampleID","collectDate","suva254")]
  if(correctFe == TRUE){
    outputTable<-combinedData[,c("domainID","siteID","sampleID","collectDate","suva254","suva254Corrected")]
  }
  
  return(outputTable)
}
