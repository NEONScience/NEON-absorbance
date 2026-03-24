################################################################################
#' @title calcSuvaX

#' @author
#' Robert Hensley \email{hensley@battelleecology.org} \cr

#' @description This function calculates the specific ultra-violet absorbance
#' (SUVA) for a user-specified wavelength from National Ecological Observatory
#' Network (NEON) water chemistry data.

#' @param absorbanceData User input of the table of NEON absorbance data [data.frame]
#' @param concentrationData User input of the table of NEON concentration data [data.frame]
#' @param wavlength User input of the wavelength to return data for. [numeric]
#' @param correctFe User input of whether correction for overlapping absorbption of Fe(III)
#' should also be included. See README for more information. Defaults to FALSE. [logical]

#' @import tidyverse

#' @return This function returns a table of SUVA254 values

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords NEON, DOC, SUVA, UV Absorbance

#' @examples
#' #Using an example file
#' #outputData <- calcSuvaX(
#' #absorbanceData=swc_externalLabAbsorbanceScan,
#' #concentrationData=swc_externalLabDataByAnalyte,
#' #wavelength=401,
#' #correctFe=FALSE)

#' @export

# changelog and author contributions / copyrights
################################################################################
calcSuvaX<-function(
    absorbanceData,
    concentrationData,
    wavelength,
    correctFe=FALSE){

  if(wavelength<220|wavelength>600){
    print("Error: Wavelength is outside 200-600nm range of absorbance scan")
    return(NULL)
  }

  # Isolates specified wavelength if even, or surrounding wavelengths if odd.
  if(wavelength %% 2 == 0){
    absorbanceData<-absorbanceData[(absorbanceData$wavelength==wavelength),]
  }
  if(wavelength %% 2 != 0){
    absorbanceData<-absorbanceData[(absorbanceData$wavelength==wavelength-1|absorbanceData$wavelength==wavelength+1),]
  }

  if(nrow(absorbanceData) < 1){
    print("Error: No absorbance data")
    return(NULL)
  }

  # Averages replicate absorbance scans
  absorbanceAveraged<-absorbanceData
  absorbanceAveraged<-absorbanceAveraged |> dplyr::group_by(sampleID) |> dplyr::summarize(
                domainID=unique(domainID),siteID=unique(siteID),collectDate=unique(collectDate),
                absorbance=mean(decadicAbsorbance))
  absorbanceAveraged$wavelength<-wavelength

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
  combinedData$suvaX<-combinedData$absorbance/combinedData$DOC*100
  # Corrects for Fe absorbance and calculates suva if input condition is true
  if(correctFe == TRUE){
    combinedData$absorbanceFe<-combinedData$Fe*((-0.00000044*(combinedData$wavelength^2))-(0.00007755*combinedData$wavelength)+0.11337671)
    combinedData$absorbanceCorrected<-combinedData$absorbance-combinedData$absorbanceFe
    combinedData$suvaXCorrected<-combinedData$absorbanceCorrected/combinedData$DOC*100
  }

  # Formats output table
  outputTable<-combinedData[,c("domainID","siteID","sampleID","collectDate","suvaX")]
  names(outputTable)<-c("domainID","siteID","sampleID","collectDate",paste0("suva",wavelength))
  if(correctFe == TRUE){
    outputTable<-combinedData[,c("domainID","siteID","sampleID","collectDate","suvaX","suvaXCorrected")]
    names(outputTable)<-c("domainID","siteID","sampleID","collectDate",paste0("suva",wavelength),paste0("suva",wavelength,"Corrected"))
  }

  return(outputTable)
}
