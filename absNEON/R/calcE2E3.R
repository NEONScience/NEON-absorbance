################################################################################
#' @title calcE2E3

#' @author
#' Robert Hensley \email{hensley@battelleecology.org} \cr

#' @description This function calculates the E2:E3 (250 nm : 365 nm) absorbance 
#' ratio from National Ecological Observatory Network (NEON) water
#' chemistry data.

#' @param absorbanceData User input of the table of NEON absorbance data [dataframe]
#' @param concentrationData User input of the table of NEON concentration data [dataframe]
#' @param correctFe User input of whether correction for overlapping absorbption of Fe(III)
#' should also be included. See README for more information. Defaults to FALSE. [boolean]

#' @import plyr

#' @return This function returns a table of E2:E3 values

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords NEON, DOC, SUVA, UV Absorbance

#' @examples
#' #Using an example file
#' #outputData <- calcE2E3(
#' absorbanceData=swc_externalLabAbsorbanceScan,
#' concentrationData=swc_externalLabDataByAnalyte,
#' correctFe=FALSE)

#' @export

# changelog and author contributions / copyrights
################################################################################
calcE2E3<-function(
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
    abs250<-combinedData[(combinedData$wavelength=="250"),]
    abs250<-abs250[,c("sampleID","absorbance")]
    names(abs250)[names(abs250) == "absorbance"] <- "abs250"
    abs365<-combinedData[(combinedData$wavelength=="364"|combinedData$wavelength=="366"),]
    abs365<-plyr::ddply(abs365,c("sampleID"),summarise,abs365=mean(absorbance)) #Estimates 365 from average of 364 and 366 
    E2E3<-merge(abs250,abs365,by.x="sampleID",by.y="sampleID",all.x=T,all.y=T)
    E2E3$E2E3<-E2E3$abs250/E2E3$abs365
  
  
  #' Performs calculation with Fe correction
  if(correctFe == TRUE){
    abs250Corrected<-combinedData[(combinedData$wavelength=="250"),]
    abs250Corrected<-abs250Corrected[,c("sampleID","absorbance")]
    names(abs250Corrected)[names(abs250Corrected) == "absorbance"] <- "abs250Corrected"
    abs365Corrected<-combinedData[(combinedData$wavelength=="364"|combinedData$wavelength=="366"),]
    abs365Corrected<-plyr::ddply(abs365Corrected,c("sampleID"),summarise,abs365Corrected=mean(absorbance)) 
    E2E3Corrected<-merge(abs250Corrected,abs365Corrected,by.x="sampleID",by.y="sampleID",all.x=T,all.y=T)
    E2E3Corrected$E2E3Corrected<-E2E3Corrected$abs250Corrected/E2E3Corrected$abs365Corrected
  }
  
  #' Formats output table
  outputTable<-unique(absorbanceData[,c("domainID","siteID","sampleID","collectDate")])
  outputTable<-merge(outputTable,E2E3,by.x="sampleID",by.y="sampleID",all.x=T,all.y=F)
  outputTable<-outputTable[,c("domainID","siteID","sampleID","collectDate","E2E3")]
  
  if(correctFe == TRUE){
    outputTable<-merge(outputTable,E2E3Corrected,by.x="sampleID",by.y="sampleID",all.x=T,all.y=F)
    outputTable<-outputTable[,c("domainID","siteID","sampleID","collectDate","E2E3","E2E3Corrected")]
  }
  
  return(outputTable)
}
