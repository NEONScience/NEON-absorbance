################################################################################
#' @title formatWide

#' @author
#' Robert Hensley \email{hensley@battelleecology.org} \cr

#' @description This function re-formats National Ecological Observatory 
#' Network (NEON) UV-Vis absorbance data from long format into wide format. 

#' @param absorbanceData User input of the table of NEON absorbance data [data.frame]

#' @return This function returns a table absorbance values in wide format.

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords NEON, DOC, SUVA, UV Absorbance

#' @examples
#' #Using an example file
#' #outputData <- formatWide(
#' #absorbanceData=swc_externalLabAbsorbanceScan)

#' @export

# changelog and author contributions / copyrights
################################################################################
formatWide<-function(
    absorbanceData){
  
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
  
  # Converts to wide format
  wideFormat <- reshape(
    absorbanceAveraged,
    idvar = c("sampleID"),        
    timevar = "wavelength",      
    v.names = "absorbance", 
    direction = "wide"     
  )
  
  # Formats output table
  outputTable<-unique(absorbanceAveraged[,c("domainID","siteID","sampleID","collectDate")])
  outputTable<-merge(outputTable,wideFormat,by.x="sampleID",by.y="sampleID",all.x=T,all.y=F)
  outputTable<-outputTable[, c(2,3,1,4:195)]

  return(outputTable)
}
