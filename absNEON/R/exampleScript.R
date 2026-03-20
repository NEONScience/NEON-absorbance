siteID<-"COMO"

#' Pulls surface water chemistry data from NEON data portal and loads tables into R environment
inputData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteID, startdate="2023-10",enddate="2025-09", 
                                      package="expanded", include.provisional=T, check.size = F)

#' Assigns data table names regardless of swc or gwc prefix
absorbanceData<-inputData[[grep("externalLabAbsorbanceScan",names(inputData))]]
concentrationData<-inputData[[grep("externalLabDataByAnalyte",names(inputData))]]

outputWide<-formatWide(absorbanceData)

outputSuva254<-calcSuva254(absorbanceData,concentrationData,correctFe=FALSE)
outputSuva350<-calcSuva254(absorbanceData,concentrationData,correctFe=FALSE)
outputE2E3<-calcE2E3(absorbanceData,concentrationData,correctFe=TRUE)
outputSR<-calcSR(absorbanceData,concentrationData,correctFe=TRUE)

