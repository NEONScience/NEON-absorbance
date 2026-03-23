


# This is an example script for using the functions in the absNEON package

siteID <- "COMO"

# Pulls surface water chemistry data from NEON data portal and loads tables into R environment
inputData <- neonUtilities::loadByProduct(
  dpID = "DP1.20093.001",
  site = siteID,
  startdate = "2023-10",
  enddate = "2025-09",
  package = "expanded",
  include.provisional = T,
  check.size = F
)

# Assigns data table names regardless of swc or gwc prefix
absorbanceData <- inputData[[grep("externalLabAbsorbanceScan", names(inputData))]]
concentrationData <- inputData[[grep("externalLabDataByAnalyte", names(inputData))]]

# Convert the absorbance scan data to wide format
outputWide <- absNEON::formatWide(absorbanceData)

# Use the functions in the package...
outputSuva254 <- absNEON::calcSuva254(absorbanceData, concentrationData, correctFe =
                                        FALSE)
outputSuva350 <- absNEON::calcSuva350(absorbanceData, concentrationData, correctFe =
                                        FALSE)
outputE2E3 <- absNEON::calcE2E3(absorbanceData, concentrationData, correctFe =
                                  TRUE)
outputSR <- absNEON::calcSR(absorbanceData, concentrationData, correctFe =
                              TRUE)

