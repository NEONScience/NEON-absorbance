NEON Absorbance Package
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- ****** Description ****** -->

The NEON Absorbance Package (absNEON) contains functions to analyzing National Ecological Observatory Network (NEON) water chemistry UV-Vis absorbance data and calculate common absorbance metrics such as specific ultra-violet absorbance (SUVA).

<!-- ****** Usage ****** -->

## Usage

To use this package, the chemical properties of surface water data product (DP1.20093.001) should first be downloaded into the local environment using the neonUtilities R package.

For all the functions, the A and B replicate absorbance scans of the same sample are averaged. Replicate samples collected on the same date are not averaged. For **calcSuva254** and **calcSuva280**, discrete wavelength samples will be included if the inputs include data collected prior to the switch to full spectrum in Fall 2023. Because the required wavelengths were not measured, the other functions cannot be run on samples collected prior to the switch to full spectrum in Fall 2023. 

The **formatWide** function reformats the absorbance data from the default long format table into a wide format table, with each column containing absorbance measurements for a respective wavelength.  This function should not be used to generate the inputs for the other functions, as they utilize the default long format tables.

The **calcSuva254** function calculates the specific ultra-violet absorbance (units L/mg/m) at the 254 nm wavelength by dividing the absorbance (units 1/cm) by the dissolved organic carbon (DOC) concentration (units mg/L) and multiplying by 100 (conversion from 1/cm to 1/m). The **calcSuva280** and **calcSuva350** functions use the same method to calculate SUVA at 280 nm and 350 nm respectively.

The **calcE2E3** function calculates the E2:E3 absorbance ratio by dividing the absorbance at 250 nm to the absorbance at 365 nm.  Because NEON absorbance is measured at 2 nm intervals, the average of the absorbance at 364 nm and 366 nm is used to estimate absorbance at 365 nm.

The **calcSR** function calculates the spectral slope ratio. Absorbance values are first transformed using the natural log, and linear models are fit to the regions from 275-295 nm and 350-400 nm. The slope from 275-295 nm is divided by the slope from 350-400 nm to calculate the spectral slope.

For all the functions, the correctFe input can be used to correct for the overlapping absorbance of Fe(III). The Fe(III) absorbance coefficients for a given wavelength are multiplied by the Fe concentration, and then subtracted from the measured absorbance values. This correction assumes that all Fe is in the 3+ oxidation state (ferric). Because Fe(II) does not significantly absorb in the UV, this option will overestimate the Fe(III) absorbance (and underestimate DOC absorbance) if some Fe is in the 2+ oxidation state (ferrous).

<!-- ****** Acknowledgements ****** -->

## References
Doane, T.A. and Horwáth, W.R., 2010. Eliminating interference from iron (III) for ultraviolet absorbance measurements of dissolved organic matter. Chemosphere, 78(11), pp.1409-1415. doi: 10.1016/j.chemosphere.2009.12.062

Hansen, A.M., Kraus, T.E., Pellerin, B.A., Fleck, J.A., Downing, B.D. and Bergamaschi, B.A., 2016. Optical properties of dissolved organic matter (DOM): Effects of biological and photolytic degradation. Limnology and oceanography, 61(3), pp.1015-1032. doi: 10.1002/lno.10270

Korak, J.A. and McKay, G., 2024. Critical review of fluorescence and absorbance measurements as surrogates for the molecular weight and aromaticity of dissolved organic matter. Environmental Science: Processes & Impacts, 26(10), pp.1663-1702. doi: 10.1039/D4EM00183D

Poulin, B.A., Ryan, J.N. and Aiken, G.R., 2014. Effects of iron on optical properties of dissolved organic matter. Environmental science & technology, 48(17), pp.10098-10106. doi: 10.1021/es502670r


<!-- ****** Acknowledgements ****** -->

## Credits & Acknowledgements

<!-- HTML tags to produce image, resize, add hyperlink. -->

<!-- ONLY WORKS WITH HTML or GITHUB documents -->

<a href="http://www.neonscience.org/">
<img src="logo.png" width="300px" /> </a>

<!-- Acknowledgements text -->

The National Ecological Observatory Network is a project solely funded by the National Science Foundation and managed under cooperative agreement by Battelle. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.

<!-- ****** License ****** -->

## License

Creative Commons Zero v1.0 Universal

<!-- ****** Disclaimer ****** -->

## Disclaimer

*Information and documents contained within this package are available as-is. Codes or documents, or their use, may not be supported or
maintained under any program or service and may not be compatible with data currently available from the NEON Data Portal.*
