
# startup massage -----------------------------------------------------------------------------

.onAttach <- function(libname, pkgname){
    packageStartupMessage("Package MrnAnnoAlgo3 version: ", packageVersion("MrnAnnoAlgo3"))
}
