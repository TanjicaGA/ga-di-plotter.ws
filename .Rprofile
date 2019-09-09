#### -- Packrat Autoloader (version 0.5.0) -- ####
source("packrat/init.R")
#### -- End Packrat Autoloader -- ####

if( requireNamespace( "drat", quietly=TRUE ) )
    drat::addRepo( "ga-local", "http://gamap/r-packages" )

update_ga_packages <- function() {
    update.packages( ask=FALSE, repos="http://gamap/r-packages" )
}
