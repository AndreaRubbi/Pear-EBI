if(Sys.info()["sysname"] == "Darwin" || Sys.info()["sysname"] == "Linux") {
  dyn.load("rTripletDist.so")
  dyn.load("rQuartetDist.so")
} else {
  dyn.load("rTripletDist.dll")
  dyn.load("rQuartetDist.dll")
}

tripletDistance <- function(filename1, filename2) {
  .Call(cTripletDistance, filename1, filename2);
}

pairsTripletDistance <- function(filename1, filename2) {
  .Call(cPairsTripletDistance, filename1, filename2);
}

allPairsTripletDistance <- function(filename) {
  .Call(cAllPairsTripletDistance, filename);
}

quartetDistance <- function(filename1, filename2) {
  .Call(cQuartetDistance, filename1, filename2);
}

pairsQuartetDistance <- function(filename1, filename2) {
  .Call(cPairsQuartetDistance, filename1, filename2);
}

allPairsQuartetDistance <- function(filename) {
  .Call(cAllPairsQuartetDistance, filename);
}


