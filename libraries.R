reqd <- c("cartogram", "sp", "fMultivar", "RANN","rgdal","rgeos")
pkgs <- rownames(installed.packages())
need <- reqd[!reqd %in% pkgs]
if(length(need)) install.packages(need)
