reqd <- c("cartogram", "sp", "fMultivar", "RANN","rgdal","rgeos",
          "inflection")
pkgs <- installed.packages()[,1]
need <- reqd[!reqd %in% pkgs]
if(length(need)) install.packages(need)
