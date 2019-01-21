#from Matthias Kohl <Matthias.Kohl@stamats.de>, posted to BioC 7/9/08
heatmapCol <- function (data, col, lim) {
   nrcol <- length(col)
   data.range <- range(data)
   if (diff(data.range) == 0)
       stop("data has range 0")
   if (lim <= 0)
       stop("lim has to be positive")
   if (lim > min(abs(data.range))) {
       warning("specified bound 'lim' is out of data range\n\n
                        hence 'min(abs(range(data)))' is used")
       lim <- min(abs(data.range))
   }
   nrcol <- length(col)
   reps1 <- ceiling(nrcol * (-lim - data.range[1])/(2 * lim))
   reps2 <- ceiling(nrcol * (data.range[2] - lim)/(2 * lim))
   col1 <- c(rep(col[1], reps1), col, rep(col[nrcol], reps2))
   return(col1)
}
