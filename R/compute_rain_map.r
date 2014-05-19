## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description:R script to compute a rain map form teh CAIRS output
##
## Andreas Scheidegger --  andreas.scheidegger@eawag.ch
## =======================================================


args <- commandArgs(trailingOnly = TRUE)

if(is.na(args[1])) stop("File name is missing!")
file.rain <- args[1]
file.sensor.pos <- args[2]
output.name <- args[3]


library(lattice)
library(latticeExtra)
library(tripack)


## -----------
## read data

data <- read.table(file.rain, sep=",")
calib.coor <- read.table(file.sensor.pos, sep=",")
colnames(data) = c("x", "y", "time", "medianR", "madR", "q.10", "q.90")
colnames(calib.coor) = c("x", "y", "time", "sensor")


## -- define colors

## Meteo Schweiz Colors: http://www.meteoschweiz.admin.ch/web/de/services/datenportal/gitterdaten/alpineprecip.html
col.for.mean <- colorRampPalette(c("#ffffea", "#ecfbc2", "#cefecc", "#9befb4",
                                   "#54bd9f", "#31a895", "#3695b5", "#076eb0",
                                   "#044f90", "#0b1e96", "#2b0246"))(100) #  , "#6a2c5b"

## --
col.for.sd  <- colorRampPalette(c("#a6d96a", "#ffffbf", "#d7191c"))(100) #


## -- reference time
REF.TIME <- strptime("1984-10-20 00:00:00", "%Y-%m-%d %H:%M:%S", tz="GMT")


## -----------
## plots

pdf(output.name, 10, 5)


rr1 <- c(0, max(data$medianR))
rr2 <- range(data$madR)

for(time in sort(unique(data$time))) {

  data.temp <- data[time==data$time,]
  calib.temp <- calib.coor

  data.temp$medianR[data.temp$medianR<0] <- 0

  date <- format(REF.TIME + time/1000 - 13, format="%Y-%m-%d %H:%M:%S")       #leap seconds...

  ## draw Thiessen polygons
  plot.mean <- tileplot(medianR ~ x + y, data=data.temp,
                        main=paste0(file.rain, ", ", date),
                        col.regions=col.for.mean,
                        cut=25,
                        at=seq(0, 1, length=30)^1.7*rr1[2],
                        aspect="iso", points=FALSE,
                        use.tripack = TRUE,
                        xlab="x-coor",
                        ylab="y-coor"
                        )

  ## add Gauges
  plot.mean <- plot.mean + layer(panel.points(calib.temp$x[calib.temp$sensor=="coor"],
                                              calib.temp$y[calib.temp$sensor=="coor"],
                                              pch=4, col='black'))

  ## add Domains

  ## this ugly hack is required because of the strange behavior of layer() in loops...
  cmd.string <- ""
  for(domain in unique(calib.temp$sensor[calib.temp$sensor!="coor"])) {
    cmd.string <- paste0(cmd.string, "
    dom <- calib.temp[calib.temp$sensor=='", domain, "',]
    dom <- dom[chull(dom),]
    plot.mean <- plot.mean + as.layer(
      xyplot(y~x, data=dom,
             panel = function(x, y) { panel.polygon(x, y, col=rgb(0,0,0,0.0), lty=2, border=gray(0.3)) }
             )
      )"
           )

  }
  eval(parse(text=cmd.string))

  plot.mad <- tileplot(madR ~ x + y, data=data.temp,
                      main=paste0(file.rain, ", ", date),
                      col.regions=col.for.sd,
                      cut=25, at=seq(rr2[1], rr2[2], length=30),
                      aspect="iso", pch=16, cex=0.25, col=gray(0.9),
                      use.tripack = TRUE,
                      xlab="x-coor",
                      ylab="y-coor"
                      )

  print(plot.mean, split=c(1,1,2,1), more=TRUE)
  print(plot.mad, split=c(2,1,2,1))

}

dev.off()
