########################################

# R code to analyse presence of spirlin in Poland

########################################

#Import data
df <- read.table(file = "data.txt", 
                 header = TRUE, 
                 dec = ".", 
                 stringsAsFactors=T)

########################################
#Load libraries
library(tidyverse)
library(lattice)
library(sp)
library(gstat)
library(ggplot2)
library(fields)
library(dismo)
library(rgeos)
library(dplyr)
library(tidyr)
library(performance)
library(see)
library(devtools)
library(outliers)
library(car)
library(effects)
library(lme4)
library(glmmTMB)
library(INLA)

#########################################
#Data exploration

#Any NAs
colSums(is.na(df))


# A multi-panel dotchart to see all variables simultaneously.
Names <- c("Year", "Dist", "Alt", "Width",
           "BasinArea", "AvT", "AvJanT", "AvJulT", "FiReach", 
           "FiArea", "Slope","Flow", "Cond", "FishSp", "pred")

# Then plot
dotplot(as.matrix(as.matrix(df[,Names])),
        groups=FALSE,
        strip = strip.custom(bg = 'white',
                             par.strip.text = list(cex = 1.2)),
        scales = list(x = list(relation = "free", draw = TRUE),
                      y = list(relation = "free", draw = FALSE)),
        col = 2, cex  = 0.5, pch = 16,
        xlab = list(label = "Value of the variable", cex = 1.2),
        ylab = list(label = "Order of the data", 
                    cex = 1.2))

# Remove outliers in BasinArea
df1 <- df[which(df$BasinArea < 150000),]

# Remove outliers in size of reach
df2 <- df1[which(df1$FiReach < 3000),]

#Drop 3 catchments with few data
df3 <- df2 [-which(df2$RBasin == "Laba"), ]
df4 <- df3 [-which(df3$RBasin == "Dunaj"), ]
df5 <- df4 [-which(df4$RBasin == "Niemen"), ]

df6 <- droplevels(df5)
table(df6$RBasin)

#Select variables of interest
df7 <-subset(df6, select=c("Lat", "Lon", "ROrder", "Geol", "RBasin", "Ecoreg", "Mon",
                          "Geomor", "FloodVal", "Source", "Tree", "Substr", "Ripar",
                          "Year", "Dist", "Alt", "Width","BasinArea", "AvT", "AvJanT", 
                          "AvJulT", "FiReach", "FiArea", "Slope","FishSp", "spirl",
                          "Flow", "FlowDis", "Cond", "Bush", "Dam", "Lake", "pred"))


#Any NAs
colSums(is.na(df7))

#Remove NAs
df8 <- na.omit(df7)

#Data lost
nrow(df) - nrow(df8)
#690 rows

#Assign factors
df8$fOrder  <- as.factor(df8$ROrder)
df8$fGeol   <- as.factor(df8$Geol)
df8$fEco    <- as.factor(df8$Ecoreg)
df8$fBasin  <- as.factor(df8$RBasin)
df8$fMon    <- as.factor(df8$Mon)
df8$fYear   <- as.factor(df8$Year)
df8$fGeom   <- as.factor(df8$Geomor)
df8$fFVal   <- as.factor(df8$FloodVal)
df8$fSource <- as.factor(df8$Source)
df8$fTree   <- as.factor(df8$Tree)
df8$fSub    <- as.factor(df8$Substr)
df8$fRip    <- as.factor(df8$Ripar)
df8$fFDis   <- as.factor(df8$FlowDis)
df8$fBush   <- as.factor(df8$Bush)
df8$fDam    <- as.factor(df8$Dam)
df8$fLake   <- as.factor(df8$Lake)

# Standardize
df8$FR.std      <- (df8$FiReach-mean(df8$FiReach))/sd(df8$FiReach)
df8$Spec.std    <- (df8$FishSp-mean(df8$FishSp))/sd(df8$FishSp)
df8$alt.std     <- (df8$Alt-mean(df8$Alt))/sd(df8$Alt)
df8$slope.std   <- (df8$Slope-mean(df8$Slope))/sd(df8$Slope)
df8$basin.std   <- (df8$BasinArea-mean(df8$BasinArea))/sd(df8$BasinArea)
df8$wid.std     <- (df8$Width-mean(df8$Width))/sd(df8$Width)
df8$flow.std    <- (df8$Flow-mean(df8$Flow))/sd(df8$Flow)
df8$AvJanT.std  <- (df8$AvJanT-mean(df8$AvJanT))/sd(df8$AvJanT)
df8$dist.std    <- (df8$Dist-mean(df8$Dist))/sd(df8$Dist)
df8$cond.std    <- (df8$Cond-mean(df8$Cond))/sd(df8$Cond)
df8$Pred.std    <- (df8$pred-mean(df8$pred))/sd(df8$pred)

# Model formulation - Bernoulli GLM with all variables
B1 <- glm(spirl ~ fOrder + fGeol + fEco + fGeom + fFVal + fSource + fBasin +
                  fTree + fSub + fRip + fFDis  + fBush + fDam + fLake + 
                  Spec.std + alt.std + wid.std + basin.std +  Pred.std + 
                  AvJanT.std + slope.std + flow.std + dist.std + cond.std +
                  offset(FR.std),
                  data = df8,
                  family = binomial(link = "logit"))

vif(B1)

# Plot residuals versus fitted values
par(mfrow = c(1,1), mar = c(5,5,1,1), cex.lab = 1.2)

Fitted <- fitted(B1)
Resid  <- resid(B1, type = "pearson")
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.2)
plot(x = Fitted, y = Resid,
     xlab = "Fitted values", 
     ylab = "Pearson Residuals")
abline(h = 0, lty = 2)

###########################################

#Evidence for spatial dependency?

# Make a variogram of the Pearson residuals
mydata <- data.frame(Resid, df8$Lat, df8$Lon)
coordinates(mydata)    <- c("df8.Lat", "df8.Lon")
max(coordinates(mydata)[,1]) - min(coordinates(mydata)[,1])

Vario <- variogram(Resid ~ 1, 
                   cutoff = 5,
                   width = 0.5,
                   mydata,
                   cressie = TRUE)
plot(Vario, 
     main = "", 
     xlab = list(label = "Lat/Lon (degrees)", cex = 1.5), 
     ylab = list(label = "Semi-variogram", cex = 1.5), 
     pch = 16, fill = 3, cex = 1.5)

# Evidence of spatial dependency in the data

###########################################

# Spatial model

#Define a mesh for the sampling locations
Loc <- cbind(df8$Lon, df8$Lat)
head(Loc)
dim(Loc)

# Distances between sampling locations
D <- dist(Loc)
length(D)
par(mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between sites (lat/lon)",
     ylab = "Frequency")

plot(x = sort(D), 
     y = (1:length(D))/length(D), 
     type = "l",
     xlab = "Distance between sites",
     ylab = "Cumulative proportion")

#Make a mesh
Bound <- inla.nonconvex.hull(Loc)
mesh  <- inla.mesh.2d(boundary = Bound, 
                      max.edge  = c(0.1,0.8), 
                      cutoff    = 0.2)
par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(mesh, asp=1, main = "")
points(Loc, 
       col = "black",  
       pch = 16, 
       cex = 0.7)

mesh$n

# Sampling locations matching points on the mesh
A2 <- inla.spde.make.A(mesh, loc = Loc)
dim(A2)  #1144 observations on a 1868 grid


# Define the Matern correlation on the mesh using spde 
# (stochastic partial differential equation)

# spde   <- inla.spde2.matern(mesh, alpha = 2)

# PC priors
spde  <- inla.spde2.pcmatern(mesh, alpha = 2,  
                             prior.range = c(1, 0.001), #low prob range
                             prior.sigma = c(2,  0.05))

w.index <- inla.spde.make.index(
           name   = 'w', 
           n.spde = spde$n.spde,
          n.group = 1,
          n.repl  = 1)

# Create a data frame
N <- nrow(df8)
X <- data.frame(Intercept = rep(1, N), 
                basin.std  = df8$basin.std,
                Spec.std = df8$Spec.std,
                FR.std  = df8$FR.std,
                alt.std  = df8$alt.std,
                slope.std  = df8$slope.std,
                basin.std  = df8$basin.std,
                wid.std  = df8$wid.std,
                flow.std  = df8$flow.std,
                AvJanT.std  = df8$AvJanT.std,
                dist.std  = df8$dist.std,
                cond.std  = df8$cond.std,
                Pred.std  = df8$Pred.std,
                fOrder = df8$fOrder,
                fEco     = df8$fEco,  
                fMon     = df8$fMon, 
                fYear   = df8$fYear,   
                fGeol    = df8$fGeol,
                fBasin  = df8$fBasin,
                fGeom     = df8$fGeom,
                fFVal     = df8$fFVal,
                fSource   = df8$fSource,
                fTree     = df8$fTree,
                fSub     = df8$fSub,
                fRip     = df8$fRip, 
                fDis    = df8$fFDis,
                fBush    = df8$fBush,
                fDam     = df8$fDam,
                fLake  = df8$fLake
                )
head(X)

# Mesh points sampled
stk2 <- inla.stack(
        tag  = "Fit",
        data = list(y = df8$spirl),  
        A    = list(A2, 1),                      
        effects = list(                 
                w = w.index,          
                X = as.data.frame(X)))

# Specify a priori models
f01 <- y ~ offset(FR.std) + f(w, model = spde)
f02 <- y ~ offset(FR.std) + flow.std + fSub + wid.std + f(w, model = spde)
f03 <- y ~ offset(FR.std) + slope.std + dist.std + alt.std + f(w, model = spde)
f04 <- y ~ offset(FR.std) + cond.std + basin.std + f(w, model = spde)
f05 <- y ~ offset(FR.std) + AvJanT.std + basin.std + f(w, model = spde)
f06 <- y ~ offset(FR.std) + fBasin + fDis + Spec.std + f(w, model = spde)
f07 <- y ~ offset(FR.std) + fDam + fLake + fDis + fRip + f(w, model = spde)
f08 <- y ~ offset(FR.std) + fGeol + fGeom + fFVal + f(w, model = spde)
f09 <- y ~ offset(FR.std) + fTree + alt.std + fRip + f(w, model = spde)
f10 <- y ~ offset(FR.std) + basin.std + fOrder + Spec.std + f(w, model = spde)
f11 <- y ~ offset(FR.std) + alt.std + basin.std + f(w, model = spde)
f12 <- y ~ offset(FR.std) + cond.std + fRip + f(w, model = spde)
f13 <- y ~ offset(FR.std) + fGeol + dist.std + f(w, model = spde)
f14 <- y ~ offset(FR.std) + fDis + flow.std + fDam + Pred.std + f(w, model = spde)
f15 <- y ~ offset(FR.std) + Spec.std + Pred.std + f(w, model = spde)

#Execute
I1 <- inla(f01,
           family = "binomial",
           data=inla.stack.data(stk2),
           control.compute = list(cpo = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(stk2)))
I2 <- inla(f02,
           family = "binomial",
           data=inla.stack.data(stk2),
           control.compute = list(cpo = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(stk2)))
I3 <- inla(f03,
           family = "binomial",
           data=inla.stack.data(stk2),
           control.compute = list(cpo = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(stk2)))
I4 <- inla(f04,
           family = "binomial",
           data=inla.stack.data(stk2),
           control.compute = list(cpo = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(stk2)))
I5 <- inla(f05,
           family = "binomial",
           data=inla.stack.data(stk2),
           control.compute = list(cpo = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(stk2)))
I6 <- inla(f06,
           family = "binomial",
           data=inla.stack.data(stk2),
           control.compute = list(cpo = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(stk2)))
I7 <- inla(f07,
           family = "binomial",
           data=inla.stack.data(stk2),
           control.compute = list(cpo = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(stk2)))
I8 <- inla(f08,
           family = "binomial",
           data=inla.stack.data(stk2),
           control.compute = list(cpo = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(stk2)))
I9 <- inla(f09,
           family = "binomial",
           data=inla.stack.data(stk2),
           control.compute = list(cpo = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(stk2)))
I10 <- inla(f10,
           family = "binomial",
           data=inla.stack.data(stk2),
           control.compute = list(cpo = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(stk2)))
I11 <- inla(f11,
           family = "binomial",
           data=inla.stack.data(stk2),
           control.compute = list(cpo = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(stk2)))
I12 <- inla(f12,
           family = "binomial",
           data=inla.stack.data(stk2),
           control.compute = list(cpo = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(stk2)))
I13 <- inla(f13,
            family = "binomial",
            data=inla.stack.data(stk2),
            control.compute = list(cpo = TRUE, waic = TRUE),
            control.predictor = list(A = inla.stack.A(stk2)))
I14 <- inla(f14,
            family = "binomial",
            data=inla.stack.data(stk2),
            control.compute = list(cpo = TRUE, waic = TRUE),
            control.predictor = list(A = inla.stack.A(stk2)))
I15 <- inla(f15,
            family = "binomial",
            data=inla.stack.data(stk2),
            control.compute = list(cpo = TRUE, waic = TRUE),
            control.predictor = list(A = inla.stack.A(stk2)))

waic2 <- c(I1$waic$waic, I2$waic$waic, I3$waic$waic,
           I4$waic$waic, I5$waic$waic, I6$waic$waic, I7$waic$waic,
           I8$waic$waic, I9$waic$waic, I10$waic$waic, I11$waic$waic,
           I12$waic$waic, I13$waic$waic,I14$waic$waic,I15$waic$waic)
WAIC     <- cbind(waic2)
rownames(WAIC) <- c("I1", "I2", "I3", "I4", "I5", "I6",
                    "I7", "I8", "I9", "I10", "I11", "I12",
                    "I13", "I14", "I15")
round(WAIC,0)


# Fitted values and Pearson residuals I6
IndexFit <- inla.stack.index(stk2, tag = "Fit")$data
mu3      <- I6$summary.fitted.values[IndexFit, "mean"]
E3       <- (df8$spirl - mu3) / sqrt(mu3)

# Variogram model I6
MyData3 <- data.frame(E3 = E3, Xkm = df8$Lon, Ykm = df8$Lat)
coordinates(MyData3) <- c("Xkm", "Ykm")
Vario3 <- variogram(E3 ~ 1, MyData3, cutoff = 5, width = 0.5)
plot(Vario3,  main = "", xlab = list(label = "Lat/Lon (degrees)", cex = 1.5), 
     ylab = list(label = "Semi-variogram", cex = 1.5), pch = 16, fill = 3, cex = 1.5)

# Model output
glm <- I6$summary.fixed[, c("mean", "0.025quant", "0.975quant")]
colnames(glm) <- c("Mean", "lowerCI", "upperCI")
round(glm, 2)
                        
#Posterior marginals of the standard deviation
post.se <- inla.tmarginal(function(x) sqrt(1 / exp(x)),
                          I6$internal.marginals.hyperpar[[1]])

par(mfrow = c(1, 1), mar = c(5, 5, 3, 1), mgp = c(2, 0.6, 0))
plot(post.se, type = "l", xlab = 'Spatial effect variance', 
     ylab = "Density", main = "")


####################################################
# Plot the spatial random field
PlotField <- function(field, mesh, ContourMap, xlim, ylim, Add=FALSE,...){
        stopifnot(length(field) == mesh$n)
        
        if (missing(xlim)) xlim <- ContourMap@bbox[1, ] 
        if (missing(ylim)) ylim <- ContourMap@bbox[2, ]
        
        proj <- inla.mesh.projector(mesh, 
                                    xlim = xlim, 
                                    ylim = ylim, 
                                    dims = c(300, 300))
       
        field.proj <- inla.mesh.project(proj, field)
        
        # Plot
        image.plot(list(x = proj$x, 
                        y = proj$y,
                        z = field.proj), 
                   xlim = xlim, 
                   ylim = ylim,
                   asp = 1,
                   add = Add,
                   ...)  
}

# Posterior mean of the spatial field at the mesh points:
w.pm <- I6$summary.random$w$mean  
w.sd <- I6$summary.random$w$sd

# Plot random spatial field
PlotField(field = w.pm, 
           mesh = mesh, 
           xlim = range(mesh$loc[,1]),
           ylim = range(mesh$loc[,2]))

# Add the sampling locations (in UTM)
points(x = Loc[,1],
       y = Loc[,2], 
     cex = 0.5, 
     col = "black", 
     pch = 16)

# Calculate the spatial field on a grid
wproj <- inla.mesh.projector(mesh)   

#####################################
# Alternative plot of the spatial random field

w.pm100_100 <- inla.mesh.project(wproj, w.pm)
w.sd100_100 <- inla.mesh.project(wproj, w.sd)

Grid <- expand.grid(Ykm = wproj$y, 
                    Xkm = wproj$x)
w.proj    <- inla.mesh.projector(mesh) 
Grid$w.pm <- as.vector(w.pm100_100)     
Grid$w.sd <- as.vector(w.sd100_100)               

col.m <- colorRampPalette(c('darkgreen', "white"))(100)
plot.wpm <- levelplot(w.pm ~ Ykm * Xkm,
                      data = Grid, 
                      aspect = "fill",
                      col.regions = col.m,
                      scales = list(draw = TRUE),
                      xlab = list("Longitude", cex = 1.2),
                      ylab = list("Latitude", cex = 1.2),
                      main = list("", cex = 1.2),
                      panel=function(...){
                              panel.levelplot(...)
                              grid.points(x = df8$Lat, 
                                          y = df8$Lon, 
                                          pch = 16,
                                          size = unit(0.4, "char"))
                              })
plot.wpm

# And do the same for the posterior standard deviation 
col.sd <- colorRampPalette(c('red3', 'white'))(100) 
plot.wsd <- levelplot(w.sd ~ Ykm * Xkm,
                      data = Grid, 
                      aspect = "fill",
                      col.regions = col.sd,
                      scales = list(draw = TRUE),
                      xlab = list("Longitude", cex = 1.2),
                      ylab = list("Latitude", cex = 1.2),
                      main = list("", cex = 1.2),
                      panel=function(...){
                              panel.levelplot(...)
                              grid.points(x = df8$Lat, 
                                          y = df8$Lon, 
                                          pch = 16,
                                          size = unit(0.4, "char"))})
plot.wsd


# Correlation function imposed on the residuals:
post.se <- inla.tmarginal(function(x) sqrt(1 / exp(x)),
                          I6$internal.marginals.hyperpar[[1]])

SpFi.w <- inla.spde2.result(inla = I6,
                            name = "w",
                            spde = spde,
                            do.transfer = TRUE)

Kappa <- inla.emarginal(function(x) x, 
                        SpFi.w$marginals.kappa[[1]] )

sigmau <- inla.emarginal(function(x) sqrt(x), 
                         SpFi.w$marginals.variance.nominal[[1]] )

r <- inla.emarginal(function(x) x, 
                    SpFi.w$marginals.range.nominal[[1]] )

c(Kappa, sigmau, r)
# Hyperparameters are κspatial = 1.39, σspatial = 3.23; range = 2.30 km.

#Plot correlation structure
D     <- as.matrix(dist(Loc[,1:2]))
d.vec <- seq(0, max(D), length = 100)      
Cor.M <- (Kappa * d.vec) * besselK(Kappa * d.vec, 1) 
Cor.M[1] <- 1

par(mfrow=c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = d.vec, 
     y = Cor.M, 
     pch = 16, 
     type = "l", 
     cex.lab = 1.5,
     xlab = "Distance (km)", 
     ylab = "Correlation",
     xlim = c(0, 4))

#Posterior marginals of the standard deviation
post.se <- inla.tmarginal(function(x) sqrt(1 / exp(x)),
                          I6$internal.marginals.hyperpar[[1]])

# Posterior distributions
par(mfrow = c(2, 3), mar = c(3, 3, 2, 2), mgp = c(2, 0.6, 0))
plot(I6$marginals.fix[[1]], type = 'l', 
     xlab = 'Intercept', ylab = 'Density', cex.lab=1.3, ylim = c(0, 0.8),
     abline (v = 0, lty = 2)) 
plot(I6$marginals.fix[[2]], type = 'l', 
     xlab = 'Oder', ylab = 'Density', cex.lab=1.3, ylim = c(0, 0.35),
     abline (v = 0, lty = 2)) 
plot(I6$marginals.fix[[3]], type = 'l', 
     xlab = 'Vistula', ylab = 'Density', cex.lab=1.3, ylim = c(0, 0.75),
     abline (v = 0, lty = 2)) 
plot(I6$marginals.fix[[4]], type = 'l', 
     xlab = 'Flow disruption', ylab = 'Density', cex.lab=1.3, ylim = c(0, 1.8),
     abline (v = 0, lty = 2)) 
plot(I6$marginals.fix[[5]], type = 'l', 
     xlab = 'Fish species', ylab = 'Density', cex.lab=1.3, ylim = c(0, 3.4),
     abline (v = 0, lty = 2)) 
plot(x = d.vec, y = Cor.M, type = "l", cex.lab = 1.3, 
     xlab = "Distance (km)", ylab = "Correlation", xlim = c(-0.1, 4),
     abline (v = 0, lty = 2))

