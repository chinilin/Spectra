# title         : spectra.R
# purpose       : preprocessing & transform raw spectra, fiting models, predict target variable
# producer      : A. Chinilin
# address       : Moscow

# import data
data <- read.table(file.choose(), header=TRUE, sep="\t") # data.frame with raw spectra
RAW.spectra <- as.data.frame(t(data))
colnames(RAW.spectra) <- RAW.spectra[1, ]
RAW.spectra <- RAW.spectra[-1, ]

library(prospectr) # https://cran.r-project.org/web/packages/prospectr/prospectr.pdf

# Moving average or runnnig mean
MA.spectra <- as.data.frame(movav(RAW.spectra, w = 11)) # window size of 11 bands
# note that the 5 first and last bands are lost in the process

# compare raw spectra and spectra with moving average
MA.spectra <- as.data.frame(t(MA.spectra))
plot(data$Wavelength[6:501], data$b1A[6:501], type = "l", xlim = c(400,900), ylim = c(0,0.3),
     xlab = "Wavelength/nm", ylab = "Reflectance")
lines(data$Wavelength[6:501], MA.spectra$b1A, type = "l", col = "red")
legend("topleft", col = c("black", "red"), lty = 1, legend = c("Raw", "Moving average"),
       text.col = c("black", "red"))
MA.spectra <- as.data.frame(t(MA.spectra))

# Savitzky-Golay filtering
SG.spectra <- as.data.frame(savitzkyGolay(RAW.spectra,
                                          m = 0, p = 2, w = 11)) # window size of 11 bands
                                                                 # m = 0 - just filtering
                                                                 # p = 2 - second polynomail order

# compare raw spectra and spectra with Savitzky-Golay filtering
SG.spectra <- as.data.frame(t(SG.spectra))
plot(data$Wavelength[6:501], data$b1A[6:501], type = "l", xlim = c(400,900), ylim = c(0,0.3),
     xlab = "Wavelength/nm", ylab = "Reflectance")
lines(data$Wavelength[6:501], SG.spectra$b1A, type = "l", col = "red")
legend("topleft", col = c("black", "red"), lty = 1, legend = c("Raw", "SG filtering"),
       text.col = c("black", "red"))
SG.spectra <- as.data.frame(t(SG.spectra))

# 1st and 2nd derivatives with "savitzkyGolay" function
FD.spectra <- as.data.frame(savitzkyGolay(RAW.spectra, m = 1, p = 2, w = 11))
SD.spectra <- as.data.frame(savitzkyGolay(RAW.spectra, m = 2, p = 2, w = 11))

FD.spectra <- as.data.frame(t(FD.spectra))
SD.spectra <- as.data.frame(t(SD.spectra))
plot(data$Wavelength[6:501], FD.spectra$b1A, type = "l", xlim = c(400,900),
     xlab = "Wavelength/nm", ylab = "")
lines(data$Wavelength[6:501], SD.spectra$b1A, type = "l", col = "red")
legend("bottomleft", legend = c("1st der", "2nd der"), lty = c(1, 1), col = 1:2, horiz = T)
FD.spectra <- as.data.frame(t(FD.spectra))
SD.spectra <- as.data.frame(t(SD.spectra))

# SG filtering + Continuum removal
CR.spectra <- as.data.frame(continuumRemoval(SG.spectra, type = "R"))
colnames(CR.spectra) <- colnames(SG.spectra)

# SG filtering + Standard Normal Variate (SNV)
SNV.spectra <- as.data.frame(standardNormalVariate(SG.spectra))

# SG filtering + SNV–Detrend
SNVD.spectra <- as.data.frame(detrend(SG.spectra, wav = as.numeric(colnames(SG.spectra))))

# SG filtering + Multiplicative Scatter Correction (MSC)
library(pls)
MSC.spectra <- as.data.frame(msc(t(RAW.spectra)))
MSC.spectra <- as.data.frame(t(MSC.spectra))

# multiplot
par(mfcol=c(4,2), mar=c(4.1,4,2.90,0.5))

# 1 plot
matplot(as.numeric(colnames(RAW.spectra)), t(RAW.spectra), type = "l",
        xlab = "Wavelength/nm", ylab = "Reflectance", main = "No preprocessing")
# 2 plot
matplot(as.numeric(colnames(MA.spectra)), t(MA.spectra), type = "l",
        xlab = "Wavelength/nm", ylab = "Reflectance", main = "Moving average filter")
# 3 plot
matplot(as.numeric(colnames(SG.spectra)), t(SG.spectra), type = "l",
        xlab = "Wavelength/nm", ylab = "Reflectance", main = "Savitzky-Golay filter")
# 4 plot
matplot(as.numeric(colnames(FD.spectra)), t(FD.spectra), type = "l",
        xlab = "Wavelength/nm", ylab = "", main = "1-st derivative")
# 5 plot
matplot(as.numeric(colnames(SD.spectra)), t(SD.spectra), type = "l",
        xlab = "Wavelength/nm", ylab = "", main = "2-nd derivative")
# 6 plot
matplot(as.numeric(colnames(SNV.spectra)), t(SNV.spectra), type = "l",
        xlab = "Wavelength/nm", ylab = "", main = "Standard Normal Variate (SNV)")
# 7 plot
matplot(as.numeric(colnames(SNVD.spectra)), t(SNVD.spectra), type = "l",
        xlab = "Wavelength/nm", ylab = "", main = "SNV-Detrend")
# 8 plot
matplot(as.numeric(colnames(MSC.spectra)), t(MSC.spectra), type = "l",
        xlab = "Wavelength/nm", ylab = "", main = "Multiplicative Scatter Correction (MSC)")
# 9 plot
matplot(as.numeric(colnames(CR.spectra)), t(CR.spectra), type = "l",
        xlab = "Wavelength/nm", ylab = "", main = "Continuum Removal (CR)")
dev.off()
# save png
# not run
png("Spectra.png", width = 4096, height = 2160, units = 'px', res = 300)

# to L8 band 2 (0,452-0,512), band 3 (0,533-0,590), band 4 (0,636-0,673) & band 5 (0,851-0,879) (Vis-NIR)
RAW.spectra.L8 <- subset(RAW.spectra, select = -c(1:52,114:133,192:236,275:451,481:506))
#-------------------------------------------------------------------------------------------
# import data
load("~/Google Drive/Ph.D. Thesis/Spectra/13_apr_2016/tr&preproc_spectra.RData")

library(pls)
library(caret)
library(doParallel)

# create new variable (on my example it`s organic carbon)
RAW.spectra$C <- c(3.21,3.71,3.55,2.67,2.09,3.16,2.87,3.40,0.74,2.07,2.96,2.93,0.98,1.81,0.86,3.47,3.35,2.67,1.81,2.45,2.03,1.87)
RAW.spectra$Kaol <- c(7.88,3.07, 3.38,2.92,6.50,24.95,26.46,7.98,2.24,3.49,0.21,5.82,8.41,8.18,9.79,6.80)

cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)
set.seed(1234)
# compile cross-validation settings
ctrl <- trainControl(method = "LOOCV", returnResamp = "final")
ctrl1 <- trainControl(method = "repeatedcv", number = 5, repeats = 10, allowParallel = TRUE)
ctrl2 <- trainControl(method = "cv", number = 5)
#-------------------------------------------------------------------------------------------
# train PLSR model
set.seed(1234)
mod1 <- train(C~., data = RAW.spectra, # change data
              method = "pls",
              metric = "RMSE",
              trControl = ctrl1,
              preProcess = c("center", "scale"))
plot(varImp(object = mod1), main = "PLSR - Variable Importance", top = 15)
# check CV profile
plot(mod1)
# or
mod1.1 <- plsr(C~., data = RAW.spectra,
               scale = T,
               validation = "CV",
               segments = 5, jackknife = T)
summary(mod1.1)
plot(RMSEP(mod1.1), legendpos = "topright")
explvar(mod1.1)
# or train PCR model
mod1.2 <- pcr(C~., data = RAW.spectra,
              scale = T,
              validation = "CV",
              segments = 5, jackknife = T)
validationplot(mod1.2)
#-------------------------------------------------------------------------------------------
# PCA-Stepwise LM
set.seed(1234)
mod2 <- train(C~., data = RAW.spectra,
              method = "lmStepAIC",
              trControl = ctrl1,
              preProcess = c("center", "scale", "pca"),
              trace = F)
# regression coeffitients
coefs <- coef(mod2$finalModel)
plot(varImp(object = mod2), main = "Stepwise LM + PCA - Variable Importance", top = 15)
#-------------------------------------------------------------------------------------------
# Ridge or lasso regression
# note, that if "alpha" is set to 0 this process runs a ridge model,
# if it’s set to 1 it runs a LASSO model and an "alpha" between 0 and 1
# results in an elastic net model
set.seed(1234)
mod3 <- train(C~., data = RAW.spectra, # change data
              method = "glmnet",
              metric = "RMSE",
              trControl = ctrl1,
              preProcess = c("center", "scale"))
plot(varImp(object = mod3), main = "Elastic Net - Variable Importance", top = 15)
#-------------------------------------------------------------------------------------------
# RF
rftg <- data.frame(mtry = seq(2, 55, by = 2)) # take a lot of time to compute
# can change parametres
# or
mtry <- as.integer(sqrt(ncol(RAW.spectra[, 1:506])))
rf.tuneGrid <- expand.grid(.mtry = mtry)
set.seed(1234)
mod4 <- train(C~., data = RAW.spectra,
              method = "rf",
              tuneGrid = rf.tuneGrid, # or rftg
              trControl = ctrl1)
#-------------------------------------------------------------------------------------------
# XGBoost
gb.tuneGrid <- expand.grid(eta=c(0.3,0.4),
                           nrounds=c(50,100),
                           max_depth=2:3, gamma=0,
                           colsample_bytree=0.8, min_child_weight=1,
                           subsample = 1)
set.seed(1234)
mod5 <- train(C~., data = RAW.spectra,
              method = "xgbTree",
              tuneGrid = gb.tuneGrid,
              trControl = ctrl1)
#-------------------------------------------------------------------------------------------
# compile models and compare perfomance
# if we use "ctrl1" or "ctrl2" in "trControl" parametres
model_list <- list(PLSR = mod1, PCA_SLM = mod2, GLMnet = mod3, RF = mod4,
                   XGBoost = mod5)
results <- resamples(model_list)
# boxplot comparing results
bwplot(results, metric = "Rsquared")
bwplot(results, metric = "RMSE")
# shut down the cluster 
stopCluster(cluster)
registerDoSEQ()
############################################################################################
# work with Landsat images to predict target variable
library(RStoolbox)
library(raster)
library(maptools)
library(ggplot2)
library(sp)
library(plotKML)

setwd("~/Google Drive/Ph.D. Thesis/Space_images/Landsat 8 OLI_TIRS/25 - Apr - 2014")
# import Landsat metadata from MTL file
metaData <- readMeta("~/Google Drive/Ph.D. Thesis/Space_images/Landsat 8 OLI_TIRS/25 - Apr - 2014/LC08_L1TP_175025_20140425_20170423_01_T1_MTL.txt")
# load Landsat bands based on metadata
lsat <- stackMeta(metaData)
# plotting of remote sensing imagery in RGB with ggplot2
ggRGB(lsat, r = 4, g = 3, b = 2,
      stretch = "lin")
# radiometric conversions and corrections
lsat.ref <- radCor(lsat, metaData = metaData,
                   method = "apref",
                   bandSet = c("B2_dn", "B3_dn", "B4_dn","B5_dn"))
# crop the data
fields <- readShapePoly("Fields_of_interest.shp")
lsat.sub <- crop(lsat.ref, extent(fields))
lsat.sub <- mask(lsat.sub, fields)
ggRGB(lsat.sub, r = 4, g = 3, b = 2, stretch = "lin")+
  ggtitle("Test fields")
# unsupervised ckustering using kmeans clustering
set.seed(21)
unC <- unsuperClass(lsat.sub,
                    nSamples = 250,
                    nClasses = 8,
                    nStarts = 3)
unC
colors <- rainbow(8)
plot(unC$map, col = colors, box = F,
     legend = F, axes = F, main = "K-means classification")
legend(567000,5593750, legend = paste0("C", 1:8), fill = colors,
       title = "Classes", horiz = F, bty = "n")
# calculates R-mode PCA for RasterBricks or RasterStacks and
# returns a RasterBrick with multiple layers of PCA scores
lsat.pca <- rasterPCA(lsat.sub,
                      spca = T,
                      nComp = 2,     # can change
                      maskCheck = T,
                      nSamples = NULL)
summary(lsat.pca$model)
# write PC rasters
PC <- writeRaster(lsat.pca$map,
                  filename = c("PC1.tiff","PC2.tiff"),
                  format = "GTiff",
                  bylayer = T,
                  datatype = "INT2S")

C.map <- predict(lsat.pca$map, mod2$finalModel,
                 progress = "text", na.rm = T)

scale <- list("SpatialPolygonsRescale", layout.scale.bar(), 
              offset = c(565300,5592250), scale = 500, fill = c("transparent","black"))
text1 <- list("sp.text", c(565300,5592310), "0")
text2 <- list("sp.text", c(565800,5592310), "500 m")
arrow <- list("SpatialPolygonsRescale", layout.north.arrow(), 
              offset = c(566750,5593650), scale = 250)
spplot(C.map, col.regions = SAGA_pal[[1]],
       #scales = list(draw = T),
       sp.layout=list(scale, text1, text2, arrow))