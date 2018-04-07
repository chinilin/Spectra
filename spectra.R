# title         : spectra.R
# purpose       : preprocessing & transform raw spectra, fiting models, predict target variable
# producer      : A. Chinilin
# address       : Moscow. RSAU-MTAA

# import data
data <- read.table(file.choose(), header=TRUE, sep="\t") # data.frame with raw spectra
RAW.spectra <- as.data.frame(t(data))
colnames(RAW.spectra) <- RAW.spectra[1, ]
RAW.spectra <- RAW.spectra[-1, ]

library(prospectr) # https://cran.r-project.org/web/packages/prospectr/prospectr.pdf

# simple moving (or running) average filter
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

# Savitzky–Golay smoothing. Savitzky-Golay filtering is a very common
# preprocessing technique. It fits a local polynomial regression
# on the signal and requires equidistant bandwidth. 
SG.spectra <- as.data.frame(savitzkyGolay(RAW.spectra,
                                          m = 0, p = 2, w = 11)) # w = 11  (must be odd)  window size of 11 bands
                                                                 # m = 0 - just filtering
                                                                 # p = 2 - second polynomial order

# compare raw spectra and spectra with Savitzky-Golay filtering
SG.spectra <- as.data.frame(t(SG.spectra))
plot(data$Wavelength[6:501], data$b1A[6:501], type = "l", xlim = c(400,900), ylim = c(0,0.3),
     xlab = "Wavelength/nm", ylab = "Reflectance")
lines(data$Wavelength[6:501], SG.spectra$b1A, type = "l", col = "red")
legend("topleft", col = c("black", "red"), lty = 1, legend = c("Raw", "SG filtering"),
       text.col = c("black", "red"))
SG.spectra <- as.data.frame(t(SG.spectra))

# 1st and 2nd derivatives with "savitzkyGolay" function. Taking (numerical)
# derivatives of the spectra can remove both additive and multiplicative effects in the
# spectra
FD.spectra <- as.data.frame(savitzkyGolay(RAW.spectra, m = 1, p = 2, w = 11)) # m = 1 - first derivative
SD.spectra <- as.data.frame(savitzkyGolay(RAW.spectra, m = 2, p = 2, w = 11)) # m = 2 - second derivative

FD.spectra <- as.data.frame(t(FD.spectra))
SD.spectra <- as.data.frame(t(SD.spectra))
plot(data$Wavelength[6:501], FD.spectra$b1A, type = "l", xlim = c(400,900),
     xlab = "Wavelength/nm", ylab = "")
lines(data$Wavelength[6:501], SD.spectra$b1A, type = "l", col = "red")
legend("bottomleft", legend = c("1st der", "2nd der"), lty = c(1, 1), col = 1:2, horiz = T)
FD.spectra <- as.data.frame(t(FD.spectra))
SD.spectra <- as.data.frame(t(SD.spectra))

# SG filtering + compute continuum–removed values. The continuum removal technique was introduced by as an effective method to highlight absorption
# features of minerals. It can be viewed as an albedo normalization technique. This technique is based on
# the computation of the continuum (or envelope) of a given spectrum. 
CR.spectra <- as.data.frame(continuumRemoval(SG.spectra, type = "R"))
colnames(CR.spectra) <- colnames(SG.spectra)

# SG filtering + Standart Normal Variate (SNV) transformation. Standard
# Normal Variate (SNV) is another simple way for normalizing spectra that intends to correct
# for light scatter.It is better to perform SNV transformation after filtering (by e.g. Savitzky–
# Golay) than the reverse.
SNV.spectra <- as.data.frame(standardNormalVariate(SG.spectra))

# SG filtering + detrend normalization. The SNV–Detrend further accounts for wavelength-dependent scattering effects (variation in curvilinearity
# between the spectra). After a SNV transformation, a 2nd–order polynomial is fit to the spectrum
# and subtracted from it.
SNVD.spectra <- as.data.frame(detrend(SG.spectra, wav = as.numeric(colnames(SG.spectra))))

# SG filtering + Multiplicative Scatter Correction (MSC)
library(pls)
MSC.spectra <- as.data.frame(msc(t(SG.spectra)))
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
        xlab = "Wavelength/nm", ylab = "", main = "Standart Normal Variate (SNV)")
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
dev.off()
#-------------------------------------------------------------------------------------------
# import data
load("~/Google Drive/Ph.D. Thesis/Spectra/13_apr_2016/tr&preproc_spectra.RData")

library(caret)
library(doParallel)

# create new variable (on my example it`s organic carbon content or kaolinite & smektite content)
colnames(RAW.spectra) <- paste(colnames(RAW.spectra), "nm")
RAW.spectra$C <- c(3.21,3.71,3.55,2.67,2.09,3.16,2.87,3.40,0.74,2.07,2.96,2.93,0.98,1.81,0.86,3.47,3.35,2.67,1.81,2.45,2.03,1.87)
RAW.spectra <- RAW.spectra[-c(5,7,8,11,13,20), ]
RAW.spectra$Kaol <- c(7.88,3.07,3.38,2.92,6.50,24.95,26.46,7.98,2.24,3.49,0.21,5.82,8.41,8.18,9.79,6.80)
RAW.spectra$Sm <- c(58.85,59.63,65.21,50.03,58.76,42.93,43.50,54.32,69.07,48.24,59.60,57.61,55.36,51.36,50.91,52.60)

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
plot(varImp(object = mod1), main = "PLSR - Variable Importance",
     top = 15, ylab = "Variable")
#-------------------------------------------------------------------------------------------
# Principal Component Analysis
pcr.grid <- expand.grid(ncomp = 1:10)
set.seed(1234)
mod2 <- train(C~., data = RAW.spectra, # change data
              method = "pcr",
              metric = "RMSE",
              trControl = ctrl1,
              tuneGrid = pcr.grid,
              preProcess = c("center", "scale"))
plot(varImp(object = mod2), main = "PCR - Variable Importance",
     top = 15, ylab = "Variable")
#-------------------------------------------------------------------------------------------
# PCA + MLR with Stepwise Selection
set.seed(1234)
mod3 <- train(C~., data = RAW.spectra,
              method = "lmStepAIC",
              trControl = ctrl1,
              preProcess = c("center", "scale", "pca"),
              trace = F)
# regression coeffitients
coefs <- coef(mod3$finalModel)
plot(varImp(object = mod3),
     main = "PCA + MLR - Variable Importance",
     top = 15, ylab = "Variable")
#-------------------------------------------------------------------------------------------
# Ridge or lasso regression
# note, that if "alpha" is set to 0 this process runs a ridge model,
# if it’s set to 1 it runs a LASSO model and an "alpha" between 0 and 1
# results in an elastic net model
set.seed(1234)
mod4 <- train(C~., data = RAW.spectra, # change data
              method = "glmnet",
              metric = "RMSE",
              trControl = ctrl1,
              preProcess = c("center", "scale"))
plot(varImp(object = mod4), main = "Lasso/Ridge - Variable Importance",
     top = 15, ylab = "Variable")
png(".png", width = 1920, height = 1080, units = 'px', res = 300)
#-------------------------------------------------------------------------------------------
# RF
rftg <- data.frame(mtry = seq(2, 55, by = 2)) # take a lot of time to compute
# can change parametres
# or
mtry <- as.integer(sqrt(ncol(RAW.spectra[, 1:496])))
rf.tuneGrid <- expand.grid(.mtry = mtry)
set.seed(1234)
mod5 <- train(C~., data = RAW.spectra,
              method = "rf",
              tuneGrid = rf.tuneGrid, # or rftg
              trControl = ctrl1,
              importance = TRUE)
plot(varImp(object = mod5), main = "Randon Forest - Variable Importance",
     top = 15, ylab = "Variable")
#-------------------------------------------------------------------------------------------
# XGBoost
gb.tuneGrid <- expand.grid(eta = c(0.3,0.4,0.5,0.6),
                           nrounds = c(50,100,150),
                           max_depth = 2:3, gamma = 0,
                           colsample_bytree = 0.8, min_child_weight = 1,
                           subsample = 1)
set.seed(1234)
mod6 <- train(C~., data = RAW.spectra,
              method = "xgbTree",
              tuneGrid = gb.tuneGrid,
              trControl = ctrl1)
plot(varImp(object = mod6), main = "XGBoost - Variable Importance",
     top = 15, ylab = "Variable")
#-------------------------------------------------------------------------------------------
# SVM
svmRadialTuneGrid <- expand.grid(sigma = c(0.05,0.0456,0.0577),
                                 C = c(1.5,1.596,1.65,1.89,1.95,2,2.2,2.44))
set.seed(1234)
mod7 <- train(C~., data = RAW.spectra,
              method = "svmRadial",
              tuneGrid = svmRadialTuneGrid,
              preProcess = c("center", "scale"),
              trControl = ctrl1)
plot(varImp(object = mod7), main = "SVM - Variable Importance",
     top = 15, ylab = "Variable")
#-------------------------------------------------------------------------------------------
# compile models and compare perfomance
# if we use "ctrl1" or "ctrl2" in "trControl" parametres
model_list <- list(PLSR = mod1, PCR = mod2, GLMNET = mod4,
                   RF = mod5, XGBoost = mod6)
results <- resamples(model_list)
summary(results)
# boxplot comparing results
bwplot(results, layout = c(3, 1)) # RMSE, MSE and R-squared
bwplot(results, metric = "Rsquared", main = "Algorithms accuracy comparing")
bwplot(results, metric = "RMSE", main = "Algorithms accuracy comparing",
       xlim = c(0,2))
stopCluster(cluster)
registerDoSEQ()
#-------------------------------------------------------------------------------------------
require(gridExtra)
grid.arrange(plot(varImp(object = mod2), main = "PCR - Variable Importance (FD spectra)",
                  top = 15, ylab = "Variable"),
             plot(varImp(object = mod5), main = "Randon Forest - Variable Importance (FD spectra)",
                  top = 15, ylab = "Variable"),
             ncol = 2, nrow = 1)
png("SOC Importance PCR&RF.png", width = 3200, height = 1800, units = 'px', res = 300)
dev.off()