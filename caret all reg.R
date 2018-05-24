# title         : caret all regression.R
# purpose       : Run a list of caret regression models in parallel and compare R2 and RMSE
# producer      : A. Chinilin
# address       : Moscow. RSAU-MTAA

models <- c("glmnet", "knn",
            "pcaNNet", "pcr", "pls", "ranger", "qrf",
            "cubist", "svmRadial", "xgbTree")

# load all libraries
library(doParallel)
library(caret)
library(dplyr)
library(DT)
cl <- makeCluster(detectCores())
registerDoParallel(cl)

# compile cross-validation settings
ctrl1 <- trainControl(method = "repeatedcv", number = 5,
                      repeats = 10, allowParallel = TRUE)

# use lapply/loop to run everything
l <- lapply(models, function(i) 
{cat("----------------------------------------------------","\n");
  set.seed(1234); cat(i," <- done\n");
  t2 <- train(C~., data = RAW.spectra, (i), trControl = ctrl1,
              preProcess = c("center", "scale"),
              metric = "RMSE")
}
)

# use lapply to print the results
results <- lapply(1:length(l), function(i) 
{cat(sprintf("%-20s",(models[i])));
  cat(round(l[[i]]$results$Rsquared[which.min(l[[i]]$results$RMSE)],4),"\t");
  cat(round(l[[i]]$results$RMSE[which.min(l[[i]]$results$RMSE)],4),"\t")
  cat(l[[i]]$times$everything[3],"\n")
}
)

# stop the parallel processing and register sequential front-end
stopCluster(cl)
registerDoSEQ()

# preallocate data types
i = 1; MAX = length(l);
x1 <- character() # Name
x2 <- numeric()   # R2
x3 <- numeric()   # RMSE
x4 <- numeric()   # time [s]
x5 <- character() # long model name

# fill data and check indexes and NA
for (i in 1:length(l)) {
  x1[i] <- l[[i]]$method
  x2[i] <- as.numeric(l[[i]]$results$Rsquared[which.min(l[[i]]$results$RMSE)])
  x3[i] <- as.numeric(l[[i]]$results$RMSE[which.min(l[[i]]$results$RMSE)])
  x4[i] <- as.numeric(l[[i]]$times$everything[3])
  x5[i] <- l[[i]]$modelInfo$label
}

# coerce to data frame
df <- data.frame(x1,x2,x3,x4,x5, stringsAsFactors = FALSE)

# print all results to R-GUI
df

# call web browser output with sortable column names
datatable(df,  options = list(
  columnDefs = list(list(className = 'dt-left', targets = c(0,1,2,3,4,5))),
  pageLength = MAX,
  order = list(list(2, 'desc'))),
  colnames = c('Num', 'Name', 'R2', 'RMSE', 'time [s]', 'Model name'),
  caption = paste('Regression results from caret models'),
  class = 'cell-border stripe')  %>%
  formatRound('x2', 3) %>%
  formatRound('x3', 3) %>%
  formatRound('x4', 3) %>%
  formatStyle(2,
              background = styleColorBar(x2, 'steelblue'),
              backgroundSize = '100% 90%',
              backgroundRepeat = 'no-repeat',
              backgroundPosition = 'center'
  )
# compile models and compare perfomance
model_list <- list(GLMNET = l[[1]], KNN = l[[2]], pcaNNet = l[[3]],
                   PCR = l[[4]], PLSR = l[[5]], RF = l[[6]],
                   QRF = l[[7]], Cubist = l[[8]], SVM = l[[9]], XGB = l[[10]])
results <- resamples(model_list)
summary(results)
# boxplot comparing results
bwplot(results, layout = c(3, 1)) # RMSE, MSE and R-squared
#-----------------------------------------------------------------------------#
# Run a list cross-validation methods with PCR method
cl <- makeCluster(detectCores())
registerDoParallel(cl)

# define all cross-validation methods
cvMethods <- c("boot", # bootstrap
               "boot632", # 0.632 bootstrap
               "LGOCV", # leave-one-group cross validation, variant of LOOCV for hierarchical data
               "LOOCV", # leave-one-out cross validation, also known as jacknife
               "cv", # cross validation
               "repeatedcv" # repeated n-fold cross validation
               )

# use R lapply function to loop through all CV methos with qrf
all <- lapply(cvMethods, function(x)
  {set.seed(1234); print(x); tc <- trainControl(method=(x))
fit1 <- train(C~., data = RAW.spectra,
              preProcess = c("center", "scale"),
              trControl = tc,
              method = "pcr") # may choose any of possible regression models
}
)

# stop cluster
stopCluster(cl)
registerDoSEQ()

# extract the used cvMethods 
myNames <- lapply(1:6, function(x) all[[x]]$control$method)

# save results
results <- sapply(all, getTrainPerf)

# change column Names to cv methods
colnames(results) <- myNames

# get the results
results