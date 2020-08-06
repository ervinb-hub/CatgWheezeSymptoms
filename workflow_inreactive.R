require(cluster)
require(NbClust)
require(clustertend)
require(caret)
require(mice)
require(dbscan)
require(fpc)
require(colorspace)
require(class)
require(tree)
require(orclus)
require(doParallel)
library(psych)

# Environment settings
rm(list = ls())

# PARAMETERS TO BE SET
MC_Simulations <- 1                     # N. of Monte Carlo Simulations
METHODS <- c("knn", "nnet", "rpart")    # Metodologies in use (only these 3)

cl <- makeCluster((detectCores() * 0.75))
registerDoParallel(cl)


#' Set categorical variables as factor data types
#' @param dframe The dataset we are working on (as data.frame)
#' @return A dataset with categorical variables as factors
factorize <- function(dframe) {
  tmp <- which(sapply(dframe, class) == "integer")
  for (i in tmp) { 
    #print(paste(names(df)[[i]], length(unique(df[,i]))))
    cnt <- length(unique(dframe[,i]))
    if (cnt <= 3) {
      dframe[,i] <- as.factor(dframe[,i])
    }
  }
  return(dframe)
}


# Try imputation with Amelia II
impWithAmelia <- function(df, catgs) {
  ridge <- c(0.01, 0.05, 0.1)
  names(ridge) <- c("1%", "5%", "10%")
  for (i in 1:length(ridge)) {
    print(paste(cat("\n"), "Trying to impute with Amelia II, ridge = ", names(ridge)[[i]]))
    a <- amelia(df, m = 1, noms = catgs, p2s = 0,  empri = ridge[[i]] * nrow(df))
    if (a$code == 1) {
      print("Chain lengths:")
      print(sapply(a$iterHist, FUN = nrow))
      break
    }
  }
  return(a)
}


# Delete predictors with more than X% missing values
remHighNaVals <- function(dframe, percentage) {
  colsToRemove <- c()
  
  for (i in names(dframe)) {
    nas <- sum(is.na(dframe[i]))
    naratio <- nas / nrow(dframe)
    if (naratio > percentage) {
      print(paste(i, nas, round(naratio,2)))
      colsToRemove <- c(colsToRemove, which(names(dframe) == i))
    }
  }
  print(colsToRemove)
  #return(dframe[, -colsToRemove])
}


#' Performs PCA on a dataset
#' @param dframe The dataset we are working on (as data.frame)
#' @param graph Whether to display the relative graph or not
#' @return A new dataset undergone the PCA transformation
doPCA <- function(dframe, graph=FALSE) {
  if (graph) {
    fa.parallel(dframe, fa="pc", n.iter = 100)
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
    abline(h = 1, col="dimgray", lty=2, lwd=2)
  }
  # pca computation
  ppc <- preProcess(dframe, method = "pca", thresh = .95)
  result <- predict(ppc, dframe)
  #set.seed(1234)
  print(paste("DF after PCA : ", hopkins(result, nrow(result) - 1)))
  return(result)
}


# Performs a set of transformations to a dataset
pre.process.data <- function(dframe, pmodel=NULL, graph=FALSE) {
  result <- list()
  catgs <- which(sapply(factorize(dframe), class) == "factor")
  # Performs a set of transformations on a dataset
  if (is.null(pmodel)) {
    # Transformation: YeoJohnson & SpatialSign
    hskew <- which(names(dframe) %in% c("grassesiu","catiu","dogiu","hdmiu","aspergillusiu"))
    ppc_yj <- preProcess(dframe[, -hskew], method = c("center","scale","YeoJohnson"))
    dframe[, -hskew] <- predict(ppc_yj, newdata = dframe[, -hskew])
    
    ppc_ss <- preProcess(dframe[, hskew], method = c("center","scale","spatialSign"))
    dframe[, hskew] <- predict(ppc_ss, dframe[, hskew])
    
    result$prep.workflow <- list(ppc_yj, ppc_ss)
    result$data <- dframe
  }
  # Designed for test sets--the tranformations applied to the training set are stored
  # in list pmodel and applied to the test set.
  else {
    if (! is.list(pmodel))
      stop("The argument pmodel needs to be a list()")
    
    tmp <- dframe
    for (i in pmodel) {
      tmp <- predict(i, newdata = tmp)
    }
    result$data <- tmp
  }
  
  if (graph) {
    labels <- paste0(substr(colnames(dframe[, -catgs]), 1, 9), ".")
    par(mfrow=c(1,1), las=2, mar=c(6,4,1,1) + 0.1)
    boxplot(dframe[, -catgs], col = "red", names = labels)
    par(las = 1)
  }
  return(result)
}


#' Graphically traces the sum of squares for different number of clusters
#' @references "R in Action, R. I. KABACOFF, Manning 2015"
#' @param data The dataset we are working on (as data.frame)
#' @param nc The maximum number of clusters to consider
#' @param seed The seed used to make results reproducible
#' @param graph Whether to display the relative graph or not
wssplot <- function(data, nc=15, seed=1234, graph=FALSE) {
  wss <- (nrow(data) - 1) * sum(apply(data, 2, var))
  for (i in 2:nc){
    #set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)
  }
  if (graph) {
    plot(1:nc, wss, type="b", xlab="Number of Clusters",
         ylab="Within groups sum of squares", col="magenta")
  }
}


#' Finds the outliers by following the boxplot criteria and assigns a rank
#' to each observation based on the number of predictors that it appears
#' as an outlier point.
#' @param dframe The dataset we are working on (as data.frame)
#' @param catg A list containing the positions of the categorical fields
#' @param threshold The min number of predictors that an object is an outlier
#' @return A list of observations found to be outliers
get.outliers <- function(dframe, catg, threshold) {
  df <- dframe[, -catg]
  tmp <- list()
  for (i in 1:length(df)) {
    tmp[[i]] <- which(df[,i] %in% boxplot.stats(as.matrix(df[,i]))$out)
  }
  result <- plyr::count(unlist(tmp))
  result <- result[order(result$freq, decreasing = TRUE ),]
  return(result[result$freq >= threshold ,]$x)
}


#' Performs the dummification of a set of categorical variables
#' @param dfame The dataset we are working on (as data.frame)
#' @param vars An array with the positions of categorical variables
#' @param return The same dataframe with dummified variables
dummify.vars <- function(dframe, vars) {
  f <- paste(vars, collapse = "+")
  f <- as.formula(paste("~ ", f))
  dummies <- dummyVars(f, data = factorize(dframe), fullRank = TRUE)
  group <- data.frame(predict(dummies, newdata = factorize(dframe)))
  dframe <- cbind.data.frame(apply(group, 2, as.integer), dframe)
  dframe[, vars] <- NULL
  return(dframe)
}


#' Performs the K-Means clustering on a dataset
#' @param dframe The dataset we are working on (as data.frame)
#' @param k A list with several numbers of clusters
#' @param graph Whether to display the relative silhouette graph or not
#' @return A list with the k-means objects for each k
do.kmeans <- function(dframe, k, graph=FALSE) {
  if (is.null(k))
    stop("k needs to be provided")
  
  par(mfrow=c(1,2), las=1)
  result <- list()
  for (i in k) {
    print(paste("k-means algorithm with k = ", i))
    km <- kmeans(dframe, i, nstart=25)
    print(table(km$cluster))
    result[[i]] <- km
  
    if (graph) {
      plot(silhouette(km$cluster, dist(dframe, method = "euclidean")), 
           col=rainbow(i), main = paste("k-means with ", i, "clusters"))
    }
  }
  return(result)
}


#' Performs PAM clustering on a dataset and estimates the optimal k
#' @param dframe The dataset we are working on (as data.frame)
#' @param k The second best number of clusters to evaluate.
#' @param graph Whether to display the relative silhouette graph or not
#' @return A list with the PAM objects for each k
do.pam <- function(dframe, k=NULL, graph=FALSE) {
  if (length(k) > 1)
    stop("k should be a scalar")
  
  #set.seed(1234)
  pamk.res <- pamk(dframe)
  print(paste("PAM algorithm with suggested k = ", pamk.res$nc))
  print(table(pamk.res$pamobject$clustering))
  result <- list("pamk" = pamk.res$pamobject)
  
  # the second best number of clusters
  if (! is.null(k)) {
    pam2 <- pam(dframe, k=k, stand = FALSE)
    print(paste("PAM algorithm with second best k = ", k))
    print(table(pam2$clustering))
    result[["pam2"]] <- pam2
  }
  
  if (graph) {
    par(mfrow=c(1,2), las=1)
    for (i in result) {
      nclusters <- max(i$clustering)
      plot(silhouette(i$clustering, dist(dframe, method = "euclidean")), 
       col=rainbow(nclusters), main = paste("PAM with ", nclusters, "clusters"))
    }
  }
  return(result)
}


#' Performs hierarchical clustering with Ward.D2 linkage on a dataset
#' @param dframe The dataset we are working on (as data.frame)
#' @param k The (scalar) number of clusters to consider
#' @param graph Whether to display the relative silhouette graph or not
#' @return The HC object containing the clusters
do.hclustering <- function(dframe, k, graph=FALSE) {
  if (is.null(k))
    stop("k parameter missing. The number of clusters needs to be specified")
  
  print(paste("HC algorithm with k = ", k))
  hc <- hclust(dist(dframe), method="ward.D2")
  result <- cutree(hc, k=k)
  print(table(result))
  
  if (graph) {
    par(mfrow=c(1,1), mar=c(5,4,1,1) + 0.1)
    plot(hc, hang=-1, cex=.6, main="Ward.D2 Linkage Clustering")
    rect.hclust(hc, k=k)
  }
  return(result)
}


#' Plots a graph to estimate the optimal eps parameter for the algorithm
#' DBSCAN. It estimates how big eps should be in order to have K MinPts
#' in its neighbourhood. It is given by an elbow point 
#' @param dframe The dataset we are working on (as data.frame)
estimate.dbscan.pars <- function(dframe) {
  color <- rainbow_hcl(5)
  leg.text <- c()
  
  for (k in seq(3, 11, length=5)) {
    kNNdist <- sort(kNNdist(dframe, k))
    if (k == 3) {
      plot(kNNdist, type = "l", xlim = c(0,1700), ylab = "kNN distance (eps)", 
           col = color[1], xlab = "Points (sample) sorted by distance")
    }
    else {
      lines(sort(kNNdist), col = head(color, 1))
    }
    leg.text <-c(leg.text, paste("k=", k))
    color <- color[-1]
  }
  
  legend("bottomright", legend = leg.text, fill=rainbow_hcl(5), inset = 0.05)
  grid()
}


#' Helper function for reading interactively user's input
#' @return The numeric value inserted in the console
read.user.input <- function() { 
  n <- readline(prompt="Enter cutoff value, [0 to exit]: ")
  n <- try(as.numeric(n))
  if (is.na(n)) {
    print("Invalid input, enter a real number")
    n <- read.user.input()
  }
  return(as.numeric(n))
}


#' Performs the OPTICS algorithm into a dataset
#' @param dframe The dataset we are working on (as data.frame)
#' @param cutoff The cutting level on a reachability plot
#' @param minPts The minimum number of points consider for a core point
#' @param graph Whether to display the relative  or reachability plot not
#' @return The OPTICS object returned
do.optics <- function(dframe, cutoff, minPts, graph=FALSE) {
  print(paste("OPTICS algorithm with ", "MinPts = ", minPts, " eps = ", cutoff))
  desc <- paste("Reachability Plot, ", "MinPts = ", minPts, " eps = ", cutoff)
  res.optcs <- optics(dframe, minPts = minPts)
  res <- extractDBSCAN(res.optcs, eps_cl = cutoff)
  print(table(res$cluster))
  if (graph) {
    plot(res, main = desc)
  }
  return(res)
}


#' Performs the ORCLUS algorithm into a dataset
#' @param dframe The dataset we are working on (as data.frame)
#' @param k The final number of clusters
#' @param s The seed for reproducible results
#' @param threshold The dimension of the final cluster-specific subspace
#' @return The ORCLUS object returned
do.orclus <- function(dframe, k, s = 246, threshold = 0.1) {
  vals <- c()
  # Find a good value for l
  for (i in seq(1, length(dframe) * 0.5, 1)) {
    #set.seed(s)
    out <- orclus(dframe, k=k, l=i, k0=10, a = 0.25, verbose = FALSE)
    vals <- c(vals, out$sparsity.coefficient)
  }
  names(vals) <- seq(1, length(dframe) * 0.5, 1)
  vals <- sort(vals)
  #print(sapply(vals, function(x) {round(x, 3)} ))
  
  idx <- max(which(vals < threshold))
  l <- as.integer(names(vals)[idx])
  print(paste("ORCLUS algorithm with ", "k = ", k, " l = ", l))
  
  # Run ORCLUS algorithm
  #set.seed(s)
  out <- orclus(dframe, k=k, l=l, k0=10, a = 0.75, verbose = FALSE)
  return(out)
}

#' @return ROC, Sens, Spec, Accuracy and Kappa
fiveClassSummary <- function(...) {
  return(c(twoClassSummary(...), defaultSummary(...)))
}

#' Calculates the Cross-Entropy between two distributions
#' @param gtruth An array representing the groud truth
#' @param pred An array representing the estimated distribution
#' @return The value of cross-entropy between the two distributions
cross.entropy <- function(gtruth, pred) {
  levels(gtruth) <- levels(pred) <- union(levels(gtruth), levels(pred))
  gtruth <- table(gtruth) / length(gtruth)
  pred <- table(pred) / length(pred)
  gtruth <- ifelse(gtruth > 0, log2(gtruth), 0)
  return(as.numeric(-(pred %*% gtruth)))
}


# List of closures with the instructions to build the datasets
make.dataset <- list()
make.dataset[[1]] <- function(x) {x <- dummify.vars(x, 'group'); return(x)}
make.dataset[[2]] <- function(x) {x$group <- as.integer(ifelse(x$group == 0, 0, 1));
                       return(x)
                     }
make.dataset[[3]] <- function(x) {return(x[, -1])}
make.dataset[[4]] <- function(x) {x <- x[, -1]; 
                       cnt <- apply(x, 2, function(y) {sum(is.na(y)) / nrow(x)});
                       b <- which(cnt > 0.25); return(x[, -b])
                     }


## MAIN The execution starts here ##
data <- read.csv("data/63_variables_pre-imputation_for_missing.csv")

# Menu --------------------------------------------------------------------
ch <- menu(c("Dataset2 - Dummified group",
              "Dataset3 - Two groups (Control/Wheeze)",
              "Dataset1 - Without group variable",
              "Dataset4 - Removes Group and variables > 25% NA"),
           graphics = FALSE, title = "Chose a dataframe")

#df <- dframes[[ch]]
df <- make.dataset[[ch]](data)

imp <- mice(df, m = 1, seed = 123, maxit = 5, printFlag = FALSE)
df <- complete(imp, 1)
print(paste("Imputed DF : ", hopkins(df, nrow(df) - 1)))

df <- pre.process.data(df, graph = TRUE)$data
df.pca <- doPCA(df, graph = FALSE)

print(paste("Scaled DF(no PCA) : ", hopkins(df, nrow(df) - 1)))

# outlier.scores <- lof(df, k=5)
# outliers <- order(outlier.scores, decreasing=TRUE)
# outliers <- outlier(df)
# outliers <- sort(outliers, decreasing = TRUE)
# df <- df[-outliers[1:8],]
# data <- data[-outliers[1:8],]

# Estimates the number of clusters
pdf(file = NULL)
sink("/dev/null")
nc <- NbClust(data=df, diss=NULL, distance="euclidean", 
              min.nc=2, max.nc=10, method="kmeans")
sink()
dev.off()

par(mfrow=c(1,2))
print("Best number of clusters")
print(table(nc$Best.nc[1,]))
barplot(table(nc$Best.nc[1,]), col = "cyan", xlab="Number of Clusters", 
        ylab="Number of Criteria")
wssplot(df, nc=10, graph=TRUE)


k <- max(nc$Best.partition)
# the second best number of clusters
tmp <- table(nc$Best.nc[1,])
tmp <- sort(tmp, decreasing = TRUE)
k2 <- as.numeric(names(tmp)[2])


# PARTITIONING APPROACHES --------------------------------------------------------
# k-means                                              
kmeans <- do.kmeans(df, c(k, k2), graph = TRUE)

# PAM
#pam <- do.pam(df, k2, graph = TRUE)

# HIERARCHICAL APPROACH
# Hierarchical Clustering
hc2 <- do.hclustering(df, k, graph = TRUE)
hc2 <- sapply(hc2, function(x) {x-1})
hc3 <- do.hclustering(df, k2, graph = TRUE)
hc3 <- sapply(hc3, function(x) {x-1})


# DENSITY BASED APPROACHES
# DBSCAN
estimate.dbscan.pars(df)

#Evaluate OPTICS parameters
n <- -1
while(n != 0) {
  n <- read.user.input()
  if (n == 0)
    break
  eps <- n
  do.optics(df, n, 7, graph = TRUE)
}
  
# OPTICS
opts <- do.optics(df, eps, 7,  graph = FALSE)

# ORCLUS
orcls <- do.orclus(df, k = 2, s = 333, threshold = 0.1)

predictions <- list(
                    kM_1st = kmeans[[k]]$cluster,
                    kM_2nd = kmeans[[k2]]$cluster,
                    HC_1st = hc2,
                    HC_2nd = hc3,
                    OPTICS = opts$cluster,
                    ORCLS = orcls$cluster
)


# Heuristic Methods --------------------------------------------------------

# Initialises S3 data structures to store simulation results
simulations <- NULL

# Loops through the classification methods
for (j in METHODS) {
  timestamp()
  # Simulations loop
  for (s in 1:MC_Simulations) {
    cat("\n")
    print(paste("Epoch: ", s))
    # Split dataset into train/test.
    set.seed(NULL)
    sample <- sample(1:nrow(data), size=0.75 * nrow(data))
    print("Resampling the data")
    print(as.integer(rownames(data))[sort(sample)])
    
    # For dataset4 only. It should remove the same columns
    if (ch == 4) {
      df.train <- data[sample, which(names(data) %in% names(df))]
      df.test <- data[-sample, which(names(data) %in% names(df))]
    }
    else {
      df.train <- make.dataset[[ch]](data[sample,])
      df.test <- make.dataset[[ch]](data[-sample,])
    }
    
    imp <- mice(df.train, m = 1, seed = 123, printFlag = FALSE)
    df.train <- complete(imp, 1)
    df <- rbind.data.frame(df.train, df.test)
    imp <- mice(df, m = 1, seed = 123, printFlag = FALSE)
    df <- complete(imp, 1)
    df.test <- tail(df, nrow(df.test))
    
    training <- pre.process.data(df.train)
    # Test set pre-processing
    df.test <- data.frame(pre.process.data(df.test, pmodel = training$prep.workflow))
    colnames(df.test) <- gsub("data.", "", colnames(df.test))
    
    for (i in 1:length(predictions)) {
      cat("\n")
      print(paste(names(predictions)[i], "model"))
      # Assignment of attribute cluster to train & test datasets
      training$data$cluster <- as.factor(predictions[[i]][sample])
      df.test$cluster <- as.factor(predictions[[i]][-sample])
      #training$data$cluster <- as.factor(ifelse(predictions[[i]][sample] == 1, "one", "two"))
      #df.test$cluster <- as.factor(ifelse(predictions[[i]][-sample] == 1, "one", "two"))
    
      # Setting parameters research area
      pGrid.knn <- expand.grid(k = c(1,3))
      pGrid.NN <- expand.grid(size = 1:5, decay = seq(0.1, 0.9, 0.1))
      pGrid.dt <- expand.grid(cp = c(0.01, 0.05, 0.1, 0.25, 0.5))
      pGrid <- list(knn = pGrid.knn, nnet = pGrid.NN, rpart = pGrid.dt)
    
      fitCtrl <- trainControl(method = "repeatedcv",
                              number = 10,
                              repeats = 5,
                              sampling = "up",
                              summaryFunction = defaultSummary, # or fiveClassSummary
                              savePredictions = TRUE,
                              #classProbs = TRUE,       # Only for ROC estimation
                              allowParallel = TRUE
      )
      
      # Preparing the list of parameters for train function
      train.call <- list(form = cluster ~ .,
                         data = training$data, 
                         method = j,
                         trControl = fitCtrl,
                         tuneGrid = pGrid[[j]], 
                         metric = "Accuracy"            # ROC for auc estimation
      )
      
      # Disables the tracing for NN due to very high verbosity
      if (j == "nnet") {
        train.call <- c(train.call, trace = FALSE)
      }
      
      if (j == "rpart") {
        train.call <- c(train.call, tuneLength = 4)
      }
      
      # The call to caret::train with the list of parameters
      model <- tryCatch(
        {
          do.call(train, train.call)
        },
        error = function(e) {
          print("Failure in CARET::TRAIN")
          print(e)
        },
        finally = {
          # Save the partial results if anything goes wrong
          save(simulations, file = "mc_simulation_results.RData")
        }
      )
  
      # Displays the best tuned results    
      print(paste("Method: ", model$modelInfo$label))
      idx <- as.numeric(rownames(model$bestTune))
      print(model$results[idx ,])
      
      pred <- predict(model, df.test, na.action=na.omit)
      #pred <- predict(model, df.test, na.action=na.omit, type="prob")
      cc <- complete.cases(df.test)
      
      #print(table(pred, df.test$cluster[cc]))
      CE <- cross.entropy(df.test$cluster[cc], pred)
      print(paste("Cross Entropy H_p(q)= ", CE))
      
      # Fills a temporary record for the simulation data.frame
      tmp.data <- model$results[idx ,]
      tmp.data$ClsAlgo <- names(predictions)[[i]]
      tmp.data$SimNum <- s
      tmp.data$CrsEntropy <- CE
      
      simulations[[j]] <- rbind.data.frame(simulations[[j]], tmp.data)
      
      # if the model is Decision Tree, print the tree representation
      if (model$method == "rpart" & length(model$finalModel$splits) > 0) {
       plot(model$finalModel, main=paste(names(predictions)[i]))
       text(model$finalModel, pretty = 0, cex = 0.75, col="blue")
      }
    } # end validation of models
  }   # end simulations loop
  
  # Reset indexes
  rownames(simulations[[j]]) <- NULL
  
  write.table(simulations[[j]], paste0(j, "_val.csv"), sep=",", row.names=FALSE)
}

# Save the simulation data
save(simulations, file = "mc_simulation_results.RData")
  
