library(sqldf) 

for (i in c("knn", "nnet", "rpart")) {
  tryCatch(
    {
      tbl <- read.csv(paste0(i, "_val.csv"))
    },
    error=function(e) {
      print(paste0("Was not possible to read from file ", i, "_val.csv"))
      stop(e)
    }
  )
  
  tmp <- sqldf("select ClsAlgo, 
             avg(Accuracy) as Accuracy, 
               avg(Kappa) as Kappa, 
               avg(AccuracySD) as AccuracySD, 
               avg(KappaSD) as KappaSD, 
               avg(CrsEntropy) as CrsEntropy
               from tbl
               group by ClsAlgo 
               order by Accuracy desc, CrsEntropy asc, Kappa desc")
  
  cat('\n')
  print(tmp)
}




