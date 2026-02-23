
print(RcppParallel::defaultNumThreads())
library(RidgeR)
dataPath <- file.path(system.file(package = "RidgeR"), "extdata/")
expr.diff <- read.table(paste0(dataPath, "Ly86-Fc_vs_Vehicle_logFC.txt"))

expr.diff <- expr.diff[,rep(1,1000000)]
#t_old <- system.time({res.old <- SecAct.inference.gsl.old(expr.diff)})
t_new <- system.time({res.new <- SecAct.inference.gsl.new(expr.diff)})
print("n = 1000000")
#print(t_old); 
print(t_new)

