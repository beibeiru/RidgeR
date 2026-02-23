
print(RcppParallel::defaultNumThreads())
library(RidgeR)
dataPath <- file.path(system.file(package = "RidgeR"), "extdata/")
expr.diff <- read.table(paste0(dataPath, "Ly86-Fc_vs_Vehicle_logFC.txt"))

t_old <- system.time({res.old <- SecAct.inference.gsl.old(expr.diff)})
t_new <- system.time({res.new <- SecAct.inference.gsl.new(expr.diff)})
print("n = 1")
print(t_old); print(t_new)

expr.diff <- expr.diff[,rep(1,10)]
t_old <- system.time({res.old <- SecAct.inference.gsl.old(expr.diff)})
t_new <- system.time({res.new <- SecAct.inference.gsl.new(expr.diff)})
print("n = 10")
print(t_old); print(t_new)

expr.diff <- expr.diff[,rep(1,100)]
t_old <- system.time({res.old <- SecAct.inference.gsl.old(expr.diff)})
t_new <- system.time({res.new <- SecAct.inference.gsl.new(expr.diff)})
print("n = 100")
print(t_old); print(t_new)

expr.diff <- expr.diff[,rep(1,1000)]
t_old <- system.time({res.old <- SecAct.inference.gsl.old(expr.diff)})
t_new <- system.time({res.new <- SecAct.inference.gsl.new(expr.diff)})
print("n = 1000")
print(t_old); print(t_new)

expr.diff <- expr.diff[,rep(1,10000)]
t_old <- system.time({res.old <- SecAct.inference.gsl.old(expr.diff)})
t_new <- system.time({res.new <- SecAct.inference.gsl.new(expr.diff)})
print("n = 10000")
print(t_old); print(t_new)

expr.diff <- expr.diff[,rep(1,100000)]
t_old <- system.time({res.old <- SecAct.inference.gsl.old(expr.diff)})
t_new <- system.time({res.new <- SecAct.inference.gsl.new(expr.diff)})
print("n = 100000")
print(t_old); print(t_new)

#expr.diff <- expr.diff[,rep(1,1000000)]
#t_old <- system.time({res.old <- SecAct.inference.gsl.old(expr.diff)})
#t_new <- system.time({res.new <- SecAct.inference.gsl.new(expr.diff)})
#print("n = 1000000")
#print(t_old); 
#print(t_new)


