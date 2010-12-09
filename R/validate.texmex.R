validate.texmex <- function(){
   check <- "package:svUnit" %in% search()
   if (!check){
      stop("You need to attach the svUnit package to validate texmex")
   }
   where <- (1:length(search()))[search() == "package:texmex"]
   res <- svSuiteList(pos=where)
   cat("Running over 100 tests, including MCMC and bootstrap implementations.\n
        This will take some time...\n")
   res <- runTest(res)
   res
}

