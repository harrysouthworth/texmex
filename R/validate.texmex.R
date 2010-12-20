validate.texmex <- function(testDir){
   check <- "package:RUnit" %in% search()
   if (!check){
      check <- try(library(RUnit))
      if (class(check) == "try-error"){
          stop("You need to attach the RUnit package to validate texmex")
      }
   }
   if (missing(testDir)){
       testDir <- paste(.path.package("texmex"), "R", sep = "/")
   }
   
   res <- defineTestSuite("texmex", dirs=testDir,
                          testFuncRegexp="^test.+",
                          testFileRegexp="*.R")
   cat("Running over 100 tests, including MCMC and bootstrap implementations.\nThis will take some time...\n\n")

   res <- runTestSuite(res)

   res
}

