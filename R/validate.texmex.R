validate.texmex <- function(){
   check <- "package:RUnit" %in% search()
   if (!check){
      check <- try(library(RUnit))
      if (class(check) == "try-error"){
          stop("You need to attach the RUnit package to validate texmex")
      }
   }

   # Create a temporary directory to store tests scripts
   tempDir <- paste(sample(c(letters, LETTERS, 0:9), size=20), collapse="")
   tempDir <- paste(".", tempDir, sep = "")
   on.exit(system(paste("rm -r", tempDir)))
   system(paste("mkdir", tempDir))

   # Get tests and dump to temporary directory
   wh <- (1:length(search()))[search() == "package:texmex"]
   tests <- objects(wh)[substring(objects(wh), 1, 5) == "test."]
   for (i in 1:length(tests)){
      dump(tests[i], file=paste(tempDir, "/", tests[i], ".R", sep = ""))
   }

   res <- defineTestSuite("texmex", dirs=tempDir,
                          testFuncRegexp="^test.+",
                          testFileRegexp="*.R")
   cat("Running over 100 tests, including MCMC and bootstrap implementations.\nThis will take some time...\n\n")

   res <- runTestSuite(res)

   res
}

