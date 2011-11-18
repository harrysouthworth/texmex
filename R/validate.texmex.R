validate.texmex <- function () {
   check <- "package:RUnit" %in% search()
   if (!check) {
       check <- try(library(RUnit))
       if (class(check) == "try-error") {
           stop("You need to attach the RUnit package to validate texmex")
       }
   }

   wh <- (1:length(search()))[search() == "package:texmex"]
   tests <- objects(wh)[substring(objects(wh), 1, 5) == "test."]

   # Create temporary directory to put tests into
   if (.Platform$OS.type == "windows"){ sep <- "\\" }
   else { sep <- "/" }

   dir <- file.path(tempdir(), "texmex.tests", fsep = sep)
#   cmd <- paste("mkdir", dir)
#   try(system(cmd), silent=TRUE)
   dir.create(dir)

   for (i in 1:length(tests)) {
       str <- paste(dir, sep, tests[i], ".R", sep = "")
       dump(tests[i], file = str)
   }
   res <- defineTestSuite("texmex", dirs = dir, testFuncRegexp = "^test.+", testFileRegexp = "*.R")
   cat("Running over 250 tests, including MCMC and bootstrap implementations.\nThis will take some time...\n\n")
   res <- runTestSuite(res)
   res
}

