validate.texmex <-
function () {
   wh <- (1:length(search()))[search() == "package:sombrero"]
   tests <- objects(wh)[substring(objects(wh), 1, 5) == "test."]

   # Create temporary directory to put tests into
   if (.Platform$OS.type == "windows"){ sep <- "\\" }
   else { sep <- "/" }

   dir <- file.path(tempdir(), "texmex.tests", fsep = sep)
   dir.create(dir)

   for (i in 1:length(tests)) {
       str <- paste(dir, sep, tests[i], ".R", sep = "")
       dump(tests[i], file = str)
   }
   res <- defineTestSuite("texmex", dirs = dir, testFuncRegexp = "^test.+", testFileRegexp = "*.R")
   cat("Running over 1000 checks, including MCMC and bootstrap implementations.\nThis will take some time...\n\n")
   res <- runTestSuite(res)
   res
}
