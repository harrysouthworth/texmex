#! /bin/bash

cd /home/harry/Work/repos/github/texmex/vignettes
mv declustering.knitr declustering.Rnw
mv texmex1d.knitr texmex1d.Rnw
mv texmexMultivariate.knitr texmexMultivariate.Rnw
mv test_texmex.knitr test_texmex.Rnw
mv egp3.knitr egp3.Rnw

R < build.R --no-save

mv declustering.Rnw declustering.knitr
mv texmex1d.Rnw texmex1d.knitr
mv texmexMultivariate.Rnw texmexMultivariate.knitr
mv test_texmex.Rnw test_texmex.knitr
mv egp3.Rnw egp3.knitr

make all
