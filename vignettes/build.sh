#! /bin/bash

cd /home/harry/Work/repos/googlecode/texmex/vignettes
mv declustering.knitr declustering.Rnw
mv texmex1d.knitr texmex1d.Rnw
mv texmexMultivariate.knitr texmexMultivariate.Rnw

R < build.R --no-save

mv declustering.Rnw declustering.knitr
mv texmex1d.Rnw texmex1d.knitr
mv texmexMultivariate.Rnw texmexMultivariate.knitr

make all
