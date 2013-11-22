#! /bin/tcsh

cd /home/kpzv097/googlecode/texmex/vignettes
mv declustering.knitr declustering.Rnw
mv texmex1d.knitr texmex1d.Rnw
mv texmexMultivariate.knitr texmexMultivariate.Rnw

R < build.R

mv declustering.Rnw declustering.knitr
mv texmex1d.Rnw texmex1d.knitr
mv texmexMultivariate.Rnw texmexMultivariate.knitr

make all
