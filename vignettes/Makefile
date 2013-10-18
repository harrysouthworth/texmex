all:
	pdflatex  texmex1d.tex
	bibtex  texmex1d
	pdflatex  texmex1d.tex
	pdflatex  texmex1d.tex
	rm *.log texmex1d.aux texmex1d.bbl texmex1d.blg texmex1d.tex
	pdflatex  texmexMultivariate.tex
	bibtex  texmexMultivariate
	pdflatex  texmexMultivariate.tex
	pdflatex  texmexMultivariate.tex
	rm *.log texmexMultivariate.aux texmexMultivariate.bbl texmexMultivariate.blg texmexMultivariate.tex	
	pdflatex  declustering.tex
	bibtex  declustering
	pdflatex  declustering.tex
	pdflatex  declustering.tex
	rm *.log declustering.aux declustering.bbl declustering.blg declustering.tex
	
texmexMultivariate:
	pdflatex  texmexMultivariate.tex
	bibtex  texmexMultivariate
	pdflatex  texmexMultivariate.tex
	pdflatex  texmexMultivariate.tex
	rm *.log texmexMultivariate.aux texmexMultivariate.bbl texmexMultivariate.blg texmexMultivariate.tex
	
texmex1d:
	pdflatex  texmex1d.tex
	bibtex  texmex1d
	pdflatex  texmex1d.tex
	pdflatex  texmex1d.tex
	rm *.log texmex1d.aux texmex1d.bbl texmex1d.blg texmex1d.tex
	
declustering:
	pdflatex  declustering.tex
	bibtex  declustering
	pdflatex  declustering.tex
	pdflatex  declustering.tex
	rm *.log declustering.aux declustering.bbl declustering.blg declustering.tex
	
clean:
	rm *.log *.aux *.bbl *.blg *.tex
  
