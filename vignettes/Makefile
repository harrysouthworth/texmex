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
	pdflatex  egp3.tex
	bibtex  egp3
	pdflatex egp3.tex
	bibtex egp3
	pdflatex egp3.tex
	rm *.log egp3.aux egp3.tex egp3.toc egp3.bbl egp3.blg
	pdflatex test_texmex.tex
	bibtex test_texmex
	pdflatex test_texmex.tex
	rm *.log test_texmex.aux test_texmex.tex test_texmex.toc test_texmex.blg test_texmex.bbl
	mv texmex1d.pdf ../inst/doc
	mv texmexMultivariate.pdf ../inst/doc
	mv declustering.pdf ../inst/doc
	mv egp3.pdf ../inst/doc
	mv test_texmex.pdf ../inst/doc
	
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
	
egp3:
	pdflatex  egp3.tex
	bibtex  egp3
	pdflatex egp3.tex
	bibtex egp3
	pdflatex egp3.tex
	pdflatex egp3.tex
	rm *.log egp3.aux egp3.toc egp3.bbl egp3.blg egp3.tex
	
test_texmex:
	pdflatex test_texmex.tex
	bibtex test_texmex
	pdflatex test_texmex.tex
	pdflatex test_texmex.tex
	rm *.log test_texmex.aux test_texmex.bbl test_texmex.blg  test_texmex.tex
	
clean:
	rm *.log *.aux *.bbl *.blg *.tex
