EPStoPDF = epstopdf

ALL: paper.ps paper.pdf 

eps_source = $(wildcard *.eps)

pdf_source = $(eps_source:.eps=.pdf)

paper.dvi: paper.tex visualization.tex $(eps_source)
	latex paper.tex < /dev/null
	latex paper.tex < /dev/null
	latex paper.tex < /dev/null
	latex paper.tex < /dev/null

paper.pdf: paper.tex visualization.tex $(pdf_source)
	pdflatex paper.tex < /dev/null
	pdflatex paper.tex < /dev/null
	pdflatex paper.tex < /dev/null
	pdflatex paper.tex < /dev/null

pdf:	paper.pdf 

%.ps: %.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) paper.dvi 
	$(RM) paper.blg paper.log paper.pdf
	$(RM) paper.aux paper.ps paper.bbl paper.toc
	$(RM) *~

.PHONY: clean
