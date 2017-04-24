Text:=RCCResearchII-201509
Bib:=MultiTissue.bib
pdf:
	cp ~/Dropbox/Literature/Library/bibs/MultiTissue.bib ./
	pdflatex $(Text).tex && bibtex $(Text).aux && pdflatex $(Text).tex && pdflatex $(Text).tex
clean:
	find . -type f ! \( -name $(Text).pdf -o -name $(Text).tex -o -name Makefile -o -name $(Bib) \) -delete
