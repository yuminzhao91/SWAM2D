NAME= manual

all: $(NAME).pdf

$(NAME).pdf:

%.pdf : %.dvi
	dvipdf ${NAME}

%.dvi : %.tex
	latex ${NAME}
	bibtex ${NAME}
	makeindex -s ${NAME}.ist ${NAME}
	latex ${NAME}
	latex ${NAME}

clean:
	rm -f *~ *.out *.aux *.dvi *.ps *.pdf *.toc *.lof *.log *.lot *.flg *.bbl *.blg *.idx *.ind *.ilg *.snm *.nav *.fig.bak
