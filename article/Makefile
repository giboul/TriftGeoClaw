.Phony: all

MAIN ?= ImpulsionAvalanches

main:
	pdflatex --shell-escape $(MAIN).tex

bib:
	bibtex $(MAIN).aux

all:
	make main
	make bib
	make main
	make main