PAPERS := paper.pdf paper_en.pdf

.PHONY: all check clean

all: check $(PAPERS)

%.pdf: %.tex
	latexmk -pdf -interaction=nonstopmode -halt-on-error $<

check:
	python3 revision_por_patas.py

clean:
	latexmk -c paper.tex paper_en.tex
