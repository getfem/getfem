


all : doc.ps


doc.ps : doc.tex
	latex doc; latex doc;
	dvips -o doc_i.ps -p1 doc
	mpage -2 doc_i.ps > doc.ps
	rm -f doc_i.ps *.toc *.aux *.log *.bak *.dvi
