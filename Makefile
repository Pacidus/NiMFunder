CXXFLAGS = -Wall -Wextra
sources = main.cpp
entetes = afficher.h
objets = $(sources:.cc=.o)
%: %.o
	$(LINK.cc) -o$@ $ˆ

Sim: $(objets)

###

main.pdf: $(sources) $(entetes) Makefile
	u2ps -o -$ˆ | ps2pdf -@<

clean:
	rm -f ∗~ ∗.o ∗.bak

mrproper: clean
	rm -f hello

depend :
	makedepend $(sources)
