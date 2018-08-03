MAIN = main.cpp
HPP = class.hpp

energy_check.out:$(MAIN)  $(HPP)
	g++ -Wall -O3 -o energy_check.out $(MAIN) -lm

clean:
	rm -f energy_check.out

distclean:
	make clean
	rm -f *~
