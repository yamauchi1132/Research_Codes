MAIN = main.cpp
HPP = user_define.hpp class.hpp

bound_check.out:$(MAIN)  $(HPP)
	g++ -Wall -O3 -o bound_check.out $(MAIN) -lm

clean:
	rm -f bound_check.out

distclean:
	make clean
	rm -f *~
