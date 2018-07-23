MAIN = main.cpp
HPP = user_define.hpp class.hpp

two_body_initial.out:$(MAIN) $(HPP)
	g++ -Wall -O3 -o two_body_initial.out $(MAIN) -lm

clean:
	rm -f two_body_initial.out

distclean:
	make clean
	rm -f *~
