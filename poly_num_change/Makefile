MAIN = main.cpp
HPP = user_define.hpp class.hpp

poly_num_change.out:$(MAIN) $(HPP)
	g++ -Wall -O3 -o poly_num_change.out $(MAIN) -lm

clean:
	rm -f poly_num_change.out

distclean:
	make clean
	rm -f *~
