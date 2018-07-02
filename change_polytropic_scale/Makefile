MAIN = main.cpp
HPP = user_define.hpp

change_scale.out:$(MAIN) $(HPP)
	g++ -Wall -O3 -o change_scale.out $(MAIN) -lm

clean:
	rm -f change_scale.out

distclean:
	make clean
	rm -f *~
