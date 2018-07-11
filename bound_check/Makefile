MAIN = main.cpp

bound_check.out:$(MAIN) 
	g++ -Wall -O3 -o bound_check.out $(MAIN) -lm

clean:
	rm -f bound_check.out

distclean:
	make clean
	rm -f *~
