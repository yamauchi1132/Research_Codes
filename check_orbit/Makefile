MAIN = main.cpp

two_body_initial.out:$(MAIN) 
	g++ -Wall -O3 -o check_orbit.out $(MAIN) -lm

clean:
	rm -f check_orbit.out

distclean:
	make clean
	rm -f *~
