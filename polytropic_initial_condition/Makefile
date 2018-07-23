MAIN = main.cpp
HPP = user_defined.hpp lane_emden.hpp make_and_output_data.hpp class.hpp

polytropic_condition.out:$(MAIN) $(HPP)
	g++ -Wall -O3 -o polytropic_condition.out $(MAIN) -lm

clean:
	rm -f polytropic_condition.out

distclean:
	make clean
	rm -f *~

