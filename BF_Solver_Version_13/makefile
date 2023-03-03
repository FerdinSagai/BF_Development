SOURCEBEAST = A.1.Fine_laden_gas_solver.cpp\

OBJECTBEAST = A.1.Fine_laden_gas_solver.o\

BEAST_COMPILE:
	g++ -c -mcmodel=large $(SOURCEBEAST)  -fopenmp
	g++ -o GF_SOLVER $(OBJECTBEAST) -fopenmp
	
clean:
	rm -f *.dat *.tmp *.o *.rep
	rm -f *.m
	rm -f GF_SOLVER

clear:
	rm -f *.dat *.tmp *.o *.rep
	rm -f *.m
	rm -f GF_SOLVER

reset:
	rm -f 2.* 3.* 4.* 5.* 10.*
