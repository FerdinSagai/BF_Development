SOURCEIP = A.0.GUI_Standin.f90\

OBJECTIP = A.0.GUI_Standin.o\

BEAST_COMPILE:
	gfortran -c $(SOURCEIP) 
	gfortran -o INPUT_GUI $(OBJECTIP)
	
clean:
	rm -f *.dat *.tmp *.o *.rep
	rm -f *.m
	rm -f INPUT_GUI

clear:
	rm -f *.dat *.tmp *.o *.rep
	rm -f *.m
	rm -rf Results Executables
	rm -f INPUT_GUI
	
reset:
	rm -f 2.* 3.* 4.* 5.* 10.*
