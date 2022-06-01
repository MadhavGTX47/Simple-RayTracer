shade : parser.o
	g++ parser.o -o shade

parser.o : parser.cpp
	g++ -I/usr/include/eigen3/ -O3 -c parser.cpp 

clean:
	rm *.o *.ppm *.exe output out
	 
