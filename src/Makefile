myProgram: main.o initial_conditions.o prim_to_cons.o solvers.o outputFile.o header.H
	g++ main.o initial_conditions.o prim_to_cons.o solvers.o outputFile.o -o myProgram -O3

main_new.o: main.cpp header.H
	g++ -c main_new.cpp -o main_new.o -O3

initial_conditions.o: initial_conditions.cpp header.H
	g++ -c initial_conditions.cpp -o initial_conditions.o -O3

prim_to_cons.o: prim_to_cons.cpp header.H
	g++ -c prim_to_cons.cpp -o prim_to_cons.o -O3

solvers.o: solvers.cpp header.H
	g++ -c solvers.cpp -o solvers.o -O3

outputFile.o: outputFile.cpp header.H
	g++ -c outputFile.cpp -o outputFile.o -O3
