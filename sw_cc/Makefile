TARGET1 = router
CC = g++

all: $(TARGET1)

router: main.o parser.o data_structure.o
	$(CC) -o $(TARGET1) -O3 -std=c++11 main.o parser.o data_structure.o 
main.o: main.cc
	$(CC) -c -O3 main.cc 
parser.o: data_structure.o
	$(CC) -c -O3 -std=c++11 parser.cc 
data_structure.o: data_structure.cc
	$(CC) -c -O3 -std=c++11 data_structure.cc 

clean:
	rm -f *.o *.gch $(TARGET1)
