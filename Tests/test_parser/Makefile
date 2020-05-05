TARGET1 = parser
CC = g++

all: $(TARGET1)

parser: parser.cc
	$(CC) -o $(TARGET1) -O3 -std=c++11 parser.cc parser.h

clean:
	rm -f *.o *.gch $(TARGET1)
