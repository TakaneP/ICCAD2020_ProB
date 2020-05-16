CXX		=	g++ -std=c++11 -O3
EXE		=   router	
OBJ_DIR	=	objs

CFLAGS = -g -O3

all: $(EXE)

$(EXE): $(OBJ_DIR)/main.o $(OBJ_DIR)/parser.o $(OBJ_DIR)/data_structure.o
	$(CXX) $^ $(CFLAGS) -o $@
	
$(OBJ_DIR)/main.o: main.cc parser.h data_structure.h | obj_dir
	$(CXX) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/parser.o: parser.cc parser.h data_structure.h | obj_dir
	$(CXX) $(CFLAGS) -c $< -o $@
	
$(OBJ_DIR)/data_structure.o: data_structure.cc data_structure.h | obj_dir
	$(CXX) $(CFLAGS) -c $< -o $@

obj_dir:
	mkdir -p $(OBJ_DIR)

clean:
	rm -rf *.gch objs/