CXX		=	g++ -std=c++11 
EXE		=   router	
OBJ_DIR	=	objs
FLUTE_DIR = flute-3.1
CFLAGS = -g -O3

FLUTE_SRC     = $(FLUTE_DIR)/dist.c $(FLUTE_DIR)/dl.c $(FLUTE_DIR)/err.c $(FLUTE_DIR)/heap.c \
				$(FLUTE_DIR)/mst2.c $(FLUTE_DIR)/neighbors.c $(FLUTE_DIR)/flute.c $(FLUTE_DIR)/flute_mst.c
FLUTE_OBJ     = $(FLUTE_SRC:.c=.o)

all: $(EXE)

$(EXE): $(OBJ_DIR)/main.o $(OBJ_DIR)/parser.o $(OBJ_DIR)/data_structure.o $(OBJ_DIR)/segment_tree.o $(OBJ_DIR)/RSMT.o $(OBJ_DIR)/output.o $(FLUTE_OBJ)
	$(CXX) $^ $(CFLAGS) -o $@ 
	
$(OBJ_DIR)/main.o: main.cc parser.h data_structure.h segment_tree.h output.h | obj_dir
	$(CXX) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/parser.o: parser.cc data_structure.h parser.h segment_tree.h | obj_dir
	$(CXX) $(CFLAGS) -c $< -o $@
	
$(OBJ_DIR)/RSMT.o: RSMT.cc data_structure.h RSMT.h | obj_dir
	$(CXX) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/data_structure.o: data_structure.cc data_structure.h | obj_dir
	$(CXX) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/segment_tree.o: segment_tree.cc segment_tree.h | obj_dir
	$(CXX) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/output.o: output.cc output.h
	$(CXX) $(CFLAGS) -c $< -o $@

obj_dir:
	mkdir -p $(OBJ_DIR)

clean:
	rm -rf *.gch objs/ && cd $(FLUTE_DIR) && make clean && cd ..
