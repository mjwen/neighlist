CC = g++
CFLAGS = -O2 -g -Wall # include file should be added here using ``-I''
LDFLAGS =             # Extra flags to give to compilers when they are supposed
                      # to invoke the linker, ‘ld’, such as -L.
LDLIBS = -lm          # Library flags or names given to compilers when they are supposed
                      # to invoke the linker  e.g. -lm


SOURCES = example.cpp neighbor_list.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# linkding
neigh: $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

#.cpp.o:
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@


.PHONY: clean
clean:
	rm *.o neigh *.xyz
