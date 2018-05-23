CC = g++
CFLAGS = -O2 -g -Wall -fPIC
LDFLAGS =
LDLIBS = -lm

SOURCES = src/neighbor_list.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# linking
libneigh.so: $(OBJECTS)
	$(CC) $(LDFLAGS) -shared -o $@ $^ $(LDLIBS)

# creating objects
src/%.o: src/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@


.PHONY: clean
clean:
	rm *.so src/*.o examples/*.o examples/graphite
