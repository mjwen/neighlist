CC = g++
CFLAGS = -O2 -g -Wall -I../src
LDFLAGS = -L../
LDLIBS = -lm -lneigh

SOURCES = graphite.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# linking
graphite: $(OBJECTS)
	$(CC) $(LDFLAGS)  -Wl,-rpath=../ -o $@ $^ $(LDLIBS)

# creating objects
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@


.PHONY: clean
clean:
	rm *.o graphite *.xyz
