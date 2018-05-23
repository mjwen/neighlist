CXX = g++
CXXFLAGS = -std=c++11 -g -Wall			# include file should be added here using ``-I''
LDFLAGS =             # Extra flags to give to compilers when they are supposed
                      # to invoke the linker, ‘ld’, such as -L.
LDLIBS = -lm          # Library flags or names given to compilers when they are supposed
                      # to invoke the linker  e.g. -lm


SOURCES = main.cpp neighbor_list.cpp padding.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# linkding
exec: $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

.o:
	$(CXX) $(CXXFLAGS) $< -o $@



.PHONY: clean
clean:
	rm *.o exec *.xyz
