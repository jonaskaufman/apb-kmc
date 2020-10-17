CXX = clang++
CXXFLAGS = -std=c++11 -O3
OBJS = apb-kmc.o grid.o simulation.o wrapper.o

apb-kmc: $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS)
	rm -f apb-kmc