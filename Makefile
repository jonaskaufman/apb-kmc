CXX = g++-9
CXXFLAGS = -std=c++17 -I submodules/nlohmann-json/include -I submodules/kmc-lotto/include -O3 -DNDEBUG
OBJS = apb-kmc.o grid.o calculator.o simulation.o wrapper.o

apb-kmc: $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS)
	rm -f apb-kmc
