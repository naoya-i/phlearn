
all: src/main.o src/ilp_loglinear.o src/storage.o src/util.o src/scw.o src/eqc.o
	LIBRARY_PATH=$$HOME/src/gurobi600/linux64/lib/; \
	g++ -g -std=c++11 -o bin/learn src/main.o src/ilp_loglinear.o src/storage.o src/util.o src/scw.o src/eqc.o \
	-L $$HOME/src/phillip/bin -l phillip -lgurobi_c++ -lgurobi60

.cpp.o:
	CPLUS_INCLUDE_PATH=$$HOME/src/phillip/src; \
	g++ -g -std=c++11 -c -o $(<:.cpp=.o) $<
