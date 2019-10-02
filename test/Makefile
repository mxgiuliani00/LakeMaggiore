

####### Compiler, tools and options

CC            = gcc
CXX           = g++
CXXFLAGS      = -c -g -Wall 

####### Compile
all: LakeMaggioreSim

LakeMaggioreSim: main_lakemaggiore.o model_lakemaggiore.o param_function.o ncRBF.o rbf.o ann.o pwLinear.o catchment.o lake.o lakemaggiore.o maggiore_panperduto.o moeaframework.o utils.o
	$(CXX) main_lakemaggiore.o model_lakemaggiore.o param_function.o ncRBF.o rbf.o ann.o pwLinear.o catchment.o lake.o lakemaggiore.o maggiore_panperduto.o moeaframework.o utils.o -o LakeMaggioreSim

main_lakemaggiore.o: main_lakemaggiore.cpp ../src/model_lakemaggiore.h \
	../src/moeaframework.h
	$(CXX) $(CXXFLAGS) main_lakemaggiore.cpp

model_lakemaggiore.o: ../src/model_lakemaggiore.cpp ../src/lake.h \
	../src/lakemaggiore.h \
	../src/catchment.h \
	../src/maggiore_panperduto.h \
	../src/param_function.h \
	../src/ncRBF.h \
	../src/rbf.h \
	../src/ann.h \
	../src/pwLinear.h
	$(CXX) $(CXXFLAGS) ../src/model_lakemaggiore.cpp

ncRBF.o: ../src/ncRBF.cpp ../src/ncRBF.h \
	../src/param_function.h
	$(CXX) $(CXXFLAGS) ../src/ncRBF.cpp

rbf.o: ../src/rbf.cpp ../src/rbf.h \
	../src/param_function.h
	$(CXX) $(CXXFLAGS) ../src/rbf.cpp

ann.o: ../src/ann.cpp ../src/ann.h \
	../src/param_function.h
	$(CXX) $(CXXFLAGS) ../src/ann.cpp

pwLinear.o: ../src/pwLinear.cpp ../src/pwLinear.h \
	../src/param_function.h
	$(CXX) $(CXXFLAGS) ../src/pwLinear.cpp

param_function.o: ../src/param_function.cpp ../src/param_function.h
	$(CXX) $(CXXFLAGS) ../src/param_function.cpp

catchment.o: ../src/catchment.cpp ../src/catchment.h
	$(CXX) $(CXXFLAGS) ../src/catchment.cpp

maggiore_panperduto.o: ../src/maggiore_panperduto.cpp ../src/maggiore_panperduto.h
	$(CXX) $(CXXFLAGS) ../src/maggiore_panperduto.cpp

lakemaggiore.o: ../src/lakemaggiore.cpp ../src/lakemaggiore.h \
	../src/lake.h
	$(CXX) $(CXXFLAGS) ../src/lakemaggiore.cpp

lake.o: ../src/lake.cpp ../src/lake.h
	$(CXX) $(CXXFLAGS) ../src/lake.cpp

moeaframework.o: ../src/moeaframework.c ../src/moeaframework.h
	$(CXX) $(CXXFLAGS) ../src/moeaframework.c

utils.o: ../src/utils.cpp ../src/utils.h
	$(CXX) $(CXXFLAGS) ../src/utils.cpp

clean:
	rm -rf *.o
	rm LakeMaggioreSim
