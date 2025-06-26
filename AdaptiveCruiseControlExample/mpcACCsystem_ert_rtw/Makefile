.phony: all clean symbols

APP = mpcACCsystem
NaN = rtGetNaN
TYPES = rt_nonfinite

all : $(APP)

$(APP):
# -g adds debug information
	g++ -g -O0 -o $@ ert_main.cpp $@.cpp $(NaN).cpp $(TYPES).cpp -I. -std=c++14

run:
	./$(APP)

symbols:
	objdump --demangle -t $(APP)

clean:
	rm -rf $(APP)
