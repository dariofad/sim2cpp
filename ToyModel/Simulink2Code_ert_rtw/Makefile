.phony: all clean symbols

TOY = Simulink2Code

all : $(TOY)

$(TOY):
# -g adds debug information
	g++ -g -O0 -o $@ main.cpp $@.cpp -I. -std=c++14

run:
	./$(TOY)

symbols:
	objdump --demangle -t $(TOY)

clean:
	rm -rf $(TOY)
