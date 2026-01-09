CXX = g++
CXXFLAGS = -std=c++14 -Wall -Wextra
TARGET = fuel_control
SOURCES = main.cpp AbstractFuelControl_M1.cpp \
	AbstractFuelControl_M1_data.cpp \
	rt_nonfinite.cpp \
	rt_logging.cpp \
	rt_zcfcn.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET) \
	-I.

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

run:
	./$(TARGET)

.PHONY: all clean
