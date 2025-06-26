.phony: all clean symbols

APP = dualACC
EGO = egoCar
NaN = rtGetNaN
TYPES = rt_nonfinite
LEAD = leadCar

all : $(APP)

$(APP):
# -g adds debug information
	g++ -g -O0 -o $@ main.cpp\
		$(EGO)_ert_rtw/$(EGO).cpp\
		$(EGO)_ert_rtw/$(NaN).cpp\
		$(EGO)_ert_rtw/$(TYPES).cpp\
		$(LEAD)_ert_rtw/$(LEAD).cpp\
		-I. -I$(EGO)_ert_rtw -I$(LEAD)_ert_rtw\
		-std=c++14

run:
	./$(APP)

symbols:
	objdump --demangle -t $(APP)

clean:
	rm -rf $(APP)
