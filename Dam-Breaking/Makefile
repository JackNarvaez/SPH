all: Main

Main: main.cpp SPHfunctions.cpp SPHfunctions.h
	g++ $< SPHfunctions.cpp -o Main.x;\
	./Main.x

Plot: Plot.py
	python $<

clean:
	rm -f *.x *.txt *.png *.out *.gif