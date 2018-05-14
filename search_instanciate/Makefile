JAVAC=javac
sources = $(wildcard *.java)
classes = $(sources:.java=.class)

all: $(classes)

%.class : %.java
	$(JAVAC) -sourcepath ./ -cp /home/blambin/CPSolver/choco-3.3.1/choco-solver-3.3.1-with-dependencies.jar $<

clean :
	rm -f *.class

