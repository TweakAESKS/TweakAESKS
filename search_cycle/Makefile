# Generic Makefile for compiling a simple executable.


CC := g++
#CC := icpc
SRCDIR := src
BUILDDIR := build
USERDEFINES :=
#CFLAGS := -fast -s -Wall -Wno-unused-variable -std=c++11 -DNDEBUG -opt-multi-version-aggressive -xHost
#CFLAGS := -fast -s -Wall -Wno-unused-variable -std=c++11 -DNDEBUG -opt-multi-version-aggressive -xHost -prof-use -prof-dir=profile/
#CFLAGS := -fast -Wall -Wno-unused-variable -std=c++11 -DNDEBUG -opt-multi-version-aggressive -xHost -prof-gen -prof-dir=profile/
CFLAGS := -O3 -s -mtune=native -march=native -Wall -Wno-unused-variable -std=c++11 -DNDEBUG -fopenmp
#CFLAGS := -O3 -s -Wall -Wno-unused-variable -std=c++11 -DNDEBUG
#CFLAGS := -O3 -march=native -Wall -Wno-unused-variable -std=c++11 -DNDEBUG
#CFLAGS := -g -march=native -Wall -Wno-unused-variable -std=c++11 -DNDEBUG
#CFLAGS := -O2 -g -march=native -Wall -Wno-unused-variable -std=c++11 -DNDEBUG
#CFLAGS := -O2 -g -pg -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls -march=native -Wall -Wno-unused-variable -std=c++11 -DNDEBUG
#CFLAGS := -O2 -g -fno-omit-frame-pointer -march=native -Wall -Wno-unused-variable -std=c++11 -DNDEBUG
#CFLAGS := -O3 -std=c++11 -DNDEBUG
#CFLAGS := -O3 -s -march=core2 -Wall -Wno-unused-variable -std=c++11 -DNDEBUG -static

LIBS := -lboost_program_options

TARGET := searchBestPerm

SOURCES := $(shell find $(SRCDIR) -type f -name *.cpp)
HEADERS := $(shell find $(SRCDIR) -type f -name *.hpp)



#OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.cpp=.o))

#DEPS := $(OBJECTS:.o=.deps)
DEPS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.cpp=.deps))


all: $(TARGET)

$(TARGET): $(OBJECTS)
	@echo " Linking ..."; $(CC) $(USERDEFINES) $(CFLAGS) $^ $(LIBS) -o $(TARGET)




$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp  $(HEADERS)
	@mkdir -p $(BUILDDIR)
	@echo " CC $<"; $(CC) $(USERDEFINES) $(CFLAGS) -MD -MF $(@:.o=.deps) -c -o $@ $<

clean:
	@echo " Cleaning..."; $(RM) -r $(BUILDDIR) $(TARGET) *~

-include $(DEPS)

#$(BUILDDIR)/GFElement.o:   .FORCE

.FORCE:

.PHONY: clean .FORCE
