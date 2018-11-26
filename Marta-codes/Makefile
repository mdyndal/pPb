F_FILES = $(wildcard source/*.f)
OBJ_DIR = obj
OBJ_FILES = $(patsubst source/%.f,$(OBJ_DIR)/%.o,$(F_FILES))
LIB_FILES = main.o $(OBJ_FILES)
VPATH = source/

###############################################################################

FC = gfortran
FFLAGS = -g -O0 -Wunused

RFLAGS = $(shell root-config --cflags)
RLIBS = $(shell root-config --libs)
RHEAD = $(shell root-config --incdir)

CC = clang++
CFLAGS = -lgfortran -Wall -I$(RHEAD) -g

RM = /bin/rm
RMFLAGS = -rf

.PHONY: all
all: generator cgenerator

#
# Make the executable
#
generator: $(LIB_FILES)
	$(FC) $(FFLAGS) -o $@ $^

nice:
	$(RM) $(RMFLAGS) *.o $(OBJ_DIR)

clean: nice
	$(RM) $(RMFLAGS) generator cgenerator main.oxx

cgenerator: main.oxx $(OBJ_FILES)
	$(CC) $(CFLAGS) -o $@ $^ $(RLIBS)
#
# Make the objects
#
$(OBJ_DIR)/%.o: %.f source/%.f
	$(FC) -c $(FFLAGS) $< -o $@

%.oxx: %.cxx
	$(CC) $(CFLAGS) -c $(RFLAGS) $< -o $@

$(OBJ_FILES): | $(OBJ_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)
