CC=gcc
CFLAGS=-g -Wall -pedantic
IDIR=src
ODIR=src

LIBS=-lfftw3 -lm -lpthread

_DEPS = pwsim.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = pwsim.o files.o geo.o utils.o angles.o fou.o tf.o wavegen.o calacv.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

pwsim: $(OBJ) $(DEPS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
