CC=gcc
CFLAGS=-g -Wall -Wextra -pedantic
IDIR=src
ODIR=src
SAN_ODIR=build/asan

LIBS=-lfftw3 -lm -lpthread

_DEPS = pwsim.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = pwsim.o files.o geo.o utils.o angles.o fou.o tf.o wavegen.o calacv.o globals.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))
SAN_OBJ = $(patsubst %,$(SAN_ODIR)/%,$(_OBJ))

SAN_FLAGS=-fsanitize=address,undefined -fno-omit-frame-pointer -g -Wall -Wextra -pedantic


$(ODIR)/%.o: $(IDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(SAN_ODIR)/%.o: $(IDIR)/%.c $(DEPS)
	@mkdir -p $(SAN_ODIR)
	$(CC) -c -o $@ $< $(SAN_FLAGS)

pwsim: $(OBJ) $(DEPS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

pwsim_san: $(SAN_OBJ) $(DEPS)
	$(CC) -o $@ $^ $(SAN_FLAGS) $(LIBS)

sanitize: pwsim_san

valgrind: pwsim
	valgrind --leak-check=full --track-origins=yes --error-exitcode=1 ./pwsim params/params_nsf_full.dat

.PHONY: clean clean_san sanitize valgrind

clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ pwsim

clean_san:
	rm -rf $(SAN_ODIR) pwsim_san
