# 
#
CC = gcc49
CFLAGS  = -Wall -Ofast -ftree-vectorizer-verbose=1 -fopenmp
LDFLAGS = -lm -lgsl -lgslcblas 
#
SRCMIB = \
geodesics.c geodesics_gsl.c geometry_utils.c harm_model.c \
harm_utils.c init_harm_data.c make_ppm.c jnu_mixed.c main.c radiation.c set_model_param.c tetrads.c
 
OBJMIB = \
geodesics.o geodesics_gsl.o geometry_utils.o harm_model.o \
harm_utils.o init_harm_data.o make_ppm.o jnu_mixed.o main.o radiation.o set_model_param.o tetrads.o

ibothros2d: $(OBJMIB) makefile 
	$(CC) $(CFLAGS) -o ibothros2d $(OBJMIB) $(LDFLAGS)

$(OBJMIB) : makefile decs.h defs.h constants.h

clean:
	rm *.o ibothros2d

cleanup:
	rm ibothros2d*.ppm ibothros2d.dat Unsorted.txt


