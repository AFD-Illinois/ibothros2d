# 
#
CC = gcc5
CFLAGS  = -Wall -Ofast 
LDFLAGS = -lm -lgsl -lgslcblas 
#

SRC = \
geodesics.c geodesics_gsl.c geometry_utils.c harm_model.c \
harm_utils.c init_iharm2dv3_data.c make_ppm.c jnu_mixed.c \
main.c radiation.c set_model_param.c tetrads.c
 
OBJ = \
geodesics.o geodesics_gsl.o geometry_utils.o harm_model.o \
harm_utils.o init_iharm2dv3_data.o make_ppm.o jnu_mixed.o \
main.o radiation.o set_model_param.o tetrads.o

ibothros2d: $(OBJ) makefile 
	$(CC) $(CFLAGS) -o ibothros2d $(OBJ) $(LDFLAGS)

$(OBJ) : makefile decs.h defs.h constants.h $(SRC)

clean:
	rm *.o ibothros2d

cleanup:
	rm ibothros2d*.ppm ibothros2d.dat Unsorted.txt


