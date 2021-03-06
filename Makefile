# Choose compiler
CC = gcc

# Optimization flags: cc is generic, gcc and icc are for PIII
CFLAGS_gcc = -g -Wall -DHAVE_INLINE #Usually gcc is needed for debugging
CFLAGS_icc = -O3 -xP  -DHAVE_INLINE

# Set linker flags: i.e. mpe libraries or options like -static
LFLAGS_gcc =

# Compilation and linking flags 
CFLAGS = $(CFLAGS_$(CC)) -c
LFLAGS = $(LFLAGS_$(CC)) -lm -lgsl -lgslcblas

C_SRC     = main.c init.c misc.c debug.c propagatorSympl.c forcesHamilton.c propagatorRK.c forcesDiscreet.c propagatorMP.c
OBJ_OSDH  = main.o init.o misc.o debug.o propagatorSympl.o forcesHamilton.o 
OBJ_OSDF  = main.o init.o misc.o debug.o propagatorSympl.o forcesDiscreet.o 
OBJ_RKDH  = main.o init.o misc.o debug.o propagatorRK.o    forcesHamilton.o 
OBJ_RKDF  = main.o init.o misc.o debug.o propagatorRK.o    forcesDiscreet.o 
OBJ_MPDH  = main.o init.o misc.o debug.o propagatorMP.o    forcesHamilton.o

.c.o:
	$(CC)  $(CFLAGS) $*.c

mpdh: $(OBJ_MPDH)
	$(CC)  $(LFLAGS) -o rod $(OBJ_MPDH)

osdh: $(OBJ_OSDH)
	$(CC)  $(LFLAGS) -o rod $(OBJ_OSDH)

osdf: $(OBJ_OSDF)
	$(CC)  $(LFLAGS) -o rod $(OBJ_OSDF)

rkdh: $(OBJ_RKDH)
	$(CC)  $(LFLAGS) -o rod $(OBJ_RKDH)

rkdf: $(OBJ_RKDF)
	$(CC)  $(LFLAGS) -o rod $(OBJ_RKDF)

clean:
		rm -f *.o rod
