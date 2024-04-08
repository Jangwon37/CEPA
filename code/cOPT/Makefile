# Makefile
CC=g++
CFLAGS=-m64 -c -Wall -O2
LDFLAGS=-m64
MAINOBJS=main.o 
OBJS=CenData.o Density.o WKM.o PHI.o OPT.o MyMath.o
TARGET=copt

all: $(OBJS) $(MAINOBJS)
	$(CC) -o $(TARGET) $(LDFLAGS) $(MAINOBJS) $(OBJS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o *.so tags $(TARGET) $(TARGET_TOOL) $(TARGET_TEST) $(OBJS)
