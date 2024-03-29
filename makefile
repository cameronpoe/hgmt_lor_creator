CC=gcc
BASEFLAGS= -Wall -lm
CFLAGS= $(BASEFLAGS) -O2
CDBGFLAGS= $(BASEFLAGS) -g
DEPS=vector_ops.h

hgmt_lor_creator: hgmt_lor_creator.o vector_ops.o
	$(CC) -o hgmt_lor_creator $^  $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

debug_hgmt_lor_creator: hgmt_lor_creator.c vector_ops.o
	$(CC) -o debug_hgmt_lor_creator $^ $(CDBGFLAGS)
