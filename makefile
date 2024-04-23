CC=gcc
BASEFLAGS= -Wall -lm
CFLAGS= $(BASEFLAGS) -O2
CDBGFLAGS= $(BASEFLAGS) -g
DEPS=vector_ops.h helper_functions.h hgmt_structs.h llist.h compton_chain_ordering.h

hgmt_lor_creator: hgmt_lor_creator.o vector_ops.o helper_functions.o llist.o compton_chain_ordering.o
	$(CC) -o hgmt_lor_creator $^  $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
