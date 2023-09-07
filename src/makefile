# Makefile for stepStone 
CC= gcc
CFLAGS= -O4
LFLAGS= -lm 

SOURCES=step_fastq.c step_shred.c step_commsPlot.c step_depthPlot.c step_freqPlot.c step_coverage.c step_construct.c step_linkStones.c step_breakSort.c step_Eulerpath.c step_chr2chr.c step_number.c step_short-number.c step_breakProcess.c step_cleanProcess.c stepBreakPoint.c stepBreakPoint-S.c step_processStones.c step_shortReads.c step_edgeStones.c step_sortStones.c stepStone.c

OBJS = $(patsubst %.c,%.o,$(SOURCES)) fast.o
EXECS = $(patsubst %.c,%,$(SOURCES))
EXECS_BIN = $(patsubst %.c,step-bin/%,$(SOURCES))
COMPILE = $(CC) $(CFLAGS) 


all:  cleanall iprint $(OBJS) executables clean oprint

executables:
	for exe in $(EXECS);  do $(COMPILE) -o $$exe $$exe.o fast.o $(LFLAGS); cp $$exe step-bin/.; done

%.o: %.c fasta.h
	$(CC) $(CFLAGS)  -c $< -o $@

iprint:
	@echo '+++ Compiling All ... '

oprint:
	@echo 'All Done '


clean: 
	@echo '+++ Cleaning Up ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS)
	@cp step-bin/stepStone .
	@cp step-bin/step_chr2chr .
	@cp step-bin/step_Eulerpath .

cleanall: 
	@echo '+++ Cleaning All ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS) fast.o
	@rm -f $(EXECS_BIN)
