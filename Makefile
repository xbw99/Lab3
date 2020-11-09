all:
	scc main -ww -vv -g -G
test:	
	./main beachbus.pgm
clean:
	rm -rf *o main
	rm -rf *.cc *.h *.o *.si
	rm -rf ba0.pgm
