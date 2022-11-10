build:
	gcc src/twave.c -lm -shared -fPIC -o lib/libtwave.so

test:
	gcc test.c lib/libtwave.so -lm -o test.out
	./test.out
clean:
	rm -f *.out