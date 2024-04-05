sift_dbg: sift_1b.cpp
	g++ -g -Ofast sift_1b.cpp -o sift_1b

sift: sift_1b.cpp
	g++ -Ofast -fopenmp sift_1b.cpp -o sift_1b

run: sift
	./sift_1b

memprof: sift
	heaptrack ./sift_1b
