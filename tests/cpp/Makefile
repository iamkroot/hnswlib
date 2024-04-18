sift_dbg: sift_1b.cpp
	g++ -g -Ofast sift_1b.cpp -o sift_1b

sift: sift_1b.cpp
	g++ -Ofast -fopenmp sift_1b.cpp -o sift_1b

msift_dbg: sift_1b_multi.cpp
	g++ -g -Ofast sift_1b_multi.cpp -o sift_1b_multi

msift: sift_1b_multi.cpp
	g++ -Ofast sift_1b_multi.cpp -o sift_1b_multi
# 	g++ -Ofast -fopenmp sift_1b_multi.cpp -o sift_1b_multi

run: msift
	./sift_1b_multi

memprof: msift
	heaptrack ./sift_1b_multi
