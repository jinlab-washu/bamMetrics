sampath = \/opt\/samtools\/bin\/samtools/

bamMetrics: bamMetrics.cpp
	sed -i.bak 's/$$SAMTOOLS/$(sampath)' bamMetrics.cpp
	$(CC) bamMetrics.cpp -O2 -o bamMetrics

clean:
	rm bamMetrics.cpp
	mv bamMetrics.cpp.bak bamMetrics.cpp
	rm bamMetrics
