all: ofdm ofdm_tdl 

ofdm: ofdm.cpp 
	g++ ofdm.cpp -litpp -o ofdm

ofdm_tdl: ofdm_tdl.cpp
	g++ ofdm_tdl.cpp -litpp -o ofdm_tdl

.PHONY: clean

clean:
	$(RM) ofdm ofdm_tdl
