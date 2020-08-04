all: mimo_ofdm mimo_ofdm_alamouti 

mimo_ofdm: mimo_ofdm.cpp
	g++ mimo_ofdm.cpp -litpp -o mimo_ofdm

mimo_ofdm_alamouti: alamouti/mimo_ofdm_alamouti.cpp
	g++ alamouti/mimo_ofdm_alamouti.cpp -litpp -o alamouti/mimo_ofdm_alamouti

.PHONY: clean

clean:
	$(RM) mimo_ofdm alamouti/mimo_ofdm_alamouti
