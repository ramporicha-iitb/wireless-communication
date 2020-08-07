all: mimo_ofdm alamouti/mimo_ofdm_alamouti precoding/unquantized precoding/vector_quantization

mimo_ofdm: mimo_ofdm.cpp
	g++ mimo_ofdm.cpp -litpp -o mimo_ofdm

alamouti/mimo_ofdm_alamouti: alamouti/mimo_ofdm_alamouti.cpp
	g++ alamouti/mimo_ofdm_alamouti.cpp -litpp -o alamouti/mimo_ofdm_alamouti

precoding/unquantized: precoding/unquantized.cpp
	g++ precoding/unquantized.cpp -litpp -o precoding/unquantized

precoding/vector_quantization: precoding/vector_quantization.cpp precoding/codebook.cpp
	g++ precoding/vector_quantization.cpp precoding/codebook.cpp -litpp -o precoding/vector_quantization

.PHONY: clean

clean:
	$(RM) mimo_ofdm alamouti/mimo_ofdm_alamouti precoding/unquantized
