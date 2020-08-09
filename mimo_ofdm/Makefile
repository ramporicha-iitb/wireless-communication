all: mimo_ofdm alamouti/mimo_ofdm_alamouti precoding/unquantized precoding/vector_quantization precoding/scalar_parameter_quantization

mimo_ofdm: mimo_ofdm.cpp
	g++ -ggdb mimo_ofdm.cpp -litpp -o mimo_ofdm

alamouti/mimo_ofdm_alamouti: alamouti/mimo_ofdm_alamouti.cpp
	g++ -ggdb alamouti/mimo_ofdm_alamouti.cpp -litpp -o alamouti/mimo_ofdm_alamouti

precoding/unquantized: precoding/unquantized.cpp
	g++ -ggdb precoding/unquantized.cpp -litpp -o precoding/unquantized

precoding/vector_quantization: precoding/vector_quantization.cpp precoding/codebook.cpp
	g++ -ggdb precoding/vector_quantization.cpp precoding/codebook.cpp -litpp -o precoding/vector_quantization

precoding/scalar_parameter_quantization: precoding/scalar_parameter_quantization.cpp precoding/codebook.cpp precoding/parameterization.cpp
	g++ -ggdb precoding/scalar_parameter_quantization.cpp precoding/codebook.cpp precoding/parameterization.cpp -litpp -o precoding/scalar_parameter_quantization

.PHONY: clean

clean:
	$(RM) mimo_ofdm alamouti/mimo_ofdm_alamouti precoding/unquantized precoding/vector_quantization precoding/scalar_parameter_quantization
