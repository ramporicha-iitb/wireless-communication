#include <cassert>
#include <itpp/itcomm.h>

using namespace std;
using namespace itpp;

#define DEBUG 0

#if DEBUG
#define DEBUG_MSG(str) do { std::cout << str << std::endl; } while( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif

#define BATCH_SIZE          1e3
#define DEBUG_SYM_COUNT     4

#define CP_LEN              16
#define NUM_SUBCARRIER      128
#define NUM_OFDM_FRAMES     1e5
#define OUTPUT_FILENAME     "ofdm_tdl.it"

int main()
{
    RNG_randomize();
  
    int batch_size;
    int Nsub, Ncp, Nframe;

    Ncp = CP_LEN;
    Nsub = NUM_SUBCARRIER;
    Nframe = NUM_OFDM_FRAMES;
    batch_size = BATCH_SIZE;

    assert(Nframe % batch_size == 0);

    // Modulators
    OFDM ofdm;
    ofdm.set_parameters(Nsub, Ncp);
    BPSK_c modulator;
    DEBUG_MSG("Symbols: " << modulator.get_symbols());
    
    // Channels
    AWGN_Channel awgn_channel;

    double Es, Eb;
    Es = 1.0;                                   // Energy per symbol
    Eb = Es / modulator.bits_per_symbol();
    
    vec EbN0dB = "0:15";                        // Eb/N0 values in dB
    vec N0 = Eb * pow(inv_dB(EbN0dB), -1.0);    // N0 : variance of the complex valued noise
     
    // Transmitted bits / symbols
    bvec tx_bits, rx_bits;
    cvec tx_symbols;
    cvec rx_symbols;
    cvec tx_ofdm_symbols;
    cvec rx_ofdm_symbols;
    cvec rx_symbols_equalized;

    // Bit error rates 
    vec ber;
    BERC berc;
    ber.set_size(EbN0dB.length());

    // Channel coefficient row numbers corresponding to the begining of the OFDM frame
    ivec channel_rows = "0:" + to_string(Nsub + Ncp) + ":" + to_string(batch_size * (Nsub + Ncp) - 1);

    Real_Timer timer;
    timer.tic();
    
    for (int snr_loop = 0; snr_loop < EbN0dB.length(); ++snr_loop) 
    {
        size_t err_count = 0;
        size_t total_bit_count = 0;
        
        TDL_Channel tdl_channel(vec("0 -3 -6 -9"), ivec("0 1 2 3"));
        tdl_channel.set_norm_doppler(1e-6);

        for (int frame_loop = 0; frame_loop < Nframe / batch_size; ++frame_loop)
        {
            /*-------------------------- TRANSMITTER -----------------------------*/
            // Generate a vector of random bits to transmit        
            int num_tx_bits = batch_size * Nsub * modulator.bits_per_symbol();
            tx_bits = randb(num_tx_bits);
            DEBUG_MSG("Number of Tx bits: " << num_tx_bits);
            DEBUG_MSG("Sample Tx bits: " << tx_bits.left(DEBUG_SYM_COUNT * modulator.bits_per_symbol()));

            // Modulate the bits to symbols
            tx_symbols = modulator.modulate_bits(tx_bits);
            DEBUG_MSG("Number of Tx symbols: " << tx_symbols.length());
            DEBUG_MSG("Sample Tx symbols: " << tx_symbols.left(DEBUG_SYM_COUNT));

            // OFDM modulation
            tx_ofdm_symbols = ofdm.modulate(tx_symbols);
            DEBUG_MSG("Number of Tx OFDM symbols: " << tx_ofdm_symbols.length());
            DEBUG_MSG("Sample Tx OFDM symbols: " << tx_ofdm_symbols.left(DEBUG_SYM_COUNT));

            /*----------------------- TDL + AWGN CHANNEL -------------------------*/
            // Set the noise variance of the AWGN channel
            awgn_channel.set_noise(N0(snr_loop));

            // Generate channel filter coefficients
            // Returns a matrix with one tap per column
            cmat tdl_channel_coeff;
            tdl_channel.generate(batch_size * (Nsub + Ncp), tdl_channel_coeff);
            DEBUG_MSG("Channel coefficient matrix size: " << tdl_channel_coeff.rows() << "x" << tdl_channel_coeff.cols());
            DEBUG_MSG("Channel coefficients, row[0]: " << tdl_channel_coeff.get_row(0));

            // Run the transmited symbols through the channel 
            cvec tdl_channel_output;
            tdl_channel.filter_known_channel(tx_ofdm_symbols, tdl_channel_output, tdl_channel_coeff);
            DEBUG_MSG("Channel output size: " << tdl_channel_output.size());
            DEBUG_MSG("Channel output: " << tdl_channel_output.left(10));
            tdl_channel_output = tdl_channel_output.left(batch_size * (Nsub + Ncp));

            rx_ofdm_symbols = awgn_channel(tdl_channel_output);
            DEBUG_MSG("Number of Rx OFDM symbols: " << rx_ofdm_symbols.length());
            DEBUG_MSG("Sample Rx OFDM symbols: " << rx_ofdm_symbols.left(DEBUG_SYM_COUNT));
        
            /*------------------------- RECEIVER --------------------------------*/
            // OFDM demodulation
            rx_symbols = ofdm.demodulate(rx_ofdm_symbols);
            DEBUG_MSG("Number of Rx symbols: " << rx_symbols.length());
            DEBUG_MSG("Sample Rx symbols: " << rx_symbols.left(DEBUG_SYM_COUNT));

            // Zero-Force equalization
            cvec hfft;
            // Channel coefficients at the beginning of the OFDM frame
            cmat relevant_channel_coeff = tdl_channel_coeff.get_rows(channel_rows); 
            for (int i = 0; i < rx_symbols.size() / Nsub; ++i)
            {
                // FFT of channel coefficient matrix
                hfft = concat(hfft, fft(relevant_channel_coeff.get_row(i), Nsub));
            }
            rx_symbols_equalized = elem_div(rx_symbols, hfft);
            
            // Demodulate the received symbols into bits
            rx_bits = modulator.demodulate_bits(rx_symbols_equalized);
            DEBUG_MSG("Number of Rx bits: " << rx_bits.length());
            DEBUG_MSG("Sample Rx bits: " << rx_bits.left(DEBUG_SYM_COUNT));
        
            // Calculate the bit error rate
            berc.clear();
            berc.count(tx_bits, rx_bits);
            err_count += berc.get_errors();
            total_bit_count += berc.get_total_bits();
            DEBUG_MSG("Error count: " << err_count << "\t\t" << "Total count: " << total_bit_count);
        }

        ber(snr_loop) = (double)err_count / total_bit_count;
        cout << "Eb/N0 (dB): " << EbN0dB[snr_loop] << "\t\t" << "BER: " << ber(snr_loop) << endl;
        DEBUG_MSG("------------------------------------------------------------------------------");
    }
    DEBUG_MSG("Elapsed time: " << timer.toc());

    // Print the results
    DEBUG_MSG("Eb/N0 (dB): " << EbN0dB);
    DEBUG_MSG("BER: " << ber);
    
    // Save the results to file
    cout << "Saving results to " << OUTPUT_FILENAME << endl;
    it_file fp;
    fp.open(OUTPUT_FILENAME);
    fp << Name("ebno_dB") << EbN0dB;
    fp << Name("ber") << ber;
    fp.close();

    return 0;
}   
