#include <itpp/itcomm.h>

using namespace std;
using namespace itpp;

#define DEBUG 0
#if DEBUG
#define DEBUG_MSG(str) do { std::cout << str << std::endl; } while( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif

#define DEBUG_SYM_COUNT     4

#define CP_LEN              16
#define NUM_SUBCARRIER      128
#define NUM_OFDM_FRAMES     1e5
#define OUTPUT_FILENAME     "ofdm.it"

int main()
{
    RNG_randomize();
    
    // Modulators
    BPSK_c modulator;
    OFDM ofdm;
    ofdm.set_parameters(NUM_SUBCARRIER, CP_LEN);
    DEBUG_MSG ("Symbols: " << modulator.get_symbols());
    
    // Channels
    AWGN_Channel awgn_channel;

    double Es, Eb;
    Es = 1.0;                               // Energy per symbol
    Eb = Es / modulator.bits_per_symbol();
    
    vec EbN0dB, N0;
    EbN0dB = linspace(0.0, 15.0, 16);       // Eb/N0 values in dB
    N0 = Eb * pow(inv_dB(EbN0dB), -1.0);    // N0 : variance of the complex valued noise
    
    // Number of bits transmitted for each Eb/N0 valaue
    int num_tx_bits;
    num_tx_bits = NUM_OFDM_FRAMES * NUM_SUBCARRIER * modulator.bits_per_symbol();
   
    // Transmitted bits / symbols
    bvec tx_bits, rx_bits;
    cvec tx_symbols;
    cvec rx_symbols;
    cvec tx_ofdm_symbols;
    cvec rx_ofdm_symbols;

    // Bit error rates 
    vec ber;
    BERC berc;
    ber.set_size(EbN0dB.length());
    
    Real_Timer timer;
    timer.tic();
    
    for (int i = 0; i < EbN0dB.length(); i++) {

        /*-------------------------- TRANSMITTER -----------------------------*/
        // Generate a vector of random bits to transmit        
        tx_bits = randb(num_tx_bits);
        DEBUG_MSG("Number of Tx bits: " << num_tx_bits);
        DEBUG_MSG("Tx bits: " << tx_bits.right(DEBUG_SYM_COUNT * modulator.bits_per_symbol()));

        // Modulate the bits to symbols
        tx_symbols = modulator.modulate_bits(tx_bits);
        DEBUG_MSG("Number of Tx symbols: " << tx_symbols.length());
        DEBUG_MSG("Tx symbols: " << tx_symbols.right(DEBUG_SYM_COUNT));

        // OFDM modulation
        tx_ofdm_symbols = ofdm.modulate(tx_symbols);
        DEBUG_MSG("Number of Tx OFDM symbols: " << tx_ofdm_symbols.length());
        DEBUG_MSG("Tx OFDM symbols: " << tx_ofdm_symbols.right(DEBUG_SYM_COUNT));

        /*----------------------- AWGN CHANNEL -------------------------------*/
        // Set the noise variance of the AWGN channel
        awgn_channel.set_noise(N0(i));

        // Run the transmited symbols through the channel 
        rx_ofdm_symbols = awgn_channel(tx_ofdm_symbols);
        DEBUG_MSG("Number of Rx OFDM symbols: " << rx_ofdm_symbols.length());
        DEBUG_MSG("Rx OFDM symbols: " << rx_ofdm_symbols.right(DEBUG_SYM_COUNT));
        
        /*------------------------- RECEIVER --------------------------------*/
        // OFDM demodulation
        rx_symbols = ofdm.demodulate(rx_ofdm_symbols);
        DEBUG_MSG("Number of Rx symbols: " << rx_symbols.length());
        DEBUG_MSG("Rx symbols: " << rx_symbols.right(DEBUG_SYM_COUNT));

        // Demodulate the received symbols into bits
        rx_bits = modulator.demodulate_bits(rx_symbols);
        DEBUG_MSG("Number of Rx bits: " << rx_bits.length());
        DEBUG_MSG("Rx bits: " << rx_bits.right(DEBUG_SYM_COUNT));

        // Calculate the bit error rate
        berc.clear();
        berc.count(tx_bits, rx_bits);
        ber(i) = berc.get_errorrate();
        cout << "Eb/N0 (dB): " << EbN0dB[i] << "\t\t" << "BER: " << ber(i) << endl;
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
