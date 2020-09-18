#include <cassert>
#include <itpp/srccode/vq.h>
#include <itpp/itcomm.h>

using namespace std;
using namespace itpp;

using namespace std;
using namespace itpp;

#define DEBUG 0

#if DEBUG
#define DEBUG_MSG(str) do { std::cout << str << std::endl; } while( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif

#define CYCLIC_PREFIX_LEN   16
#define NUM_SUBCARRIER      128
#define NUM_OFDM_FRAMES     1e4

#define FRAME_FB_INTV       8
#define SUBCARRIER_FB_INTV  8
#define CODEBOOK_SIZE       16

#define NUM_TX_ANTENNA      2
#define NUM_RX_ANTENNA      2

#define OUTPUT_FILENAME     "scalar_parameter_quantization.it"

void get_scalar_parameters(cmat V, 
        Array<vec>& phases, 
        Array<vec>& rotation_angles);

void generate_scalar_parameter_codebook(
        int t, int n, int cb_size,
        vec& phase_cb, 
        Array<vec>& rotation_angle_cb);

void construct_precoding_matrix(
        int t, int n,
        const Array<vec>& phases, 
        const Array<vec>& rotation_angles,
        cmat& VQmat);

int main(int argc, char *argv[])
{
    RNG_randomize();

    int Nsub, Ncp, Nframe;

    Ncp = CYCLIC_PREFIX_LEN;
    Nsub = NUM_SUBCARRIER;
    Nframe = NUM_OFDM_FRAMES;

    // MIMO model
    int Nt, Nr;
    Nt = NUM_TX_ANTENNA;
    Nr = NUM_RX_ANTENNA;

    assert(Nr >= Nt);
    assert(Nframe >= FRAME_FB_INTV && Nframe % FRAME_FB_INTV == 0);
    assert(Nsub >= SUBCARRIER_FB_INTV && Nsub % SUBCARRIER_FB_INTV == 0);
    assert(ceil(log2(CODEBOOK_SIZE)) == floor(log2(CODEBOOK_SIZE)));        // Integer multiple

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

    // Bit error rates
    vec ber;
    BERC berc;
    ber.set_size(EbN0dB.length());

    // SVD matrices
    size_t fb_count = Nsub / SUBCARRIER_FB_INTV;
    vec S[fb_count];
    cmat U[fb_count], V[fb_count], Uh[fb_count], VQ[fb_count];;

    // Codebook
    vec phase_codebook;                     
    Array<vec> rot_angle_codebook(Nt);
    generate_scalar_parameter_codebook(Nr, Nt, CODEBOOK_SIZE, 
            phase_codebook, rot_angle_codebook);

    // Quantizer
    Scalar_Quantizer phase_quantizer;
    Array<Scalar_Quantizer> rot_angle_quantizer(Nt);

    phase_quantizer.set_levels(phase_codebook);
    for (int l = 0; l < Nt; ++l)
    {
        rot_angle_quantizer(l).set_levels(rot_angle_codebook(l));
    }
    DEBUG_MSG("Phase-angle codebook: \n" << phase_codebook);
    DEBUG_MSG("Rotation-angle codebooks: \n" << rot_angle_codebook);

    Real_Timer timer;
    timer.tic();

    for (int snr_loop = 0; snr_loop < EbN0dB.length(); ++snr_loop)
    {
        size_t err_count = 0;
        size_t total_bit_count = 0;

        TDL_Channel *tdl_channel[Nr][Nt]; // Nr x Nt independent channels
        for (int i = 0; i < Nr; ++i) {
            for (int j = 0; j < Nt; ++j) {
                tdl_channel[i][j] = new TDL_Channel(vec("0 -3 -6 -9"), ivec("0 1 2 3"));
                tdl_channel[i][j]->set_norm_doppler(1e-6);
            }
        }

        for (int frame_loop = 0; frame_loop < Nframe; ++frame_loop)
        {
            /*-------------------------- FEEDBACK ----------------------------*/
            if (frame_loop % FRAME_FB_INTV == 0)
            {
                // Generate channel filter coefficients
                // Returns a matrix with one tap per column for each independent channel

                // tdl_channel_coeff[i][j]:
                // Channel coefficient matrix between i-th rx antenna and j-th tx antenna
                cmat tdl_channel_coeff[Nr][Nt];
                for (int i = 0; i < Nr; ++i)
                {
                    for (int j = 0; j < Nt; ++j)
                    {
                        tdl_channel[i][j]->generate(Nsub + Ncp, tdl_channel_coeff[i][j]);
                    }
                }

                // FFT of channel coefficient vector at the begining OFDM frame
                cvec hfft[Nr][Nt];
                cvec relevant_channel_coeff[Nr][Nt];
                for (int i = 0; i < Nr; ++i)
                {
                    for (int j = 0; j < Nt; ++j)
                    {
                        relevant_channel_coeff[i][j] = tdl_channel_coeff[i][j].get_row(0);
                        hfft[i][j] = fft(relevant_channel_coeff[i][j], Nsub);
                    }
                }
          
                cmat h_mat(Nr,Nt);      //Channel matrix at each time index
                size_t fb_mat_count = 0;
                for (int sc = 0; sc < Nsub; ++sc)
                {
                    if (sc % SUBCARRIER_FB_INTV == 0)
                    {
                        for (int i = 0; i < Nr; ++i)
                        {
                            for (int j = 0; j < Nt; ++j)
                            {
                                h_mat(i,j) = hfft[i][j](sc);
                            }
                        }

                        // Channel matrix -> SVD -> Precoding matrix
                        svd(h_mat, U[fb_mat_count], S[fb_mat_count], V[fb_mat_count]);

                        // Paramterization: Get scalar parameters for given precoding matrix using Givens rotation
                        cmat& VQref = VQ[fb_mat_count];
                        const cmat& Vref = V[fb_mat_count];
                        Array<vec> phases(Vref.cols()), rot_angles(Vref.cols()-1);
                        get_scalar_parameters(Vref, phases, rot_angles);
                       
                        // Precoding matrix quantization (Receiver)
                        Array<vec> quantized_phases(Vref.cols()), quantized_rot_angles(Vref.cols()-1);
                        for (int k = 0; k < Vref.cols(); ++k)
                        {
                            // In practise we would use encode/decode function
                            quantized_phases(k) = phase_quantizer.Q(phases(k));

                            for (int l = 0; l < Nt-k-1; ++l)
                            {
                                quantized_rot_angles(k) = concat(quantized_rot_angles(k), rot_angle_quantizer(l).Q(rot_angles(k)(l)));
                            }
                        }
                        hermitian_transpose(U[fb_mat_count], Uh[fb_mat_count]);

                        DEBUG_MSG("Phases: \n" << phases);
                        DEBUG_MSG("Quantized phases: \n" << quantized_phases);
                        DEBUG_MSG("Rotation angles: \n" << rot_angles);
                        DEBUG_MSG("Quantized rotation angles: \n" << quantized_rot_angles);
                        
                        ++fb_mat_count;
                        
                        // Reconstruct precoding matrix from quantized scalar parameter feedback (Transmitter)
                        construct_precoding_matrix(Vref.rows(), Vref.cols(), 
                                quantized_phases, quantized_rot_angles, VQref);
                        DEBUG_MSG("Precoding matrix (V): \n" << Vref);
                        DEBUG_MSG("Quantized precoding matrix (VQ): \n" << VQref);
                    }
                }
                continue;
            }
 
            /*-------------------------- TRANSMITTER -----------------------------*/
            // Generate a vector of random bits to transmit

            // row[i]: data transmitted on i-th tx antenna
            int num_tx_bits = Nsub * modulator.bits_per_symbol();
            bmat tx_bits = randb(Nt, num_tx_bits);
            DEBUG_MSG("Tx bits size: " << tx_bits.rows() << "x" << tx_bits.cols());

            // Modulate the bits to symbols
            cmat tx_symbols(Nt, Nsub);
            for (int i = 0; i < Nt; ++i)
            {
                tx_symbols.set_row(i, modulator.modulate_bits(tx_bits.get_row(i)));
            }
            tx_symbols *= 1.0/sqrt(Nt);   // Normalized power
            DEBUG_MSG("Tx symbols size: " << tx_symbols.rows() << "x" << tx_symbols.cols());

            // Precoding
            // 0 to (SUBCARRIER_FB_INTV - 1) : index 0
            // SUBCARRIER_FB_INTV to (2*SUBCARRIER_FB_INTV - 1) : index 1 and so on.
            size_t fb_mat_index;
            for (int i = 0; i < Nsub; ++i)
            {
                fb_mat_index = i/SUBCARRIER_FB_INTV;
                tx_symbols.set_col(i, VQ[fb_mat_index] * tx_symbols.get_col(i));
            }

            // OFDM modulation
            cmat tx_ofdm_symbols(Nt, Nsub + Ncp);
            for (int i = 0; i < Nt; ++i)
            {
                fb_mat_index = i/SUBCARRIER_FB_INTV;
                tx_ofdm_symbols.set_row(i, ofdm.modulate(tx_symbols.get_row(i)));
            }
            DEBUG_MSG("Tx OFDM symbols size: " << tx_ofdm_symbols.rows() << "x" << tx_ofdm_symbols.cols());

            /*----------------------- TDL + AWGN CHANNEL -------------------------*/
            // Generate channel filter coefficients
            // Returns a matrix with one tap per column for each independent channel

            // tdl_channel_coeff[i][j]: 
            // Channel coefficient matrix between i-th rx antenna and j-th tx antenna
            cmat tdl_channel_coeff[Nr][Nt];
            for (int i = 0; i < Nr; ++i)
            {
                for (int j = 0; j < Nt; ++j)
                {
                    tdl_channel[i][j]->generate(Nsub + Ncp, tdl_channel_coeff[i][j]);
                }
            }

            // Run the transmited symbols through the channel
            cvec channel_output_vec;
            cmat tdl_channel_output(Nr, (Nsub + Ncp)); // Channel output at each Rx antenna
            tdl_channel_output.clear();

            for (int i = 0; i < Nr; ++i)
            {
                for (int j = 0; j < Nt; ++j)
                {
                    tdl_channel[i][j]->filter_known_channel(
                            tx_ofdm_symbols.get_row(j), channel_output_vec, tdl_channel_coeff[i][j]);

                    // Channel output at i-th Rx antenna = 
                    // Sum total of channel outputs due to each Tx antenna at i-th Rx antenna
                    tdl_channel_output.set_row(i, tdl_channel_output.get_row(i) +
                            channel_output_vec.left(Nsub + Ncp));
                }
            }

            // Set the noise variance of the AWGN channel
            awgn_channel.set_noise(N0(snr_loop));

            cmat rx_ofdm_symbols(Nr, (Nsub + Ncp));
            for (int i = 0; i < Nr; ++i)
            {
                rx_ofdm_symbols.set_row(i, awgn_channel(tdl_channel_output.get_row(i)));
            }

            /*------------------------- RECEIVER --------------------------------*/
            // OFDM demodulation
            cmat rx_symbols(Nr, Nsub);
            for (int i = 0; i < Nr; ++i)
            {
                rx_symbols.set_row(i, ofdm.demodulate(rx_ofdm_symbols.get_row(i)));
            }
            DEBUG_MSG("Rx symbols size: " << rx_symbols.rows() << "x" << rx_symbols.cols());

            // Postcoding and Zero-Force equalization
            cmat rx_symbols_equalized(Nr, Nsub);
            mat sigma_inv = zeros(Nt, Nr);

            for (size_t i = 0; i < Nsub; ++i)
            {
                fb_mat_index = i/SUBCARRIER_FB_INTV;
                for (size_t i = 0; i < Nt; ++i)
                {
                    sigma_inv(i,i) = 1/S[fb_mat_index](i);
                }
                // sigma_inv * Uh[fb_mat_index] * rx_symbols.get_col(i):
                // returns matrix with a single column
                rx_symbols_equalized.set_col(i,
                        cvectorize(sigma_inv * (Uh[fb_mat_index] * rx_symbols.get_col(i))));
            }

            bmat rx_bits;
            for (int i = 0; i < Nt; ++i)
            {
                rx_bits.append_row(modulator.demodulate_bits(rx_symbols_equalized.get_row(i)));
            }

            // Calculate the bit error rate
            berc.clear();
            berc.count(rvectorize(tx_bits), rvectorize(rx_bits));
            err_count += berc.get_errors();
            total_bit_count += berc.get_total_bits();
            DEBUG_MSG("Error count: " << err_count << "\t\t" << "Total count: " << total_bit_count);
        }
        
        ber(snr_loop) = (double)err_count / total_bit_count;
        cout << "Eb/N0 (dB): " << EbN0dB[snr_loop] << "\t\t" << "BER: " << ber(snr_loop) << endl;
        DEBUG_MSG("------------------------------------------------------------------------------");

        // Free dynamic memory
        for (int i = 0; i < Nr; ++i) {
            for (int j = 0; j < Nt; ++j) {
                delete tdl_channel[i][j];
            }
        }
    }
    cout << endl << "Elapsed time: " << timer.toc() << endl;

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
