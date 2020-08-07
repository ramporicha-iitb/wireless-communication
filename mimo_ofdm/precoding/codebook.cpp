#include <cassert>
#include <itpp/itbase.h>
#include <itpp/srccode/vqtrain.h>
#include <itpp/srccode/vq.h>

using namespace std;
using namespace itpp;

#define NUM_TRAIN_SAMPLES     1e6

/*
 * size of channel matrix (H) is Nr x Nt
 * size of precoding vector (V) is Nt x Nt
 * Nt column orthonormal vectors
 */

void generate_vector_codebook(
        int Nr, 
        int Nt, 
        int cb_size, 
        Array<mat>& codebook) 
{
    assert(Nt == codebook.size());

    vec s;
    cmat H, U, V;

    Array< Array<vec> > training_data(Nt);

    for (int i = 0; i < Nt; ++i)
    {
        training_data(i).set_size(NUM_TRAIN_SAMPLES);
    }

    for (int i = 0; i < NUM_TRAIN_SAMPLES; ++i)
    {
        H = randn_c(Nr,Nt);
        svd(H, U, s, V);

        for (int j = 0; j < V.cols(); ++j)
        {
            // j-th column orthonormal vector
            Array<vec>& tmp = training_data(j);

            // i-th training data of j-th column orthonormal vector
            tmp(i) = concat(real(V.get_col(j)), imag(V.get_col(j))); 
        }
    }

    // Train using LBG algorithm
    // codebook(i): codebook for i-th column orthonormal vector
    for (int i = 0; i < Nt; ++i)
    {
        codebook(i) = lbg(training_data(i), cb_size, 9999, false);
    }
}
