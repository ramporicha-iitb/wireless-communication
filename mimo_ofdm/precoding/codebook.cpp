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

void generate_scalar_parameter_codebook (
        int t, int n, int cb_size,
        vec& phase_cb, Array<vec>& rotation_angle_cb)
{
    assert(rotation_angle_cb.size() == t);

    /* Precoding matrix V is of size t x n, t >= n
     *
     * Phases Phi(k,j) is uniformly disributed over (-pi, pi] for all k and j
     *
     * PDF of rotation angle theta(k,l) - 2l * sin(theta(k,l)) ^ (2l-1) * cos(theta(k,l)),
     * 0 <= theta(k,l) < pi/2, for 1 <= k <= n, 1 <= l <= t-k
     */

    // Train - Phase angle
    // Uniform random numbers between -pi and pi
    vec rand_phases = 2*pi*randu(NUM_TRAIN_SAMPLES) - pi;
    phase_cb = sqtrain(rand_phases, cb_size);

    // Train - Rotation angles
    vec rand_theta, rand_rotation_angles;

    for (int l=1; l <= t; ++l)
    {
        rand_theta = pi/2*randu(NUM_TRAIN_SAMPLES);

        // random numbers for rotation angles
        rand_rotation_angles = 2*l*elem_mult(pow(sin(rand_theta), 2*l - 1), cos(rand_theta));

        // rotation_angle_cb(l): codebook for theta(k,l)
        rotation_angle_cb(l-1) = sqtrain(rand_rotation_angles, cb_size);
    }
}
