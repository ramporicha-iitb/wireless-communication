// https://en.wikipedia.org/wiki/Givens_rotation

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

void get_scalar_parameters(cmat V, Array<vec>& phases, Array<vec>& rotation_angles) 
{
    /*
     * V : t x n matrix (t >= n)
     */

    int t,n;
    double theta;
    
    t = V.rows();
    n = V.cols();

    // Code donot handle non-square matrix parameterization as of now
    assert(t == n);

    cvec D(n);  // vector of size n
    cmat G;     // Givens matrix of size n x n
   
    phases.set_size(n);
    rotation_angles.set_size(n);

    DEBUG_MSG("V matrix: \n" << round_to_zero(V));

    for (int k = 0; k < n; ++k)
    {
        D.ones();
        D.set_subvector(k, exp(1j * angle(V.get_row(k).right(n-k))));
        
        // Store phase-angle parameters
        phases(k) = angle(V.get_row(k).right(n-k)); 
        
        V = V * diag(conj(D));

        vec rot_angles;
        for (int l = 0; l < n-k-1; ++l)
        {
            G = itpp::eye_c(n);
            
            theta = atan2(real(V(k,n-l-2)), real(V(k,n-l-1)));
    
            // Store rotation-angle parameters 
            rot_angles = concat(rot_angles, theta);

            /*
             * [sin -cos]
             * [cos  sin]
             */
            G(n-l-1,n-l-1) = sin(theta);
            G(n-l-2,n-l-1) = -cos(theta);
            G(n-l-1,n-l-2) = cos(theta);
            G(n-l-2,n-l-2) = sin(theta);

            V = V * G;

        }
        rotation_angles(k) = rot_angles;
        DEBUG_MSG("V matrix: \n" << round_to_zero(V));
    }
}

// Reconstruction of precoding matrix (V) from phase-angle and rotation-angle parametes
void construct_precoding_matrix(int t, int n,
        const Array<vec>& phases, const Array<vec>& rotation_angles,
        cmat& VQmat)
{
    double theta;
    int phase_index = 0;
    int angle_index = 0;

    cvec D(n);  // vector of size n
    cmat G;     // Givens matrix of size n x n
    
    VQmat = eye_c(n);
    DEBUG_MSG("V matrix: \n" << round_to_zero(VQmat));

    for (int k = 0; k < n; ++k)
    {
        D.ones();
        D.set_subvector(k, exp(1j * phases(k)));

        phase_index += (n-k);

        VQmat = diag(D) * VQmat;

        for (int l = 0; l < t-k-1; ++l)
        {
            theta = rotation_angles(k)(l);
            G = itpp::eye_c(n);

            G(n-l-1,n-l-1) = sin(theta);
            G(n-l-2,n-l-1) = -cos(theta);
            G(n-l-1,n-l-2) = cos(theta);
            G(n-l-2,n-l-2) = sin(theta);

            VQmat = transpose(G) * VQmat;
        }
        DEBUG_MSG("V matrix: \n" << round_to_zero(VQmat));
    }

    cmat eye_tilde = zeros_c(t,n);
    eye_tilde.set_submatrix(0,0,eye_c(n));
    VQmat = eye_tilde * VQmat;
    DEBUG_MSG("V matrix: \n" << round_to_zero(VQmat));
}
