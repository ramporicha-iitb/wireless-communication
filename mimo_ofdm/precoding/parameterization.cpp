#include <cassert>
#include <itpp/itcomm.h>

using namespace std;
using namespace itpp;

void get_scalar_parameters(cmat V, Array<vec>& phases, Array<vec>& rotation_angles) 
{
    /*
     * V : t x n matrix (t >= n)
     */

    assert(phases.size() == V.cols());
    assert(rotation_angles.size() == V.cols()-1);

    int t,n;
    double theta;
    
    t = V.rows();
    n = V.cols();

    cvec D(t);  // vector of size t
    cmat G;     // Givens matrix of size t x t
   
    for (int k = 0; k < n; ++k)
    {
        D.ones();
        D.set_subvector(k, exp(1j * angle(V.get_col(k).right(t-k))));
        
        // Store phase-angle parameters
        phases(k) = angle(V.get_col(k).right(t-k)); 
        
        V = diag(conj(D)) * V;
        
        for (int l = 0; l < t-k-1; ++l)
        {
            G = itpp::eye_c(t);
            
            theta = atan2(real(V(t-l-1,k)), real(V(t-l-2,k)));

            // Store rotation-angle parameters 
            rotation_angles(k) = concat(rotation_angles(k), theta);
            
            G(t-l-1,t-l-1) = cos(theta);
            G(t-l-2,t-l-1) = -sin(theta);
            G(t-l-1,t-l-2) = sin(theta);
            G(t-l-2,t-l-2) = cos(theta);

            V = transpose(G) * V;
        }
    }
}
