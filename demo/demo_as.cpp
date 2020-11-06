#include "castor/matrix.hpp"
#include "AnalyticalScattering.hpp"

using namespace castor;

int main()
{
    double k = 5; // wave number
    double R = 1; // sphere radius

    matrix<double> Xobs({-10.,0.,0.}); // outside the sphere

    AnalyticalScattering<double> as;

    auto Uobs = as.sphereHelmholtz(Xobs, R, k, neumann, domain);
    disp(Uobs,2);
}