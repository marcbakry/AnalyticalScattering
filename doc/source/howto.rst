
.. _label-howto:

How-to
======

The analytical computations assumes that the wave is propagating along the :math:`-\vec e_z` axis.

Using ``AnalyticalScattering`` is easy. First, include the header file and the ``matrix`` header from the **castor** project.

.. code:: c++

    #include "castor/matrix.hpp"
    #include "analyticalscattering.hpp"

**Remark:** ``AnalyticalScattering`` does not belong to the ``castor`` namespace.

Then, set the wavenumber, the radius of the sphere, and some observation point.

.. code:: c++

    double k = 5.;  // wave number
    double R = 1.;  // sphere radius

    matrix<double> Xobs({-10., 0., 0.});    // outside the sphere

Now assuming a Neumann boundary condition, we ask for a *domain* computation (``Xobs`` is not on the sphere so it is not a *boundary* computation).

.. code:: c++

    AnalyticalScattering<double> as;    // create the object
    
    auto Uobs = as.sphereHelmholtz(Xobs, R, k, neumann, domain);
    disp(Uobs);     // display the value of the scattered field at Xobs


.. code:: text

    Matrix 1x1 of type 'St7complexIdE' (16 B):
       (0.01207,-0.03683)