
.. _label-howto:

How-to
======

The analytical computations assumes that the wave is propagating along the :math:`-\vec e_z` axis.

Using ``AnalyticalScattering`` is easy. First, include the header file and the ``matrix`` header from the `castor <http://leprojetcastor.gitlab.labos.polytechnique.fr/castor>`_ project.

.. code:: c++

    #include "castor/matrix.hpp"
    #include "analyticalscattering.hpp"

**Remark:** ``AnalyticalScattering``  *belongs* to the ``castor::`` namespace so one can use the following in the *preamble*

.. code:: c++

    using namespace castor;

Then, set the wavenumber, the radius of the sphere, and some observation point.

.. code:: c++

    double k = 5.;  // wave number
    double R = 1.;  // sphere radius

    matrix<double> Xobs({-10., 0., 0.});    // outside the sphere

Now assuming a Neumann boundary condition, we ask for a *domain* computation (``Xobs`` is not on the sphere so it is not a *boundary* computation).

.. code:: c++

    AnalyticalScattering<double> as;    // create the object
    
    auto Uobs = as.sphereHelmholtz(Xobs, R, k, neumann, domain);
    disp(Uobs,2);     // display the value of the scattered field at Xobs


.. code:: text

    Matrix 1x1 of type 'St7complexIdE' (16 B):
       (0.01207,-0.03683)

In order to compile the code, it is possible to use one of the command given in the :ref:`label-requirements` section.

The corresponding file ``demo_as.cpp`` can be found in the ``demo/`` subfolder of this project.
