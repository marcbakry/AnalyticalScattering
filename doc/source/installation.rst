.. _label-install:

Installation
++++++++++++

Clone the following repository

.. code:: text

    git clone https://github.com/marcbakry/AnalyticalScattering.git

The header ``AnalyticalScattering.hpp`` can be found in the ``include/`` subfolder.


.. _label-requirements:

Requirements
++++++++++++

``AnalyticalScattering`` may, or may not use an external library. If you use a **fully** ``c++-17`` **-compliant compiler** (more precisely, the ``<cmath>`` Bessel-function-related part of the ``c++-17`` standard), then you should define ``CXX17_COMPLIANT`` (by passing the ``-DCXX17_COMPLIANT`` option to the compiler see below). **For your information, the** ``g++`` **compiler (which is an alias for** ``clang`` **) from the developer tools in MacOS is NOT compliant**. The typical compiler command is (assuming the compiler is ``g++`` and the main is ``overview.cpp``)

.. code:: text

    g++ -std=c++17 -I/path/to/analyticalscattering_hpp -I/path/to/matrix_hpp -DCXX17_COMPLIANT overview.cpp -o testAnayticalScattering

**Remark:** if you are using ``cmake``, simply add ``add_definitions(-DCXX17_COMPLIANT)`` in the ``CMakeLists.txt`` file.

In the other case, one must link with the ``GNU Scientific Library`` (see `GSL <https://www.gnu.org/software/gsl/>`_ ) which is a C-library for scientific computing. The compilation and installation are straightforward as it does not hae dependencies. Once the library has been installed the command line above becomes

.. code:: text

    g++ -std=c++14 -I/path/to/analyticalscattering_hpp -I/path/to/matrix_hpp overview.cpp -lgsl -lgslcblas -o testAnayticalScattering

Note that if for one reason or the other one already links to a ``cblas`` library, one can omit ``-lgslcblas``.