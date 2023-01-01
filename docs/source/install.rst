###############
Install Femwell
###############

As femwell depends on the complex version of slepc4py, it needs to be installed first.
This can either be done using

.. code::

    conda/mamba install slepc4py=*=complex*
    
or using pip using

.. code::

    export PETSC_CONFIGURE_OPTIONS="--with-scalar-type=complex"
    pip install petsc
    pip install petsc4py
    pip install slepc
    pip install slepc4py

As only the mode solving depends on slepc4py, it is for now left out of the requirements to avoid a failing installation via pip.
But that's only a temporary solution, hopefully there will be an easier way to install slepc4py via pip.

After installing slepc4py, femwell can be installed via 

.. code::
    
    pip install femwell
