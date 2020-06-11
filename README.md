# multiple-scenario
Nested Schur decomposition method for parameter fitting with multiple experiments

This module takes a dict of pyomo models with the specified set of global parameters.
See the example for a detailed look at how to use this. The example requires the latest version of Kipet.


# Example

First, import the class:

    from NestedSchurDecomposition import NestedSchurDecomposition as NSD  

At the moment, an initial guess for the parameter values is needed. Each parameter considered needs to be provided:

    d_init_guess = {'k1' : (2.5, (0.0, 5.0)),
                    'k2' : (0.5, (0.0, 1.0)),
                   }

You also need to add the name of the parameters you are fitting in your pyomo models to the options dict:

    options = {'parameter_var_name': 'P'}
    
Next, create an instance of the NestedSchurDecomposition class with your models, an initial guess for the parameter values, and the options dict.
    
    nsd = NSD(models, d_init_guess, options)
   
To start the process, simply run the following:
    
    results = nsd.nested_schur_decomposition()
