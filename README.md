# multiple-scenario
Nested Schur decomposition method for parameter fitting with multiple experiments modeled using Pyomo

# NOTE: Update coming soon!

## Usage
To use the module in your own code, you first need to import the module (here as NSD, for brevity).

    from NestedSchurDecomposition import NestedSchurDecomposition as NSD

When creating an instance of the NSD class, there are several arguments required.
In the description of the problem, it was mentioned that the user needs to supply some information about the parameters in addition to the Pyomo models. In the declaration of the NSD instance, a dictionary of the models (models), a dictionary of the parameter information (d_init_guess), and the name of the parameters (parameter_var_name) within the Pyomo models need to be included:

    nsd = NSD(models, d_init_guess, parameter_var_name, options)

The dictionary of models is pretty straightforward; the keys are not used in the model and should be used to designate the model (such as a name or experiment number). In the example problem included, this is a dictionary of the form:

    models = {0: model1, 1: model2, 2: model3, 3 :model4}

Each model should contain a subset of the global parameters that are to be fit using the NSD algorithm.
At this moment, an initial guess and bounds for the global parameters are required. This is also passed as a dictionary where each parameter is a key with a tuple that contains both the initial guess and a tuple of the lower and upper bounds.


In the example problem, the global parameter set is {k1, k2}. Three of the models contain the full set of global parameters while the fourth only uses k1. This does not present an issue with the NSD algorithm. The first parameter k1 has an initial value of 2.5 and is within [0.0, 5.0]. The second parameter k2 has an initial value of 0.5 and is within [0.0, 1.0]. The second argument to the NSD instance should then be a dictionary of the form (using our example):

    d_init_guess = {'k1' : (2.5, (0.0, 5.0)),
                    'k2' : (0.5, (0.0, 1.0)),
                   }

The last necessary parameter to include for the NSD is the name of the global parameter set in the models, which should (necessarily!) be the same. For example Kipet uses the designation 'P' for the kinetic parameters (and those that aren't) to be fit. Since your own models may not be using Kipet, this is a required argument in the declaration of the NSD instance.

Finally, there is an options dictionary that can be passed as well. The reader should take a look at the code to see which options are currently included. These have to do mainly with tolerances, the method used (such as trust-constr), etc.

Once the models, the initial parameter values, and the parameter set name are provided, the NSD method can be invoked using the following:

    results = nsd.nested_schur_decomposition()

where results is a SciPy Results object. The final parameter values can be taken from either the results output or through the nsd.parameter_opt attribute:
    
    >>> nsd.parameters_opt
    >>> {'k1': 2.499957616139173, 'k2': 0.8000022861374527}
