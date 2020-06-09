"""
This module contains the nested Schur method for decomposing multiple
experiments into a bi-level optimization problem.

This version is adapted from the module used within Kipet - this version works
with all pyomo models and does not require Kipet.

Author: Kevin McBride 2020

"""
import copy
from pathlib import Path
import pickle
import shutil
from string import Template

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
from pyomo.core.base.PyomoModel import ConcreteModel

from pyomo.environ import (
    Constraint, 
    Objective,
    Param, 
    Set,
    SolverFactory,
    Suffix,
    )

from scipy.optimize import (
    Bounds,
    minimize,
    )

global opt_dict
opt_dict = dict()
global opt_count
opt_count = -1

global global_param_name
global global_constraint_name
global parameter_var_name
global global_set_name

# If for some reason your model has these attributes, you will have a problem
global_set_name = 'global_parameter_set'
global_param_name = 'd_params_nsd_globals'
global_constraint_name = 'fix_params_to_global_nsd_constraint'

# Header template
iteration_spacer = Template('\n' + '#'*30 + ' $iter ' + '#'*30 + '\n')

class NestedSchurDecomposition():
    
    """Nested Schur Decomposition approach to parameter fitting using multiple
    experiments
    
    """
    def __init__(self, models, d_info, kwargs=None):
        
        """Args:
            models (dict): A dict of pyomo Concrete models
            
            d_info: A dict of the global parameters including the iniital
                value and a tuple of the bounds, i.e. {'p1' : 2.0, (0.0, 4.0)}
        
            kwargs (dict): Optional arguments for the algorithm (incomplete)
        
        """
        # The models should be entered in as a dict (for now)
        self.models_dict = copy.copy(models)
        
        # The global parameter information is needed, especially the bounds
        self.d_info = d_info
        self.d_init = {k: v[0] for k, v in d_info.items()}
        
        # Arrange the kwargs
        self._kwargs = {} if kwargs is None else copy.copy(kwargs)
        
        # Options - inner problem optimization
        self.ncp = self._kwargs.pop('ncp', 3)
        self.nfe = self._kwargs.pop('nfe', 50)
        self.verbose = self._kwargs.pop('verbose', False)
        self.sens = self._kwargs.pop('use_k_aug', True)
        self.parameter_var_name = self._kwargs.pop('parameter_var_name', None)
        self.objective_name = self._kwargs.pop('objective_name', None)
        self.gtol = self._kwargs.pop('gtol', 1e-12)
        
        global parameter_var_name
        parameter_var_name = self.parameter_var_name
        
        if self.parameter_var_name is None:
            raise ValueError('NSD requires that the parameter attribute be provided')
        
        # Options - nested Schur decomposition
        self.remove_files = self._kwargs.pop('remove_files', True)
        self.d_init_user = self._kwargs.pop('d_init', True)
        self.method = self._kwargs.pop('method', 'trust-constr')
        
        # Run various assertions that the model is correctly structured
        self._test_models()
        
        # Add the global constraints to the model
        self._add_global_constraints()
        self._prep_models()
        
    def _test_models(self):
        """Sanity check on the input models"""
        
        for model in self.models_dict.values():
            
            # Check if the models are even models
            assert(isinstance(model, ConcreteModel) == True)
        
    def _add_global_constraints(self):
        """This adds the dummy constraints to the model forcing the local
        parameters to equal the current global parameter values
        
        """
        global global_param_name
        global global_constraint_name
        global global_set_name
        
        for model in self.models_dict.values():
            param_dict = {}
            for param in self.d_info.keys():
                if param in getattr(model, self.parameter_var_name):
                    param_dict.update({param: self.d_info[param][0]})

            setattr(model, global_set_name, Set(initialize=param_dict.keys()))

            setattr(model, global_param_name, Param(getattr(model, global_set_name),
                                  initialize=param_dict,
                                  mutable=True,
                                  ))
            
            def rule_fix_global_parameters(m, k):
                
                return getattr(m, parameter_var_name)[k] - getattr(m, global_param_name)[k] == 0
                
            setattr(model, global_constraint_name, 
            Constraint(getattr(model, global_set_name), rule=rule_fix_global_parameters))
        
    def _prep_models(self):
        """Prepare the model for NSD algorithm. Checks discretization,
        checks for an objective named "objective", and adds appropriate
        suffixes for sensitivity analysis.
        
        """
        for model in self.models_dict.values():        
            model.dual = Suffix(direction=Suffix.IMPORT_EXPORT)
            model.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
            model.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
            model.ipopt_zL_in = Suffix(direction=Suffix.EXPORT)
            model.ipopt_zU_in = Suffix(direction=Suffix.EXPORT)
        
        return None
            
    def _generate_bounds_object(self):
        """Creates the Bounds object needed by SciPy for minimization
        
        """
        lower_bounds = []
        upper_bounds = []
        
        for k, v in self.d_info.items():
            lower_bounds.append(v[1][0])
            upper_bounds.append(v[1][1])
        
        return Bounds(lower_bounds, upper_bounds, True)
        
    def nested_schur_decomposition(self):
        """Here is where the magic happens. This is the outer problem controlled
        by a trust region solver running on scipy. This is the only function that
        the user needs to call to run this thing.
        
        Returns:
            res (scipy.optimize.optimize.OptimizeResult): The results from the 
                trust region optimation (outer problem)
                
        """    
        print(iteration_spacer.substitute(iter='NSD Start'))
    
        d_init = self.d_init #self._generate_initial_d()
        d_bounds = self._generate_bounds_object()
    
        self.d_iter = list()
        def callback(x, *args):
            self.d_iter.append(x)
    
    
        if self.method in ['trust-exact', 'trust-constr']:
        # The args for scipy.optimize.minimize
            fun = _inner_problem
            x0 = list(d_init.values()) #list(d_init.values()) if isinstance(d_init, dict) else d_init
            args = (self.models_dict,)
            jac = _calculate_m
            hess = _calculate_M
            
            callback(x0)
            results = minimize(fun, x0, args=args, method=self.method,
                           jac=jac,
                           hess=hess,
                           callback=callback,
                           bounds=d_bounds,
                           options=dict(gtol=self.gtol,
                                      #  initial_tr_radius=0.1,
                                      #  max_tr_radius=0.1
                                        ),
                           )
            self.parameters_opt = {k: results.x[i] for i, k in enumerate(self.d_init.keys())}
            
            
        if self.method in ['newton']:
            x0 = list(d_init.values())
            results = self._run_newton_step(x0, self.models_dict)
            self.parameters_opt = {k: results[i] for i, k in enumerate(self.d_init.keys())}
        
        return results, opt_dict
    
    def _run_newton_step(self, d_init, models):
        """This runs a basic Newton step algorithm - use a decent alpha!
        
        UNDER CONSTRUCTION!
        """
        tol = 1e-6
        alpha = 0.4
        max_iter = 40
        counter = 0
        self.d_iter.append(d_init)
        
        while True:   
        
            _inner_problem(d_init, models, generate_gradients=False)
            M = opt_dict[opt_count]['M']
            m = opt_dict[opt_count]['m']
            d_step = np.linalg.inv(M).dot(-m)
            d_init = [d_init[i] + alpha*d_step[i] for i, v in enumerate(d_init)]
            self.d_iter.append(d_init)
            
            if max(d_step) <= tol:
                
                print(f'Terminating sequence: minimum tolerance in step size reached ({tol}).')
                break
            
            if counter == max_iter:
                print(f'Terminating sequence: maximum number of iterations reached ({max_iter})')
                break
            
            counter += 1
            
        return d_init
                 
def _optimize(model, d_vals, verbose=False):
    """Solves the optimization problem with optional k_aug sensitivity
    calculations (needed for Nested Schur Decomposition)
    
    Args:
        model (pyomo ConcreteModel): The current model used in parameter
            fitting
            
        d_vals (dict): The dict of global parameter values
        
        verbose (bool): Defaults to false, option to see solver output
        
    Returns:
        model (pyomo ConcreteModel): The model after optimization
        
    """
    global global_param_name
    global parameter_var_name
    
    delta = 1e-12
    ipopt = SolverFactory('ipopt')
    tmpfile_i = "ipopt_output"

    for k, v in getattr(model, parameter_var_name).items():
        getattr(model, parameter_var_name)[k].unfix()

    for k, v in getattr(model, global_param_name).items():
        getattr(model, global_param_name)[k] = d_vals[k]
    
    results = ipopt.solve(model,
                          symbolic_solver_labels=True,
                          keepfiles=True, 
                          tee=verbose, 
                          logfile=tmpfile_i)
    
    return  model

def _inner_problem(d_init_list, models, generate_gradients=False, initialize=False):
    """Calculates the inner problem using the scenario info and the global
    parameters d
    
    Args:
        d_init_last (list): list of the parameter values
        
        models (dict): the dict of pyomo models used as supplemental args
        
        generate_gradients (bool): If the d values do not line up with the 
            current iteration (if M or m is calculated before the inner prob),
            then the inner problem is solved to generate the corrent info
            
        initialize (bool): Option only used in the initial optimization before 
            starting the NSD routine
    
    Returns:
        Either returns None (generating gradients) or the scalar value of the 
        sum of objective function values from the inner problems
        
    """    
    global opt_count
    global opt_dict
    global global_constraint_name
    global parameter_var_name
    global global_param_name
     
    opt_count += 1   
        
    options = {'verbose' : False}
    _models = copy.copy(models) 
  
    m = 0
    M = 0
    Si = []
    Ki = []
    Ei = []
    
    objective_value = 0
    
    if opt_count == 0:
        print(iteration_spacer.substitute(iter=f'Inner Problem Initialization'))
        print(f'Initial parameter set: {d_init_list}')
    else:
        print(iteration_spacer.substitute(iter=f'Inner Problem {opt_count}'))
        print(f'Current parameter set: {d_init_list}')
    
    for k, model in _models.items():
        
        valid_parameters = dict(getattr(model, parameter_var_name).items()).keys()
        
        if isinstance(d_init_list, dict):
            d_init = {k: d_init_list[k] for k in valid_parameters}
        else:
            d_init = {param: d_init_list[i] for i, param in enumerate(valid_parameters)}
        
        # Optimize the inner problem
        model_opt = _optimize(model, d_init)
        
        kkt_df, var_ind, con_ind_new = _get_kkt_matrix(model_opt)
        duals = [model_opt.dual[getattr(model_opt, global_constraint_name)[key]] for key, val in getattr(model_opt, global_param_name).items()]
        col_ind  = [var_ind.loc[var_ind[0] == f'{parameter_var_name}[{v}]'].index[0] for v in valid_parameters]
        dummy_constraints = _get_dummy_constraints(model_opt)
        dc = [d for d in dummy_constraints]
        
        # Perform the calculations to get M and m
        K = kkt_df.drop(index=dc, columns=dc)
        E = np.zeros((len(dummy_constraints), K.shape[1]))
        K_i_inv = np.linalg.inv(K.values)
         
        for i, indx in enumerate(col_ind):
            E[i, indx] = 1
            
        S = E.dot(K_i_inv).dot(E.T)
        M += np.linalg.inv(S)
        m += np.array(duals)
        objective_value += model_opt.objective.expr()

        Si.append(S)
        Ki.append(K_i_inv)
        Ei.append(E)

    # Save the results in opt_dict - needed for further iterations
    opt_dict[opt_count] = { 'd': d_init_list,
                            'obj': objective_value,
                            'M': M,
                            'm': m,
                            'S': Si,
                            'K_inv': Ki,
                            'E': Ei,
                            } 
    
    if not generate_gradients:
        return objective_value
    else:
        return None

def _get_kkt_matrix(model):
    
    nlp = PyomoNLP(model)
    varList = nlp.get_pyomo_variables()
    conList = nlp.get_pyomo_constraints()
    duals = nlp.get_duals()
    
    J = nlp.extract_submatrix_jacobian(pyomo_variables=varList, pyomo_constraints=conList)
    H = nlp.extract_submatrix_hessian_lag(pyomo_variables_rows=varList, pyomo_variables_cols=varList)
    
    var_index_names = [v.name for v in varList]
    con_index_names = [v.name for v in conList]

    J_df = pd.DataFrame(J.todense(), columns=var_index_names, index=con_index_names)
    H_df = pd.DataFrame(H.todense(), columns=var_index_names, index=var_index_names)
    
    var_index_names = pd.DataFrame(var_index_names)
    
    KKT_up = pd.merge(H_df, J_df.transpose(), left_index=True, right_index=True)
    KKT = pd.concat((KKT_up, J_df))
    KKT = KKT.fillna(0)
    
    return KKT, var_index_names, con_index_names

def _get_dummy_constraints(model):
    """Get the locations of the contraints for the local and global parameters
    
    Args:
        models (pyomo ConcreteModel): The current model in the inner problem
        
    Returns:
        dummy_constraints (str): the names of the dummy contraints for the 
            parameters
    
    """
    global global_constraint_name
    global parameter_var_name
    global global_param_name
    
    dummy_constraint_name = global_constraint_name
    dummy_constraint_template = Template(f'{dummy_constraint_name}[$param]')
    parameters = getattr(model, global_param_name).keys()
    dummy_constraints = [dummy_constraint_template.substitute(param=k) for k in parameters]
    
    return dummy_constraints

def _calculate_M(x, scenarios):
    """Calculates the Hessian, M
    This is scipy.optimize.minimize conform
    Checks that the correct data is retrieved
    Needs the global dict to get information from the inner optimization
    
    Args:
        x (list): current parameter values
        
        scenarios (dict): The dict of scenario models
        
    Returns:
        M (np.ndarray): The M matrix from the NSD method
    
    """
    global opt_dict
    global opt_count
    
    if opt_count == 0 or any(opt_dict[opt_count]['d'] != x):
        _inner_problem(x, scenarios, generate_gradients=True)

    M = opt_dict[opt_count]['M']
    return M
    
def _calculate_m(x, scenarios):
    """Calculates the jacobian, m
    This is scipy.optimize.minimize conform
    Checks that the correct data is retrieved
    Needs the global dict to get information from the inner optimization
    
    Args:
        x (list): current parameter values
        
        scenarios (dict): The dict of scenario models
        
    Returns:
        m (np.ndarray): The m matrix from the NSD method
    
    """
    global opt_dict
    global opt_count
    
    if opt_count == 0 or any(opt_dict[opt_count]['d'] != x):
        _inner_problem(x, scenarios, generate_gradients=True)
    
    m = opt_dict[opt_count]['m']
    return m