# set some minimizer specific options # 
# Useable: L-BFGS-B(fast), Powell(slow), Nelder-Mead(veeery slow) 
# other minimizers show bad performance, were not optimized and are highly experimental ... 
known_minimizers = ("L-BFGS-B", "Nelder-Mead", "SLSQP", "TNC", "Powell") 
minimizer_options = {} 
minimizer_options.update({ "L-BFGS-B":
                             {   'method'   : "L-BFGS-B", 
                                 'args'     : (), 
                                 'jac'      : None, 
                                 'tol'      : None, 
                                 'callback' : None, 
                                 'options'  : 
                                     {'disp': None, 
                                      'iprint': -1, 
                                      'gtol': 1e-20, 
                                      'eps': 5e-08, #'eps':  5e-08, 
                                      'maxiter': 150000, 
                                      'ftol': 1e-12, 
                                      'maxcor': 100, # ftol: 1e-12 
                                      'maxfun': 150000} 
                             }   
                          }) 
minimizer_options.update({ 'Nelder-Mead':   
                             {   'method': 'Nelder-Mead',  
                                 'tol' : None, 
                                 'callback' : None, 
                                 'options': 
                                     {'disp': True, 
                                      'maxiter': None,
                                      'return_all': False, 
                                      'maxfev': None, 
                                      'xtol': 0.0001, 
                                      'ftol': 0.001} 
                             }   
                          }) 
minimizer_options.update({ 'SLSQP':
                             {   'method': 'SLSQP', 
                                 'args': (), 
                                 'tol' : None, 
                                 'callback' : None, 
                                 'jac':None, 
                                 'constraints':(), 
                                 'options': 
                                     {'disp': False, 
                                      'iprint': 1, 
                                      'eps': 1e-6, #'eps': 1.4901161193847656e-08, 
                                      'maxiter': 100, 
                                      'ftol': 1e-06} 
                             }   
                          }) 
minimizer_options.update({ 'TNC':
                             {   'method': 'TNC', 
                                 'args': (), 
                                 'tol' : None, 
                                 'callback' : None, 
                                 'jac':None, 
                                 'options': 
                                     {'disp': False, 
                                      'minfev': 0, 
                                      'scale': None, 
                                      'rescale': -1, 
                                      'offset': None, 
                                      'gtol': -1, 
                                      'eps': 1e-8, 
                                      'eta': -1, 
                                      'maxiter': None, 
                                      'maxCGit': -1, 
                                      'mesg_num': None, #'eps': 1e-08,  
                                      'ftol': -1, 
                                      'xtol': -1, 
                                      'stepmx': 0, 
                                      'accuracy': 0} 
                             }   
                          })  
minimizer_options.update({ 'Powell':        
                             {   'method': 'Powell', 
                                 'args': (), 
                                 'callback' : None, 
                                 'jac':None, 
                                 'options': 
                                     {'maxiter': None, 
                                      'maxfun': None, 
                                      'ftol' : 1e-8, 
                                      'xtol' : 1e-9, 
                                      'full_output': 0,  
                                      'disp': True, 
                                      'retall': 0 } 
                             }   
                          }) # works well, but slow:  'ftol' : 1e-8, 'xtol' : 1e-9, 
