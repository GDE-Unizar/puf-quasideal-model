"""
This module contains the functions used as console-scripts.
"""
from multiprocessing import Pool as _Pool
from time import time as _time
from argparse import ArgumentParser as _ArgumentParser
from datetime import datetime as _datetime
import puf_qideal_model as pqm


def run_in_parallel(func, args):
    """
    This function takes a 'func' method and a list of arguments 'args', and iterates the execution of 'func' in parallel. The method 'func' might receive any number and types of arguments, but these will be passed in order as they appear in 'args'.

    Parameters
    ----------
    func : method
        Method that can take any number and types of position-based arguments.

    args : list
        List that contains the arguments to pass to 'func'.
        
    Returns
    -------
    list
        List with the results of each process.
    """
    with _Pool() as pool:
        return pool.starmap(func, args)

    return results
    
                                        
def master(ninst=10, nrep=10, nbits=100, pintra=0.01, pinter=0.5, nsim=10000, ncores=1,verbose=False):
    """
    Helper method to parallelize the execution of 'pqm.montecarlo_quasideal_Dks_distribution' function.
    """
    ## Calculamos el número de puntos por proceso
    nsim_per_proc = nsim//ncores
    if verbose and nsim%ncores!=0:
        print(f"El número de simulaciones no puede repartirse equitativamente entre núcleos, lo cual es una condición necesaria para la ejecución del programa. Por ello, se realizarán en su lugar {nsim_per_proc*ncores} simulaciones")
        print()
        
    ## Generamos la matriz de argumentos para la función 'run_in_parallel'
    args = [[nsim_per_proc,ninst,nrep,nbits,pintra,pinter,verbose] for i in range(ncores)]
    
    ## Ejecutamos la simulación
    start = _time()
    montecarlo_batch = run_in_parallel(func=pqm.montecarlo_quasideal_Dks_distribution, args=args)
    end = _time()
    intra=[dks_value for montecarlo_set in montecarlo_batch for dks_value in montecarlo_set.intra]
    inter=[dks_value for montecarlo_set in montecarlo_batch for dks_value in montecarlo_set.inter]
       
    datalog = f"""execdate: {_datetime.now()}
exectime: {end - start:.3f} s
ninst: {ninst}
nrep: {nrep}
nbits: {nbits}
pintra: {pintra}
pinter: {pinter}
nsim: {nsim} 
ncores: {ncores},
nsim_per_proc: {nsim_per_proc}
              """
    return pqm.DksDistribution(ninst=ninst,nrep=nrep,nbits=nbits,pintra=pintra,pinter=pinter,intra=intra,inter=inter),datalog
    
    
def pqm_simulation():
    """
    Console-script to simulate a quasi-ideal PUF.
    """
    # Create the parser
    parser = _ArgumentParser(description='Quasideal PUF intra/inter-hamming distance simulation.')

    # Optional argument with default value
    parser.add_argument('-i', '--ninst', type=int, default=10, help='Number of PUF instances in each simulated experiment (default: 10)')
    parser.add_argument('-r', '--nrep', type=int, default=10, help='Number of repetitions in each simulated experiment (default: 10)')
    parser.add_argument('-b', '--nbits', type=int, default=100, help='Number of bits in the PUFs responses (default: 100)')
    parser.add_argument('-p', '--pintra', type=float, default=0.01, help='average intra-distance as measured (default: 0.01)')
    parser.add_argument('-u', '--pinter', type=float, default=0.5, help='average inter-distance as measured (default: 0.5)')
    parser.add_argument('-o', '--out', type=str, default='', help='The name of the output file. Note that regardless of the name chosen, the Hamming distance distribution file is always preceded by "hamming". By default, each file is named after the parameters selected.')

    # Parse the arguments
    args = parser.parse_args()

    # Main logic of your script
    try:           
        if args.out=='':
            hamming_name=f'-i{args.ninst}-r{args.nrep}-b{args.nbits}-p{args.pintra:.3f}-u{args.pinter:.3f}'
        else:
            hamming_name=args.out
        pqm.sim_quasideal_hamming_distribution(ninst=args.ninst,nrep=args.nrep,nbits=args.nbits,pintra=args.pintra,pinter=args.pinter,return_data_only=False).save(f'hamming{hamming_name}.txt')
        
    except Exception as e:
        print(f'An error occurred: {e}')
        

def pqm_dks_distribution():
    """
    Console-script to simulate a quasi-ideal PUF.
    """
    # Create the parser
    parser = _ArgumentParser(description='Kolmogorov-Smirnov distribution for a quasideal PUF model.')

    # Optional argument with default value
    parser.add_argument('-i', '--ninst', type=int, default=10, help='Number of PUF instances in each simulated experiment (default: 10)')
    parser.add_argument('-r', '--nrep', type=int, default=10, help='Number of repetitions in each simulated experiment (default: 10)')
    parser.add_argument('-b', '--nbits', type=int, default=100, help='Number of bits in the PUFs responses (default: 100)')
    parser.add_argument('-p', '--pintra', type=float, default=0.01, help='average intra-distance as measured (default: 0.01)')
    parser.add_argument('-u', '--pinter', type=float, default=0.5, help='average inter-distance as measured (default: 0.5)')
    parser.add_argument('-s', '--nsim', type=int, default=10000, help='Number of experiments simulated in order to make Dks distribution statistics (default: 10000)')
    parser.add_argument('-c', '--ncores', type=int, default=1, help='Number of cores used for parallelization (default: 1)')
    parser.add_argument('-o', '--out', type=str, default='', help='The name of the output file. Note that regardless of the name chosen, the Dks distribution file is always preceded by "Dks". By default, each file is named afeter the parameters selected.')
    
    # Optional flag (store_true implies it's a flag and stores True if present)
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose mode')

    # Parse the arguments
    args = parser.parse_args()

    # Main logic of your script
    if args.verbose:
        print('Verbose mode is on. Evolution and logs will be shown here.')
        print()

    try:
        Dks_distribution,datalog = master(args.ninst,args.nrep,args.nbits,args.pintra,args.pinter,args.nsim,args.ncores,args.verbose)
        
        if args.verbose:
            print(datalog)
            
        if args.out=='': 
            Dks_name=f'-i{args.ninst}-r{args.nrep}-b{args.nbits}-p{args.pintra:.3f}-u{args.pinter:.3f}'
        else:
            Dks_name=args.out
        Dks_distribution.save(f'Dks{Dks_name}.txt')

    except Exception as e:
        print(f'An error occurred: {e}')
          