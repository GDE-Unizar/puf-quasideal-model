"""
This package provides some useful classes and functions for working with a quasi-ideal PUF model, which allows hypothesis testing with experimental data to check whether the reliability and uniqueness distributions fit the binomial distribution well. This check in turn allows to reliably calculate the identification error rates for the proposed PUF in question.
"""


from numpy import sqrt as _sqrt,\
                  histogram as _histogram,\
                  linspace as _linspace,\
                  array as _array
                  
from numpy.random import choice as _choice

from matplotlib.pyplot import hist as _hist,\
                              plot as _plot,\
                              legend as _legend,\
                              show as _show
                  
from scipy.stats import fit as _fit,\
                        binom as _binom,\
                        johnsonsb as _johnsonsb
                        
def _read_dictionary_from_file(file_path):
    """
    Auxiliary method to extract a dictionary from a TXT file with possibly commented lines.
    """
    dictionary = {}

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            if ':' in line:
                key, value = line.split(':', 1)
                dictionary[key.strip()] = value.strip()

    return dictionary


## Module's classses
class HammingDistribution:
    """
    This class contains the intra / inter Hamming distance distributions for a PUF experiment simulated on the basis of a PUF quaideal model.
    
    Parameters
    ----------
    ninst : int, optional, 'None' by default
        Number of instances simulated

    nrep : int, optional, 'None' by default
        Number of repetitions

    nbits : int, optional, 'None' by default
        Number of PUF's responses bits

    intra : list of floats, optional, 'None' by default
        Dataset of intra-distance values

    inter : list of floats, optional, 'None' by default
        Dataset of inter-distance values    
    
    Methods
    -------
    load(filename)
        This function allows you to load a Hamming distribution in text format. Caution: this function overrides any parameters defined when creating the class object.
        
    save(filename)
        Saves the Hamming distribution object to a file.
        
    plot(intra=True, intra_fit=True, inter=True, inter_fit=True)
        Plots the Hamming distribution. If intra or intra_fit is True, it plots the intra-distribution and its fit. If inter or inter_fit is True, plots the inter-distribution and its fit.
    """
    def __init__(self, ninst=None, nrep=None, nbits=None, intra=None, inter=None):
        """
        Constructor method.
        """
        self.ninst=ninst
        """int: Number of instances simulated"""
        self.nrep=nrep
        """int: Number of repetitions"""
        self.nbits=nbits
        """int: Number of PUF's responses bits"""
        self.intra=intra
        """list of floats: Dataset of intra-distance values"""
        self.inter=inter
        """list of floats: Dataset of inter-distance values"""
        if not(isinstance(intra, type(None)) or isinstance(nbits, type(None))):
            self.pintra=_array(self.intra).mean()/self.nbits
            """Average Hamming intra-distance"""
        if not(isinstance(inter, type(None)) or isinstance(nbits, type(None))):
            self.pinter=_array(self.inter).mean()/self.nbits
            """Average Hamming intra-distance"""
        if not(isinstance(intra, type(None)) or isinstance(nbits, type(None))):
            self.dks_intra=abs(_histogram(self.intra, bins=self.nbits+1, range=(0,self.nbits+1), density=True)[0].cumsum()-_binom.pmf(k=range(self.nbits+1),n=self.nbits,p=self.pintra).cumsum()).max()
            """Dks value for the introduced Hamming intra-distance distribution against the theoretical binomial of paramaters n=nbits, p=pintra"""
        if not(isinstance(inter, type(None)) or isinstance(nbits, type(None))):
            self.dks_inter=abs(_histogram(self.inter, bins=self.nbits+1, range=(0,self.nbits+1), density=True)[0].cumsum()-_binom.pmf(k=range(self.nbits+1),n=self.nbits,p=self.pinter).cumsum()).max()
            """Dks value for the introduced Hamming inter-distance distribution against the theoretical binomial of paramaters n=nbits, p=pinter"""        
        
    def load(self, filename):
        """
        This function allows you to load a Hamming distribution in text format. Caution: this function overrides any parameters defined when creating the class object.
        
        Parameters
        ----------
        filename : str, optional, 'None' by default
            Path to the file that contains the intra/inter Hamming distributios which is used to feed this class' object. The format of this file is:
                # Intra- / inter-hamming distributions 
                ninst: [Number of instances simulated]
                nrep: [Number of repetitions]
                nbits: [Number of PUF's responses bits]
                pintra: [Average intra-distance simulated]
                pinter: [Average inter-distance simulated]
                intra: [List of intra-distance values]
                inter: [List of inter-distance values]        
        """
        data_read = _read_dictionary_from_file(filename)
        self.ninst=int(data_read['ninst'])
        self.nrep=int(data_read['nrep'])
        self.nbits=int(data_read['nbits'])
        self.intra=list(map(float, data_read['intra'].strip(' []\n').split(',')))
        self.inter=list(map(float, data_read['inter'].strip(' []\n').split(',')))
        self.pintra=_array(self.intra).mean()/self.nbits
        self.pinter=_array(self.inter).mean()/self.nbits
        self.dks_intra=abs(_histogram(self.intra, bins=self.nbits+1, range=(0,self.nbits+1), density=True)[0].cumsum()-_binom.pmf(k=range(self.nbits+1),n=self.nbits,p=self.pintra).cumsum()).max()
        self.dks_inter=abs(_histogram(self.inter, bins=self.nbits+1, range=(0,self.nbits+1), density=True)[0].cumsum()-_binom.pmf(k=range(self.nbits+1),n=self.nbits,p=self.pinter).cumsum()).max()
            
    def save(self, filename):
        """
        Saves the Hamming distribution object to a file.
        
        Parameters
        ----------
        filename : str
            The name of the file to which the Hamming distribution should be saved.
        
        Returns
        -------
        None
        """
        data_to_write = f'# Intra- / inter-hamming distributions\n# Average intra-distance = {100*self.pintra:.3f}%\n# Average inter-distance = {100*self.pinter:.3f}%\nninst: {self.ninst}\nnrep: {self.nrep}\nnbits: {self.nbits}\nintra: {[int(i) for i in self.intra]}\ninter: {[int(i) for i in self.inter]}'
        with open(filename, "w") as f:
            f.write(data_to_write)
            
    def intra_fit(self, x):
        """
        Returns the binomial mass function value for the Hamming intra-distribution at x.
        
        Parameters
        ----------
        x : int or list of int
            The value for which to calculate the probability value.
        
        Returns
        -------
        float
            The probability value at x.
        """
        return _binom.pmf(k=x, n=self.nbits,p=self.pintra)
        
    def inter_fit(self, x):
        """
        Returns the binomial mass function value for the Hamming inter-distribution at x.
        
        Parameters
        ----------
        x : int or list of int
            The value for which to calculate the probability value.
        
        Returns
        -------
        float
            The probability value at x.
        """
        return _binom.pmf(k=x, n=self.nbits,p=self.pinter)        

    def plot(self, intra=True, intra_fit=True, inter=True, inter_fit=True):
        """
        Plots the Hamming distribution. If intra or intra_fit is True, it plots the intra-distribution and its fit. If inter or inter_fit is True, plots the inter-distribution and its fit.
        
        Parameters
        ----------
        intra : bool, optional, 'True' by default
            Whether to plot the intra-distribution.
            
        intra_fit : bool, optional, 'True' by default
            Whether to plot the fit of the intra-distribution.
            
        inter : bool, optional, 'True' by default
            Whether to plot the inter-distribution.
            
        inter_fit : bool, optional, 'True' by default
            Whether to plot the fit of the inter-distribution.
        
        Returns
        -------
            None
        """
        if intra:
            _,intra_bin_edges,_ = _hist(self.intra, bins=self.nbits+1, range=(0,self.nbits+1), density=True, align='left', rwidth=0.9)
        if intra_fit:
            _plot(range(self.nbits+1), _binom.pmf(k=range(self.nbits+1), n=self.nbits, p=self.pintra), label=f'n={self.nbits:.0f}\np={self.pintra:.3f}')
        if inter:
            _,inter_bin_edges,_ = _hist(self.inter, bins=self.nbits+1, range=(0,self.nbits+1), density=True, align='left', rwidth=0.9)
        if inter_fit:
            _plot(range(self.nbits+1), _binom.pmf(k=range(self.nbits+1), n=self.nbits, p=self.pinter), label=f'n={self.nbits:.0f}\np={self.pinter:.3f}')
        _legend()
        _show()
        
        
class DksDistribution:
    """
    This class contains the Kolmogorov-Smirnov (Dks) distance distributions between the ideal and simulated binomial distributions for the intra / inter Hamming distance based on a quasideal PUF model.
    
    ninst : int, optiona, 'None' by default
        Number of instances simulated

    nrep : int, optiona, 'None' by default
        Number of repetitions

    nbits : int, optiona, 'None' by default
        Number of PUF's responses bits

    pintra : float, optiona, 'None' by default
        Average intra-distance simulated

    pinter : float, optiona, 'None' by default
        Average inter-distance simulated

    intra : list of floats, optional, 'None' by default
        Dataset of Dks intra-distance values

    inter : list of floats, optional, 'None' by default
        Dataset of Dks inter-distance values
    
    Methods
    -------
    load(filename)
        his function allows to load a Dks distribution in text format.
        
    save(filename)
        Saves the Dks distribution object to a file.
        
    plot(intra=True, intra_fit=True, inter=True, inter_fit=True)
        Plots the Hamming distribution. If intra or intra_fit is True, it plots the intra-distribution and its fit. If inter or inter_fit is True, plots the inter-distribution and its fit.
        
    intra_fit(x)
        Returns the probability density function (pdf) value for the Dks intra-distribution at x.
        
    intra_alpha_value(alpha)
        This function returns the Dks coordinate corresponding to a given 'alpha' value for the intra-Dks distribution.
        
    intra_p_value(observed_intra_Dks):
        Given the Johnson SB fit to the Dks intra distribution, this function computes the probability of observing a Dks value as extreme, or more extreme, as the value observed (i.e., computes the area below the Johnson SB density from the observed data to infinity). If the value returned by this method is smaller than a significance level alpha, the hypothesis that the observed Hamming intra-distance distribution fits a binomial must be rejected.
        
    inter_fit(x):
        Returns the probability density function (pdf) value for the Dks inter-distribution at x.
        
    inter_alpha_value(alpha)
        This function returns the Dks coordinate corresponding to a given 'alpha' value for the inter-Dks distribution.        
        
    inter_p_value(observed_inter_Dks):
        Given the Johnson SB fit to the Dks inter distribution, this function computes the probability of observing a Dks value as extreme, or more extreme, as the value observed (i.e., computes the area below the Johnson SB density from the observed data to infinity). If the value returned by this method is smaller than a significance level alpha, the hypothesis that the observed Hamming inter-distance distribution fits a binomial must be rejected.

    plot(self, intra=True, intra_fit=True, inter=True, inter_fit=True, intra_bins=10, inter_bins=10):
        Plots the Hamming distribution. If intra or intra_fit is True, plots the intra-distribution and its fit. If inter or inter_fit is True, plots the inter-distribution and its fit.
    """    
    def __init__(self, ninst=None, nrep=None, nbits=None, pintra=None, pinter=None, intra=None, inter=None):
        """
        Constructor method.
        """
        self.ninst=ninst
        """int: Number of instances simulated"""
        self.nrep=nrep
        """int: Number of repetitions"""
        self.nbits=nbits
        """int: Number of PUF's responses bits"""
        self.pintra=pintra
        """float: Average intra-distance simulated"""
        self.pinter=pinter
        """float: Average inter-distance simulated"""
        self.intra=intra
        """list of floats: Dataset of Dks intra-distance values"""
        self.inter=inter
        """list of floats: Dataset of Dks inter-distance values"""
        if not isinstance(intra, type(None)):
            self.intra_fit_params=_johnsonsb.fit(self.intra)
            """list of floats: Parameters of the Johnson SB distribution that fits best the dataset of Dks intra-distance provided. If len(intra)<2 this attribute is not created"""
        if not isinstance(inter, type(None)):
            self.inter_fit_params=_johnsonsb.fit(self.inter)
            """list of floats: Parameters of the Johnson SB distribution that fits best the dataset of Dks inter-distance provided. If len(inter)<2 this attribute is not created"""            

    def load(self, filename):
        """
        This function allows to load a Dks distribution in text format.
        
        Parameters
        ----------
        filename : str, optional, 'None' by default
            Path to the file that contains the Dks intra/inter distributios which is used to feed this class' object. The format of this file is:
                # Dks Intra / inter distributions
                ninst: [Number of instances simulated]
                nrep: [Number of repetitions]
                nbits: [Number of PUF's responses bits]
                pintra: [Average intra-distance simulated]
                pinter: [Average inter-distance simulated]
                intra: [List of Dks intra-distance values]
                inter: [List of Dks inter-distance values]        
        """
        data_read = _read_dictionary_from_file(filename)
        self.ninst=int(data_read['ninst'])
        self.nrep=int(data_read['nrep'])
        self.nbits=int(data_read['nbits'])
        self.pintra=float(data_read['pintra'])
        self.pinter=float(data_read['pinter'])
        self.intra=list(map(float, data_read['intra'].strip(' []\n').split(',')))
        self.inter=list(map(float, data_read['inter'].strip(' []\n').split(','))) 
        if len(self.intra)>1:
            self.intra_fit_params=_johnsonsb.fit(self.intra)
        if len(self.inter)>1:
            self.inter_fit_params=_johnsonsb.fit(self.inter)
    
    def save(self, filename):
        """
        Saves the Dks distribution object to a file.
        
        Parameters
        ----------
        filename : str
            The name of the file to which the Hamming distribution should be saved.
        
        Returns
        -------
        None
        """
        data_to_write = f'# Dks Intra / inter distributions\nninst: {self.ninst}\nnrep: {self.nrep}\nnbits: {self.nbits}\npintra: {self.pintra}\npinter: {self.pinter}\nintra: {[float(i) for i in self.intra]}\ninter: {[float(i) for i in self.inter]}'
        with open(filename, "w") as f:
            f.write(data_to_write)

    def intra_fit(self, x):
        """
        Returns the probability density function (pdf) value for the Dks intra-distribution at x.
        
        Parameters
        ----------
        x : float
            The value for which to calculate the pdf value.
        
        Returns
        -------
        float
            The pdf value at x.
        """
        return _johnsonsb.pdf(x, *self.intra_fit_params)
        
    def intra_alpha_value(self, alpha):
        """
        This function returns the Dks coordinate corresponding to a given 'alpha' value for the intra-Dks distribution.

        Parameters
        ----------
        alpha : float
            Value (per one) of the cumulative area on the right tail whose corresponding Dks value is to be returned.
            
        Returns
        -------
        float
            Dks value such that the cumulative probabilioty from it to +infinity equals 'alpha'.
        """
        return _johnsonsb.isf(alpha, *self.intra_fit_params)
        
    def intra_p_value(self, observed_intra_Dks):
        """
        Given the Johnson SB fit to the Dks intra distribution, this function computes the probability of observing a Dks value as extreme, or more extreme, as the value observed (i.e., computes the area below the Johnson SB density from the observed data to infinity). If the value returned by this method is smaller than a significance level alpha, the hypothesis that the observed Hamming intra-distance distribution fits a binomial must be rejected.
        
        Parameters
        ----------
        observed_intra_Dks : float
            The observed value of Dks in an experimental Hamming intra-distance distribution.

        Returns
        -------
        float
            The p-value corresponding to the observed Dks.
        """
        return _johnsonsb.sf(observed_intra_Dks, *self.intra_fit_params)

    def inter_fit(self, x):
        """
        Returns the probability density function (pdf) value for the Dks inter-distribution at x.
        
        Parameters
        ----------
        x : float
            The value for which to calculate the pdf value.
        
        Returns
        -------
        float
            The pdf value at x.
        """
        return _johnsonsb.pdf(x, *self.inter_fit_params)
        
    def inter_alpha_value(self, alpha):
        """
        This function returns the Dks coordinate corresponding to a given 'alpha' value for the inter-Dks distribution.

        Parameters
        ----------
        alpha : float
            Value (per one) of the cumulative area on the right tail whose corresponding Dks value is to be returned.
            
        Returns
        -------
        float
            Dks value such that the cumulative probabilioty from it to +infinity equals 'alpha'.
        """
        return _johnsonsb.isf(alpha, *self.inter_fit_params)        
        
    def inter_p_value(self, observed_inter_Dks):
        """
        Given the Johnson SB fit to the Dks inter distribution, this function computes the probability of observing a Dks value as extreme, or more extreme, as the value observed (i.e., computes the area below the Johnson SB density from the observed data to infinity). If the value returned by this method is smaller than a significance level alpha, the hypothesis that the observed Hamming inter-distance distribution fits a binomial must be rejected.
        
        Parameters
        ----------
        observed_inter_Dks : float
            The observed value of Dks in an experimental Hamming inter-distance distribution.

        Returns
        -------
        float
            The p-value corresponding to the observed Dks.
        """
        return _johnsonsb.sf(observed_inter_Dks, *self.inter_fit_params)        

    def plot(self, intra=True, intra_fit=True, inter=True, inter_fit=True, intra_bins=10, inter_bins=10):
        """
        Plots the Hamming distribution. If intra or intra_fit is True, plots the intra-distribution and its fit. If inter or inter_fit is True, plots the inter-distribution and its fit.
        
        Parameters
        ----------
        intra : bool, optional, 'True' by default
            Whether to plot the Dks intra-distribution.
            
        intra_fit : bool, optional, 'True' by default
            Whether to plot the fit of the Dks intra-distribution.
            
        inter : bool, optional, 'True' by default
            Whether to plot the Dks inter-distribution.
            
        inter_fit : bool, optional, 'True' by default
            Whether to plot the fit of the Dks inter-distribution.
            
        intra_bins : int, optional, '10' by default
            The number of bins for the Dks intra-distribution plot.
            
        inter_bins : int, optional, '10' by default
            The number of bins for the Dks inter-distribution plot.
        
        Returns
        -------
            None
        """    
        if intra:
            _,intra_bin_edges,_ = _hist(self.intra, bins=intra_bins, density=True, rwidth=0.9)
        if intra_fit:
            _plot(_linspace(0,intra_bin_edges[-1],50), [self.intra_fit(x) for x in _linspace(0,intra_bin_edges[-1],50)], label=f'a={self.intra_fit_params[0]:.3f}\nb={self.intra_fit_params[1]:.3f}\nloc={self.intra_fit_params[2]:.3f}\nscale={self.intra_fit_params[3]:.3f}')
        if inter:
            _,inter_bin_edges,_ = _hist(self.inter, bins=inter_bins, density=True, rwidth=0.9)
        if inter_fit:
            _plot(_linspace(0,inter_bin_edges[-1],50), [self.inter_fit(x) for x in _linspace(0,inter_bin_edges[-1],50)], label=f'a={self.inter_fit_params[0]:.3f}\nb={self.inter_fit_params[1]:.3f}\nloc={self.inter_fit_params[2]:.3f}\nscale={self.inter_fit_params[3]:.3f}')
        _legend()
        _show()
        
 
## Module's functions
def sim_quasideal_hamming_distribution(ninst, nrep, nbits, pintra, pinter, return_data_only=True, verbose=False):
    """
    This function simulates the result of intra/inter-Hamming distances for a PUF experiment using a quasi-ideal model. The model uses the experimentally observed quantities `pintra` and `pinter` to calculate the parameters
    `p_inst` and `p_rep`. `p_inst` determines the probability that a bit is '1' in an instance, and `p_rep`
    determines the probability that a bit is '1' in a repetition (technically, that the difference from a 
    'golden key' is '1'). This ensures that the distributions of intra/inter-Hamming distances follow _binomials 
    with parameters `pintra` and `pinter` respectively.

    Parameters
    ----------
    ninst : int
        Number of instances simulated
    nrep : int
        Number of repetitions
    nbits : int
        Number of PUF's responses bits
    pintra : float
        Average intra-distance simulated
    pinter : float
        Average inter-distance simulated
    return_data_only : bool, optional, 'True' by default
        If 'True', this function only returns a double-column numpyarray with the intra / inter hamming distributed values respectively. On the other hand, if 'False', this function return a whole 'HammingDistribution' object.
    verbose : bool, optional, 'False' by default
        If True, prints the calculated `p_rep` and `p_inst` values. Default is False.

    Returns
    -------
    tuple of list of float or `HammingDistribution` object
        A tuple containing `intra_set` and `inter_set`, which are lists of intra and inter-distances respectively, or a `HammingDistribution` object with the same content.

    """
    p_rep = 1/2-_sqrt(1-2*pintra)/2
    p_inst = 1/2-1/2*_sqrt((2*pinter-1)/(2*pintra-1))
    if verbose:
        print(f'p_rep: {p_rep}\np_inst: {p_inst}')
    
    # quasi-ideal data generation:
    key=[]
    data=[]
    for i in range(ninst):
        key.append(_choice([1,0],p=[p_inst,1-p_inst],size=nbits))
        for j in range(nrep):
            data.append(key[i]^_choice([1,0],p=[p_rep,1-p_rep],size=nbits))
    data = _array(data).reshape(ninst,nrep,nbits)
    
    # Intra-distance computation:
    intra_set=[]
    for i in range(ninst):
        for j in range(nrep):
            for k in range(j+1,nrep,1):
                intra_set.append((data[i][j]^data[i][k]).sum())

    # Inter-distance computation:
    inter_set=[]
    for i in range(nrep):
        for j in range(ninst):
            for k in range(j+1,ninst,1):
                inter_set.append((data[j][i]^data[k][i]).sum())
    
    if return_data_only:
        return [intra_set, inter_set]
    else:
        return HammingDistribution(ninst=ninst,nrep=nrep,nbits=nbits,intra=intra_set,inter=inter_set)
    
    
def montecarlo_quasideal_Dks_distribution(N, ninst, nrep, nbits, pintra, pinter, verbose=False):
    """
    This function computes the distributions of the Kolmogorov-Smirnov statistic of a quasi-ideal model against binomial distributions of parameters n=nbits, p=pintra, p=pinter. The function returns a 2D list of size `N` of these KS indices.

    Parameters
    ----------
    ninst : int
        Number of instances simulated
    nrep : int
        Number of repetitions
    nbits : int
        Number of PUF's responses bits
    pintra : float
        Average intra-distance simulated
    pinter : float
        Average inter-distance simulated
    verbose : bool, optional, 'False' by default
        If True, prints the calculated `p_rep` and `p_inst` values. Default is False.

    Returns
    -------
    `DksDistribution` object
        `DksDistribution` object containing the set of intra / inter Dks computed.

    """    
    fit_intra = _binom.pmf(k=range(nbits+1),n=nbits,p=pintra)
    fit_inter = _binom.pmf(k=range(nbits+1),n=nbits,p=pinter)
    
    intra_set=[] 
    inter_set=[]
    for i in range(N):
        hamming_distributed_values = sim_quasideal_hamming_distribution(ninst=ninst,nrep=nrep,nbits=nbits,pintra=pintra,pinter=pinter)
        
        hist_intra,_ = _histogram(hamming_distributed_values[0], bins=nbits+1, range=(0,nbits+1), density=True)
        hist_inter,_ = _histogram(hamming_distributed_values[1], bins=nbits+1, range=(0,nbits+1), density=True) 
    
        intra_set.append(abs(hist_intra.cumsum()-fit_intra.cumsum()).max())
        inter_set.append(abs(hist_inter.cumsum()-fit_inter.cumsum()).max())
        
        if verbose:
            if (i+1)%100==0:
                print(f"{(i+1)/100:.0f}/{N/100:.0f}",end='\r')
                
    return DksDistribution(ninst=ninst,nrep=nrep,nbits=nbits,pintra=pintra,pinter=pinter,intra=intra_set,inter=inter_set)
