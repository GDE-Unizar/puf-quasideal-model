from setuptools import setup, find_packages

setup(
    name="PUF quasideal model",  
    version="0.1",  
    author="gds",  
    author_email="gds@unizar.es", 
    description="Python package to run and process a PUF quasideal model",  
    packages=find_packages(),  
    entry_points={'console_scripts':[
        'pqm-simulation=puf_qideal_model.console_scripts:pqm_simulation',
        'pqm-dks-distribution=puf_qideal_model.console_scripts:pqm_dks_distribution']}
)
