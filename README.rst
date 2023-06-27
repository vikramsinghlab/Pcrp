===============================
pcrp
===============================

pcrp (plant circadian rhythm protein predictor) uses random walk with restart and subgraph automorphism to learn representations for vertices in a protein interaction network and return nodes having similar representations to the seed nodes.

Usage
-----

**Example Usage**
    ``$pcrp data/TulsiPIN_Wt.txt data/Tulsi_genes_present.txt data/TulsiPIN.out data/OrcaOut_OteEDG data/TulsiPIN_Wt_Nodes.txt data/Ote2GO.txt data/Ote2KO.txt``

**--output**: *output_filename*

    The output is in a file named plantCR_linkage.txt inside the results directory created on runtime. The output file is tab separated where the first column contains the predicted circadian clock-associated proteins and the second column contains the algorithms supporting the prediction, i.e. `RWR` (random walk with restart) or `GDV` (graphlet degree vector) or `RWR and GDV` (when both are supporting the prediction). Intermediate results are saved in the `temp` directory created at runtime in the current working directory.
        
        Ote100040090151 GDV
        Ote100214260021 RWR and GDV
        Ote100203530131 RWR
        ...

**Full Command List**
    The complete list of command line options is available with ``$pcrp --help``


Requirements
------------
* wheel>=0.23.0
* Cython>=0.20.2
* argparse>=1.2.1
* networkx>=2.6.3
* numpy>=1.21.2
* pandas>=1.3.5
* scikit_learn>=1.0.2
* scipy>=1.7.3
* sympy>=1.9
* tqdm>=4.65.0

(may have to be independently installed) 

Data
----

The data used in this study can be accessed `here. <https://drive.google.com/file/d/1gOr3o86V2y-g470U3FK-ue4N3PN9itjJ/view?usp=sharing>`__

Installation
------------
#. cd pcrp
#. pip install -r requirements.txt 
#. python setup.py install


Citing
------
If you find `pcrp` useful in your research, we ask that you cite the following paper::


@article{Singh2022.03.02.482599,
    author = {Vikram Singh and Vikram Singh},
    title = {Characterizing circadian connectome of *O. tenuiflorum* using an integrated network theoretic framework},
    elocation-id = {2022.03.02.482599},
    year = {2022},
    doi = {10.1101/2022.03.02.482599},
    publisher = {Cold Spring Harbor Laboratory},
    URL = {https://www.biorxiv.org/content/early/2022/03/02/2022.03.02.482599},
    eprint = {https://www.biorxiv.org/content/early/2022/03/02/2022.03.02.482599.full.pdf},
    journal = {bioRxiv}
}
