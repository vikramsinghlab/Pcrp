===============================
pcrp
===============================

pcrp (plant circaidan rhythm protein predictor) uses random walk with restart and subgraph automorphism to learn representations for vertices in a protein interaction netowrk and return nodes having similar representations to the seed nodes.

Usage
-----

**Example Usage**
    ``$pcrp data/TulsiPIN_Wt.txt data/Tulsi_genes_present.txt data/TulsiPIN.out data/OrcaOut_OteEDG data/TulsiPIN_Wt_Nodes.txt data/Ote2GO.txt data/Ote2KO.txt``

**--output**: *output_filename*

    The output is present in a file named plantCR_linkage.txt inside the results directory created on runtime. The output file is tab separated where first column contains the predicted circadian clock associated proteins and the second column conatins the lagorithm supporting the prediction i.e. `RWR` (random walk with restart) or `GDV` (graphlet degree vector) or `RWR and GDV` (when both are supporing the prediction). Intermediate results are saved in `temp` directory created at runtime in the current working directory.
        
        Ote100040090151 GDV
        Ote100214260021 RWR and GDV
        Ote100203530131 RWR
        ...

**Full Command List**
    The full list of command line options is available with ``$pcrp --help``


Requirements
------------
* numpy
* pandas

(may have to be independently installed) 



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
    title = {Characterizing circadian connectome of O. tenuiflorum using an integrated network theoretic framework},
    elocation-id = {2022.03.02.482599},
    year = {2022},
    doi = {10.1101/2022.03.02.482599},
    publisher = {Cold Spring Harbor Laboratory},
    URL = {https://www.biorxiv.org/content/early/2022/03/02/2022.03.02.482599},
    eprint = {https://www.biorxiv.org/content/early/2022/03/02/2022.03.02.482599.full.pdf},
    journal = {bioRxiv}
}
