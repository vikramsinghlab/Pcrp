import pandas as pd
import sys

from typing import Tuple, Dict

def generate_seed_list(seed_file: str) -> list:
    """Read a single column seed file into a python list
    
    Parameters
    ----------
    
    seed_file (str):    File containing list of proteins, one protein id per 
                        line, to be used as seed for RWR algorithm.
                        
    Return
    ------
    
    list    A list object containing the seed ids
    """

    seed_list = []

    try:
        fp = open(seed_file, "r")
    except IOError:
        sys.exit("Error opening file {}".format(seed_file))
        
    for line in fp.readlines():
        seed_list.append( line.rstrip() )

    fp.close()
    return seed_list
    

def annot2gene(gene2annot: str) -> Tuple[dict, dict]:
    """Helper script to get annotation to proteins mapping
    
    Parameters
    ----------
    
    gene2annot (str):   A space separated file with protein id in first column
                        followed by annotaions corresponding to that id.  
    
    Return
    ------
    
    Tuple(dict, dict)   Annotations mapping dictionaries
    """
    A2G = {}; G2A = {}
    with open(gene2annot, 'r') as f1:
        for line in f1:
            Ary = line.strip().split()
            Id = Ary.pop(0)

            if Id in G2A:
                raise Exception( "Duplicate genes" )
            
            G2A[Id] = Ary
            
            for i in Ary:
                if i not in A2G:
                    A2G[i] = []
                A2G[i].append(Id)
    return A2G, G2A
    
def filter_genes(AssoGen: dict, criterion: list = [0.5, 0.5] ) -> dict:
    """Filter the ids as per the criterion specified by the user
    
    Parameters
    ----------
    
    AssoGen (dict):    A dictionary object with gene ids as keys and
                       FIS_GO	FIS_KO	MIS columns
    
    criterion (list):   A list object containing the minimum cutoff vallues 
                        (float) for FIS_GO and FIS_KO values
    
    Return
    ------
    """
    
    AssoGen = {id: AssoGen[id] for id in AssoGen 
                   if AssoGen[id][1] >= criterion[0] and AssoGen[id][2] >= criterion[1] } 
    return AssoGen
    
def my_print(AssoGen: dict, FName: str) -> None:
    """My custom print function that formats the output as per my requirements
        
        Parameters
        ----------
        
        AssoGen (dict):    A dictionary object with gene ids as keys and
                           FIS_GO	FIS_KO	MIS columns
        FName (str):    Name of the output file
    """
    
    fout = open( FName, 'w' )
    print( "Gene Id", 'FIS_GO', 'FIS_KO', 'MIS', sep = "\t", file = fout )
    for i in AssoGen:
        print(i, AssoGen[i][1], AssoGen[i][2], AssoGen[i][3], sep = "\t", file = fout )
    fout.close()
    
def read_dod(fname: str) -> Tuple[Dict[str, Dict[int, float]], list, list]:
    """Read the index file produced by sgdv function as a dictionary of 
    dictionary
    
    Parameters
    ----------
    
    fname (str): Output file of sgdv function
    """
    outdict = {}; xrow_indices = set(); yrow_indices = set()
    try:
        with open(fname, 'r') as fIn:
            for line in fIn:
                Ary = line.strip("\n").split("\t")
                
                yrow_indices.add(int(Ary[0]))
                protinds = [int(i) for i in Ary[2].split(" ")]
                xrow_indices.update(protinds)
                protsims = [float(i) for i in Ary[3].split(" ")] 

                outdict[int(Ary[0])] = dict(zip(protinds, protsims))
        
        return(outdict, list(xrow_indices), list(yrow_indices))
    except OSError as e:
        print("I/O error({0}): {1}".format(e.errno, e.strerror))
        sys.exit()

def read_col(Infile: str, 
    col_indices: list,
    sep: str = "\t", 
    header: bool = None, 
    index_col: list = False
    ) -> pd.DataFrame:
    """Read a space separated file as a pandas dataframe and return selected
    columns as data frame

    Parameters
    -----------
        FName (str):    Name of input file containing ids
        col_indices (list):    list of column indices to be extracted
        header (bool):    if the file contains a header (False)
     
    Return
    ------
        A padas data frame having as many columns as specified in col_indices 
        list
    """
    
    try:
        return pd.read_csv(Infile, usecols = col_indices, sep=sep, 
            header=header, index_col= index_col)
    except OSError as e:
        print("I/O error({0}): {1}".format(e.errno, e.strerror))
        sys.exit()

def read_list(FName: str, header: bool = False) -> list:
    """Read a file containing ids into a list

    Parameters
    -----------
        FName (str):    Name of input file containing ids
        header (bool):    if the file contains a header (False)
     
    Return
    ------
        A list of ids
    """
    
    try:
        with open(FName, 'r') as fIn:
            if header: next(fIn)
            return(fIn.read().splitlines())
    except OSError as e:
        print("I/O error({0}): {1}".format(e.errno, e.strerror))
        sys.exit()

def read_dict(FName: str, header: bool = False) -> dict:
    """Read a file containing ids into a dictionary

    Parameters
    -----------
        FName (str):    Name of input file containing ids
        header (bool):    if the file contains a header (False)
     
    Return
    ------
        A dict of ids and their corresponding values
    """
    
    try:
        with open(FName, 'r') as fIn:
            if header: next(fIn)
            dct = dict(l.split() for l in fIn.read().splitlines())
            return(dct)
    except OSError as e:
        print("I/O error({0}): {1}".format(e.errno, e.strerror))
        sys.exit()

def filter_ids(fname: str, cutoff: float) -> set:
    """Filters the ids on the basis of GO and KO scores

    Parameters
    -----------
        FName (str):    Name of input file containing ids
        cutoff (float):    Threshold value to filter the ids
     
    Return
    ------
        Ids (set) A set of ids matching the specified criterion
    """
    try:
        f1 = open(fname, 'r')
    except OSError as e:
        print("I/O error({0}): {1}".format(e.errno, e.strerror))
        sys.exit()
    
    Ids = set()    
    with f1:
        for line in f1:
            Ary = line.strip("\n").split("\t")
            
            if float(Ary[1]) > cutoff and float(Ary[2]) > cutoff:
                Ids.add(Ary[0])
    return Ids
