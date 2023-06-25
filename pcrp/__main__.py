import sys
import os
import argparse

from pathlib import Path

from pcrp.walker import Walker
from pcrp.utils import generate_seed_list, annot2gene, filter_ids, read_col
from pcrp.gdv import sgdv
from pcrp.stattests import ptest_rwr, ptest_gdv, linkage_test


def cli(argv, Inpath):
    parser = argparse.ArgumentParser(
        prog = 'pcrp',
        description = """pcrp- identify circadian rhythm genes using the 
        positional information of core clock genes in the protein inteaction
        network"""
        )

    parser.add_argument('network',
        type = str,
        help = """Edge-list formatted Weighted protein interaction network (PIN) 
        file. If weights are absent every interaction will be given a weight
        of 1"""
    )

    parser.add_argument('seed',
        type = str,
        help = """Seed file containing core circadian clock protein ids,
        one id per line"""
    )

    parser.add_argument('orcaout',
        type = str,
        help = """Output ORCA file correpsonding to PIN"""
    )

    parser.add_argument('orcaoutpath',
        type = str,
        help = """ORCA ouputdirectory path for the random ensemble"""
    )
    
    parser.add_argument('nodesmap',
        type = str,
        help = """Node indices to map protein ids in the PIN"""
    )

    parser.add_argument('protein2go',
        type = str,
        help = """Gene Ontlogy annotations correpsonding to PIN nodes"""
    )
    
    parser.add_argument('protein2ko',
        type = str,
        help = """KEGG Ontlogy annotations correpsonding to PIN nodes"""
    )
    args = parser.parse_args(argv)
    
    return args
  
def run_process(args, Outpath = Path('temp')):
    seed_list = generate_seed_list(args.seed) # Seed genes/proteins
    
    # RWR
    wk = Walker(args.network, None)  # Network File
    ProbsL, OG = wk.rwr(seed_list, 0.7, 0.1)
    
    with open(Outpath / 'PlantCR_RWR.txt', 'w') as fOut:
        for proteinid in ProbsL:
            if ProbsL[proteinid] > 1e-05 and proteinid not in seed_list:
                print(proteinid, "\t", ProbsL[proteinid], sep = "", file = fOut)
                
    # permutation test on RWR output
    AssosGene_rwr = set()
    gfit = ptest_rwr(args.network, ProbsL, len(seed_list), 10000, 1e-5)

    fPtest = open(Outpath / "PlantCR_RWR_0.05.ptest", 'w')
    for i in sorted(gfit):
        if gfit[i] < 0.05 and i not in seed_list:
            AssosGene_rwr.add(i)
            print(i, gfit[i], sep = "\t", file = fPtest)
    
    # permutation test on GDV output
    gdv_ind, Ind2Ids = sgdv(args.nodesmap, args.seed, args.orcaout, Outpath)
    distscore = ptest_gdv(args.orcaoutpath, Outpath / gdv_ind, Outpath)
    
    AssosGene_gdv = set()
    fOut = open(Outpath / "PlantCR_GDV_0.05.ptest", 'w')
    for k in distscore:
        keys = distscore[k].keys();
        values = distscore[k].values()
        for key, value in distscore[k].items():
            pvalue = value / 1000
            if pvalue < 0.05 and Ind2Ids[str(key)] not in seed_list:
                AssosGene_gdv.add(Ind2Ids[str(key)])
                print(Ind2Ids[str(key)], "\t", pvalue, sep = "", file = fOut)
    fOut.close()

    # Linkage test
    GO2G, G2GO = annot2gene(args.protein2go)
    KO2G, G2KO = annot2gene(args.protein2ko)
     
    algo = 'RWR'
    
    sys.stdout.write("Running linkage test on candidate RWR genes ...\n")
    # Gene ontology based
    GO_SimMat, GoMax, GoMax_id, cgenes_go, sgenes_go = linkage_test(OG, GO2G, G2GO, seed_list, AssosGene_rwr)

    OutDict_go = {}; OutDict_ko = {}; OutDict_lnk = {}
    fOut = open(Outpath / f'PlantCR_Linkage_GO_{algo}.txt', 'w')
    for i in range(len(cgenes_go)):
        print(cgenes_go[i], sgenes_go[GoMax_id[i][1]], GoMax[i], sep = "\t", file = fOut)

        if cgenes_go[i] not in OutDict_lnk: OutDict_lnk[cgenes_go[i]] = [0, 0]
        OutDict_lnk[cgenes_go[i]][0] = GoMax[i]
    fOut.close()
  
    # KEGG ontology based
    KO_SimMat, KoMax, KoMax_id, cgenes_ko, sgenes_ko = linkage_test(OG, KO2G, G2KO, seed_list, AssosGene_rwr)

    fOut_ko = open(Outpath / f'PlantCR_Linkage_KO_{algo}.txt', 'w')
    for i in range(len(cgenes_ko)):
        print(cgenes_ko[i], sgenes_ko[KoMax_id[i][1]], KoMax[i], sep = "\t", file = fOut_ko)

        if cgenes_ko[i] not in OutDict_lnk: OutDict_lnk[cgenes_ko[i]] = [0, 0]
        OutDict_lnk[cgenes_ko[i]][1] = KoMax[i]
    fOut_ko.close()

    # verify the output

    fOut_lnk = open(Outpath / f"PlantCR_Linkage_{algo}.txt", 'w')
    for k in OutDict_lnk:
        print(k, *OutDict_lnk[k], sep = "\t", file = fOut_lnk)
    fOut_lnk.close()

    algo = 'GDV'

    sys.stdout.write("Running linkage test on candidate GDV genes ...\n")
    # Gene ontology based
    GO_SimMat, GoMax, GoMax_id, cgenes_go, sgenes_go = linkage_test(OG, GO2G, G2GO, seed_list, AssosGene_gdv)

    OutDict_go = {}; OutDict_ko = {}; OutDict_lnk = {}
    fOut = open(Outpath / f'PlantCR_Linkage_GO_{algo}.txt', 'w')
    for i in range(len(cgenes_go)):
        print(cgenes_go[i], sgenes_go[GoMax_id[i][1]], GoMax[i], sep = "\t", file = fOut)

        if cgenes_go[i] not in OutDict_lnk: OutDict_lnk[cgenes_go[i]] = [0, 0]
        OutDict_lnk[cgenes_go[i]][0] = GoMax[i]
    fOut.close()
    
    # KEGG ontology based    
    KO_SimMat, KoMax, KoMax_id, cgenes_ko, sgenes_ko = linkage_test(OG, KO2G, G2KO, seed_list, AssosGene_gdv)

    fOut_ko = open(Outpath / f'PlantCR_Linkage_KO_{algo}.txt', 'w')
    for i in range(len(cgenes_ko)):
        print(cgenes_ko[i], sgenes_ko[KoMax_id[i][1]], KoMax[i], sep = "\t", file = fOut_ko)

        if cgenes_ko[i] not in OutDict_lnk: OutDict_lnk[cgenes_ko[i]] = [0, 0]
        OutDict_lnk[cgenes_ko[i]][1] = KoMax[i]
    fOut_ko.close()

    # verify the output

    fOut_lnk = open(Outpath / f"PlantCR_Linkage_{algo}.txt", 'w')
    for k in OutDict_lnk:
        print(k, *OutDict_lnk[k], sep = "\t", file = fOut_lnk)
    fOut_lnk.close()
    
    s_rwr = filter_ids(Outpath / "PlantCR_Linkage_RWR.txt", 0.75)
    s_gdv = filter_ids(Outpath / "PlantCR_Linkage_GDV.txt", 0.75)
    
    union_ids = s_rwr.union(s_gdv)
    c_gdv = union_ids.intersection(s_gdv)
    c_rwr = union_ids.intersection(s_rwr)
    
    Algo = ['RWR', 'GDV', 'RWR and GDV']
    Results = Path('results')
    Results.mkdir(exist_ok = True)
    fOut = open(Results / "PlantCR_linkage.txt", 'w')
    fOut_ids = open(Outpath / "PlantCR_linkage_ids.txt", 'w')
    for id in union_ids:
        Flag = None
        if id in c_rwr and id in c_gdv:
            Flag = 2
        elif id in c_rwr:
            Flag = 0
        else:
            Flag = 1
        
        print(id, Algo[Flag], file = fOut)
        print(id, file = fOut_ids)
    fOut.close()
    fOut_ids.close()
    
def main():
    Outpath = Path('temp')
    Outpath.mkdir(exist_ok = True)
    args = cli(sys.argv[1:], Path('data'))
    run_process(args, Outpath)

if __name__ == '__main__':
    sys.exit(main())
