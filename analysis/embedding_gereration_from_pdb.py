import esm
import Bio.PDB
import os, re
import glob
import numpy as np
import torch
import argparse


def calc_residue_dist(residue_one, residue_two,ca='CA') :
    """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one[ca].coord - residue_two[ca].coord
    return np.sqrt(np.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two,ca='CA') :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = np.zeros((len(chain_one), len(chain_two)),'float')
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist(residue_one, residue_two, ca)
    return answer

def data_gener_from_pdb(outdir,mode='esmif',indir=None):
    ## output directory
    try:
        os.makedirs(outdir)
    except:
        print('Output directory existed: %s' % outdir)
    
    ## input pdbs 
    if indir==None:
        pdb_list = glob.glob('./dataset_5/Structure/ESMfold/cdhit90'+'/*.pdb')+glob.glob('./dataset_5/Structure/ESMfold/casPedia'+'/*.pdb')
    else:
        pdb_list = glob.glob(indir+'/*.pdb')
    
    ## esmfold_if
    if mode == 'esmif':
        ## load esmfold_if
        model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()   ## convert pdb structure to embedding (dim=512)
        model = model.eval().cuda()
        with torch.no_grad():
            for this_pdb in pdb_list:
                name = re.sub('\.pdb$','',os.path.basename(this_pdb))
                ## converting
                chain_id = 'A'
                structure = esm.inverse_folding.util.load_structure(this_pdb, chain_id)
                coords, native_seq = esm.inverse_folding.util.extract_coords_from_structure(structure)
                rep = esm.inverse_folding.util.get_encoder_output(model, alphabet, coords)
                ## save
                np.save(os.path.join(outdir,name),rep.cpu().detach().numpy())
    ## distance map
    elif mode == 'distance':
        for this_pdb in pdb_list:
            name = re.sub('\.pdb$','',os.path.basename(this_pdb))
            structure = Bio.PDB.PDBParser().get_structure(name, this_pdb)
            model = structure[0]
            dist_matrix = calc_dist_matrix(model["A"], model["A"],ca='CA')
            np.save(os.path.join(outdir,name),dist_matrix)
            
def main():

    parser = argparse.ArgumentParser(description='coverting pdb to representations')
    parser.add_argument('--outdir', type=str, default='./dataset_5/DistanceMap', help='output directiory')
    parser.add_argument('--mode', type=str, default='distance', help='"esmif" or "distance".')
    parser.add_argument('--indir', type=str, default=None, help='input data directory. default None.')
    args = parser.parse_args()
    
    data_gener_from_pdb(args.outdir,args.mode,args.indir)

if __name__ == '__main__':

    main()
