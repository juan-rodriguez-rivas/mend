#!/usr/bin/env python
from Bio import AlignIO
from Bio import SeqIO
import subprocess
import random
import os.path
import inspect
import sys
import numpy


# list of natural aminoacids
alphabet = [
'-',
'A','C','D','E',
'F','G','H','I',
'K','L','M','N',
'P','Q','R','S',
'T','V','W','Y',
]

# add 1 to the indices (gaps will have index 1)
AAMap = {alphabet[i]: i+1 for i in range(len(alphabet))}

def get_arguments():
    import argparse
    
    # Define command line input
    command_parser = argparse.ArgumentParser(description='Given a paired alignment (two FASTA files), it computes and print out the MEND score for each pair in inter-domain (or inter-protein) positions')
    command_parser.add_argument('-ali1', type=str, required=True, help='Input file of the first domain (or protein) alignment in FASTA format')
    command_parser.add_argument('-ali2', type=str, required=True, help='Input file of the second domain (or protein) alignment in FASTA format')
    command_parser.add_argument('-o', type=str, required=True, help='File output path. Default standard output')
    command_parser.add_argument('-num_rands', type=int, default=10, help='Optonal. The number of randomization to apply. 10, by default')
    command_parser.add_argument('-seed', type=int, default=1, help='Optional. Seed to be use in the randomizations. By default, seed=1')
    command_parser.add_argument('-temp_dir', type=str, default='.', help='Optional. Directory where temporary files will be stored. By default, the working directory is used')
    command_parser.add_argument('-keep_temp', type=bool, default=False, help='Optional, if true temporary files are kept, if false they are removed. By default, false')

    
    # Load command line arguments
    args = command_parser.parse_args()
    
    return args


def main(file_ali1, file_ali2, num_rands, seed, file_out, temp_dir, keep_temp):
    
    # TODO: Check input parameters
    
    # Randomize pairings
    out_filename = os.path.basename(file_out)
    randomized_alis_root = temp_dir + '/'+ out_filename + '.ali2.shuf_'
    file_out_seed = temp_dir + '/' + out_filename + '.seed'
    randomize_pairings(file_ali2, randomized_alis_root, file_out_seed, num_rands, seed)
    
    # Compute DCA models
    computed_models = temp_dir + '/' + out_filename
    num_pos = compute_dca_models(file_ali1, file_ali2, randomized_alis_root, num_rands, computed_models)

    # Apply MEND correction
    mean_maximum = compute_mean_maximum(num_rands, computed_models)
    mend_scores = compute_mend_scores(computed_models, mean_maximum, num_pos)
    print_mend_scores(mend_scores, num_pos, file_out)
    
    # Remove temporary files
    if(not keep_temp):
        clean_temporary_files(num_rands, randomized_alis_root, file_out_seed, computed_models)
        
    
def clean_temporary_files(num_rands, randomized_alis_root, file_out_seed, computed_models):
    
    os.remove(file_out_seed)
    root = computed_models + '.correct_pairing'
    remove_files_root(root)
    
    for i in range(1, num_rands+1):
        
        os.remove(randomized_alis_root + str(i))
        root = computed_models + '.shuf_' + str(i)
        remove_files_root(root)

    
        
def remove_files_root(root):
    
    os.remove(root)
    os.remove(root + '.scores')
    os.remove(root + '.dim')
    os.remove(root + '.out')
    os.remove(root + '.err')
    os.remove(root + '.apc.scores.inter')
    os.remove(root + '.apc.scores.intra1')
    os.remove(root + '.apc.scores.intra2')
    os.remove(root + '.raw.scores.inter')
    os.remove(root + '.raw.scores.intra1')
    os.remove(root + '.raw.scores.intra2')
    


def randomize_pairings(ali_file, out_root, file_out_seed, num_rands=10, seed=1):
    """ Given a FASTA alignment, write num_rands FASTA files where the order of sequences is randomized 
    """
    
    # Load alignemnt
    ali = AlignIO.read(open(ali_file), "fasta")
    filename = os.path.basename(ali_file)
    
    # Set seed and save it in an output file
    random.seed(seed)
    
    seed_file = open(file_out_seed, 'w')
    seed_file.write(str(seed) + "\n")
    seed_file.close()
     
    # Generate output alignments files with randomized pairings
    for i in range(1,num_rands+1):
        out_root_current = out_root + str(i)
        f = open(out_root_current, 'w')
        random.shuffle(ali._records)
        AlignIO.write(ali, f, 'fasta')
        f.close()
        

def compute_dca_models(file_ali1, file_ali2, randomized_alis_root, num_rands, out_root):
    
    # Find mpl executable path
    filename_this = inspect.getframeinfo(inspect.currentframe()).filename
    path_this = os.path.dirname(os.path.abspath(filename_this))
    path_mpl_exec = path_this + '/mpl/mpl'
    
    # Compute model with the correct pairing
    # Prepare input mpl file
    out_root_current = out_root + '.correct_pairing'
    make_mpl_input(file_ali1, file_ali2, out_root_current)
    
    # Run mpl
    cmd = [path_mpl_exec, "-i", out_root_current, "-l", "0.01", "-g"]
    file_out = open(out_root_current + ".out", 'w')
    file_err = open(out_root_current + ".err", 'w')
    subprocess.run(cmd, stdout = file_out, stderr = file_err)
    
    # Renumber scores and apply APC
    file_scores = out_root_current + '.scores'
    renumber_scores(out_root_current + ".dim", file_scores, out_root_current)
    
    # Computed models for MSAs with random pairings
    for i in range(1, num_rands+1):
    
        # Prepare input file for mpl
        out_root_current = out_root + '.shuf_' + str(i)
        file_rand_current = randomized_alis_root + str(i)
        make_mpl_input(file_ali1, file_rand_current, out_root_current)
        
        # Run mpl
        cmd = [path_mpl_exec, "-i", out_root_current, "-l", "0.01", "-g"]
        file_out = open(out_root_current + ".out", 'w')
        file_err = open(out_root_current + ".err", 'w')
        subprocess.run(cmd, stdout = file_out, stderr = file_err)
            
        # Renumber scores and apply APC
        file_scores = out_root_current + '.scores'
        num_pos = renumber_scores(out_root_current + ".dim", file_scores, out_root_current)
    
    return num_pos

def make_mpl_input(file_ali1, file_ali2, out_file):
    
    # Parse the MSA files
    try: 
        fasta_seqs1 = list(SeqIO.parse(open(file_ali1),'fasta'))
    except:
        sys.stderr.write("check FASTA file: "+file_ali1+"\n")
        exit()
    
    try: 
        fasta_seqs2 = list(SeqIO.parse(open(file_ali2),'fasta'))
    except:
        sys.stderr.write("check FASTA file: "+file_ali2+"\n")
        exit()
    
    if len(fasta_seqs1) != len(fasta_seqs2): 
        sys.stderr.write("Number of sequences in the two alignments differ\n")
        exit()
    else: 
        ns0 = len(fasta_seqs1)
        
    # MSA is a list of sequences
    msa1        = []
    msa2        = []
    msa         = []
    # Headers is a list of header
    headers1    = []
    headers2    = []
    
    # Loop over the paired alignment
    for i in range(ns0):
        name1, sequence1 = fasta_seqs1[i].id, str(fasta_seqs1[i].seq)
        name2, sequence2 = fasta_seqs2[i].id, str(fasta_seqs2[i].seq)
        if i == 0: 
            n1      = len(sequence1)
            n2      = len(sequence2)
            ntot    = n1+n2
        
        # Filter out unnatural aminoacids
        if numpy.sum([1 for aa in sequence1 if not aa in alphabet]) > 0: continue
        if numpy.sum([1 for aa in sequence2 if not aa in alphabet]) > 0: continue
        headers1.append(name1)
        headers2.append(name2)
        msa1.append([ x for x in sequence1])
        msa2.append([ x for x in sequence2])
        msa.append([ x for x in sequence1+sequence2])
    
    # Switch to a numpy matrix for slicing
    mat1    = numpy.asarray(msa1)
    mat2    = numpy.asarray(msa2)
    mat     = numpy.asarray(msa)
    
    # Number of sequences
    ns  = len(mat[:,0])
    # Number of positions in domain1
    nv1 = len(mat1[0,:])
    # Number of positions in domain2
    nv2 = len(mat2[0,:])
    # Total number of positions 
    nv  = len(mat[0,:])
    
    # Dump a temp file
    fprm=open(out_file + '.dim', 'w')
    fprm.write("%i %i %i #np1 np2 nseq \n" % (nv1, nv2, ns))
    
    
    with open(out_file,'w') as f: 
        # Dump input file(s) for mpl
        for s in range(ns): 
            f.write(' '.join('%2i' % AAMap[x] for x in mat[s,:])+'\n')


def renumber_scores(file_dim, file_in, out_root):

    # Read from dat
    for line in open(file_dim,'r'):
        nv1, nv2, ns = map(int, line.split()[:3]) # read the first 3 fields only
        break # read the first line only


    # Read scores
    with open(file_in) as f: lines=[line for line in f]
    nlines=sum([1 for line in lines if len(line.split())==3])
    if int(numpy.rint(0.5*(numpy.sqrt(8.0*nlines)+1.0))) != nv1+nv2:
        sys.stderr.write("ERROR: check num. of positions")
        exit()
    mat=numpy.zeros((nv1+nv2,nv1+nv2))
    for line in lines:
        a,b,s=line.split()
        a=int(a)
        b=int(b)
        mat[a-1,b-1]=float(s)
        mat[b-1,a-1]=float(s)
        
  
    # Compute averages over rows/columns
    sm1     = numpy.mean(mat[:nv1,:nv1])
    cm1     = numpy.mean(mat[:nv1,:nv1],axis=0)
    sm2     = numpy.mean(mat[nv1:,nv1:])
    cm2     = numpy.mean(mat[nv1:,nv1:],axis=0)
    sm12    = numpy.mean(mat[:nv1,nv1:])
    cm12    = numpy.mean(mat[:nv1,nv1:],axis=0)
    rm12    = numpy.mean(mat[:nv1,nv1:],axis=1)

    # Correct for APC as Baker's lab
    mat_apc = numpy.copy(mat)
    for a in range(nv1):
        for b in range(nv1):
            apc                     = cm1[a]*cm1[b]/sm1
            mat_apc[a,b]            = "{0:.5f}".format(mat[a,b]-apc)
    for a in range(nv2):
        for b in range(nv2):
            apc                     = cm2[a]*cm2[b]/sm2
            mat_apc[nv1+a,nv1+b]    = "{0:.5f}".format(mat[nv1+a,nv1+b]-apc)
    for a in range(nv1):
        for b in range(nv2):
            apc                     = rm12[a]*cm12[b]/sm12
            mat_apc[a,nv1+b]        = "{0:.5f}".format(mat[a,nv1+b]-apc)
            mat_apc[nv1+b,a]        = "{0:.5f}".format(mat[nv1+b,a]-apc)


    # Write files with raw scores
    u1      = open(out_root + '.raw.scores.intra1', 'w')
    u2      = open(out_root + '.raw.scores.intra2', 'w')
    uint    = open(out_root + '.raw.scores.inter', 'w')

    for a in range(nv1+nv2-1):
        for b in range(a+1,nv1+nv2):
            if b <= a: continue
            inda=a+1
            indb=b+1
            if inda <= nv1 and indb <= nv1:
                # intra-domain1
                sp=mat[a,b]
                u1.write(' '.join([str(x) for x in [inda,indb,sp]])+'\n')
            elif inda <= nv1 and indb > nv1:
                # inter-domain
                sp=mat[a,b]
                uint.write(' '.join([str(x) for x in [inda,indb-nv1,sp]])+'\n')
            elif inda > nv1 and indb > nv1:
                # intra-domain2
                sp=mat[a,b]
                u2.write(' '.join([str(x) for x in [inda-nv1,indb-nv1,sp]])+'\n')

    u1.close()
    u2.close()
    uint.close()

    # Write files with APC scores
    u1      = open(out_root + '.apc.scores.intra1', 'w')
    u2      = open(out_root + '.apc.scores.intra2', 'w')
    uint    = open(out_root + '.apc.scores.inter', 'w')

    for a in range(nv1+nv2-1):
        for b in range(a+1,nv1+nv2):
            if b <= a: continue
            inda=a+1
            indb=b+1
            if inda <= nv1 and indb <= nv1:
                # intra-domain1
                sp=mat_apc[a,b]
                u1.write(' '.join([str(x) for x in [inda,indb,sp]])+'\n')
            elif inda <= nv1 and indb > nv1:
                # inter-domain
                sp=mat_apc[a,b]
                uint.write(' '.join([str(x) for x in [inda,indb-nv1,sp]])+'\n')
            elif inda > nv1 and indb > nv1:
                # intra-domain2
                sp=mat[a,b]
                u2.write(' '.join([str(x) for x in [inda-nv1,indb-nv1,sp]])+'\n')

    u1.close()
    u2.close()
    uint.close()
    
    return [nv1,nv2]
    

def compute_mean_maximum(num_rands, randomized_models):
    
    maxima = list()
    
    # Find the maximum for each shuffling
    for i in range(1, num_rands+1):
         
        file_in_current = randomized_models + ".shuf_" + str(i) + '.apc.scores.inter'
        fh = open(file_in_current, 'r')
        int_list = list()
         
        for line in fh:
            fields = line.strip().split()
            int_list.append(round(float(fields[2]), 4))
        fh.close()
         
        maxima.append(max(int_list))
    
    # Return mean maximum
    array = numpy.array(maxima).astype(numpy.float)
    return round(numpy.mean(array), 5)


def compute_mend_scores(computed_models, mean_maximum, num_pos):
    
    mend_scores = numpy.zeros((num_pos[0],num_pos[1]))
    
    # Load scores with the correct pairing and compute mend score
    file_in = computed_models + '.correct_pairing.apc.scores.inter'
    fh = open(file_in, 'r')
    for line in fh:
        fields = line.strip().split()
        pos1 = int(fields[0])
        pos2 = int(fields[1])
        score = float(fields[2])
        
        mend_scores[pos1-1,pos2-1] = score/mean_maximum
    fh.close()
    
    return mend_scores


def print_mend_scores(mend_scores, num_pos, file_out):
    
    num_pos1 = num_pos[0]
    num_pos2 = num_pos[1]
    
    fh = open(file_out, 'w')
    for i in range(0,num_pos1):
        pos1 = str(i + 1)
        for j in range(0,num_pos2):
            pos2 =  str(j + 1)
            mend_score = str(round(mend_scores[i,j], 5))
            fh.write(pos1 + "\t" + pos2 + "\t" + mend_score + "\n")
    fh.close()



if __name__ == '__main__':
    
    args = get_arguments()
     
    file_ali1 = args.ali1
    file_ali2 = args.ali2
    file_out = args.o
    num_rands = args.num_rands
    seed = args.seed
    temp_dir = args.temp_dir
    keep_temp = args.keep_temp
    
#     Testing
#     file_ali1 = 'test/PF01333.14-PF08802.5.domain1.fasta'
#     file_ali2 = 'test/PF01333.14-PF08802.5.domain2.fasta'
#     file_out = 'mend_scores_new'
#     num_rands = 2
#     seed = 1
#     temp_dir = 'test'
#     keep_temp = False
   
    main(file_ali1, file_ali2, num_rands, seed, file_out, temp_dir, keep_temp)


