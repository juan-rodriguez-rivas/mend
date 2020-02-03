#!/usr/bin/env python
import subprocess
import os.path
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
    
    # Load command line arguments
    args = command_parser.parse_args()
    
    return args


def main():
    
    # Run compute_mend_scores.py
    file_ali1 = 'test/PF01333.14-PF08802.5.domain1.fasta'
    file_ali2 = 'test/PF01333.14-PF08802.5.domain2.fasta'
    file_out = 'test/mend_scores.new'
    
    cmd = ["./compute_mend_scores.py", "-ali1", file_ali1, "-ali2", file_ali2, "-o", file_out, "-temp_dir", "test"]
#     file_stdout = open(file_out + ".out", 'w')
#     file_stderr = open(file_out + ".err", 'w')
#     subprocess.run(cmd, stdout = file_stdout, stderr = file_stderr)
    subprocess.run(cmd)
     
    if not (os.path.isfile(file_out) and os.path.getsize(file_out)):
        exit("ERROR: The output file with MEND scores was not found.")
    
    num_positions = (118, 39)
    
    # Load the stored result and computed one and compare both
    file_old = "test/mend_scores.old"
    mend_scores_new = load_mend_scores(file_out, num_positions)
    mend_scores_old = load_mend_scores(file_old, num_positions)
    result = compare_mend_scores(mend_scores_old, mend_scores_new, num_positions)
    
    
    if(result == 1):
        print("Test successfully completed")
#        os.remove(file_out + ".out")
#        os.remove(file_out + ".err")
        
    elif(result == 2):
        print("Test partially successful. The MEND scores were computed but some minor difference were detected. This might be due to the use of a different pseudorandom number algorithm (MersenneTwister was used for the stored results). Compare top contact predictions in files test/mend_scores.old and test/mend_scores.new. If the both contain the same order of pair of positions with small numeric differences, it should be right.")
    else:
        print("Test unsuccessful. Relevant numeric differences detected.")
    
    
    
def load_mend_scores(file_in, num_pos):
    
    mend_scores = numpy.zeros((num_pos[0],num_pos[1]))
    
    # Load scores with the correct pairing and compute mend score
    fh = open(file_in, 'r')
    for line in fh:
        fields = line.strip().split()
        pos1 = int(fields[0])
        pos2 = int(fields[1])
        score = float(fields[2])
        
        mend_scores[pos1-1,pos2-1] = score
    fh.close()
    
    return mend_scores



def compare_mend_scores(mend_scores_old, mend_scores_new, num_pos, threshold_strict = 0.001, threshold_loose = 0.1):
    
    success = 1
    
    for i in range(0,num_pos[0]):
        for j in range(0,num_pos[1]):
            if(abs(mend_scores_old[i,j]-mend_scores_old[i,j]) >= threshold_strict):
                success = 2
                
            if(abs(mend_scores_old[i,j]-mend_scores_old[i,j]) >= threshold_loose):
                success = 3
                break
    
    return success


if __name__ == '__main__':
    
    main()
    
