import pandas as pd
import h2o
import numpy as np
import os
import sys

this_dir = os.path.abspath(os.path.dirname(__file__))

def generate_potential_peptides(seq,max_intervening_distance = 20,model='combined'):
    n_cleavage_residues = 2 # REMOVE THIS - FIX IT
    output = []
    for start in range(1,len(seq)+1): # start point
        for length in range(1,9): # length of first precursor
            for distance in range(0,max_intervening_distance+1):
                for reversed in (True,False):
                    if (distance==0 and reversed == False): 
                        continue # don't allow spliced peptides that look like normal nonspliced peptides
                    if start + distance + 9 - 1 > len(seq):
                        continue
                    else:
                        start1 = start
                        length1 = length
                        start2 = start1 + distance + length1
                        end1 = start1 + length1 - 1
                        end2 = start1 + distance + 9 - 1
                        length2 = 9 - length1
                        if reversed == True:
                            start1, start2 = start2, start1
                            end1,end2 = end2,end1
                            length1,length2=length2,length1
                        peptide1 = seq[start1-1:end1]
                        peptide2 = seq[start2-1:end2]
                        peptide = peptide1+peptide2
                        binding1_1 = peptide[length1 - 1]
                        if length1 > 1:
                            binding1_2 = peptide[length1 - 2]
                        else:
                            binding1_2 = None
                        binding2_1 = peptide[length1]
                        if length2 > 1:
                            binding2_2 = peptide[length1 + 1]
                        else:
                            binding2_2 = None
                        try:
                            cleavage1 = seq[max(0, start1 - (n_cleavage_residues + 1)):start1 - 1]
                            cleavage2 = seq[end1:min(end1 + n_cleavage_residues, len(seq))]
                            cleavage3 = seq[max(0, start2 - (n_cleavage_residues + 1)):start2 - 1]
                            cleavage4 = seq[end2:min(end2 + n_cleavage_residues, len(seq))]
                        except:
                            continue
                        output.append(
                            [peptide,start1,end1,start2,end2,length1, binding1_2, binding1_1, binding2_1, binding2_2, reversed, distance,
                             cleavage1,
                             cleavage2, cleavage3, cleavage4])
    df = pd.DataFrame(output,
                             columns=['peptide','start1','end1','start2','end2', 'length1', 'binding1_2', 'binding1_1', 'binding2_1',
                                      'binding2_2',
                                      'reversed', 'distance', 'cleavage1', 'cleavage2', 'cleavage3',
                                      'cleavage4'])
    df[['seq_1', 'seq_2', 'seq_3', 'seq_4', 'seq_5', 'seq_6', 'seq_7', 'seq_8',
               'seq_9']] = pd.DataFrame([list(x) for x in df['peptide']])
    for cleavage in ['cleavage1', 'cleavage2', 'cleavage3', 'cleavage4']:
        cols = []
        for i in range(1, n_cleavage_residues + 1):
            cols.append(cleavage + '_' + str(i))
        #print(cols)
        #print([list(x) for x in df_output[cleavage]])
        df[cols] = pd.DataFrame([list(x) for x in df[cleavage]])
    df_original = df.copy(deep=True) # Change this to something more elegant
    df.drop('peptide',inplace=True,axis=1)
    df.drop(['start1','end1','start2','end2'],inplace=True,axis=1)
    df.drop(['cleavage1','cleavage2','cleavage3','cleavage4'],inplace=True,axis=1)
    
    df[['length1', 'distance']] = df[['length1', 'distance']].applymap(str)
    df.drop(['distance'], axis=1, inplace=True)
    df.fillna('',inplace=True)
    h2o.init(max_mem_size="2G")
    h2o.remove_all()
    df_h2o = h2o.H2OFrame(python_obj=df)
    clf = h2o.load_model(os.path.join(this_dir,'models',model))
    results = clf.predict(test_data = df_h2o).as_data_frame()
    percentiles = np.fromfile(os.path.join(this_dir,'percentiles'))
    results['percentile'] = pd.cut(results['p1'], percentiles, labels=[round(x * 0.1,1) for x in range(0, 1000)])
    results[['start1','end1','start2','end2','peptide']] = df_original[['start1','end1','start2','end2','peptide']]
    results = results.sort_values(by='p1',ascending=False)[['peptide','start1','end1','start2','end2','p1','percentile']]
    return(results)

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print('No input argument given')
        quit()
    if len(sys.argv) >= 4:
        print('Too many arguments given')
        quit()
    input = sys.argv[1]
    if os.path.isfile(os.path.join(this_dir, input)):
        with open(os.path.join(this_dir, input), 'r') as f:
            input_seq = ''.join(f.read().split())
    else:
        input_seq = input
    if len(sys.argv) == 3:
        output = sys.argv[2]
    else:
        output = None
    results = generate_potential_peptides(seq=input_seq)   
    if output == None:
        print(results)
    else:
        with open(os.path.join(this_dir,output),'w') as f:
            results.to_csv(f,index=None)
    
                        