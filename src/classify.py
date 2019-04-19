import pandas as pd
from subprocess import call
from os import remove
def build_identity_matrix(queries,subjects,blast_df):
    q_s_ident={}
    for row in blast_df.itertuples():
        q=row[1].strip()
        s=row[2].strip()
        i=row[3]
        if q not in q_s_ident:
            q_s_ident[q]={}
        q_s_ident[q][s]=i
            

    dm = []
    for el1 in queries:
        temp=[el1]
        for el2 in subjects:
            if el1 not in q_s_ident:
                if el2 not in q_s_ident:
                    temp.append(0)
                else:
                    try:
                        temp.append(q_s_ident[el2][el1])
                    except KeyError:
                        temp.append(0)
                
            else:
                try:
                    temp.append(q_s_ident[el1][el2])
                except KeyError:
                    temp.append(0)
        dm.append(temp)
    dm = pd.DataFrame(dm)
    dm.set_index(0,inplace=True)
    dm.columns = subjects
    return(dm)

def run_blast(queryfile,subjectfile,resultfile,makeDB=''):
    blastn_args=['blastn','-query',queryfile,
             '-outfmt','6','-out',resultfile]
    if makeDB:
        makeblastdb_args = ['makeblastdb','-in',subjectfile,'-dbtype','nucl','-out',makeDB]
        blastn_args.extend(['-db',makeDB])
        print("Launching: {}".format(' '.join(makeblastdb_args)))
        call(makeblastdb_args)
    else:
        blastn_args.extend(['-subject',subjectfile])
    print("Launching: {}".format(' '.join(blastn_args)))
    call(blastn_args)
    results = pd.read_csv(resultfile,sep='\t',header=None)
    results.columns = ['qseqid','sseqid','pident',
           'length','mismatch','gapopen',
           'qstart','qend','sstart','send','evalue','bitscore']
    return(results)

def classify(identity_matrix, cutoff):
    OTU_class = {}
    refs = list(identity_matrix.columns)
    for row in identity_matrix.itertuples():
        classification=[]
        for i,ident in enumerate(row[1:]):
            if ident >= cutoff:
                classification.append(refs[i])
        OTU_class[row[0]] = ';'.join(classification)
    return(OTU_class)

def classify_with_references(reference_fasta,MiSeq_fasta,feature_table,outfile,cutoff=95.0):
    blst_results='tempfile'
    results=run_blast(reference_fasta,MiSeq_fasta,resultfile=blst_results)
    refs = list(set(results['qseqid']))
    miseqASV = list(set(results['sseqid']))
    class_mtrx = build_identity_matrix(miseqASV,refs,results)
    otu_class=classify(class_mtrx,cutoff)

    feature = pd.read_csv(feature_table,sep='\t',header = 1)#,skiprows=2)
    ids = list(feature.columns)[0]
    otu_align = []
    for el in feature[ids]:
        try:
            otu_align.append(otu_class[el])
        except KeyError:
            otu_align.append('')
    feature['classification'] = otu_align
    remove(blst_results)
    feature.to_csv(outfile,index=False)

