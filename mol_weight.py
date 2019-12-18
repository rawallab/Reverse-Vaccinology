import subprocess,os,pandas as pd,time,glob,re,sys,shutil

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, _verify_alphabet





if __name__=="__main__":

  cudir=os.getcwd()
  infl=sys.argv[1] 
  start_time=time.time()
  filename=infl.split("/")[-1]
  name=filename.split(".")[0]
  inpath = os.path.join(cudir+"/input_files")
  
  if os.path.exists(inpath) and os.path.isdir(inpath):
    shutil.rmtree(inpath)
  os.mkdir(cudir+"/input_files/")
  dirpath = os.path.join(cudir+"/"+name)
  
  if os.path.exists(dirpath) and os.path.isdir(dirpath):
    shutil.rmtree(dirpath)
  os.mkdir(cudir+"/"+name)
   
  fasta_sequence = SeqIO.parse(open(cudir+"/input_files/"+filename),"fasta")
  for fasta in fasta_sequence:
           
        name1, sequence = fasta.id, str(fasta.seq)
        
        my_prot = Seq(sequence, IUPAC.protein)
        if _verify_alphabet(my_prot)==False:
          with open(cudir+"/false_sequence444.fasta", "a") as handle:
                 count = SeqIO.write(fasta, handle, "fasta")
          continue 
        else:
         
         with open(cudir+"/"+name+"/"+name+".fasta", "a") as handle:
                 count = SeqIO.write(fasta, handle, "fasta")
  fasta_sequence = SeqIO.parse(open(cudir+"/"+name+"/"+name+".fasta"),"fasta")
  with open(cudir+"/"+name+"_molecular_weight.xlsx","w")as fleout:
   with open(cudir+"/"+name+"_mol_weight.xlsx","w") as molfl:
     fleout.write("Sequence_id\tDescription\tmolecular_weight\tGRAVY\tIsoelectric_point\tAromaticity\n")
     molfl.write("Sequence_id\tmolecular_weight\n")
     for fasta in fasta_sequence:
        molwei=0
        name1, sequence = fasta.id, str(fasta.seq)
        description=fasta.description
        analysed_seq = ProteinAnalysis(sequence)
        fleout.write(name1+"\t"+description+"\t"+str(round((analysed_seq.molecular_weight()/1000),2))+"\t"+str(round(analysed_seq.gravy(),2))+"\t"+str(round(analysed_seq.isoelectric_point(),2))+"\t"+str(round(analysed_seq.aromaticity(),2))+"\n")  
        if ((analysed_seq.molecular_weight())/1000)<110:
            molwei=1
        else:
            molwei=0 
        molfl.write(name1+"\t"+str(molwei)+"\n")
   # print time.time()-start_time
            
        
