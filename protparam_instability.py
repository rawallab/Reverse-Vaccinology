import subprocess,os,pandas as pd
import time,glob
import sys,shutil
from multiprocessing import Process
from Bio import SeqIO
import os,re
import threading
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import os
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, _verify_alphabet
from Bio import SeqIO


def protparm(cudir,filename,name):
   
   
   fasta_sequence = SeqIO.parse(open(cudir+"/"+name+"/"+filename),"fasta")
   for fasta in fasta_sequence:
        name1, sequence = fasta.id, str(fasta.seq)
        ##print sequence
        X = ProteinAnalysis(sequence)
        ##print name1+"\t"+str(X.instability_index())
        if float(round(X.instability_index(),2))<40:
             ii=(round(X.instability_index(),2))
             stab="stable"
             stab_coff=1
        else:
            ii=(round(X.instability_index(),2))
            stab="unstable"   
            stab_coff=0  
   return ii,stab,stab_coff  




if __name__=="__main__":
  
  cudir=os.getcwd()
  infl=input("Enter the file name") 
  inpath = os.path.join(cudir+"/input_file")
  
  if os.path.exists(inpath) and os.path.isdir(inpath):
    shutil.rmtree(inpath)
  os.mkdir("input_file")
  filename=infl.split("\t")[-1]
  name=filename.split(".")[0]
  dirpath = os.path.join(cudir+"/"+name)
  
  if os.path.exists(dirpath) and os.path.isdir(dirpath):
    shutil.rmtree(dirpath)
  os.mkdir(cudir+"/"+name)
  subprocess.call(["cp",infl,cudir+"/input_file/"]) 
  fasta_sequence = SeqIO.parse(open(cudir+"/input_file/"+filename),"fasta")
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
  
  for fasta in fasta_sequence:
     with open(name+"_Stability_Binary_result.txt","a") as final_fl: 
      with open(name+"_Stability_value.txt","a") as inter_fl:
       
        name1, sequence = fasta.id, str(fasta.seq)
        #print name1,sequence 
        description=fasta.description
        with open(cudir+"/"+name+"/"+"sequence.fasta", "w") as handle:
                 count = SeqIO.write(fasta, handle, "fasta")
      
        iidx,stable,stable_coff=protparm(cudir,"sequence.fasta",name)
       
                   
        inter_fl.write(name1+"\t"+description+"\t"+str(iidx)+"\t"+stable+"\n")
        final_fl.write(name1+"\t"+description+"\t"+str(stable_coff)+"\n")   
         
