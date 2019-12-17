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



if __name__=="__main__":
  
  cudir=os.getcwd()
  infl=sys.argv[1] 
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
         ##print(my_prot, _verify_alphabet(my_prot))
         with open(cudir+"/"+name+"/"+name+".fasta", "a") as handle:
                 count = SeqIO.write(fasta, handle, "fasta")
  fasta_sequence = SeqIO.parse(open(cudir+"/"+name+"/"+name+".fasta"),"fasta")
  hel=0
  for fasta in fasta_sequence:
     with open(name+"_tmpred_binary.txt","a") as final_fl: 
      with open(name+"_tmpred_values.txt","a") as inter_fl:
        #print fasta 
        name1, sequence = fasta.id, str(fasta.seq)
        #print (name1,sequence) 
        description=fasta.description
        #print description  
        #print name1 
        #input("press key") 
        with open(cudir+"/"+name+"/"+"sequence.fasta", "w") as handle:
                 count = SeqIO.write(fasta, handle, "fasta")
        subprocess.call(["cp","-R",cudir+"/"+name+"/sequence.fasta","/sadika/soft/TMPred-master"])
        subprocess.call(["./TMPred","sequence.fasta"],cwd="/sadika/soft/TMPred-master/")
        subprocess.call(["cp","/sadika/soft/TMPred-master/sequence.fasta.tmpred",cudir+"/"+name+"/sequence.fasta.tmpred"])
        helices=[]
        ttlscr=[]
        tmpr=0 
        files=[f for f in os.listdir(cudir+"/"+name+"/") if f.endswith("tmpred")]
        for fl in files:
          #print (fl)
          with open(cudir+"/"+name+"/"+fl,"r") as tool_fl:
            for tuples in tool_fl:
                tuples=tuples.strip()
                if "strong transmembrane helices" in tuples:
                          tuples=list(filter(None,tuples.split(" ")))
                          helices.append(tuples[0])
                          ttlscr.append(tuples[-1]) 
                elif "no probable" in tuples:
                          helices.append(0)
                          ttlscr.appennd(0)
                  
            if len(helices)>0:
               hel=int(max(helices))   
               if hel<=1:
                  tmpr=1
               else:
                  tmpr=0
                        
            else:
                         hel=0
                         tmpr=0
    
                          
        inter_fl.write(name1+"\t"+description+"\t"+str(hel)+"\n")
        final_fl.write(name1+"\t"+description+"\t"+str(tmpr)+"\n")   
         
