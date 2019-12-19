import time,os,sys
import subprocess
from Bio import SeqIO

def blast_sequences(seqfl1,seqfl2,outfilenm):
 """
   Blast the two protomes sequences
 """ 
 subprocess.call(["./blastp","-query",seqfl1,"-subject",seqfl2,"-evalue","0.00001","-out",outfilenm,"-outfmt","7 qseqid sseqid pident evalue bitscore length qlen slen"],cwd="/sadika/ncbi-blast-2.7.1+/bin")

if __name__=='__main__':
  cudir=os.getcwd()
  print cudir
  start_time=time.time()
  file1=sys.argv[1]
  type(file1)
  fastafl=(file1.split("/")[-1])
  fastafl=(fastafl.split(".")[0])
   
  os.mkdir(fastafl)
  fls=[f for f in os.listdir(cudir+"/dataDeg") if f.endswith(".faa")]
  print fls

  for q,files in enumerate(fls):
     print files 
     file2=cudir+"/dataDeg/"+files
     file3=cudir+"/"+files+"_blast.txt"
     blast_sequences(file1,file2,file3)
     subprocess.call(["python","try.py",file1,file2,file3,str(q+1),fastafl],cwd="/home/sadika/Desktop/Random_seq/deg1")
  
  subprocess.call(["python","parse_fl.py",file1,fastafl],cwd=cudir)  
 
  with open(cudir+"/"+fastafl+"/"+fastafl+"_Human_result.txt","w") as finalgut:
   fasta_sequence = SeqIO.parse(open(file1),"fasta")
   for fasta in fasta_sequence:
           
        name1, sequence = fasta.id, str(fasta.seq)
        f=0
        with open(cudir+"/"+fastafl+"/"+fastafl+"_final_fl.txt","r") as gutfl:
         for lines in gutfl:
            lines=lines.strip()
            if name1 in lines:
              
                lines=lines.split("\t")
                print lines
                finalgut.write(name1+"\t"+str(1)+"\n")
                #raw_input()
                f=1
                break
         if f==0:       
               finalgut.write(name1+"\t"+str(0)+"\n")
             
            
           
            
  
       
  print time.time()-start_time
#/home/sadika/Desktop/database_tcruzi/proteins/deg-np-15.2/degaa-np.fasta
