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

def batch_iterator(iterator, batch_size):

    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch





  


def blast_seq(cudir,filename,name):
  fileobj=open(cudir+"/"+name+"/tcruzi_allergy.csv","w")
  fileobj1=open(cudir+"/"+name+"/tcruzi_allergy_output.csv","w")
  fileobj1.write("qseqid\tsseqid\tqlen\tslen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\talign_length\t%_identity\tsstrand\n")
  subprocess.call(["./blastp","-query",cudir+"/"+name+"/"+filename,"-db","/sadika/ncbi-blast-2.7.1+/allergy_protein/Allergen_protein_sequence_1-2088.fasta","-evalue","0.00001","-outfmt","7 qseqid  sseqid qlen  slen qstart qend sstart send evalue bitscore length pident sstrand "],cwd="/sadika/ncbi-blast-2.7.1+/bin",stdout=fileobj)
  fileobj.close()
  with open(cudir+"/"+name+"/tcruzi_allergy.csv","r") as flrd:
     for lines in flrd:
       if lines.startswith("#"):
          continue
       else:
          fileobj1.write(lines)
  fileobj1.close()


def allergentool(cudir,filename,name):
   '''
    run protein for finding out to be allergenic or not
   '''
   blast_seq(cudir,filename,name)
   fileop=open(cudir+"/"+name+"/"+name+".fasta").readlines()
   file_result=open(cudir+"/"+name+"/"+name+"_allergy_out.txt","w")
   protein_ls=[]
   for lines in fileop:
      if lines.startswith(">"):
        lines=lines.strip()
        lines=lines.strip(">")
        protein_ls.append(lines)
   fileop1=open(cudir+"/"+name+"/tcruzi_allergy.csv").readlines()
   for prot in protein_ls:
     prot=prot.split(" ")
     alevalue=[]
     perid=[]
     for i,lines in enumerate(fileop1):
       if i==0:
         continue
       else:
        if lines.startswith(prot[0]):
         lines=lines.split("\t")
         ###print float(lines[8])
         alevalue.append(float(lines[8]))
         perid.append(float(lines[-2]))  
     if len(perid)>0:   
       if max(perid)>=70:
         file_result.write(prot[0]+"\tStrongly Allergen\n")
        
       elif max(perid)>=50 and max(perid)<70:
        if min(alevalue)<float(0.0000001):
           file_result.write(prot[0]+"\tnot likely to be Allergen\n")
           
        else:
           file_result.write(prot[0]+"\t likely to be Allergen\n")
           
       elif max(perid)<50: 
         file_result.write(prot[0]+"\tNot likely to be Allergen\n") 
          
     else:
       file_result.write(prot[0]+"\tNot Allergen\n")   
         
    

if __name__=="__main__":
    cudir=os.getcwd()
    fls=[f for f in os.listdir("/home/sadika/Desktop/tcruzi_files_strain/") if f.endswith(".fasta")]
    for files in fls:
       print (files)
       infl="/home/sadika/Desktop/tcruzi_files_strain/"+files
       inpath = os.path.join(cudir+"/input_file")
       if os.path.exists(inpath) and os.path.isdir(inpath):
           shutil.rmtree(inpath)
       os.mkdir("input_file")
       filename=infl.split("/")[-1]
       name=filename.split(".")[0]
       dirpath = os.path.join(cudir+"/"+name)
       if os.path.exists(dirpath) and os.path.isdir(dirpath):
          shutil.rmtree(dirpath)
       os.mkdir(name)
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
       with open(name+"_Final_Resultant_values.txt","a") as inter_fl: 
                inter_fl.write("Sequence_id\tDescription\tallergenval\n")
                
       for fasta in fasta_sequence:
          with open(name+"_Binary_File.txt","a") as final_fl:
              with open(name+"_values_file.txt","a") as inter_fl:
                 name1, sequence = fasta.id, str(fasta.seq)
                 description=fasta.description
                 with open(cudir+"/"+name+"/"+"sequence.fasta", "w") as handle:
                    count = SeqIO.write(fasta, handle, "fasta")
                
                
                 p8=threading.Thread(target=allergentool,args=(cudir,"sequence.fasta",name))
                
                 p8.start()
                
                 p8.join()  
                
                 files=cudir+"/"+name+"/sequence.fasta"
        
                
                 allergen=0 
                 allergenval="" 
                
                 files=[f for f in os.listdir(cudir+"/"+name+"/") if f.endswith(".txt") or f.endswith(".tmpred")]
                 for fl in files:
                     with open(cudir+"/"+name+"/"+fl,"r") as tool_fl:
                            if fl.endswith("allergy_out.txt"):

                                for tuples in tool_fl:
                                    tuples=tuples.strip()
                                    if tuples.startswith(name1):
                                        tuples=tuples.split("\t")
                                        if tuples[-1]=="Not likely to be Allergen":
                                            allergen=1
                                            allergenval="Not likely to be Allergen"

                                        elif tuples[-1]=="Not Allergen":
                                            allergen=1
                                            allergenval="Not Allergen"

                                        elif tuples[-1]=="Allergen":
                                            allergen=0
                                            allergenval="Allergen"

                                        elif tuples[-1]=="Strongly likely to be Allergen": 
                                            allergen=0 
                                            allergenval="Strongly likely to be Allergen" 

                           
                
                 inter_fl.write(name1+"\t"+description+"\t"+allergenval+"\n")
                 final_fl.write(name1+"\t"+description+"\t"+str(allergen)+"\n")
                 fls=glob.glob(cudir+'/'+name+"/*")
                 for fl in fls:
                    if fl.endswith(".fasta"):
                        continue
                    else:  
                        os.remove(fl)  
      
