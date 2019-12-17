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
    infl=input("Enter the file name") 
    filename=infl.split("\t")[-1]
    name=filename.split(".")[0]
    '''
    dirpath = os.path.join(cudir+"/"+name)
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
    os.mkdir(cudir+"/"+name)
    subprocess.call(["cp",infl,cudir+"/input_file/"]) 
    '''
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
    #subprocess.call(["cp","-R",cudir+"/"+name+"/"+name+".fasta","/sadika/soft/psortb_commandline_docker"])
    #subprocess.call(["perl","psortb","-p","-i",cudir+"/"+name+"/"+name+".fasta","-o","long","-r",cudir+"/"+name+"/"],cwd="/sadika/soft/psortb_commandline_docker")
    subprocess.call(["perl","psortb","-n","-i",cudir+"/"+name+"/"+name+".fasta","-o","long","-r",cudir+"/"+name+"/"],cwd="/sadika/soft/psortb_commandline_docker")
    
        
