'''
[kamal@localhost Typanosoma_Species]$ python try.py 
Enter the filename:TbruceiTREU927.fasta
Enter the 2nd filename:TbruceiLister427.fasta
Enter the blast resultfile:Tbrucei_927_427.txt
'''
import pandas as pd,sys,os
#/home/rawal-lab/Desktop/ncbi_tcruzi_species/raw_data/clbrener_seq.fasta
def parse_fasta(file1,fastafl):
    """
     Parsing fasta sequences
    """
    print "fast file"+fastafl
    print file1
    with open(file1,"r") as filrd:
      with open(file1+"_proteinid_name.txt","w") as filewr: 
       line=[lines.strip() for lines in filrd if lines.startswith(">")]
       
       for row in line:
        if "is_pseudo" in row:
           row=row.split("|")
           prot=row[5].split("=")
           proteinid=row[0].strip()
           protein_name=prot[1]
           filewr.write((proteinid.strip(">").strip())+"\t"+protein_name+"\n")
        else:
          #print(row)
          proteinid=""
          protein_name=""
          if "RecName:" in row:
            row=row.split(";")
            #print row
            proteinid=(row[0].lstrip(">")).split("RecName: Full=")[0]
            protein_name=(row[0].lstrip(">")).split("RecName: Full=")[1]
            filewr.write(proteinid.strip(">")+"\t"+protein_name+"\n")
          else: 
            
           
           row=(row.lstrip(">")).split("[Trypanosoma")
           prot=row[0].split()
           proteinid=prot[0]
           protein_name=""
           for i in range(1,len(prot)):
             protein_name+=prot[i].strip(",")+" "
           if "[Trypanosoma" in protein_name:
                  #print 1
                  protein_name=protein_name.split("[Trypanosoma")
                  filewr.write(proteinid.strip(">")+"\t"+protein_name+"\n")
           else:
                  #print 2
                  filewr.write(proteinid.strip(">")+"\t"+protein_name+"\n")
          
           
def extract_blast_Result(blstfile,file1,file2,fastafl):
    df=pd.read_csv(file1+"_proteinid_name.txt",sep="\t",names=("Queryid","Query_Protein_name"))
    #print df
    #raw_input()
    qrcvls=[]
    df1=pd.read_csv(file2+"_proteinid_name.txt",sep="\t",names=("SSeqid","Subject_Protein_name"))
    #print df1
    #raw_input()
    with open(file2+"_blast_intermediate.txt","w") as blst:
     blst.write("Queryid\tSSeqid\tperc_identity\tevalue\tbit_score\talignment_length\tquery_length\tsubject_length\n")
     with open(blstfile,"r") as blstresult:
        for lines in blstresult:
           if lines.startswith("#"):
              continue
           else:
              blst.write(lines)
    df3=pd.read_csv(file2+"_blast_intermediate.txt",sep="\t")
    
    queylen=df3["query_length"].tolist()
    alignlen=df3["alignment_length"].tolist()
    for i in range(len(queylen)):
      qry_cov=(float(alignlen[i])/float(queylen[i]))
      #print qry_cov
      qrcvls.append(qry_cov)
    df3["Query_Coverage"]=qrcvls
    df3=df3.loc[df3['perc_identity']>30]
    df3=df3.loc[df3['Query_Coverage']>=0.70]
    print df3
   
    #print df3
    #raw_input()
    df3=df3.merge(df,how='left',left_on='Queryid',right_on='Queryid')
    df3=df3.merge(df1,how='left',left_on='SSeqid',right_on='SSeqid')          
    #print df3  
    df3.to_csv(file2+"_blast_result.txt",sep="\t",index=False)
    df4=pd.DataFrame()
    #df4=df3.groupby(['Queryid',"Query_Protein_name",'SSeqid',"Subject_Protein_name","evalue","bit_score"],as_index=False)['perc_identity'].max()
    idx=df3.groupby(['Queryid'])['perc_identity'].transform(max)==df3["perc_identity"]
    df4=df3[idx]
    df4=df4.drop_duplicates(subset='Queryid',keep='last')
    
    df4.to_csv(cudir+"/"+fastafl+"/"+fastafl+"_Final_max_result.txt",sep="\t",index=False)
    
if __name__=='__main__':
     
    cudir=os.getcwd()
    file1=sys.argv[1]
    file2=sys.argv[2]
    blstfile=sys.argv[3]
    fastafl=sys.argv[5]
    print fastafl
    parse_fasta(file1.strip(),fastafl)
    
    parse_fasta(file2.strip(),fastafl)
    filename=str(file1.split("/")).split(".")[0]
    fl1=filename.replace("[\'","")
    filename1=str(file2.split("/")).split(".")[0]
    fl2=filename1.replace("[\'","")
    
    
   
    extract_blast_Result(blstfile,file1,file2,fastafl)
    
               
