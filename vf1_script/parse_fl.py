import os,pandas as pd,sys
  



if "__main__"==__name__:
   file1=sys.argv[1]
   fastafl=sys.argv[2]
   cudir=os.getcwd() 
   fls=[f for f in os.listdir(cudir+"/"+fastafl) if f.endswith("Final_max_result.txt")] 
   print fls
   for i,f in enumerate(fls):
      if i ==0:
        df=pd.read_csv(cudir+"/"+fastafl+"/"+f,sep="\t")
      else:
        df1=pd.read_csv(cudir+"/"+fastafl+"/"+f,sep="\t")
        df=df.append(df1)
   df.to_csv(cudir+"/"+fastafl+"/"+fastafl+"_final_fl.txt",sep="\t") 
   
