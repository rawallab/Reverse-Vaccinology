import os,pandas as pd,sys
  



if "__main__"==__name__:
   file1=sys.argv[1]
   cudir=os.getcwd() 
   fls=[f for f in os.listdir(cudir+"/dataHuman") if f.endswith("Final_max_result.txt")] 
   for i,f in enumerate(fls):
      if i ==0:
        df=pd.read_csv(cudir+"/dataHuman/"+f,sep="\t")
      else:
        df1=pd.read_csv(cudir+"/dataHuman/"+f,sep="\t")
        df=df.append(df1)
   df.to_csv(cudir+"/human_final_fl.txt",sep="\t") 
   
