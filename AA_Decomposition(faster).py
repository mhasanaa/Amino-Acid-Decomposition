# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 21:13:08 2023

@author: Hasan
"""
import pandas as pd
import time
from datetime import datetime
import sys
sys.path.insert(1, 'C:/Users/Hasan/HKUST/Haibin SU - group - covid19/SingleSiteAnalysis')
import ImportantFunc as Imp
start_time = time.process_time()
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt

def CorrectMut(df,col):
    def Fix144(x):
        if ('144-' and '146-') not in x:
            x = x.replace('Y145-', 'Y144-')
        return x
    df[col] = df[col].apply(Fix144)
    Initial = ['E156G;F157-;R158-','G142D;V143-;Y144-;Y145-']
    Final = ['E156-;F157-;R158G','G142-;V143-;Y144-;Y145D']
    for i in range(len(Initial)):
        df[col] = df[col].str.replace(Initial[i], Final[i])
    return df

## Reading the input file
start, end = 12, 39
Thres = 5
StrMatch = '.xlsx'
path = 'C:/Users/Hasan/OneDrive - HKUST Connect/Summer Project 2023/Consensus seq/Correct Input/InputVar/Mainfocus/'
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
onlyTSV = [f for f in onlyfiles if StrMatch in f]
## Dictionary for each AA and domains
AA = ['A','I','L','V','F','W','Y','C','M','D','E','N','Q','S','T','H','K','R','G','P','-']
dict_AA = {}
for aa in AA:
    dict_AA[aa] = []
dfAA = pd.DataFrame({'AA':AA})

DomName = ['NTD', 'N1', 'N2', 'N3',
           'RBD', 'RBM', 'CTD1-2', 'S2']
DomRange = [[1,324],[14,20],[140,158],[245,264],
            [325,540],[437,508],[541,699],[700,1273]]
Domain_columns =['MonthIndex']
Domain_columns.extend(AA)
# for Variant in onlyTSV:
## Uncomment this line to do calculation for single strain
Variant ='B.1.1.7(Alpha)'+StrMatch
Variants = ['B.1.1.7(Alpha)','BA.2']
MonthRange = [[12, 20],[25, 39]]
for count,Var0 in enumerate(Variants): # 
    Variant = Var0+StrMatch
    ## Reading the raw input and the variant reference
    dfInp = pd.read_excel(path+Variant,sheet_name='True')
    dfVarRef = pd.read_excel("InputVariantsLMCM.xlsx",sheet_name='VOCI')
    dfVarSel = dfVarRef.loc[dfVarRef['class']==Variant.replace(StrMatch, "")]
    dfVarSel['Var_Marker'] = [1]*len(dfVarSel)
    dfVarSel = dfVarSel[['Mutation','Var_Marker']]
    ## Preparing 
    dfInp = dfInp[['MonthIndex','mutation info|insertion info','class']].copy(deep=True)
    dfInp = dfInp.loc[dfInp['class']!='others'].reset_index(drop=True)
    dfInp[['mutation info','insertion info']] = dfInp['mutation info|insertion info'].str.split('|', expand=True)
    dfInp = dfInp.drop(['class','insertion info'], axis=1)
    ## Correcting Mutation
    dfInp0 = CorrectMut(dfInp,'mutation info')
    ## Iteration over month
    dict_Domain = {}
    for domain in DomName:
        dict_Domain[domain] = pd.DataFrame(columns = Domain_columns,dtype='float64')
        
    Month = MonthRange[count]
    start = Month[0]
    end = Month[1]
    for i,z in enumerate(range(start,end+1)):
        dfMonthly = dfInp0.loc[dfInp0['MonthIndex']==z].reset_index(drop=True)
        ## Selecting Mutation above the threshold
        MutList = [mut  for i in range(len(dfMonthly)) for mut in dfMonthly.at[i,'mutation info'].split(';')]
        print('Variant =',Var0,'Month =',z,' # Seqs:',len(dfMonthly))
        ## Creating DataFrame for Mutation Distribution
        MutList = Imp.count_dups(MutList)
        dfRef = pd.DataFrame({'Mutation':MutList[0],'Count':MutList[1]})
        ## Selecting the important mutation for each months
        dfRef = dfRef.loc[dfRef['Count']>=Thres].reset_index(drop=True)
        dfRef = dfRef[['Mutation','Count']].copy(deep=True)
        dfRef['AA'] = [x[-1] for x in dfRef['Mutation']]
        dfRef['Pos'] = dfRef['Mutation'].astype(str).str.extractall('(\d+)').unstack().fillna('').sum(axis=1).astype(int)
        dfRef = pd.merge(dfRef, dfVarSel,how='left',on='Mutation')
        dfRef = dfRef.loc[dfRef['Var_Marker']!=1].reset_index(drop=True)
        dfRef = dfRef.drop(['Var_Marker','Mutation'], axis=1)
        dfRef.sort_values('Pos',inplace = True, ascending=[True])
        ## Start the decomposition process: Mutation Profiling, splitting  into each domains
        ## Iteration over domain
        for j in range(len(DomRange)):
            dfMutDomain = dfRef.loc[(dfRef['Pos'] >= DomRange[j][0])&(dfRef['Pos'] <= DomRange[j][1])].reset_index(drop=True)
            dfMutDomain = dfMutDomain.drop(['Pos'], axis=1)
            dfMutDomain = dfMutDomain.groupby('AA').sum().reset_index(level=['AA'])
            dfAA_Mut = pd.merge(dfAA, dfMutDomain,how='left',on='AA').fillna(0)
            AA_Count = dfAA_Mut['Count'].tolist()
            ## Append to each domains df
            row = [z]
            row.extend(AA_Count)
            dict_Domain[DomName[j]].loc[i] = row
        print("Current time =", datetime.now().strftime("%H:%M:%S"))
        print('TIME TAKEN: ' + str(time.process_time() - start_time) + 's\n')
    ## Summing up all the same month data
    for count_dom,domain in enumerate(DomName):
        df = dict_Domain[domain]
        df['H_Aliphatic'] = df['A']+df['I']+df['L']+df['V']
        df['H_Aromatic'] = +df['F']+df['W']+df['Y']
        df['Polar'] = df['N']+df['Q']+df['S']+df['T']
        df['Non_Polar'] = df['C']+df['M']
        df['Positive'] = df['H']+df['K']+df['R']
        df['Negative'] = df['D']+df['E']
        df['Special'] = df['G']+df['P']
        ## Plotting 100% stacked plot for each domain
        dfPlot = df[['H_Aliphatic','H_Aromatic','Polar','Non_Polar',
                     'Positive','Negative','Special','-']]
        dfPlot = dfPlot.divide(dfPlot.sum(axis=1), axis=0).fillna(0)
        dfPlot['MonthIndex'] = df['MonthIndex']
        dfPlot = dfPlot.set_index('MonthIndex').rename(columns={'-':'Deletion'})
        ax = dfPlot.plot.area(color={'H_Aliphatic':'darkorange','H_Aromatic':'lime','Polar':'magenta','Non_Polar':'darkviolet',
                                     'Positive':'red','Negative':'blue','Special':'sienna','Deletion':'black'})
        ax=plt.gca()
        ax.set_xlim(start,end)
        ax.set_ylim(0,1)
        Title = '{}_{}_{}-{}'.format(Var0,domain,DomRange[count_dom][0],DomRange[count_dom][1])
        plt.title(Title)
        plt.savefig(Title+'.png',dpi=600, bbox_inches = 'tight')
    ## Print out all the results
    # FileProcessing = 'AA-Decomposition_Thres-%s-%s.xlsx'%(str(Thres),Variant)
    # with pd.ExcelWriter(FileProcessing) as writer:
    #     for domain in DomName:
    #         dict_Domain[domain].to_excel(writer, sheet_name=domain, index=None)
    # print("Current time =", datetime.now().strftime("%H:%M:%S"))
    # print('TIME TAKEN: ' + str(time.process_time() - start_time) + 's\n')