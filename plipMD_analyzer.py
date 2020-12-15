import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from itertools import groupby
import math
import sys


def analyzer(file=None):
	print ('Loading results')
	table=pd.read_excel(file,index_col=[0])

	interaction_types=list(set(table['Type']))

	heatmap=pd.DataFrame(index=sorted(list(set(table['Residue']))), columns=sorted(interaction_types))

	print ('Initializing analysis and ploting results')

	for i in list(set(table['Residue'])):
		residue=sorted(list(table[table['Residue']==i]['Type']))
		groups=groupby(residue)
		for x in groups:
			heatmap.loc[i,x[0]]=len(list(x[1]))
	heatmap.to_excel('Heatmap_summary.xlsx')

	plot_heatmap(data=heatmap)
	get_frequency(table=table)

def plot_heatmap(data=None):

	fig, ax = plt.subplots(figsize=(10,2))
	
	ax=sns.heatmap(data.transpose().fillna(0),cmap='Blues',cbar_kws=dict(label='Frequency',shrink=1,orientation='vertical',spacing='uniform',pad=0.02))

	plt.title('Interactions Summary',size='22',weight='bold')
	plt.xlabel('Residue',fontsize=20,fontweight='bold')
	plt.xticks (rotation=90,fontsize=12)
	plt.yticks (fontsize=12)
	plt.tick_params ('both',width=2,labelsize=12)
	plt.savefig('HeatMap.png',dpi=300,format='png', bbox_inches = 'tight')

def get_frequency(table=None):
	time=table['Time'].drop_duplicates()
	index=table.index.drop_duplicates()
	data=pd.DataFrame()
	hbond=[]
	saltbridge=[]
	pication=[]
	hydroph=[]
	pistack=[]
	wb=[]
	halogen=[]

	for x in index: 
		data.loc[x,'hbond']=len([i for i in table.loc[x,'Type'] if i=='hbond'])
		data.loc[x,'waterbridge']=len([i for i in table.loc[x,'Type'] if i=='waterbridge'])
		data.loc[x,'saltbridge']=len([i for i in table.loc[x,'Type'] if i=='saltbridge'])
		data.loc[x,'pication']=len([i for i in table.loc[x,'Type'] if i=='pication'])
		data.loc[x,'hydroph']=len([i for i in table.loc[x,'Type'] if i=='hydroph_interaction'])
		data.loc[x,'pistack']=len([i for i in table.loc[x,'Type'] if i=='pistack'])
		data.loc[x,'halogenbond']=len([i for i in table.loc[x,'Type'] if i=='halogenbond'])

	data.to_excel('Frequency.xlsx')
	
	return plot_frequency(data=data)

def plot_frequency(data=None):
	
	for column in data.columns:
		
		plt.rcParams['axes.linewidth'] = 1.5
		fig, ax = plt.subplots(figsize=(10,5))

		y_values=[float('nan') if x==0 else x for x in data[column]]

		ax.scatter(data.index,y_values,c='navy',s=1, marker='.')

		plt.title (column,fontsize=24,fontweight='bold')
		plt.xlabel ('Frame',fontsize=20,fontweight='bold')
		plt.ylabel ('Frequency',fontsize=20,fontweight='bold')

		yint = range(math.floor(min(data[column])), math.ceil(max(data[column]))+1)

		plt.yticks(yint)

		plt.tick_params ('both',width=2,labelsize=14)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		plt.tight_layout()
		plt.savefig(column+'.png',dpi=300,format='png',bbox_inches='tight')

if __name__ == "__main__":
	print ('USAGE: python plipMD_analizer.py Results.xlsx')
	analyzer(sys.argv[1])