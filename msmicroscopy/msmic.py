import pandas as pd
import os
import sys
import matplotlib.ticker as ticker
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from math import log
matplotlib.use("Agg")
from matplotlib import rc
rc("pdf", fonttype=42)
  
# Running instructions, if sys.argv too short.
if len(sys.argv) < 6:
    print("Correct form: python msmic.py training_set_identifiers training_set_datafile abbreviations_file query_data output_directory")
    exit(0)

output_folder = sys.argv[5]
inputdatafile = sys.argv[4]
trainingdata = sys.argv[2]
training_identifiers = sys.argv[1]
abbfile = sys.argv[3]

database = pd.read_csv(trainingdata,sep='\t')
identifiers = pd.read_csv(training_identifiers,sep='\t')
checkdata = pd.read_csv(inputdatafile)
checkdata = checkdata[['Bait','Prey','AvgSpec']]
database['Prey'] = database['Prey'].str.upper()
checkdata['Prey'] = checkdata['Prey'].str.upper()
if not os.path.isdir(output_folder): 
    os.makedirs(output_folder)
    
database['loc'] = database.apply(lambda x: identifiers[identifiers['Bait'] == x['Bait']].iloc[0]['Organelle '].lower().strip().capitalize(),axis=1)
database['uniqueforloc'] = database.apply(lambda x: len(database[database['Prey']==x['Prey']]['loc'].unique())==1 ,axis=1)
database = database[database['uniqueforloc']]

# The code below removes duplicate preys from localizations, and substitutes them with one entry with an average AvgSpec value from the two (or more) prior entries. 
nrows = []
done = set()
droprows = []
for index,row in database.iterrows(): 
    pdb = database[(database['loc']==row['loc']) & (database['Prey']==row['Prey'])]
    if pdb.shape[0] > 1: 
        donestring = pdb.iloc[0]['loc'] + '-' + pdb.iloc[0]['Prey']
        if donestring in done: continue
        done.add(donestring)
        
        ndic = {}
        for c in pdb.columns:
            if type(row[c])==float or type(row[c]) == int: 
                ndic.update({c: sum(pdb[c])/pdb.shape[0]})
            else: 
                ndic.update({c: pdb.iloc[0][c]})
        nrows.append(ndic)
        
        droprows.extend(list(pdb.index))
database = database.drop(droprows).append(nrows,ignore_index=True)


database['locnorm'] = database.apply(lambda row: 
                                    row['AvgSpec'] / database[(database['loc']==row['loc'])]['AvgSpec'].sum(),
                                    axis = 1
                                    )    
database.to_csv(os.path.join(output_folder, 'unique localization sets.csv'))
database.rename(columns={'loc': 'Localization'})['Bait Prey AvgSpec Localization'.split()].to_csv(os.path.join(output_folder, 'database_for_online.csv'),index=False)
with open(os.path.join(output_folder,'locsums.csv'),'w') as fil: 
    fil.write('localization,localization sum\n')
    for loc in database['loc'].unique(): 
        fil.write(loc + ',' + str(database[database['uniqueforloc'] & (database['loc']==loc)]['AvgSpec'].sum())+'\n')
        
locmat_unique = []
loc_index = sorted(database['loc'].unique())
locheads = sorted(checkdata['Bait'].unique())
loc_origins = ['localization','bait','prey','avgpsm']
loc_origin_data = []
for loc in loc_index: 
    nrow_uniq = []
    for bait in locheads: 
        bdf = checkdata[checkdata['Bait']==bait]
        bdf = bdf[bdf['Prey'].isin(database['Prey'].values)]
        locbdf = bdf[bdf['Prey'].isin(database[database['loc']==loc]['Prey'].values)]
        
        for _,row in locbdf.iterrows(): 
            loc_origin_data.append([loc, bait,row['Prey'],row['AvgSpec']])
        
        score_left_uni = locbdf['AvgSpec'].sum()
        score_left_uni = score_left_uni / bdf['AvgSpec'].sum()
        
        score_right_uni = database[database['loc']==loc]['AvgSpec'].sum()
        score_right_uni = score_right_uni / database['AvgSpec'].sum()
        score_uni = score_left_uni * score_right_uni
        
        nrow_uniq.append(score_uni)
    locmat_unique.append(nrow_uniq)
loc_raw = pd.DataFrame(data=locmat_unique,index = loc_index, columns = locheads)
loc_fin = pd.DataFrame(index=loc_index)
for c in loc_raw.columns: 
    loc_fin[c] = loc_raw[c]/loc_raw[c].max()
pd.DataFrame(columns = loc_origins,data=loc_origin_data).to_csv(os.path.join(output_folder, 'bait localization preys.csv'))
loc_fin.to_csv(os.path.join(output_folder,'Results.csv'))

"""
Here we just load abbreviations for localizations. 
"""
abbreviations = {}
with open(abbfile) as fil:
    for line in fil:
        line = line.strip("\n").split("\t")
        abbreviations.update({line[0].strip().lower().capitalize(): {"abb":line[1], "color": line[2]}})
        
"""
The last part is for generating figures. A multipage pdf file will be made, with one polar plot on each page.
"""

highcol = "#329932"
medcol = "#FFFF4C"
lowcol = "#FF9999"

pageindex = 0
figs = []
with PdfPages(os.path.join(output_folder, "Figures_" + output_folder + ".pdf")) as pdf:
    matplotlib.rcParams['pdf.fonttype'] = 42
    for bait in loc_fin.columns: 
        pageindex += 1
        fig = plt.figure()
        title = bait
        fig.suptitle(title,y=1.0)
        ax = fig.add_subplot(111, polar=True)
        
        keyorder = loc_fin.index
        data = loc_fin[bait].values
        N = len(data)
        theta = np.linspace(0.0+(np.pi/N), 2*np.pi+(np.pi/N), N, endpoint = False)
        radii = data
        width = 2*np.pi / N
        offset_angle = (360/N)/-2
        
        bars = ax.bar(theta, radii, width=width, bottom=0.0)#, alpha=0.3)
        for i, b in enumerate(bars):
            if radii[i] >= 0.75:
                b.set_color(highcol)
            elif radii[i] >= 0.5: 
                b.set_color(medcol)
            else:
                b.set_color(lowcol)
        temp = ax.set_rgrids([0.5, 0.75], angle=0)
        for i, t in enumerate(temp[0]):
            if i % 2 == 0:
                t.set_linewidth(0.5)
                t.set_linestyle("dotted")
            else:
                t.set_linewidth(0.75)
                t.set_linestyle("dashed")
        plt.setp(ax.spines["polar"], linewidth=1, linestyle="solid")
            
        ax.xaxis.set_major_locator(ticker.MultipleLocator((2*np.pi)/N))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator((2*np.pi)/(N*2)))
        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        
        thetaticks = np.arange(0, 360, 360/N)
        
        
        # TODO: implement this line style change again.
        #temp = ax.set_thetagrids(thetaticks)
        #print(temp)
        #for t in temp:
            #t.linestyle = "dotted"
        
        for tick in ax.xaxis.get_minor_ticks():
            tick.tick1line.set_markersize(0)
            tick.tick2line.set_markersize(0)
            tick.label1.set_horizontalalignment("center")
        minorticklabels = [""]
        for l in keyorder:
            #minorticklabels.append("")
            minorticklabels.append(abbreviations[l]["abb"])
        
        ax.set_xticklabels(minorticklabels, minor=True)
        ax.set_ylim([0,1])

        fig.tight_layout()
        pdf.savefig(fig)
        plt.cla()
        plt.clf()
