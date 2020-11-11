import os
import sys
from collections.abc import Iterable
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
    print("Correct form: python msmic.py training_set_identifiers training_set_datafile query_data abbreviations_file output_directory")
    exit(0)


training_set_file = sys.argv[1]
data_set_file = sys.argv[2]
data_only = sys.argv[3]
abbrevfile = sys.argv[4]
output_folder = sys.argv[5]

output_localization_info = True

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

training_set = {}
with open(training_set_file) as fil:
    next(fil)
    for line in fil:
        if len(line) > 0:
            line = line.strip("\n").split("\t")
            name = line[0].strip().lower()
            uniprot_id = line[1].strip().lower()
            organelle = line[2].strip().lower()
            size = line[3].strip().lower()
            crosstalk = line[4].strip().lower().split(", ")
            training_set.update({name: {"name": name, "uniprot": uniprot_id, "organelle": organelle, "crosstalk": crosstalk, "size": size}})
            
# Read in the database file
datafilelines = []
with open(data_set_file) as fil:
    datafilelines.append(next(fil))
    for line in fil:
        tline = line.lower().strip("\n").split("\t")
        if not tline[0] in training_set:
            continue
        else:
            datafilelines.append(line)

sepchr = '\t'
headers = datafilelines[0].lower().strip("\n").split(sepchr)
needed_headers = set(["bait", "prey", "avgspec"])    
dataset = {}     
for line in datafilelines[1:]:
    line = line.lower().strip("\n").split(sepchr)
    lindic = {}
    for i, h in enumerate(headers):
        if len(line) > i:
            if h in needed_headers:
                lindic.update({h: line[i]})
    if line[0] in dataset:
        dataset[line[0]].update({line[1]: lindic})
    else:
        dataset.update({line[0]: {line[1]: lindic}})
        
# Read in the experiment data file
datafilelines = []
if 'csv' == data_only.split('.')[-1]:
    sepchr = ',' 
with open(data_only) as fil:
    datafilelines.append(next(fil))
    for line in fil:
        datafilelines.append(line)
        
headers = datafilelines[0].lower().strip("\n").split(sepchr)
needed_headers = set(["bait", "prey", "avgspec"])    
        
for line in datafilelines[1:]:
    line = line.lower().strip("\n").split(sepchr)
    lindic = {}
    for i, h in enumerate(headers):
        if len(line) > i:
            if h in needed_headers:
                lindic.update({h: line[i]})
    if line[0] in dataset:
        dataset[line[0]].update({line[1]: lindic})
    else:
        dataset.update({line[0]: {line[1]: lindic}})
            
localizations = {}
testset = {}



"""
Here we just load abbreviations for localizations. 
"""
abbreviations = {}
with open(abbrevfile) as fil:
    for line in fil:
        line = line.strip("\n").split("\t")
        abbreviations.update({line[0].strip().lower(): {"abb":line[1], "color": line[2]}})




"""
Right now crosstalk is still implemented as a potential factor to consider, but it is not used afterwards anywhere. 

for each bait in dataset: check if bait is in training set. 
    If bait is in training set: 
        1. find the localization the bait is in the reference set for.
        2. for each prey protein of the bait: 
            1. if localization does not exist yet: create it
            2. else: 
                1. if prey is already in the localization database: add psm values, and average later. 
        3. add bait to localization baits
    else: 
        add bait and its dataset to testset
"""

for b_id in dataset.keys():
    if b_id in training_set.keys():
        org = training_set[b_id]["organelle"]
        cross = training_set[b_id]["crosstalk"]
        prey_ids = dataset[b_id].keys()
        for prey_id in prey_ids:
            if isinstance(dataset[b_id][prey_id], Iterable):
                if org not in localizations:
                    localizations.update({org: {"original_psm": 0.0, "preys": {}, "size": 0, "crosstalk": cross}})
                if prey_id not in localizations[org]["preys"]: 
                    localizations[org]["preys"].update({prey_id: []})
                localizations[org]["preys"][prey_id].append(float(dataset[b_id][prey_id]["avgspec"]))

                localizations[org]["size"] += 1
                localizations[org]["original_psm"] += float(dataset[b_id][prey_id]["avgspec"])
    else:
        testset.update({b_id: dataset[b_id]})
"""
keepmode defines, how duplicate preys in the SAME localization from different localization specific baits. 
options are:
1. avg: average the psm value
2. max: keep highest
3. min: keep lowest

"""        
keepmode = "avg" 
for org, orgdic in localizations.items(): 
    new_preydic = {}
    for prey_id, preyval_list in orgdic["preys"].items():   
        if keepmode == "avg": 
            new_preydic.update({prey_id: sum(preyval_list)/len(preyval_list)})
        elif keepmode == "max": 
            new_preydic.update({prey_id: max(preyval_list)})
        elif keepmode == "min": 
            new_preydic.update({prey_id: min(preyval_list)})
    orgdic["preys"] = new_preydic

# lockey will be used later to generate output data files for additional insights. Usually unnecessary, but useful for debugging.
lockey = []
for i in range(0, len(localizations.keys())*2):
    lockey.append("")

# Preytally is used to ensure each prey is unique, and not a part of multiple localizations. 
preytally = {}
for lkey, locdic in localizations.items():
    for pkey, prey in locdic["preys"].items():
        if pkey in preytally:
            preytally[pkey].append(lkey)
        else:
            preytally.update({pkey: [lkey]})

# here, we go over preytally, and add unique preys to authorized_preys.
authorized_preys = set()
for pid, lis in preytally.items():
    if len(lis) == 1:
        authorized_preys.add(pid)
        


# Here we update the localization prey dictionary, to account for possibly non-unique preys.      
i = 0  
for lkey, locdic in localizations.items():
    lockey[i*2] = (lkey + ": " + str(locdic["original_psm"]))
    i += 1
    newpreys = {}
    for pkey, prey in locdic["preys"].items():
        if pkey in authorized_preys:
            newpreys.update({pkey: prey})
    locdic["preys"] = newpreys
    locdic["size"] = len(newpreys.keys())
    locdic.update({"newpsm": 0.0})
    for pkey, prey in newpreys.items():
        locdic["newpsm"] += prey


"""
This is the main score calculator of the program. 
For each bait_id in testset: 
    lochits -dictionary is for storing localization hits
    Testset loc total psm value is used later to normalize localizations.
    bait_localization_preys will be later used to output debugging and detailed data.
"""
for bait_id in testset.keys():
    lochits = {}
    #testset[bait_id].update({"totalpsm": 0.0})
    testset[bait_id].update({"loc_total_psm": 0.0})
    bait_localization_preys = {}
    
    for prey_id in testset[bait_id].keys():
        if isinstance(testset[bait_id][prey_id], Iterable):
            # Check if prey is unique for localization
            if prey_id in authorized_preys:
                testset[bait_id]["loc_total_psm"] += float(testset[bait_id][prey_id]["avgspec"])
            for org in localizations.keys():
                # if prey is in reference set, it will be used for calculation
                if prey_id in localizations[org]["preys"]:
                    if org in lochits:
                        lochits[org]["psmhits"] += float(localizations[org]["preys"][prey_id])
                        lochits[org]["psmhits_test"] += float(testset[bait_id][prey_id]["avgspec"])
                    else:
                        lochits.update({org: {"psmhits": float(localizations[org]["preys"][prey_id]), "psmhits_test": float(testset[bait_id][prey_id]["avgspec"])}})
                    if org in bait_localization_preys:
                        bait_localization_preys[org].append(prey_id + "\t" + str(localizations[org]["preys"][prey_id]))
                    else:
                        bait_localization_preys.update({org: [prey_id + "\t" + str(localizations[org]["preys"][prey_id])]})

    # This loop calculates total psm values for each localization. 
    for o in lochits.keys():
        totalpsm = 0.0
        for f in localizations[o]["preys"].keys():
            totalpsm += float(localizations[o]["preys"][f])
        lochits[o].update({"totalpsm": totalpsm})
    
    # Here we calculate the scores for the bait for each localization
    locscores = []
    maxscore = 0.0
    for o in localizations.keys():
        if o in lochits:
            if testset[bait_id]["loc_total_psm"] > 0.0 and lochits[o]["totalpsm"] > 0.0:
                newscore = ((lochits[o]["psmhits_test"])/testset[bait_id]["loc_total_psm"]) * (lochits[o]["psmhits"]/lochits[o]["totalpsm"])
            else:
                newscore = 0.0
        else:
            newscore = 0.0
        newscore = newscore*100
        if newscore > maxscore:
            maxscore = newscore
        locscores.append([o, newscore])
    
    if maxscore > 0.0:
        for ls in locscores:
            ls[1] = ls[1]/maxscore
    testset[bait_id].update({"new_loc_score": locscores, "lochits": bait_localization_preys})#, "originals": original_locscores})

"""
Here we generate the main result output. 
"""
lines = [["Localization"]]
orkeys = list(localizations.keys())
orkeys.sort()
for o in orkeys:
    lines.append([o])

tkeys = list(testset.keys())
tkeys.sort()
for p in tkeys:
    nam = p
    lines[0].append(nam)
    for o in testset[p]["new_loc_score"]:
        found = False
        for l in lines:
            if l[0] == o[0]:
                found = True
                l.append("{:4.2f}".format(o[1]))
        if not found:
            l.append("0")
            l.append("0")
# The below is used to output localization information for use in debugging and generating reference database for R shiny server.
if output_localization_info:
    with open(os.path.join(output_folder, "Loc_info"), "w") as fil:
        fil.write("Size of training set: ")
        for i, o in enumerate(localizations.keys()):
            lockey[i*2+1] = (o + " unique: " + str(localizations[o]["newpsm"]))
        for l in localizations.keys():
            loc = localizations[l]["preys"]
            with open(os.path.join(output_folder, "unique_set_" + l + ".tsv"), "w") as fil2:
                for p in loc.keys():
                    fil2.write(p + "\t" + str(loc[p]) + "\n")
        fil.write("\n".join(lockey))
    with open(os.path.join(output_folder, "Localization_preys"), "w") as fil:
        prey_loc_fil = ["Bait\tLocalization\tPrey\tAvgSpec\n"]
        for p in tkeys:
            for locprey_key in testset[p]["lochits"].keys():
                for pr in testset[p]["lochits"][locprey_key]:
                    prey_loc_fil.append(p + "\t" + locprey_key + "\t" + pr + "\n")
        for line in prey_loc_fil:
            fil.write(line)
    if not os.path.isdir(os.path.join(output_folder, "for_R_server")): 
        os.makedirs(os.path.join(output_folder, "for_R_server"))
       
    reference_data_csv = ["Bait,Prey,AvgSpec,Localization"]
    bait_to_loc = {}
    for baitname, bdic in training_set.items(): 
        bait_to_loc.update({baitname: bdic["organelle"]})
    with open(data_set_file) as fil: 
        next(fil)
        for line in fil: 
            if len(line.strip()) < 2: continue
            bait, prey, psm = line.strip().split("\t")[:3]
            if bait.lower() not in bait_to_loc: continue
            if prey.lower() in authorized_preys: 
                reference_data_csv.append(",".join([bait, prey, psm, abbreviations[bait_to_loc[bait.lower()]]["abb"]]))
    with open(os.path.join(output_folder, "for_R_server", "reference_data.csv"), "w") as fil: 
        fil.write("\n".join(reference_data_csv))
        
        
with open(os.path.join(output_folder, "results" + ".tsv"), "w") as fil:
    for line in lines:
        fil.write("\t".join(line) + "\n")

"""
This part makes:
1. dataentries, which stores values for each localization for each bait
2. localizationNames, which stores localization names.

"""
dataentries = {}
localizationNames = set()
headers = "\t".join(lines[0]).upper().strip().split("\t")[1:]
for h in headers:
    dataentries.update({h:{}})
for line in lines[1:]:
    loc = line[0]
    line = line[1:]
    localizationNames.add(loc)
    for i, val in enumerate(line):
        val = float(val)
        dataentries[headers[i]].update({loc: val})


"""
This part is to sort localization keys
"""
temporary = sorted([(localization, abbreviations[localization]["abb"]) for localization in list(localizationNames)], key=lambda x: x[1], reverse=True) 
keyorder = [x[0] for x in temporary]
dkeys = list(dataentries.keys())
dkeys.sort()


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
    for indx, dk in enumerate(dkeys):
        pageindex += 1
        fig = plt.figure()
        title = dk
        fig.suptitle(title, y=1.0)
        d = dataentries[dk]
        data = []
        for k in keyorder:
            data.append(d[k])

        data = np.array(data)
        N = len(data)
        theta = np.linspace(0.0+(np.pi/N), 2*np.pi+(np.pi/N), N, endpoint = False)
        radii = data
        width = 2*np.pi / N

        offset_angle = (360/N)/-2

        ax = fig.add_subplot(111, polar=True)
        
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


