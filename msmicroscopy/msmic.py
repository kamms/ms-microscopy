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
	
# Running: 
if len(sys.argv) < 6:
	print("Correct form: python msmic.py training_set_identifiers training_set_datafile query_data abbreviations_file output_directory")
	exit(0)


training_set_file = sys.argv[1]
data_set_file = sys.argv[2]
data_only = sys.argv[3]
abbrevfile = sys.argv[4]
#uniprot_genename_file = sys.argv[5]
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

datafilelines = []
with open(data_set_file) as fil:
	datafilelines.append(next(fil))
	for line in fil:
		tline = line.lower().strip("\n").split("\t")
		if not tline[0] in training_set:
			continue
		else:
			datafilelines.append(line)
			
with open(data_only) as fil:
	next(fil)
	for line in fil:
		datafilelines.append(line)

dataset = {}
headers = datafilelines[0].lower().strip("\n").split("\t")
needed_headers = set(["bait", "prey", "avgspec"])	
for line in datafilelines[1:]:
	line = line.lower().strip("\n").split("\t")
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
Right now crosstalk is still implemented as a potential factor to consider, but it is not used afterwards anywhere. 
"""

for b_id in dataset.keys():
	if b_id in training_set.keys():
		org = training_set[b_id]["organelle"]
		cross = training_set[b_id]["crosstalk"]
		size = training_set[b_id]["size"]
		prey_ids = dataset[b_id].keys()
		for prey_id in prey_ids:
			if isinstance(dataset[b_id][prey_id], Iterable):
				if org in localizations:
					localizations[org]["preys"].update({prey_id: dataset[b_id][prey_id]})
					localizations[org]["size"] += 1
					localizations[org]["original_psm"] += float(dataset[b_id][prey_id]["avgspec"])
				else:
					localizations.update({org: {"original_psm": float(dataset[b_id][prey_id]["avgspec"]), 
					"preys": {prey_id: dataset[b_id][prey_id]}, "size": 1, "crosstalk": cross}})
		if org in localizations:
			if "baits" in localizations[org]:
				localizations[org]["baits"].append(b_id)
			else:
				localizations[org].update({"baits": [b_id]})
		else:
			localizations.update({org: {"baits": [b_id], "preys": {}}})
		if b_id not in localizations[org]["baits"]:
			localizations[org]["baits"].append(b_id)
			
	else:
		testset.update({b_id: dataset[b_id]})

scaledloc=False
if scaledloc:
	for l in localizations.keys():
		maxnum = 0.0
		for preyid, preyi in localizations[l]["preys"].items():
			if float(preyi["avgspec"]) > maxnum:
				maxnum=float(preyi["avgspec"])
		for preyid, preyi in localizations[l]["preys"].items():
			preyi["avgspec"] = str(float(preyi["avgspec"])/maxnum)

lockey = []
for i in range(0, len(localizations.keys())*2):
	lockey.append("")

newloc = {}
i = 0
eliminated = 0
tally = 0

preytally = {}
for lkey, locdic in localizations.items():
	for pkey, prey in locdic["preys"].items():
		if pkey in preytally:
			preytally[pkey].append(lkey)
		else:
			preytally.update({pkey: [lkey]})

authorized_preys = set()
for pid, lis in preytally.items():
	if len(lis) == 1:
		authorized_preys.add(pid)
for lkey, locdic in localizations.items():
	lockey[i*2] = (lkey + ": " + str(locdic["original_psm"]))
	i += 1
	newpreys = {}
	for pkey, prey in locdic["preys"].items():
		if pkey in authorized_preys:
			newpreys.update({pkey: prey})
	
	eliminated += len(locdic["preys"].keys()) - len(newpreys.keys())
	tally += len(newpreys.keys())
	locdic["preys"] = newpreys
	locdic["size"] = len(newpreys.keys())
	locdic.update({"newpsm": 0.0})
	for pkey, prey in newpreys.items():
		locdic["newpsm"] += float(prey["avgspec"])



for bait_id in testset.keys():
	lochits = {}
	testset[bait_id].update({"totalpsm": 0.0})
	testset[bait_id].update({"loc_total_psm": 0.0})
	bait_localization_preys = {}
	
	for prey_id in testset[bait_id].keys():
		if isinstance(testset[bait_id][prey_id], Iterable):
			testset[bait_id]["totalpsm"] += float(testset[bait_id][prey_id]["avgspec"])
			if prey_id in authorized_preys:
				testset[bait_id]["loc_total_psm"] += float(testset[bait_id][prey_id]["avgspec"])
			for org in localizations.keys():
				if prey_id in localizations[org]["preys"]:
					if org in lochits:
						lochits[org]["psmhits"] += float(localizations[org]["preys"][prey_id]["avgspec"])
						lochits[org]["psmhits_test"] += float(testset[bait_id][prey_id]["avgspec"])
					else:
						lochits.update({org: {"psmhits": float(localizations[org]["preys"][prey_id]["avgspec"]), "psmhits_test": float(testset[bait_id][prey_id]["avgspec"])}})
					if org in bait_localization_preys:
						bait_localization_preys[org].append(prey_id + "\t" + localizations[org]["preys"][prey_id]["avgspec"])
					else:
						bait_localization_preys.update({org: [prey_id + "\t" + localizations[org]["preys"][prey_id]["avgspec"]]})
				if prey_id in localizations[org]["baits"]:
					if org not in lochits:
						lochits.update({org: {"locmarkers": 1, "psmhits": 0.0, "psmhits_test": 0.0}})
					elif "locmarkers" not in lochits[org]:
						lochits[org].update({"locmarkers": 1})
					else:
						lochits[org]["locmarkers"] += 1
	for o in lochits.keys():
		totalpsm = 0.0
		for f in localizations[o]["preys"].keys():
			totalpsm += float(localizations[o]["preys"][f]["avgspec"])
		lochits[o].update({"totalpsm": totalpsm})
		
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
	original_locscores = []
	all_locscores = []
	
	if maxscore > 0.0:
		for ls in locscores:
			original_locscores.append([ls[0], ls[1]])
			all_locscores.append(ls[1])
			ls[1] = ls[1]/maxscore
	testset[bait_id].update({"new_loc_score": locscores, "originals": original_locscores, "lochits": bait_localization_preys})

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
with open(os.path.join(output_folder, "results" + ".tsv"), "w") as fil:
	for line in lines:
		fil.write("\t".join(line) + "\n")

dataentries = {}
localizationmaxes = {}
headers = "\t".join(lines[0]).upper().strip().split("\t")[1:]
for h in headers:
	dataentries.update({h:{}})
for line in lines[1:]:
	loc = line[0]
	line = line[1:]
	localizationmaxes.update({loc: 0.0})
	for i, val in enumerate(line):
		val = float(val)
		if val > localizationmaxes[loc]:
			localizationmaxes[loc] = val
		dataentries[headers[i]].update({loc: val})

abbreviations = {}
with open(abbrevfile) as fil:
	for line in fil:
		line = line.strip("\n").split("\t")
		abbreviations.update({line[0].strip().lower(): {"abb":line[1], "color": line[2]}})

        
keyorder = list(localizationmaxes.keys())
neworder = []
for k in keyorder:
	neworder.append([k, abbreviations[k]["abb"]])
neworder.sort(key=lambda x: x[1], reverse=True)
keyorder = []
for k in neworder:
	keyorder.append(k[0])
axes = []
for k in keyorder:
	axes.append(localizationmaxes[k])
dkeys = list(dataentries.keys())
dkeys.sort()


"""
The last part is for generating figures
"""

highcol = "#329932"
medcol = "#FFFF4C"
lowcol = "#FF9999"

pageindex = 0
figs = []
with PdfPages(os.path.join(output_folder, "Figures.pdf")) as pdf:
	matplotlib.rcParams['pdf.fonttype'] = 42
	for indx, dk in enumerate(dkeys):
		pageindex += 1
		fig = plt.figure()
		title = dk
		# ~ if dk in unip:
			# ~ title = unip[dk]
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
		temp = ax.set_thetagrids(thetaticks)# this was deprecated: , frac=1.1)
		for t in temp:
			t.linestyle = "dotted"
		
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


