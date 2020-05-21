import numpy as np
from bioscrape.lineage import LineageModel
from bioscrape.lineage import LineageVolumeSplitter
from bioscrape.lineage import py_SimulateInteractingCellLineage
from bioscrape.lineage import py_SimulateCellLineage
from bioscrape.lineage import py_SimulateSingleCell
from bioscrape.lineage import LineageCSimInterface
from bioscrape.lineage import py_PropagateCells
from bioscrape.lineage import py_SingleCellLineage, py_PropagateInteractingCells
import time as pytime

k1 = 1.1111
kgrow = 33.33
Kgrow = 20.2020
kdeath = 6.06
Kdeath = 5.5050
g = .020202
d = 80.808
species = ["S", "X"]
rxn1 = [[], ["S"], "massaction", {"k":k1}]
rxn2 = [["S"], [], "massaction", {"k":d}]
rxns = [rxn1, rxn2]
x0 = {"S": 0, "X":100}

#Instantiate Model
print("Instantiating Model")
M = LineageModel(species = species, reactions = rxns, initial_condition_dict = x0)
vsplit = LineageVolumeSplitter(M)
#M.create_division_event("division", {}, "massaction", {"k":.1, "species":"S"}, vsplit)
M.create_division_rule("deltaV", {"threshold":1.0}, vsplit)
#M.create_death_event("death", {}, "hillpositive", {"k":kdeath, "s1":"S", "n":2, "K":Kdeath})
M.create_volume_event("linear volume", {"growth_rate":g}, 
	"hillnegative", {"k":kgrow, "s1":"S", "n":2, "K":Kgrow})
M.py_initialize()


global_sync_period = .5
N = 1
sum_list = []
#lineage = 

lineage = None
cell_state_samples = None
single_cell_states = None
result = None
lineage_list = None
cell_state_sample_list = None
for i in [1]:
	maxtime = i*10
	dt = 0.01
	timepoints = np.arange(0, maxtime+dt, dt)
	print("Beginning Simulation", i, "for", maxtime)
	#interface = LineageCSimInterface(M)
	
	s = pytime.clock()
	#lineage = py_SimulateCellLineage(timepoints, Model = M, initial_cell_states = N)
	#cell_state_samples, sample_times = py_PropagateCells(timepoints, Model = M, initial_cell_states = N, sample_times = 10)
	#single_cell_states = py_SingleCellLineage(timepoints, Model = M)
	#lineage_list, global_results, simulator = py_SimulateInteractingCellLineage(timepoints, global_sync_period, model_list = [M],initial_cell_states = [N], global_species = ["S"], global_volume = 100, average_dist_threshold = 10.0)
	cell_state_sample_list, global_species_results, sample_times, simulator = py_PropagateInteractingCells(timepoints, global_sync_period, sample_times = 5, model_list = [M],initial_cell_states = [N], global_species = ["S"], global_volume = 0, average_dist_threshold = 10.0)
	#result = py_SimulateSingleCell(timepoints[10:], Model = M)

	e = pytime.clock()
	print("Simulation", i, "complete in", e-s, "s")

	if i > 0:
		if lineage_list is not None:
			lineage = lineage_list[0]
		if lineage is not None:
			print("Building Tree")
			sch_tree = [[]]
			sch_tree_length = 1
			for i in range(lineage.py_size()):
				sch = lineage.py_get_schnitz(i)
				if sch.py_get_parent() == None:
					sch_tree[0].append(sch)
				else:
					for j in range(len(sch_tree)):
						parent = sch.py_get_parent()
						if parent in sch_tree[j]:
							if len(sch_tree)<= j+1:
								sch_tree.append([])
								sch_tree_length += 1
							sch_tree[j+1].append(sch)
			counts = [len(sch_list) for sch_list in sch_tree]
			print("counts", counts)
			sum_list.append((maxtime, sum(counts), e-s))

		if cell_state_sample_list is not None:
			cell_state_samples = [s[0] for s in cell_state_sample_list]
		if cell_state_samples is not None:
			sum_list.append((maxtime, sum([len(l) for l in cell_state_samples]), e-s))
#raise ValueError()
import pylab as plt
print("plotting")
if len(sum_list) > 1:
	plt.figure()
	plt.subplot(121)
	plt.plot([e[0] for e in sum_list], [e[1] for e in sum_list])
	plt.xlabel("max simulation time")
	plt.ylabel("Total Cells Simulated")

	plt.subplot(122)
	plt.plot([e[1] for e in sum_list], [e[2] for e in sum_list])
	plt.xlabel("Total Cells Simulated (Lineage)\nOR\nFinal Cells Returned (Propogate)")
	plt.ylabel("CPU runtime (s)")

if result is not None:
	single_cell_states = result

if single_cell_states is not None:
	plt.figure()
	
	plt.subplot(131)
	plt.title("volume")
	plt.plot(single_cell_states["time"], single_cell_states["volume"])
	
	import pandas as pd
	#with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
	#	print("**volume**\n", single_cell_states["volume"])
	#	print("**time**\n", single_cell_states["time"])

	plt.subplot(132)
	plt.title("X")
	plt.plot(single_cell_states["time"], single_cell_states["X"])

	plt.subplot(133)
	plt.title("S")
	plt.plot(single_cell_states["time"], single_cell_states["S"])


if cell_state_samples is not None:
	plt.figure()
	print("len cell_state_samples", len(cell_state_samples), [len(s) for s in cell_state_samples])
	ax1, ax2, ax3 = plt.subplot(131), plt.subplot(132), plt.subplot(133)

	for i in range(len(sample_times)-1, -1, -1):

		volume = [cs.py_get_volume() for cs in cell_state_samples[i]]
		S = [cs.py_get_state()[0] for cs in cell_state_samples[i]]
		X = [cs.py_get_state()[1] for cs in cell_state_samples[i]]

		time = sample_times[i]

		plt.sca(ax1)
		plt.title("volume histogram")
		plt.hist(volume, log = True, label = "Sample Time = "+str(i)+" Cells="+str(len(volume)))
		if i == 0: plt.legend()

		plt.sca(ax2)
		plt.title("X histogram")
		plt.hist(X, log = True, label = "Sample Time = "+str(i)+" Cells="+str(len(X)))
		if i == 0: plt.legend()

		plt.sca(ax3)
		plt.title("S histogram")
		plt.hist(S, log = True, label = "Sample Time = "+str(i)+" Cells="+str(len(S)))
		if i == 0: plt.legend()

if lineage is not None:
	color_list = []
	for i in range(sch_tree_length):
		color_list.append((i/sch_tree_length, 0, 1.-i/sch_tree_length))

	import pylab as plt
	plt.figure(figsize = (10, 10))
	

	plt.subplot(411)
	plt.title(r"$\emptyset \leftrightarrow S$    $P(Grow) = k \frac{1}{S^2+400}$")

	plt.plot(range(len(counts)), counts)
	plt.ylabel("Cell Count (total ="+str(sum(counts))+")")
	plt.xlabel("Generation")

	print("sch_tree_length", sch_tree_length)
	plt.subplot(412)
	plt.ylabel("S per Cell")
	for i in range(sch_tree_length):
		for sch in sch_tree[i]:
			df = sch.py_get_dataframe(Model = M)
			plt.plot(df["time"], df["S"], color = color_list[i])


	plt.subplot(413)

	plt.ylabel("X per cell")
	totalX = np.zeros(len(timepoints))


	for i in range(sch_tree_length):
		for sch in sch_tree[i]:
			df = sch.py_get_dataframe(Model = M)
			start_ind = np.where(timepoints >= df["time"][0])
			start_ind = start_ind[0][0]
			end_ind = np.where(timepoints >= df["time"][len(df["time"])-1])[0][0]

			plt.plot(df["time"], df["X"], color = color_list[i])
			plt.plot(df["time"][len(df["time"])-1], df["X"][len(df["time"])-1], "x", color = color_list[i])
			plt.plot(df["time"][0], df["X"][0], "o", color = color_list[i])


			#totalX[start_ind:end_ind+1] += df["X"][:len(df["X"])]

	#plt.plot(timepoints, totalX, "--", color = "black", label = "total X")

	plt.subplot(414)
	plt.ylabel("Volume (of each cell)")
	for i in range(sch_tree_length):
	    for sch in sch_tree[i]:
	        df = sch.py_get_dataframe(Model = M)
	        plt.plot(df["time"], df["volume"], color = color_list[i])

plt.show()
