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




lineage_list = None
cell_state_sample_list = None

#M1 create X1 which is global
#M2 create X2 at a rate proprotional to X1
#As X2 accumulate, M1 stops growing

species = ["X1", "X2", "Y"]
x0 = {"X1": 0, "X2":0, "Y":100}


k1 = 11.1111
g = .3333
kgrow = 3.3333
Kgrow = 33.333
rxn1 = [[], ["X1"], "massaction", {"k":k1}]
M1 = LineageModel(species = species, reactions = [rxn1], initial_condition_dict = x0)
vsplit1 = LineageVolumeSplitter(M1)
M1.create_division_rule("deltaV", {"threshold":1.0}, vsplit1)
M1.create_volume_event("linear volume", {"growth_rate":g}, 
	"hillnegative", {"k":kgrow, "s1":"X2", "n":2, "K":Kgrow})
M1.py_initialize()

k2 = 22.2222
K2 = 22.222
rxn2 = [[], ["X2"], "hillpositive", {"k":k2, "s1":"X1","n":2, "K":K2}]
M2 = LineageModel(species = species, reactions = [rxn2], initial_condition_dict = x0)
vsplit2 = LineageVolumeSplitter(M2)
M2.create_division_rule("deltaV", {"threshold":1.0}, vsplit2)
M2.create_volume_event("linear volume", {"growth_rate":g/10}, "massaction", {"k":kgrow, "species":""})
M2.py_initialize()


maxtime = 10
dt = 0.01
global_sync_period = .5
global_volume = 1000

timepoints = np.arange(0, maxtime+dt, dt)

average_dist_threshold = 2.
global_species = ["X1", "X2"]
model_list = [M1, M2]
initial_cell_counts = [5, 5]
print("starting simulation")

#cell_state_sample_list, global_species_results, sample_times, simulator = py_PropagateInteractingCells(timepoints, global_sync_period, sample_times = 5, model_list = model_list,initial_cell_states = initial_cell_counts, global_species = global_species, global_volume = global_volume, average_dist_threshold = average_dist_threshold)
lineage_list, global_results, simulator = py_SimulateInteractingCellLineage(timepoints, global_sync_period, model_list = model_list,initial_cell_states = initial_cell_counts, global_species = global_species, global_volume = global_volume, average_dist_threshold = average_dist_threshold)
	#print("lineage_list", lineage_list)

print("simulation complete")

print("plotting")
import pylab as plt
if cell_state_sample_list is not None:
	print("cell state sample list", cell_state_sample_list)
	ax1, ax2, ax3 = plt.subplot(231), plt.subplot(232), plt.subplot(233)
	ax4, ax5, ax6 = plt.subplot(234), plt.subplot(235), plt.subplot(236)
	axes = [ax1, ax2, ax3, ax4, ax5, ax6]
	for i in range(len(sample_times)-1, -1, -1):
		for cell_type in range(2):
			print("i, cell_type = ", i, cell_type)
			#df = cell_state_sample_list[i][cell_type]
			volume = [cs.py_get_volume() for cs in cell_state_sample_list[i][cell_type]]
			x1 = [cs.py_get_state()[1] for cs in cell_state_sample_list[i][cell_type]]
			x2 = [cs.py_get_state()[2] for cs in cell_state_sample_list[i][cell_type]]

			time = sample_times[i]

			plt.sca(axes[cell_type*3])
			plt.ylabel("Cell Type "+str(cell_type+1))
			plt.title("volume histogram")
			plt.hist(volume, log = True, label = "Sample Time = "+str(i)+" Cells="+str(len(volume)))
			if i == 0: plt.legend()

			plt.sca(axes[cell_type*3+1])
			plt.title("X1 histogram")
			plt.hist(x1, log = True, label = "Sample Time = "+str(i)+" Cells="+str(len(x1)))
			if i == 0: plt.legend()

			plt.sca(axes[cell_type*3+2])
			plt.title("X2 histogram")
			plt.hist(x2, log = True, label = "Sample Time = "+str(i)+" Cells="+str(len(x2)))
			if i == 0: plt.legend()

if lineage_list is not None:
	print("lineage_list", lineage_list, [l.py_size() for l in lineage_list])
	plt.figure(figsize = (10, 10))
	axes = [plt.subplot(4, 2, 1), plt.subplot(4, 2, 3), plt.subplot(4, 2, 5), plt.subplot(4, 2, 7), 
		plt.subplot(4, 2, 2), plt.subplot(4, 2, 4), plt.subplot(4, 2, 6), plt.subplot(4, 2, 8)]

	for lin_ind in range(len(lineage_list)):

		lineage = lineage_list[lin_ind]

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
		color_list = []
		for i in range(sch_tree_length):
			color_list.append((i/sch_tree_length, 0, 1.-i/sch_tree_length))


		plt.sca(axes[4*lin_ind])
		plt.ylabel("Y per Cell")
		plt.title("Cell Type "+str(lin_ind))
		for i in range(sch_tree_length):
			for sch in sch_tree[i]:
				df = sch.py_get_dataframe(Model = model_list[lin_ind])
				plt.plot(df["time"], df["Y"], color = color_list[i])

		plt.sca(axes[4*lin_ind+1])
		plt.ylabel("X1 per Cell")
		for i in range(sch_tree_length):
			for sch in sch_tree[i]:
				df = sch.py_get_dataframe(Model = model_list[lin_ind])
				plt.plot(df["time"], df["X1"], color = color_list[i])

		plt.sca(axes[4*lin_ind+2])
		plt.ylabel("X2 per cell")
		totalX = np.zeros(len(timepoints))
		for i in range(sch_tree_length):
			for sch in sch_tree[i]:
				df = sch.py_get_dataframe(Model = model_list[lin_ind])
				start_ind = np.where(timepoints >= df["time"][0])
				start_ind = start_ind[0][0]
				end_ind = np.where(timepoints >= df["time"][len(df["time"])-1])[0][0]
				plt.plot(df["time"], df["X2"], color = color_list[i])

		plt.sca(axes[4*lin_ind+3])
		plt.ylabel("Volume (of each cell)")
		for i in range(sch_tree_length):
		    for sch in sch_tree[i]:
		        df = sch.py_get_dataframe(Model = model_list[lin_ind])
		        plt.plot(df["time"], df["volume"], color = color_list[i])
	axes = []

plt.show()