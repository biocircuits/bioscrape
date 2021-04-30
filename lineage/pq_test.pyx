# Scratchspace for testing priority_queue
from libcpp.queue cimport priority_queue
from libcpp.utility cimport pair
from cython.operator import dereference
from libcpp cimport bool
import numpy as np

from bioscrape.random cimport choose



cdef class LineageVolumeCellState():
	cdef double initial_volume #Stores the birth Volume
	cdef double initial_time #Stores the time the Cell was "born"
	#divided = -1: Cell Not divided
	#divided E [0, num_division_rules): DivisionRule divided caused the cell to divide
	#divided E [num_division_rules, num_division_rules + num_division_events]: Division Event divided-num_division_rules caused the cell to divide
	cdef int divided
	#dead = -1: Cell Not dead
	#dead E [0, num_death_rules): DeathRule divided caused the cell to die
	#dead E [num_death_rules, num_death_rules + num_death_events]: DeathEvent dead-num_death_rules caused the cell to die
	cdef int dead
	cdef state_set
	cdef double volume

	def __init__(self, v0 = 0, t0 = 0, state = []):
		self.set_initial_vars(v0, t0)
		self.set_volume(v0)
		self.set_time(t0)
		self.divided = -1
		self.dead = -1
		if len(state) == 0:
			self.state_set = 0
		else:
			self.state_set = 1
			self.py_set_state(state)

	def get_state_set(self):
		return self.state_set

	def py_set_state(self, state):
		self.state_set = 1
		return super().py_set_state(np.asarray(state))

	cdef void set_initial_vars(self, double volume, double time):
		self.initial_volume = volume
		self.initial_time = time

	cdef double get_initial_volume(self):
		return self.initial_volume

	cdef void set_state_comp(self, double val, unsigned comp_ind):
		self.state[comp_ind] = val

	cdef double get_state_comp(self, unsigned comp_ind):
		return self.state[comp_ind]

	def py_get_initial_volume(self):
		return self.get_initial_volume()

	cdef double get_initial_time(self):
		return self.initial_time

	def py_get_initial_time(self):
		return self.get_initial_time()

	cdef void set_divided(self, divided):
		self.divided = divided

	cdef int get_divided(self):
		return self.divided

	cdef void set_dead(self, dead):
		self.dead = dead
	cdef int get_dead(self):
		return self.dead

	cdef void set_volume(self, vol):
		self.initial_volume = vol

	cdef void set_time(self, time):
		self.initial_time = time







cdef class DummyClass():
	'''
	Convenience class for development. Can safely remove.
	'''
	cdef int status
	cdef int x1, x3, x5

	def __init__(self, status):
		self.status = status
		self.x1 = 1
		self.x3 = 1
		self.x5 = 1

	def get_status(self):
		return self.status

	def set_status(self, new_status):
		self.status = new_status

def py_test_capped_queue():
	return test_capped_queue()

cdef test_capped_queue():
	times = [51,22,63,34,25,46,67,38,59, 100, 101, 102, 103]
	events = [LineageVolumeCellState(t0 = i) for i in times]
	q = CappedStateQueue(9)
	for e in events:
		q.push_event(e)
		print("Still heap? " + str(q.is_still_heap()) + ".")
		print([(<LineageVolumeCellState>e).initial_time for e in q.get_underlying_array()])

	print([(<LineageVolumeCellState>e).initial_time for e in q.get_underlying_array()])

	while not q.empty():
		print("________________________________")
		print("(pop_size is " + str(q.pop_size) + ")")
		print((<LineageVolumeCellState>q.pop_event()).initial_time)
		print([(<LineageVolumeCellState>e).initial_time for e in q.get_underlying_array()])




# ctypedef int* int_pointer 
# ctypedef pair[pair[float, bool], void *] event_pair
# ctypedef pair[void *, void *] test_pair

# cdef class DummyClass():
# 	cdef int status
# 	cdef int x1, x3, x5

# 	def __init__(self, status):
# 		self.status = status
# 		self.x1 = 1
# 		self.x3 = 1
# 		self.x5 = 1

# 	def get_status(self):
# 		return self.status

# 	def set_status(self, new_status):
# # 		self.status = new_status


# def py_pq_test():
# 	return pq_test()

# cdef pq_test():
# 	cdef int int1, int2, int3
# 	cdef priority_queue[event_pair] pq 
# 	cdef event_pair pair1, pair2, pair3
# 	cdef pair[float, bool] pair1a, pair2a, pair3a
# 	cdef test_pair tp
# 	cdef DummyClass d1, d2, d3
# 	cdef list obj_list

# 	d1 = DummyClass(55)
# 	d2 = DummyClass(66)
# 	d3 = DummyClass(11)

# 	tp.first = <void *>d1
# 	tp.second = <void *>("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
# 	print(<long>(&tp.first))
# 	print(<long>(&tp.second))
# 	(<DummyClass>(tp.first)).x5 = <int>'B'
# 	(<DummyClass>(tp.first)).x1 = <int>'B'
# 	(<DummyClass>(tp.first)).x3 = <int>'B'
# 	(<DummyClass>(tp.first)).status = <int>'B'
# 	print(<str>tp.second)

# 	obj_list = [d1, d2]
# 	obj_list.append(d3)

# 	int1 = 1
# 	int2 = -5
# 	int3 = 10

# 	pair1a.first = -0.5
# 	pair1a.second = True
# 	pair1.first = pair1a
# 	pair1.second = <void *> d1

# 	pair2a.first = 1
# 	pair2a.second = False
# 	pair2.first = pair2a
# 	pair2.second = <void *> d2

# 	pair3a.first = 10
# 	pair3a.second = True
# 	pair3.first = pair3a
# 	pair3.second = <void *> d3
	

# 	pq.push(pair1)
# 	pq.push(pair2)
# 	pq.push(pair3)

# 	obj_list[0].set_status(100)
# 	del obj_list[-1]
# 	print(len(obj_list))

# 	while not pq.empty():
# 		apair = pq.top()
# 		i = <DummyClass> apair.second
# 		print(f"{apair.first}, {i.status}")
# 		pq.pop()
# 	return 0