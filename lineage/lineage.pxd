from bioscrape.simulator cimport DelayVolumeCellState

cdef class LineageVolumeCellState(DelayVolumeCellState):
	cdef double initial_volume #Stores the birth Volume
	cdef double initial_time #Stores the time the Cell was "born"
	cdef int divided
	cdef int dead
	cdef state_set

	cdef void set_initial_vars(self, double volume, double time)

	cdef double get_initial_volume(self)

	cdef void set_state_comp(self, double val, unsigned comp_ind)

	cdef double get_state_comp(self, unsigned comp_ind)
	
	cdef double get_initial_time(self)

	cdef void set_divided(self, divided)

	cdef int get_divided(self)

	cdef void set_dead(self, dead)

	cdef int get_dead(self)