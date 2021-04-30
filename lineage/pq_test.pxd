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

	cdef void set_initial_vars(self, double volume, double time)

	cdef double get_initial_volume(self)

	cdef void set_state_comp(self, double val, unsigned comp_ind)

	cdef double get_state_comp(self, unsigned comp_ind)

	cdef double get_initial_time(self)

	cdef void set_divided(self, divided)

	cdef int get_divided(self)

	cdef void set_dead(self, dead)
	cdef int get_dead(self)

	cdef void set_volume(self, vol)
	cdef void set_time(self, time)