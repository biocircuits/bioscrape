# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=True


import numpy as np
cimport numpy as np
cimport random as cyrandom
from vector cimport vector
from libc.math cimport fabs
from types cimport Model, Delay, Propensity, Rule
from scipy.integrate import odeint, ode
import sys
import warnings
import logging


##################################################                ####################################################
######################################              DELAY QUEUE TYPES                   ##############################
#################################################                     ################################################

cdef class DelayQueue:
    cdef void add_reaction(self, double time, unsigned rxn_id, double amount):
        """
        Add a reaction to the queue.
        :param time: (double) The time at which the reaction will occur
        :param rxn_id: (unsigned) The ID of the reaction, i.e. its column index in the stoich matrix
        :param amount: (double) How many of the reaction occurs?
        :return: None
        """

        pass

    def py_add_reaction(self,double time, unsigned rxn_id, double amount):
        self.add_reaction(time, rxn_id, amount)


    cdef double get_next_queue_time(self):
        """
        Find the nearest queued time. Note that it's possible no reaction occurs at the next queue time.
        This is possible in the case where the queue internally updates with some time resolution even if no reactions
        happen.
        :return: (double) The next queue time.
        """

        return 0.0

    def py_get_next_queue_time(self):
        return self.get_next_queue_time()

    cdef void get_next_reactions(self, double *rxn_array):
        """
        Find the next reaction time. rxn_array must have at least enough room available for the number of reactions.

        :param rxn_array: (double *) A place to store how many of each reaction occurs at the next queued time.
        :return: None
        """
        pass

    def py_get_next_reactions(self, np.ndarray[np.double_t, ndim=1] rxn_array):
        self.get_next_reactions(<double*> rxn_array.data)


    cdef void advance_time(self):
        """
        Advance the queue to the next queue relevant time and perform whatever internal updates are necessary. Make
        sure to call get_next_reactions() before advance_time() or you will never know what reactions occurred.
        :return: None
        """
        pass

    def py_advance_time(self):
        self.advance_time()


    cdef DelayQueue copy(self):
        """
        Cope the DelayQueue and return a new totally independent but duplicate one.
        :return: (DelayQueue) The copied DelayQueue
        """
        return None

    def py_copy(self):
        return self.copy()


    cdef void set_current_time(self, double t):
        """
        Set the current time for the queue.
        :param t: (double) the time.
        :return: None
        """
        pass

    def py_set_current_time(self, double t):
        self.set_current_time(t)

    cdef DelayQueue clear_copy(self):
        """
        Copy the DelayQueue and return a new one with the same config, but with no reactions contained in it.
        :return: (DelayQueue) A new and clear DelayQueue.
        """
        return None

    def py_clear_copy(self):
        return self.clear_copy()

    cdef np.ndarray binomial_partition(self, double p):
        """
        Partition the delay queue into two delay queues with reactions switching according to probability p
        :param p: (double) The binomial parameter 0 < p < 1
        :return: (np.ndarray) A length 2 array of objects. Each of these has to be casted back to a DelayQueue type
        """
        return None
    def py_binomial_partition(self, double p):
        return self.binomial_partition(p)

cdef class ArrayDelayQueue(DelayQueue):
    def __init__(self, np.ndarray queue, double dt, double current_time):
        """
        Initialize with a queue, dt resolution, and current time. The queue should have one row for each reaction and
        the max future time that can be handled is dt*(number of queue columns).

        :param queue: (np.ndarray) 2-D array containing the current time.
        :param dt: (double) The time resolution dt
        :param current_time: (double) The current time.
        """
        self.num_reactions = queue.shape[0]
        self.num_cols = queue.shape[1]
        self.queue = queue
        self.next_queue_time = current_time + dt
        self.dt = dt
        self.start_index = 0

    @staticmethod
    def setup_queue(unsigned num_reactions, unsigned queue_length, double dt):
        """
        Static method to create an empty ArrayDelayQueue given a desired length, number of reactions, and dt.

        :param num_reactions: (unsigned) number of reactions in the system
        :param queue_length: (unsigned) length of the queue
        :param dt: (double) time step
        :return: The created array delay queue object
        """
        return ArrayDelayQueue(np.zeros((num_reactions,queue_length)), dt, 0.0)


    cdef DelayQueue copy(self):
        cdef ArrayDelayQueue a = ArrayDelayQueue(self.queue, self.dt, 0.0)
        a.num_reactions = self.num_reactions
        a.num_cols = self.num_cols
        a.next_queue_time = self.next_queue_time
        a.dt = self.dt
        a.start_index = self.start_index
        a.queue = self.queue.copy()
        return a

    cdef DelayQueue clear_copy(self):
        """
        Copy the DelayQueue and return a new one with the same config, but with no reactions contained in it.
        :return: (DelayQueue) A new and clear DelayQueue.
        """
        cdef ArrayDelayQueue a = ArrayDelayQueue(self.queue, self.dt, 0.0)
        a.num_reactions = self.num_reactions
        a.num_cols = self.num_cols
        a.next_queue_time = self.next_queue_time
        a.dt = self.dt
        a.start_index = self.start_index
        a.queue = np.zeros((self.num_reactions,self.num_cols))
        return a


    cdef void set_current_time(self, double t):
        """
        Set the current time
        :param t: (double) the time
        :return: None
        """
        self.next_queue_time = t + self.dt

    cdef void add_reaction(self, double time, unsigned rxn_id, double amount):
        """
        Add a reaction to the queue. If the reaction time is past the max time supported by the queue length, then
        truncate to the maximum queue time. Round to the nearest dt grid point as well when inserting.
        :param time: (double) the time at which the reaction occurs
        :param rxn_id: (unsigned) the id of the reaction
        :param amount: (double) how many of the reaction occurs (typically 1.0)
        :return:
        """
        # Round to the nearest entry in the delay queue.
        cdef int index = int( (time - self.next_queue_time) / self.dt + 0.5 )
        # Don't let the index get too small or too big, truncate to fit into the queue
        if index < 0:
            index = 0
        elif index >= int(self.num_cols):
            index = self.num_cols-1

        # Shift by the start index offset
        index = (index + self.start_index) % self.num_cols

        self.queue[rxn_id,index] += amount

    cdef double get_next_queue_time(self):
        return self.next_queue_time

    cdef void get_next_reactions(self, double *rxn_array):
        cdef unsigned i
        for i in range(self.num_reactions):
            rxn_array[i] = self.queue[i,self.start_index]



    cdef np.ndarray binomial_partition(self, double p):
        """
        Partition the delay queue into two delay queues with reactions switching according to probability p
        :param p: (double) The binomial parameter 0 < p < 1
        :return: (np.ndarray) A length 2 array of objects. Each of these has to be casted back to a DelayQueue type
        """
        cdef ArrayDelayQueue q1 = self.clear_copy()
        cdef ArrayDelayQueue q2 = self.clear_copy()

        cdef unsigned time_points = q1.queue.shape[1]
        cdef unsigned num_reactions = q1.queue.shape[0]

        cdef unsigned time_index = 0
        cdef unsigned reaction_index = 0

        for time_index in range(time_points):
            for reaction_index in range(num_reactions):
                q1.queue[reaction_index,time_index] = cyrandom.binom_rnd_f(self.queue[reaction_index,time_index],p)
                q2.queue[reaction_index,time_index] = self.queue[reaction_index,time_index] - q1.queue[reaction_index,time_index]

        cdef np.ndarray a = np.empty(2,dtype=object)

        a[0] = q1
        a[1] = q2

        return a


    cdef void advance_time(self):
        # advance time by dt
        self.next_queue_time += self.dt
        # clear the current next queued time fully
        cdef unsigned i
        for i in range(self.num_reactions):
            self.queue[i,self.start_index] = 0
        # advanced the start index by 1 cycling around the end.
        self.start_index = (self.start_index + 1) % self.num_cols

##################################################                ####################################################
######################################              SIMULATION INTERFACES               ##############################
#################################################                     ################################################

cdef class CSimInterface:
    cdef np.ndarray get_update_array(self):
        return self.update_array

    def py_get_update_array(self):
        return self.get_update_array()


    cdef np.ndarray get_delay_update_array(self):
        return self.delay_update_array
    def py_get_delay_update_array(self):
        return self.get_delay_update_array()

    #Checks model or interface is valid. Meant to be overriden by the subclass
    cdef void check_interface(self):
        logging.info("No interface Checking Implemented")
    # meant to be overriden by the subclass
    cdef double compute_delay(self, double *state, unsigned rxn_index):
        return 0.0

    # must be overriden by subclass
    cdef void compute_propensities(self, double *state, double *propensity_destination, double time):
        pass
    cdef void compute_volume_propensities(self, double *state, double *propensity_destination, double volume, double time):
        pass

    # by default stochastic propensities are assumed to be the same as normal propensities. This may be overwritten by the subclass, however.
    cdef void compute_stochastic_propensities(self, double *state, double *propensity_destination, double time):
        self.compute_propensities(state, propensity_destination, time)

    # by default stochastic propensities are assumed to be the same as normal propensities. This may be overwritten by the subclass, however.
    cdef void compute_stochastic_volume_propensities(self, double *state, double *propensity_destination, double volume, double time):
        self.compute_volume_propensities(state, propensity_destination, volume, time)

    cdef unsigned requires_delay(self):
        return self.delay_flag

    cdef np.ndarray get_initial_state(self):
        return self.initial_state

    def py_get_initial_state(self):
        return self.get_initial_state()

    cdef void set_initial_state(self, np.ndarray a):
        self.initial_state = a

    def py_set_initial_state(self, np.ndarray a):
        self.set_initial_state(a)

    cdef unsigned get_num_reactions(self):
        return self.num_reactions

    def py_get_num_reactions(self):
        return self.get_num_reactions()

    cdef unsigned get_num_species(self):
        return self.num_species

    def py_get_num_species(self):
        return self.get_num_species()

    cdef double get_initial_time(self):
        return self.initial_time

    def py_get_initial_time(self):
        return self.get_initial_time()

    cdef void set_initial_time(self, double t):
        self.initial_time = t

    def py_set_initial_time(self, double t):
        self.set_initial_time(t)

    cdef void set_dt(self, double dt):
        self.dt = dt


    def py_set_dt(self, double dt):
        self.set_dt(dt)


    cdef double get_dt(self):
        return self.dt

    def py_get_dt(self):
        return self.get_dt()

    cdef double* get_param_values(self):
        return <double*> 0

    def py_get_param_values(self):
        return None

    cdef unsigned get_num_parameters(self):
        return 0

    def py_get_num_parameters(self):
        return self.get_num_parameters()

    cdef unsigned get_number_of_rules(self):
        return 0

    def py_get_number_of_rules(self):
        return self.get_number_of_rules()

    cdef void apply_repeated_rules(self, double *state, double time):
        pass

    cdef void apply_repeated_volume_rules(self, double *state, double volume, double time):
        pass

    def py_apply_repeated_rules(self, np.ndarray[np.double_t, ndim=1] state, double time=0.0):
        self.apply_repeated_rules(<double*> state.data,time)


    # Prepare for determinsitic simulation by creating propensity buffer and also doing the compressed stoich matrix
    cdef void prep_deterministic_simulation(self):
        # Clear out the vectors
        self.S_indices.clear()
        self.S_values.clear()
        cdef unsigned r
        cdef unsigned s
        # keep track of nonzero indices and the coefficients as well
        for s in range(self.num_species):
            # Add vectors for that row
            self.S_indices.push_back(vector[int]())
            self.S_values.push_back(vector[int]())
            for r in range(self.num_reactions):
                if self.update_array[s,r]+self.delay_update_array[s,r] != 0:
                    self.S_indices[s].push_back(r)
                    self.S_values[s].push_back(self.update_array[s,r]+self.delay_update_array[s,r])
        # Create proper size propensity buffer.
        self.propensity_buffer = np.zeros(self.num_reactions,)

        # Set the global simulation object to this model
        global global_sim
        global_sim = self

    def py_prep_deterministic_simulation(self):
        self.prep_deterministic_simulation()

    # Compute deterministic derivative
    cdef void calculate_deterministic_derivative(self, double *x, double *dxdt, double t):
        # Get propensities before doing anything else.
        cdef double *prop = <double*> (self.propensity_buffer.data)
        self.compute_propensities(x,  prop, t)

        cdef unsigned s
        cdef unsigned j
        for s in range(self.num_species):
            dxdt[s] = 0;
            for j in range(self.S_indices[s].size()):
                dxdt[s] += prop[ self.S_indices[s][j]  ] * self.S_values[s][j]


    def py_calculate_deterministic_derivative(self, np.ndarray[np.double_t,ndim=1] x, np.ndarray[np.double_t,ndim=1] dx,
                                              double t):
        self.calculate_deterministic_derivative(<double*> x.data, <double*> dx.data, t)



cdef class ModelCSimInterface(CSimInterface):
    def __init__(self, external_model):
        self.model = external_model
        #Check Model and initialization
        if not self.model.initialized:
            self.model.py_initialize()
            logging.info("Uninitialized Model Passed into ModelCSimInterface. Model.py_initialize() called automatically.")
        self.check_interface()
        self.c_propensities = self.model.get_c_propensities()
        self.c_delays = self.model.get_c_delays()
        self.c_repeat_rules = self.model.get_c_repeat_rules()
        self.update_array = self.model.get_update_array()
        self.delay_update_array = self.model.get_delay_update_array()
        self.initial_state = self.model.get_species_values()
        self.np_param_values = (self.model.get_params_values())
        self.c_param_values = <double*>(self.np_param_values.data)
        self.num_reactions = self.update_array.shape[1]
        self.num_species = self.update_array.shape[0]
        self.dt = 0.01

    cdef unsigned get_number_of_species(self):
        return self.num_species

    cdef unsigned get_number_of_reactions(self):
        return self.num_reactions
        
    cdef void check_interface(self):
        if not self.model.initialized:
            raise RuntimeError("Model has been changed since CSimInterface instantiation. CSimInterface no longer valid.")

    cdef double compute_delay(self, double *state, unsigned rxn_index):
        return  (<Delay> (self.c_delays[0][rxn_index])).get_delay(state, self.c_param_values)

    cdef void compute_propensities(self, double *state, double *propensity_destination, double time):
        cdef unsigned rxn
        for rxn in range(self.num_reactions):
            propensity_destination[rxn] = (<Propensity> (self.c_propensities[0][rxn]) ).get_propensity(state, self.c_param_values, time)

    cdef void compute_volume_propensities(self, double *state, double *propensity_destination, double volume, double time):
        cdef unsigned rxn
        for rxn in range(self.num_reactions):
            propensity_destination[rxn] = (<Propensity> (self.c_propensities[0][rxn]) ).get_volume_propensity(state, self.c_param_values,
                                                                                                              volume, time)
    cdef void compute_stochastic_propensities(self, double *state, double *propensity_destination, double time):
        cdef unsigned rxn
        for rxn in range(self.num_reactions):
            propensity_destination[rxn] = (<Propensity> (self.c_propensities[0][rxn]) ).get_stochastic_propensity(state,
                                                                                                       self.c_param_values, time)
    cdef void compute_stochastic_volume_propensities(self, double *state, double *propensity_destination, double volume, double time):
        cdef unsigned rxn
        for rxn in range(self.num_reactions):
            propensity_destination[rxn] = (<Propensity> (self.c_propensities[0][rxn]) ).get_stochastic_volume_propensity(state, self.c_param_values, volume, time)

    cdef unsigned get_number_of_rules(self):
        return self.c_repeat_rules[0].size()

    cdef void apply_repeated_rules(self, double *state, double time):
        cdef unsigned rule_number
        for rule_number in range(self.c_repeat_rules[0].size()):
            (<Rule> (self.c_repeat_rules[0][rule_number])).execute_rule(state, self.c_param_values, time)

    cdef void apply_repeated_volume_rules(self, double *state, double volume, double time):
        cdef unsigned rule_number
        for rule_number in range(self.c_repeat_rules[0].size()):
            (<Rule> (self.c_repeat_rules[0][rule_number])).execute_volume_rule(state, self.c_param_values, volume, time)

    cdef np.ndarray get_initial_state(self):
        return self.initial_state

    cdef void set_initial_state(self, np.ndarray a):
        np.copyto(self.initial_state,a)

    cdef double* get_param_values(self):
        return self.c_param_values

    def py_get_param_values(self):
        return self.np_param_values

    cdef unsigned get_num_parameters(self):
        return self.np_param_values.shape[0]

cdef class SafeModelCSimInterface(ModelCSimInterface):
    def __init__(self, external_model, max_volume = 1000, max_species_count = 1000):
        self.max_volume = max_volume
        self.max_species_count = max_species_count
        super().__init__(external_model)
        self.initialize_reaction_inputs()
        
    cdef void initialize_reaction_inputs(self):
        self.rxn_ind = 0
        self.s_ind = 0
        cdef unsigned ind = 0

        #Stores a list of the species index that are inputs to reaction r in a numpy array. List is over when -1 is reached
        empty_array = -np.ones((self.num_reactions, self.num_species, 2), dtype = np.int32)
        self.reaction_input_indices = empty_array.data
        #Parallel array to the one above
        for self.rxn_ind in range(self.num_reactions):
            ind = 0
            for self.s_ind in range(self.num_species):
                #IF a species s is consumed either by the reaction or delay reaction
                if (self.update_array[self.s_ind, self.rxn_ind] < 0) or (self.delay_update_array[self.s_ind, self.rxn_ind] < 0):
                    #add s to the reaction_update_indices[rxn_ind, :, 0] vector
                    self.reaction_input_indices[self.rxn_ind, ind, 0] = self.s_ind

                    #set self.reaction_input_indices[self.rxn_ind, ind, 1] vector the maximum amount of the species that could be consumed
                    if (self.update_array[self.s_ind, self.rxn_ind] < 0) and (self.delay_update_array[self.s_ind, self.rxn_ind] < 0):
                        self.reaction_input_indices[self.rxn_ind, ind, 1] = -(self.update_array[self.s_ind, self.rxn_ind]+self.delay_update_array[self.s_ind, self.rxn_ind])
                    else:
                         self.reaction_input_indices[self.rxn_ind, ind, 1] = -min(self.update_array[self.s_ind, self.rxn_ind], self.delay_update_array[self.s_ind, self.rxn_ind])
                    ind += 1

    cdef void compute_stochastic_propensities(self, double *state, double *propensity_destination, double time):
        self.check_count_function(state, 1)
        self.rxn_ind = 0
        for self.rxn_ind in range(self.num_reactions):
            self.prop_is_0 = 0
            self.s_ind = 0
            while self.reaction_input_indices[self.rxn_ind, self.s_ind, 0] > -1 and self.prop_is_0 == 0:
                if state[self.reaction_input_indices[self.rxn_ind, self.s_ind, 0]] < self.reaction_input_indices[self.rxn_ind, self.s_ind, 1]:
                    propensity_destination[self.rxn_ind] = 0
                    self.prop_is_0 = 1
                self.s_ind+=1
            if self.prop_is_0 == 0:
                propensity_destination[self.rxn_ind] = (<Propensity> (self.c_propensities[0][self.rxn_ind]) ).get_stochastic_propensity(state, self.c_param_values, time)

    cdef void compute_stochastic_volume_propensities(self, double *state, double *propensity_destination, double volume, double time):
        self.check_count_function(state, volume)
        self.rxn_ind = 0
        for self.rxn_ind in range(self.num_reactions):
            self.prop_is_0 = 0
            self.s_ind = 0
            while self.reaction_input_indices[self.rxn_ind, self.s_ind, 0] != -1 and self.prop_is_0 == 0:
                if state[self.reaction_input_indices[self.rxn_ind, self.s_ind, 0]] < self.reaction_input_indices[self.rxn_ind, self.s_ind, 1]:
                    propensity_destination[self.rxn_ind] = 0
                    self.prop_is_0 = 1
                self.s_ind+=1
            if self.prop_is_0 == 0:
                propensity_destination[self.rxn_ind] = (<Propensity> (self.c_propensities[0][self.rxn_ind]) ).get_stochastic_volume_propensity(state, self.c_param_values, volume, time)

    cdef void check_count_function(self, double *state, double volume):
        self.s_ind = 0

        for self.s_ind in range(self.num_species):
            if state[self.s_ind] > self.max_species_count:
                warnings.warn("Species #"+str(self.s_ind)+"="+str(state[self.s_ind])+" > Max Count="+str(self.max_species_count))
            elif state[self.s_ind] < 0:
                 warnings.warn("Species #"+str(self.s_ind)+"="+str(state[self.s_ind])+" < 0")
        if volume > self.max_volume:
            warnings.warn("Volume="+str(volume)+" > Max Volume="+str(self.max_volume))
        elif volume <= 0:
            warnings.warn("Volume="+str(volume)+" > Max Volume="+str(self.max_volume))

    # Compute deterministic derivative
    # This version safegaurds against species counts going negative 
    # by not allowing reactions to fire if they consume a species of 0 concentration.
    cdef void calculate_deterministic_derivative(self, double *x, double *dxdt, double t):
        # Get propensities before doing anything else.
        cdef double *prop = <double*> (self.propensity_buffer.data)
        self.compute_propensities(x,  prop, t)

        cdef unsigned s
        cdef unsigned j
        for s in range(self.num_species):

            #Reset negative species concentrations to 0
            if x[s] < 0:
                x[s] = 0

            dxdt[s] = 0
            for j in range(self.S_indices[s].size()):
                if self.S_values[s][j] <= 0 and x[s] <= 0:
                    pass #Skip reactions that consume species of 0 concentration
                else:
                    dxdt[s] += prop[ self.S_indices[s][j]  ] * self.S_values[s][j]

            #Verify that species do not go negative.
            if x[s] <= 0 and dxdt[s] < 0:
                dxdt[s] = 0

cdef class SSAResult:
    def __init__(self, np.ndarray timepoints, np.ndarray result):
        self.timepoints = timepoints
        self.simulation_result = result

    def py_get_timepoints(self):
        return self.get_timepoints()

    def py_get_result(self):
        return self.get_result()

    #Returns a Pandas Data Frame, if it is installed. If not, a Numpy array is returned.
    def py_get_dataframe(self, Model = None):
        try:
            import pandas
            if Model == None:
                warnings.warn("No Model passed into py_get_dataframe. No species names will be attached to the data frame.")
                df = pandas.DataFrame(data = self.get_result())
            else:
                columns = Model.get_species_list()
                df = pandas.DataFrame(data = self.get_result(), columns = columns)
            df['time'] = self.timepoints
            return df

        except ModuleNotFoundError:
            warnings.warn("py_get_dataframe requires the pandas Module to return a Pandas Dataframe object. Numpy array being returned instead.")
            return self.py_get_result()



cdef class DelaySSAResult(SSAResult):
    def __init__(self, np.ndarray timepoints, np.ndarray result, DelayQueue queue):
        self.timepoints = timepoints
        self.final_delay_queue = queue
        self.simulation_result = result

    def py_get_delay_queue(self):
        return self.get_delay_queue()


cdef class VolumeSSAResult(SSAResult):
    def __init__(self,np.ndarray timepoints,np.ndarray result,np.ndarray volume,unsigned divided):
        super().__init__(timepoints, result)
        self.volume = volume
        self.cell_divided_flag = divided

    def py_cell_divided(self):
        return self.cell_divided()
    def py_get_volume(self):
        return self.get_volume()

    def py_get_dataframe(self, Model = None):
        df = super().py_get_dataframe(Model = Model)
        try:
            import pandas
            df["volume"] = self.volume
            return df
        except ModuleNotFoundError:
            warnings.warn("py_get_dataframe requires the pandas Module to return a Pandas Dataframe object. Numpy array being returned instead. It is highly recommended that you install Pandas")
            return df

    cdef VolumeCellState get_final_cell_state(self):
        """
        Get the final cell state from a simulation result.
        :return: (VolumeCellState) The final cell state.
        """
        cdef VolumeCellState cs = VolumeCellState()
        cdef unsigned final_index = (<np.ndarray[np.double_t,ndim=1]> self.timepoints).shape[0]-1
        cs.set_time(self.timepoints[final_index])
        cs.set_volume(self.volume[final_index])
        cs.set_state(self.simulation_result[final_index,:])
        cs.set_volume_object(self.volume_object)
        return cs

    def py_get_final_cell_state(self):
        return self.get_final_cell_state()

    cdef Schnitz get_schnitz(self):
        return Schnitz(self.timepoints,self.simulation_result,self.volume)

    def py_get_schnitz(self):
        return self.get_schnitz()


cdef class DelayVolumeSSAResult(VolumeSSAResult):
    def __init__(self,np.ndarray timepoints, np.ndarray result, np.ndarray volume, DelayQueue queue, unsigned divided):
        super().__init__(timepoints, result, volume, divided)
        self.final_delay_queue = queue

    def py_get_delay_queue(self):
        return self.get_delay_queue()

    cdef VolumeCellState get_final_cell_state(self):
        cdef DelayVolumeCellState cs = DelayVolumeCellState()
        cdef unsigned final_index = (<np.ndarray[np.double_t,ndim=1]> self.timepoints).shape[0]-1
        cs.set_time(self.timepoints[final_index])
        cs.set_volume(self.volume[final_index])
        cs.set_state(self.simulation_result[final_index,:])
        cs.set_delay_queue(self.final_delay_queue)
        return cs


# Cell state classes

cdef class CellState:
    def py_set_state(self, np.ndarray state):
        self.state = state

    def py_get_state(self):
        return self.state

    def py_set_time(self, double t):
        self.time = t

    def py_get_time(self):
        return self.time

    def py_set_species(self, model, specie, value):
        ind = model.get_species_index(specie)
        self.state[ind] = value

    def py_get_dataframe(self, Model = None):
        try:
            import pandas
            if Model == None:
                warnings.warn("No Model passed into py_get_dataframe. No species names will be attached to the data frame.")
                df = pandas.DataFrame(data = np.expand_dims(self.state, 0))
            else:
                columns = Model.get_species_list()
                df = pandas.DataFrame(data = np.expand_dims(self.state, 0), columns = columns)
                
            df['time'] = self.time
            return df

        except ModuleNotFoundError:
            warnings.warn("py_get_dataframe requires the pandas Module to return a Pandas Dataframe object. Numpy array being returned instead. It is highly recommended that you install Pandas")
            return self.state

cdef class DelayCellState(CellState):
    def py_get_delay_queue(self):
        return self.delay_queue

    def py_set_delay_queue(self, DelayQueue q):
        self.delay_queue = q

cdef class VolumeCellState(CellState):
    def __init__(self):
        self.volume_object = None

    def py_set_volume(self, double volume):
        self.volume = volume

    def py_get_volume(self):
        return self.volume

    def py_set_volume_object(self, Volume v):
        self.volume_object = v

    def py_get_volume_object(self):
        return self.volume_object

    def __setstate__(self, state):
        self.time = state[0]
        self.volume = state[1]
        self.state = state[2].copy()

    def __getstate__(self):
        return (self.time, self.volume, self.state)

    def py_get_dataframe(self, Model = None):
        df = super().py_get_dataframe(Model = Model)
        try:
            import pandas
            df["volume"] = self.volume
        except ModuleNotFoundError:
            warnings.warn("py_get_dataframe requires the pandas Module to return a Pandas Dataframe object. Numpy array being returned instead. It is highly recommended that you install Pandas")
            pass
        return df

cdef class DelayVolumeCellState(VolumeCellState):
    def py_get_delay_queue(self):
        return self.delay_queue

    def py_set_delay_queue(self, DelayQueue q):
        self.delay_queue = q

##################################################                ####################################################
######################################              LINEAGE                           ################################
#################################################                     ################################################

import sys
cdef Lineage simulate_cell_lineage(CSimInterface sim, Volume v, np.ndarray timepoints,
                                   VolumeSimulator vsim, VolumeSplitter vsplit):

    """
    Simulate a cell lineage given the reactions, the volume and partitioning models, and the desired simulator over
    a given set of time points.
    :param sim: (CSimInterface) The reaction system.
    :param v: (Volume) The volume model.
    :param timepoints: (np.ndarray)
    :param vsim: (VolumeSimulator) The simulator to use.
    :param vsplit: (VolumeSplitter) The partitioning model to use.
    :return: (Lineage) Simulated cell lineage output.
    """

    # Prepare a lineage structure to save the data output.
    cdef Lineage l = Lineage()
    cdef np.ndarray[np.double_t, ndim=1] c_timepoints = timepoints
    cdef np.ndarray[np.double_t, ndim=1] c_truncated_timepoints
    cdef double final_time = c_timepoints[c_timepoints.shape[0]-1]

    # Simulate the first cell.
    cdef VolumeSSAResult r = vsim.volume_simulate(sim,v, timepoints)
    cdef Schnitz s = r.get_schnitz()
    cdef Schnitz daughter_schnitz1, daughter_schnitz2
    cdef VolumeCellState cs = r.get_final_cell_state()
    cdef VolumeCellState d1, d2, d1final, d2final
    cdef np.ndarray daughter_cells
    s.set_parent(None)

    l.add_schnitz(s)

    cdef unsigned list_index = 0
    cdef list old_schnitzes = []
    cdef list old_cell_states = []

    if cs.get_time() < final_time:
        old_schnitzes.append(s)
        old_cell_states.append(cs)

    cdef unsigned debug_index
    cdef VolumeCellState dcs

    while list_index < len(old_cell_states):
        #Debugging purely
        # print("list index:", list_index)
        # for debug_index in range(len(old_cell_states)):
        #     dcs = old_cell_states[debug_index]
        #     print("final time: ", dcs.get_time())
        # sys.stdout.flush()

        cs = old_cell_states[list_index]
        s = old_schnitzes[list_index]
        list_index += 1

        # Partition the cell state
        daughter_cells = vsplit.partition(cs)
        d1 = <VolumeCellState>(daughter_cells[0])
        d2 = <VolumeCellState>(daughter_cells[1])

        #Create a new timepoint array and simulate the first daughter and queue if it doesn't reach final time.
        c_truncated_timepoints = c_timepoints[c_timepoints > cs.get_time()]

        sim.set_initial_time(d1.get_time())
        sim.set_initial_state(d1.get_state())
        v.initialize(<double*> d1.get_state().data,<double*> 0, d1.get_time(), d1.get_volume())

        r = vsim.volume_simulate(sim,v,c_truncated_timepoints)
        daughter_schnitz1 = r.get_schnitz()
        d1final = r.get_final_cell_state()

        # Add on the new daughter if final time wasn't reached.
        if d1final.get_time() < final_time - 1E-9:
            #print("daughter 1 final time: ", d1final.get_time())
            old_schnitzes.append(daughter_schnitz1)
            old_cell_states.append(d1final)

        #Do the exact same thing for the other daughter.
        sim.set_initial_time(d2.get_time())
        sim.set_initial_state(d2.get_state())
        v.initialize(<double*> d2.get_state().data,<double*> 0, d2.get_time(), d2.get_volume())

        r = vsim.volume_simulate(sim,v,c_truncated_timepoints)

        daughter_schnitz2 = r.get_schnitz()
        d2final = r.get_final_cell_state()


        if d2final.get_time() < final_time - 1E-9:
            #print("daughter 2 final time: ", d2final.get_time())
            old_schnitzes.append(daughter_schnitz2)
            old_cell_states.append(d2final)

        # Set up daughters and parent appropriately.
        daughter_schnitz1.set_parent(s)
        daughter_schnitz2.set_parent(s)
        s.set_daughters(daughter_schnitz1,daughter_schnitz2)

        # Add daughters to the lineage
        l.add_schnitz(daughter_schnitz1)
        l.add_schnitz(daughter_schnitz2)

    return l

def py_simulate_cell_lineage(CSimInterface sim, Volume v, np.ndarray timepoints,
                                   VolumeSimulator vsim, VolumeSplitter vsplit):
    return simulate_cell_lineage(sim,v,timepoints,vsim,vsplit)



cdef Lineage simulate_delay_cell_lineage(CSimInterface sim, DelayQueue q, Volume v, np.ndarray timepoints,
                                   DelayVolumeSimulator dvsim, DelayVolumeSplitter dvsplit):
    """
    Simulate a cell lineage given the reactions, the volume, delay and partitioning models, and the desired simulator
    over a given set of time points.
    :param sim: (CSimInterface) The reaction system.
    :param q: (DelayQueue) The DelayQueue object containing the initial set of queued reactions.
    :param v: (Volume) The volume model.
    :param timepoints: (np.ndarray)
    :param vsim: (VolumeSimulator) The simulator to use.
    :param vsplit: (VolumeSplitter) The partitioning model to use.
    :return: (Lineage) Simulated cell lineage output.
    """

    # Prepare a lineage structure to save the data output.
    cdef Lineage l = Lineage()
    cdef np.ndarray[np.double_t, ndim=1] c_timepoints = timepoints
    cdef np.ndarray[np.double_t, ndim=1] c_truncated_timepoints
    cdef double final_time = c_timepoints[c_timepoints.shape[0]-1]

    # Simulate the first cell.
    cdef DelayVolumeSSAResult r = dvsim.delay_volume_simulate(sim,q,v,timepoints)
    cdef Schnitz s = r.get_schnitz()
    cdef Schnitz daughter_schnitz1, daughter_schnitz2
    cdef DelayVolumeCellState cs = r.get_final_cell_state()
    cdef DelayVolumeCellState d1, d2, d1final, d2final
    cdef np.ndarray daughter_cells
    s.set_parent(None)

    l.add_schnitz(s)

    cdef unsigned list_index = 0
    cdef list old_schnitzes = []
    cdef list old_cell_states = []

    if cs.get_time() < final_time:
        old_schnitzes.append(s)
        old_cell_states.append(cs)

    cdef unsigned debug_index
    cdef DelayVolumeCellState dcs

    while list_index < len(old_cell_states):
        #Debugging purely
        # print("list index:", list_index)
        # for debug_index in range(len(old_cell_states)):
        #     dcs = old_cell_states[debug_index]
        #     print("final time: ", dcs.get_time())
        # sys.stdout.flush()

        cs = old_cell_states[list_index]
        s = old_schnitzes[list_index]
        list_index += 1

        # Partition the cell state
        daughter_cells = dvsplit.partition(cs)
        d1 = <DelayVolumeCellState>(daughter_cells[0])
        d2 = <DelayVolumeCellState>(daughter_cells[1])

        #Create a new timepoint array and simulate the first daughter and queue if it doesn't reach final time.
        c_truncated_timepoints = c_timepoints[c_timepoints > cs.get_time()]

        sim.set_initial_time(d1.get_time())
        sim.set_initial_state(d1.get_state())
        v.initialize(<double*> d1.get_state().data,<double*> 0, d1.get_time(), d1.get_volume())

        r = dvsim.delay_volume_simulate(sim,d1.get_delay_queue(),v,c_truncated_timepoints)
        daughter_schnitz1 = r.get_schnitz()
        d1final = r.get_final_cell_state()

        # Add on the new daughter if final time wasn't reached.
        if d1final.get_time() < final_time - 1E-9:
            #print("daughter 1 final time: ", d1final.get_time())
            old_schnitzes.append(daughter_schnitz1)
            old_cell_states.append(d1final)

        #Do the exact same thing for the other daughter.
        sim.set_initial_time(d2.get_time())
        sim.set_initial_state(d2.get_state())
        v.initialize(<double*> d2.get_state().data,<double*> 0, d2.get_time(), d2.get_volume())

        r = dvsim.delay_volume_simulate(sim,d2.get_delay_queue(),v,c_truncated_timepoints)
        daughter_schnitz2 = r.get_schnitz()
        d2final = r.get_final_cell_state()


        if d2final.get_time() < final_time - 1E-9:
            #print("daughter 2 final time: ", d2final.get_time())
            old_schnitzes.append(daughter_schnitz2)
            old_cell_states.append(d2final)

        # Set up daughters and parent appropriately.
        daughter_schnitz1.set_parent(s)
        daughter_schnitz2.set_parent(s)
        s.set_daughters(daughter_schnitz1,daughter_schnitz2)

        # Add daughters to the lineage
        l.add_schnitz(daughter_schnitz1)
        l.add_schnitz(daughter_schnitz2)

    return l

def py_simulate_delay_cell_lineage(CSimInterface sim, DelayQueue q, Volume v, np.ndarray timepoints,
                                   DelayVolumeSimulator dvsim, DelayVolumeSplitter dvsplit):
    return simulate_delay_cell_lineage(sim,q,v,timepoints,dvsim,dvsplit)



##################################################                ####################################################
######################################              VOLUME SPLITTERS                  ################################
#################################################                     ################################################

cdef class VolumeSplitter:
    cdef np.ndarray partition(self, VolumeCellState parent):
        """
        Split the parent state into two VolumeCellState's.

        MUST BE SUBCLASSED.
        :param parent: (VolumeCellState)
        :return: (np.ndarray) Length 2 array containing daughter objects. These must be casted back to VolumeCellState.
        """
        raise NotImplementedError('partition() not implemented for VolumeSplitter')

    def py_partition(self, VolumeCellState parent):
        cdef np.ndarray ans = self.partition(parent)
        return <VolumeCellState> (ans[0]), <VolumeCellState> (ans[1])

cdef class DelayVolumeSplitter:
    cdef np.ndarray partition(self, DelayVolumeCellState parent):
        raise NotImplementedError('partition() not implemented for DelayVolumeSplitter')

    def py_partition(self, DelayVolumeCellState parent):
        cdef np.ndarray ans = self.partition(parent)
        return <DelayVolumeCellState> (ans[0]), <DelayVolumeCellState> (ans[1])


cdef class PerfectBinomialVolumeSplitter(VolumeSplitter):
    cdef np.ndarray partition(self, VolumeCellState parent):
        # set up daughters d and e
        cdef VolumeCellState d = VolumeCellState()
        cdef VolumeCellState e = VolumeCellState()

        # set times
        d.set_time(parent.get_time())
        e.set_time(parent.get_time())

        # set volume
        d.set_volume(parent.get_volume() / 2.0)
        e.set_volume(parent.get_volume() / 2.0)

        # partition the states
        cdef np.ndarray[np.double_t,ndim=1] dstate = parent.get_state().copy()
        cdef np.ndarray[np.double_t,ndim=1] estate = parent.get_state().copy()
        cdef unsigned length = dstate.shape[0]

        cdef unsigned i = 0
        cdef unsigned amount = 0

        for i in range(length):
            amount = cyrandom.binom_rnd_f(dstate[i],0.5)
            dstate[i] = <double> amount
            estate[i] -= dstate[i]

        d.set_state(dstate)
        e.set_state(estate)

        # create return structure
        cdef np.ndarray ans = np.empty(2, dtype=np.object)
        ans[0] = d
        ans[1] = e

        return ans


cdef class GeneralVolumeSplitter(VolumeSplitter):
    """
    A volume splitting class that splits the cell into two cells and can split species
    binomially, perfectly, or by duplication.
    """
    def __init__(self):
        self.partition_noise = 0

    def py_set_partitioning(self, dict options, Model m):
        """
        Set how each species is partitioned.
        :param options(dict: str -> [str]): A dictionary where keys specify partioning mode and values are lists of
         species partitioned using that mode. "perfect" means partitioned perfectly +/- 0.5 to maintain round numbers.
         "binomial" means partitioned binomially. "duplicate" means the species is duplicated with the same number into
         each daughter cell.
        :return: None
        """

        # clear vectors
        self.binomial_indices.clear()
        self.perfect_indices.clear()
        self.duplicate_indices.clear()

        # Create a set to keep track of used up indices.
        a = set(range(m.get_number_of_species()))

        if 'perfect' in options:
            L = options['perfect']
            for species in L:
                index = m.get_species_index(species)
                if index >= 0:
                    self.perfect_indices.push_back(index)
                    a.discard(index)
        if 'duplicate' in options:
            L = options['duplicate']
            for species in L:
                index = m.get_species_index(species)
                if index >= 0:
                    self.duplicate_indices.push_back(index)
                    a.discard(index)

        # Add everything else to binomial. NO need to check "binomial" entry since it is the default anyway.
        L = list(a)
        for index in L:
            self.binomial_indices.push_back(index)


    def py_set_partition_noise(self, double noise):
        """
        Set the noise, pick a uniform volume that is anywhere from 0.5 - noise to 0.5 of the original volume.
        :param noise: the noise param, should be < 0.5, ideally probably something like 0.05 to 0.1
        :return: none
        """
        self.partition_noise = noise

    cdef np.ndarray partition(self, VolumeCellState parent):
        # set up daughters d and e
        cdef VolumeCellState d = VolumeCellState()
        cdef VolumeCellState e = VolumeCellState()

        # set times
        d.set_time(parent.get_time())
        e.set_time(parent.get_time())

        # simulate partitioning noise
        cdef double p = 0.5 - cyrandom.uniform_rv()*self.partition_noise
        cdef double q = 1 - p

        # set volume
        d.set_volume(parent.get_volume() * p )
        e.set_volume(parent.get_volume() * q)

        # partition the states, copying already takes care of duplication replications.
        cdef np.ndarray[np.double_t,ndim=1] dstate = parent.get_state().copy()
        cdef np.ndarray[np.double_t,ndim=1] estate = parent.get_state().copy()
        cdef unsigned length = dstate.shape[0]

        cdef unsigned loop_index = 0
        cdef unsigned species_index = 0
        cdef unsigned amount = 0

        cdef double d_value = 0.0

        # take care of perfect splitting
        for loop_index in range(self.perfect_indices.size()):
            species_index = self.perfect_indices[loop_index]
            d_value = p * dstate[species_index]
            amount = <int> (d_value+0.5)
            if fabs(d_value-amount) <= 1E-8:
                dstate[species_index] = <double> amount
            else:
                if cyrandom.uniform_rv() <= p:
                    dstate[species_index] = <int> d_value + 1
                else:
                    dstate[species_index] = <int> d_value

            estate[species_index] -= dstate[species_index]


        # take care of binomial splitting
        for loop_index in range(self.binomial_indices.size()):
            species_index = self.binomial_indices[loop_index]
            amount = cyrandom.binom_rnd_f(dstate[species_index],p)
            dstate[species_index] = <double> amount
            estate[species_index] -= dstate[species_index]


        # set states and return
        d.set_state(dstate)
        e.set_state(estate)

        # create return structure
        cdef np.ndarray ans = np.empty(2, dtype=np.object)
        ans[0] = d
        ans[1] = e

        return ans




cdef class PerfectBinomialDelayVolumeSplitter(DelayVolumeSplitter):
    cdef np.ndarray partition(self, DelayVolumeCellState parent):
        # set up daughters d and e
        cdef DelayVolumeCellState d = DelayVolumeCellState()
        cdef DelayVolumeCellState e = DelayVolumeCellState()

        # set times
        d.set_time(parent.get_time())
        e.set_time(parent.get_time())

        # set volume
        d.set_volume(parent.get_volume() / 2.0)
        e.set_volume(parent.get_volume() / 2.0)

        # partition the states
        cdef np.ndarray[np.double_t,ndim=1] dstate = parent.get_state().copy()
        cdef np.ndarray[np.double_t,ndim=1] estate = parent.get_state().copy()
        cdef unsigned length = dstate.shape[0]

        cdef unsigned i = 0
        cdef unsigned amount = 0

        for i in range(length):
            amount = cyrandom.binom_rnd_f(dstate[i],0.5)
            dstate[i] = <double> amount
            estate[i] -= dstate[i]

        d.set_state(dstate)
        e.set_state(estate)

        # Partition the delay queue as well.

        cdef np.ndarray[object, ndim=1] array = parent.get_delay_queue().binomial_partition(0.5)

        cdef DelayQueue q1 = <DelayQueue> array[0]
        cdef DelayQueue q2 = <DelayQueue> array[1]

        d.set_delay_queue(q1)
        e.set_delay_queue(q2)

        # create return structure
        cdef np.ndarray ans = np.empty(2, dtype=np.object)
        ans[0] = d
        ans[1] = e

        return ans


cdef class CustomSplitter(VolumeSplitter):
    def __init__(self, split_function):
        self.split_function = split_function

    cdef np.ndarray partition(self, VolumeCellState parent):
        """
        Split the parent state into two VolumeCellState's according to a
        user-supplied function.

        :param parent: (VolumeCellState)
        :return: (np.ndarray) Length 2 array of VolumeCellStates.
        """
        return self.split_function(parent)

##################################################                ####################################################
######################################              CS                        ################################
#################################################                     ################################################

# Regular simulations with no volume or delay involved.
cdef class RegularSimulator:
    cdef SSAResult simulate(self, CSimInterface sim, np.ndarray timepoints):
        """
        Perform a simple regular stochastic simulation with no delay or volume involved.

        MUST BE SUBCLASSED.
        :param sim: (CSimInterface) The reaction system. Must have time initialized.
        :param timepoints: (np.ndarray) The time points (must be greater than the initial time).
        :return: (SSAResult) The simulation result.
        """
        raise NotImplementedError("simulate function not implemented for RegularSimulator")

    def py_simulate(self, CSimInterface sim, np.ndarray timepoints):
        #suggested that interfaces do some error checking on themselves to prevent kernel crashes.
        sim.check_interface()
        return self.simulate(sim,timepoints)


cdef void* global_simulator
cdef np.ndarray global_derivative_buffer

def py_global_derivative_buffer():
    global global_derivative_buffer
    return global_derivative_buffer

def py_set_globals(CSimInterface sim):
    global global_simulator
    global global_derivative_buffer
    global_simulator = <void*> sim
    global_derivative_buffer = np.empty(sim.py_get_num_species(),)

def rhs_global(np.ndarray[np.double_t,ndim=1] state, double t):
    global global_simulator
    global global_derivative_buffer
    (<CSimInterface>global_simulator).apply_repeated_rules(<double*> state.data,t)
    (<CSimInterface>global_simulator).calculate_deterministic_derivative( <double*> state.data,
                                                                         <double*> global_derivative_buffer.data, t)
    return global_derivative_buffer

def rhs_ode(double t, np.ndarray[np.double_t, ndim=1] state):
    return rhs_global(state,t)

cdef class DeterministicSimulator(RegularSimulator):
    """
    A class for implementing a deterministic simulator.
    """

    def __init__(self):
        self.atol = 1E-8
        self.rtol = 1E-8
        self.mxstep = 500000

    def py_set_tolerance(self, double atol, double rtol):
        self.set_tolerance(atol, rtol)

    def py_set_mxstep(self, unsigned mxstep):
        self.mxstep = mxstep

    def _helper_simulate(self, CSimInterface sim, np.ndarray timepoints):
        """
        This function allows for definition of the rhs function inside the simulation function. Otherwise, Cython
        does not allow closures inside cdef functions. Having a separate rhs function is also impossible because
        then the first argument becomes self.
        """
        cdef np.ndarray S = sim.get_update_array() + sim.get_delay_update_array()
        cdef np.ndarray x0 = sim.get_initial_state().copy()

        cdef unsigned num_species = S.shape[0]
        cdef unsigned num_reactions = S.shape[1]

        global global_simulator
        global global_derivative_buffer

        global_simulator = <void*> sim
        global_derivative_buffer = np.empty(num_species,)

        cdef unsigned index = 0
        cdef unsigned steps_allowed = 500
        cdef np.ndarray[np.double_t, ndim=2] results

        while True:
            results, full_output = odeint(rhs_global, x0, timepoints,atol=self.atol, rtol=self.rtol,
                                         mxstep=steps_allowed, full_output=True)

            if full_output['message'] == 'Integration successful.':
                if sim.get_number_of_rules() > 0:
                    for index in range(timepoints.shape[0]):
                        sim.apply_repeated_rules( &(results[index,0]),timepoints[index] )
                return SSAResult(timepoints,results)

            logging.info('odeint failed with mxstep=%d...' % (steps_allowed))

            # make the mxstep bigger if the user specified a bigger max
            if steps_allowed >= self.mxstep:
                break
            else:
                steps_allowed *= 10
                if steps_allowed > self.mxstep:
                    steps_allowed = self.mxstep

        sys.stderr.write('odeint failed entirely\n')

        return SSAResult(timepoints,results * np.nan)


    cdef SSAResult simulate(self, CSimInterface sim, np.ndarray timepoints):
        return self._helper_simulate(sim,timepoints)

cdef class DeterministicDilutionSimulator(RegularSimulator):
    """
    A class for implementing a deterministic simulator with dilution.
    """

    def __init__(self):
        self.atol = 1E-8
        self.rtol = 1E-8
        self.mxstep = 500000
        self.dilution_rate = 0.0

    def py_set_tolerance(self, double atol, double rtol):
        self.set_tolerance(atol, rtol)

    def py_set_mxstep(self, unsigned mxstep):
        self.mxstep = mxstep

    def py_set_dilution_rate(self, double rate):
        self.dilution_rate = rate

    def _helper_simulate(self, CSimInterface sim, np.ndarray timepoints):
        """
        This function allows for definition of the rhs function inside the simulation function. Otherwise, Cython
        does not allow closures inside cdef functions. Having a separate rhs function is also impossible because
        then the first argument becomes self.
        """
        cdef np.ndarray S = sim.get_update_array() + sim.get_delay_update_array()
        cdef np.ndarray x0 = sim.get_initial_state().copy()

        cdef unsigned num_species = S.shape[0]
        cdef unsigned num_reactions = S.shape[1]

        global global_simulator
        global global_derivative_buffer

        global_simulator = <void*> sim
        global_derivative_buffer = np.empty(num_species,)

        cdef unsigned index = 0
        cdef unsigned steps_allowed = 500
        cdef np.ndarray[np.double_t, ndim=2] results

        def rhs_dilution(np.ndarray[np.double_t,ndim=1] state, double t):
            return rhs_global(state, t) - self.dilution_rate*state


        while True:
            results, full_output = odeint(rhs_dilution, x0, timepoints,atol=self.atol, rtol=self.rtol,
                                         mxstep=steps_allowed, full_output=True)

            if full_output['message'] == 'Integration successful.':
                if sim.get_number_of_rules() > 0:
                    for index in range(timepoints.shape[0]):
                        sim.apply_repeated_rules( &(results[index,0]), timepoints[index] )
                return SSAResult(timepoints,results)

            logging.info('odeint failed with mxstep=%d...' % (steps_allowed))

            # make the mxstep bigger if the user specified a bigger max
            if steps_allowed >= self.mxstep:
                break
            else:
                steps_allowed *= 10
                if steps_allowed > self.mxstep:
                    steps_allowed = self.mxstep

        sys.stderr.write('odeint failed entirely\n')

        return SSAResult(timepoints,results * np.nan)


    cdef SSAResult simulate(self, CSimInterface sim, np.ndarray timepoints):
        return self._helper_simulate(sim,timepoints)


def perr(string):
    sys.stderr.write(string + '\n')
    sys.stderr.flush()

cdef class SSASimulator(RegularSimulator):
    """
    A class for implementing a stochastic SSA simulator.
    """
    cdef SSAResult simulate(self, CSimInterface sim, np.ndarray timepoints):
        cdef np.ndarray[np.double_t,ndim=1] c_timepoints = timepoints
        cdef np.ndarray[np.double_t,ndim=1] c_current_state = sim.get_initial_state().copy()
        cdef np.ndarray[np.double_t,ndim=2] c_stoich = sim.get_update_array() + sim.get_delay_update_array()
        cdef np.ndarray[np.double_t,ndim=2] c_delay_stoich = sim.get_delay_update_array()

        cdef unsigned num_species = c_stoich.shape[0]
        cdef unsigned num_reactions = c_stoich.shape[1]
        cdef unsigned num_timepoints = len(timepoints)

        cdef double final_time = timepoints[num_timepoints-1]

        cdef double current_time = sim.get_initial_time()
        cdef double dt = sim.get_dt()
        cdef double proposed_time = 0.0
        cdef double Lambda = 0.0

        cdef np.ndarray[np.double_t,ndim=2] c_results = np.zeros((num_timepoints, num_species),dtype=np.double)
        cdef np.ndarray[np.double_t,ndim=1] c_propensity = np.zeros(num_reactions)

        # Now do the SSA part
        cdef unsigned current_index = 0
        cdef unsigned reaction_choice = 0
        cdef unsigned species_index = 0


        while current_index < num_timepoints:
            # Compute propensity in place
            sim.apply_repeated_rules(<double*> c_current_state.data,current_time)
            sim.compute_stochastic_propensities(<double*> c_current_state.data, <double*> c_propensity.data,current_time)
            # Sample the next reaction time and update
            Lambda = cyrandom.array_sum(<double*> c_propensity.data,num_reactions)

            if Lambda == 0:
                current_time = current_time + dt
            else:
                current_time = current_time + cyrandom.exponential_rv(Lambda)

            # Update previous states
            while current_index < num_timepoints and c_timepoints[current_index] < current_time:
                for species_index in range(num_species):
                    c_results[current_index,species_index] = c_current_state[species_index]
                current_index += 1

            # Choose a reaction and update the state accordingly.
            if Lambda > 0:
                reaction_choice = cyrandom.sample_discrete(num_reactions, <double*> c_propensity.data , Lambda)

                for species_index in range(num_species):
                    c_current_state[species_index] += c_stoich[species_index,reaction_choice]

        return SSAResult(timepoints,c_results)

cdef class SafeModeSSASimulator(RegularSimulator):
    """
    A class for implementing a stochastic SSA simulator.
    """
    cdef SSAResult simulate(self, CSimInterface sim, np.ndarray timepoints):
        cdef np.ndarray[np.double_t,ndim=1] c_timepoints = timepoints
        cdef np.ndarray[np.double_t,ndim=1] c_current_state = sim.get_initial_state().copy()
        cdef np.ndarray[np.double_t,ndim=2] c_stoich = sim.get_update_array() + sim.get_delay_update_array()
        cdef np.ndarray[np.double_t,ndim=2] c_delay_stoich = sim.get_delay_update_array()

        cdef unsigned num_species = c_stoich.shape[0]
        cdef unsigned num_reactions = c_stoich.shape[1]
        cdef unsigned num_timepoints = len(timepoints)

        cdef double final_time = timepoints[num_timepoints-1]

        cdef double current_time = sim.get_initial_time()
        cdef double dt = sim.get_dt()
        cdef double proposed_time = 0.0
        cdef double Lambda = 0.0

        cdef np.ndarray[np.double_t,ndim=2] c_results = np.zeros((num_timepoints, num_species),dtype=np.double)
        cdef np.ndarray[np.double_t,ndim=1] c_propensity = np.zeros(num_reactions)

        # Now do the SSA part
        cdef unsigned current_index = 0
        cdef unsigned reaction_choice = 0
        cdef unsigned species_index = 0
        cdef unsigned reaction_index = 0


        while current_index < num_timepoints:
            # check for negative species
            for species_index in range(num_species):
                if c_current_state[species_index] < 0:
                    raise RuntimeError('Species %d has a negative value of %f, possible last reaction: %d' %
                            (species_index, c_current_state[species_index], reaction_choice))

            # Compute propensity in place
            sim.apply_repeated_rules(<double*> c_current_state.data,current_time)
            sim.compute_stochastic_propensities(<double*> c_current_state.data, <double*> c_propensity.data,current_time)
            # Sample the next reaction time and update
            Lambda = cyrandom.array_sum(<double*> c_propensity.data,num_reactions)

            for reaction_index in range(num_reactions):
                if c_propensity[reaction_index] < 0:
                    raise RuntimeError('Reaction %d has a negative propensity of %f' %
                          (reaction_index, c_propensity[reaction_index]))

            if Lambda == 0:
                current_time = current_time + dt
            else:
                current_time = current_time + cyrandom.exponential_rv(Lambda)

            # Update previous states
            while current_index < num_timepoints and c_timepoints[current_index] < current_time:
                for species_index in range(num_species):
                    c_results[current_index,species_index] = c_current_state[species_index]
                current_index += 1

            # Choose a reaction and update the state accordingly.
            if Lambda > 0:
                reaction_choice = cyrandom.sample_discrete(num_reactions, <double*> c_propensity.data , Lambda)

                for species_index in range(num_species):
                    c_current_state[species_index] += c_stoich[species_index,reaction_choice]

        return SSAResult(timepoints,c_results)




cdef class TimeDependentSSASimulator(RegularSimulator):
    """
    A class for implementing a stochastic SSA simulator.
    """
    cdef SSAResult simulate(self, CSimInterface sim, np.ndarray timepoints):
        cdef np.ndarray[np.double_t,ndim=1] c_timepoints = timepoints
        cdef np.ndarray[np.double_t,ndim=1] c_current_state = sim.get_initial_state().copy()
        cdef np.ndarray[np.double_t,ndim=2] c_stoich = sim.get_update_array() + sim.get_delay_update_array()
        cdef np.ndarray[np.double_t,ndim=2] c_delay_stoich = sim.get_delay_update_array()

        cdef unsigned num_species = c_stoich.shape[0]
        cdef unsigned num_reactions = c_stoich.shape[1]
        cdef unsigned num_timepoints = len(timepoints)

        cdef double final_time = timepoints[num_timepoints-1]

        cdef double current_time = sim.get_initial_time()
        cdef double proposed_time = 0.0
        cdef double Lambda = 0.0
        cdef double dt = sim.get_dt()

        cdef np.ndarray[np.double_t,ndim=2] c_results = np.zeros((num_timepoints, num_species),dtype=np.double)
        cdef np.ndarray[np.double_t,ndim=1] c_propensity = np.zeros(num_reactions)

        # Now do the SSA part
        cdef unsigned current_index = 0
        cdef unsigned reaction_choice = 0
        cdef unsigned species_index = 0

        cdef unsigned apply_reaction_flag = 0

        while current_index < num_timepoints:
            # Compute propensity in place
            sim.apply_repeated_rules(<double*> c_current_state.data,current_time)
            sim.compute_stochastic_propensities(<double*> c_current_state.data, <double*> c_propensity.data,current_time)
            # Sample the next reaction time and update
            Lambda = cyrandom.array_sum(<double*> c_propensity.data,num_reactions)

            if Lambda == 0:
                proposed_time = final_time + 1
            else:
                proposed_time = current_time + cyrandom.exponential_rv(Lambda)

            # update current time to either proposed time or the current time + dt.
            # in the first case, apply a reaction. in the latter case, do nothing.

            if proposed_time <= current_time + dt:
                current_time = proposed_time
                apply_reaction_flag = 1
            else:
                current_time += dt
                apply_reaction_flag = 0


            # Update previous states
            while current_index < num_timepoints and c_timepoints[current_index] < current_time:
                for species_index in range(num_species):
                    c_results[current_index,species_index] = c_current_state[species_index]
                current_index += 1

            # Choose a reaction and update the state accordingly.
            if apply_reaction_flag:
                reaction_choice = cyrandom.sample_discrete(num_reactions, <double*> c_propensity.data , Lambda)

                for species_index in range(num_species):
                    c_current_state[species_index] += c_stoich[species_index,reaction_choice]

        return SSAResult(timepoints,c_results)







cdef class DelaySimulator:
    """
    Interface class for defining a simulator with delay.
    """
    cdef DelaySSAResult delay_simulate(self, CSimInterface sim, DelayQueue dq, np.ndarray timepoints):
        """
        Perform a delay SSA given a reaction system and an initial set of queued reactions.

        MUST BE SUBCLASSED.
        :param sim: (CSimInterface) The reaction interface.
        :param dq: (DelayQueue) The delay queue.
        :param timepoints: (np.ndarray) The set of time points
        :return: (DelaySSAResult) The result.
        """

        raise NotImplementedError("Did not implement simulate function for DelaySimulator")

    def py_delay_simulate(self, CSimInterface sim, DelayQueue dq, np.ndarray timepoints):
        sim.check_interface()
        return self.delay_simulate(sim,dq,timepoints)


cdef class DelaySSASimulator(DelaySimulator):
    """
    A class for doing delay simulations using the stochastic simulation algorithm.
    """
    cdef DelaySSAResult delay_simulate(self, CSimInterface sim, DelayQueue q, np.ndarray timepoints):
        # Set up the needed variables in C

        cdef np.ndarray[np.double_t,ndim=1] c_timepoints = timepoints
        cdef np.ndarray[np.double_t,ndim=1] c_current_state = sim.get_initial_state().copy()
        cdef np.ndarray[np.double_t,ndim=2] c_stoich = sim.get_update_array()
        cdef np.ndarray[np.double_t,ndim=2] c_delay_stoich = sim.get_delay_update_array()

        cdef unsigned num_species = c_stoich.shape[0]
        cdef unsigned num_reactions = c_stoich.shape[1]
        cdef unsigned num_timepoints = c_timepoints.shape[0]


        cdef double current_time = sim.get_initial_time()
        q.set_current_time(current_time)
        cdef double final_time = c_timepoints[num_timepoints-1]
        cdef double proposed_time = 0.0
        cdef double future_reaction_time = 0.0
        cdef double Lambda = 0.0
        cdef np.ndarray[np.double_t,ndim=2] c_results = np.zeros((num_timepoints,num_species))
        cdef np.ndarray[np.double_t,ndim=1] c_propensity = np.zeros(num_reactions)
        cdef np.ndarray[np.double_t,ndim=1] c_q_rxn_amt = np.zeros(num_reactions)

        cdef unsigned current_index = 0
        cdef unsigned reaction_choice = 4294967295
        cdef unsigned species_index = 4294967295
        cdef unsigned reaction_index = 4294967295

        cdef unsigned move_to_queued_time = 0
        cdef next_queue_time = -10.0

        # Do the SSA part now

        while current_index < num_timepoints:
            # Compute the propensity in place
            sim.apply_repeated_rules(<double*> c_current_state.data, current_time)
            sim.compute_stochastic_propensities(<double*> (c_current_state.data), <double*> (c_propensity.data), current_time)
            Lambda = cyrandom.array_sum(<double*> (c_propensity.data), num_reactions)

            # Either we are going to move to the next queued time, or we move to the next reaction time.
            if Lambda == 0:
                proposed_time = final_time+1
            else:
                proposed_time = current_time + cyrandom.exponential_rv(Lambda)

            next_queue_time = q.get_next_queue_time()
            if next_queue_time < proposed_time:
                current_time = next_queue_time
                move_to_queued_time = 1
            else:
                current_time = proposed_time
                move_to_queued_time = 0

            # Update the results array with the state for the time period that we just jumped through.
            while current_index < num_timepoints and c_timepoints[current_index] < current_time:
                for species_index in range(num_species):
                    c_results[current_index,species_index] = c_current_state[species_index]
                current_index += 1

            # Now update the state accordingly.

            # IF the queue won, then add to the queue by doing the queue size thing again.
            if move_to_queued_time == 1:
                # find out how much of each reaction happened from the delay queue
                q.get_next_reactions(<double*> (c_q_rxn_amt.data))

                # update the state
                for reaction_index in range(num_reactions):
                    for species_index in range(num_species):
                        c_current_state[species_index] += c_q_rxn_amt[reaction_index]*c_delay_stoich[species_index,reaction_index]

                # advance the queue in time.
                q.advance_time()


            # if an actual reaction happened, do the reaction and maybe update the queue as well.
            else:
                # select a reaction
                reaction_choice = cyrandom.sample_discrete(num_reactions, <double*> c_propensity.data, Lambda )
                # Compute the delay for the reaction.
                computed_delay = sim.compute_delay(<double *> c_current_state.data, reaction_choice)
                # Do the reaction's initial stoichiometry.
                for species_index in range(num_species):
                    c_current_state[species_index] += c_stoich[species_index,reaction_choice]

                # Check if the delay is real, and then if needed add it to the queue.
                if computed_delay > 0.0:
                    q.add_reaction(current_time+computed_delay,reaction_choice,1.0)
                else:
                    for species_index in range(num_species):
                        c_current_state[species_index] += c_delay_stoich[species_index,reaction_choice]

        # Now need to re-align the final delay queue properly so that the first time comes first etc.
        return DelaySSAResult(c_timepoints, c_results,q)



cdef class VolumeSimulator:
    """
    Interface class for doing volume simulations.
    """
    cdef VolumeSSAResult volume_simulate(self, CSimInterface sim, Volume v, np.ndarray timepoints):
        """
        This function takes in a reaction system in sim, and a volume, and simulates over timepoints.

        MUST BE SUBCLASSED.
        :param sim: (CSimInterface) The reaction system to simulate. Must be initialized to same time as v.
        :param v: (Volume) A volume model for the system. This must have been initialized to same initial time as sim.
        :param timepoints: (np.ndarray) Set of time points. First point should be >= initial time.
        :return: The simulation result.
        """
        raise NotImplementedError("Did not implement simulation function for Volume Simulator")

    def py_volume_simulate(self, CSimInterface sim, Volume v, np.ndarray timepoints):
        sim.check_interface()
        return self.volume_simulate(sim,v,timepoints)


cdef class VolumeSSASimulator(VolumeSimulator):
    """
    Volume SSA implementation.
    """
    cdef VolumeSSAResult volume_simulate(self, CSimInterface sim, Volume v, np.ndarray timepoints):
        """
        Implements a volume simulation using a regular SSA with volume updates as a dt interval set by calling
        sim.set_dt()
        """
        # Set up the needed variables in C

        cdef np.ndarray[np.double_t,ndim=1] c_timepoints = timepoints.copy()
        cdef np.ndarray[np.double_t,ndim=1] c_current_state = sim.get_initial_state().copy()
        cdef np.ndarray[np.double_t,ndim=2] c_stoich = sim.get_update_array() + sim.get_delay_update_array()

        cdef unsigned num_species = c_stoich.shape[0]
        cdef unsigned num_reactions = c_stoich.shape[1]
        cdef unsigned num_timepoints = len(timepoints)

        cdef double current_time = sim.get_initial_time()
        cdef double final_time = c_timepoints[num_timepoints-1]
        cdef double proposed_time = 0.0
        cdef double Lambda = 0.0
        cdef np.ndarray[np.double_t,ndim=2] c_results = np.zeros((num_timepoints,num_species))
        cdef np.ndarray[np.double_t,ndim=1] c_propensity = np.zeros(num_reactions)
        cdef np.ndarray[np.double_t,ndim=1] c_volume_trace = np.zeros(num_timepoints,)

        cdef unsigned current_index = 0
        cdef unsigned reaction_choice = 4294967295
        cdef unsigned species_index = 4294967295

        cdef double delta_t = sim.get_dt()
        cdef double next_queue_time = delta_t+current_time
        cdef unsigned move_to_queued_time = 0

        cdef double current_volume = v.get_volume()
        cdef unsigned cell_divided = 0


        # Do the SSA part now

        while current_index < num_timepoints:
            # Compute the propensity in place
            sim.apply_repeated_volume_rules(<double*> c_current_state.data, current_volume, current_time)
            sim.compute_stochastic_volume_propensities(<double*> (c_current_state.data), <double*> (c_propensity.data),
                                            current_volume, current_time)
            Lambda = cyrandom.array_sum(<double*> (c_propensity.data), num_reactions)

            # Either we are going to move to the next queued time, or we move to the next reaction time.
            if Lambda == 0:
                proposed_time = final_time+1
            else:
                proposed_time = current_time + cyrandom.exponential_rv(Lambda)
            if next_queue_time < proposed_time:
                current_time = next_queue_time
                next_queue_time += delta_t
                move_to_queued_time = 1
            else:
                current_time = proposed_time
                move_to_queued_time = 0

            # Update the results array with the state for the time period that we just jumped through.
            while current_index < num_timepoints and c_timepoints[current_index] < current_time:
                for species_index in range(num_species):
                    c_results[current_index,species_index] = c_current_state[species_index]
                c_volume_trace[current_index] = current_volume
                current_index += 1

            # Now update the state accordingly.

            # IF the queue won, then update the volume and continue on or stop if the cell divided.
            if move_to_queued_time == 1:
                # Update the volume
                current_volume += v.get_volume_step(<double*>(c_current_state.data), <double*> sim.get_param_values(),
                                                    current_time, current_volume, delta_t)
                v.set_volume(current_volume)

                # IF the cell divided, just return and bounce from here!!!!
                if v.cell_divided(<double*>(c_current_state.data), <double*> sim.get_param_values(),
                                  current_time,current_volume,delta_t):
                    cell_divided = True
                    break

            # if an actual reaction happened, do the reaction and maybe update the queue as well.
            else:
                # select a reaction
                reaction_choice = cyrandom.sample_discrete(num_reactions, <double*> c_propensity.data, Lambda )

                # Do the reaction's initial stoichiometry.
                for species_index in range(num_species):
                    c_current_state[species_index] += c_stoich[species_index,reaction_choice]

        # Now need to re-align the final delay queue properly so that the first time comes first etc.
        # return DelaySSAResult(c_results,np.roll(c_delay_queue,delay_queue_length-queue_index,axis=1))

        if cell_divided:
            c_timepoints = c_timepoints[:(current_index)]
            c_volume_trace = c_volume_trace[:(current_index)]
            c_results = c_results[:current_index,:]

        cdef VolumeSSAResult vsr = VolumeSSAResult(c_timepoints,c_results,c_volume_trace,cell_divided)
        vsr.set_volume_object(v)

        return vsr


cdef class DelayVolumeSimulator:
    """
    Interface class for doing simulations with delay and volume.
    """
    cdef DelayVolumeSSAResult delay_volume_simulate(self, CSimInterface sim, DelayQueue q,
                                                    Volume v, np.ndarray timepoints):
        """
        This function performs a simulation with delay and volume present in the system. The parameters must pre
         initialized to the same time with the timepoints being greater than or equal to the initial time.

         MUST BE SUBCLASSED BY IMPLEMENTATIONS
        :param sim: (CSimInterface) Reaction interface.
        :param q: (DelayQueue) The delay queue.
        :param v: (Volume) Volume model.
        :param timepoints: (np.ndarray) time points
        :return: The simulation result.
        """
        raise NotImplementedError("Did not implement simulation function for delay/volume simulator.")

    def py_delay_volume_simulate(self,CSimInterface sim, DelayQueue q, Volume v, np.ndarray timepoints):
        sim.check_interface()
        return self.delay_volume_simulate(sim,q,v,timepoints)

cdef class DelayVolumeSSASimulator(DelayVolumeSimulator):
    """
    SSA implementation for doing simulations with delay and volume.
    """
    cdef DelayVolumeSSAResult delay_volume_simulate(self, CSimInterface sim, DelayQueue q,
                                                    Volume v, np.ndarray timepoints):

        # Set up the needed variables in C

        cdef np.ndarray[np.double_t,ndim=1] c_timepoints = timepoints.copy()
        cdef np.ndarray[np.double_t,ndim=1] c_current_state = sim.get_initial_state().copy()
        cdef np.ndarray[np.double_t,ndim=2] c_stoich = sim.get_update_array()
        cdef np.ndarray[np.double_t,ndim=2] c_delay_stoich = sim.get_delay_update_array()

        cdef unsigned num_species = c_stoich.shape[0]
        cdef unsigned num_reactions = c_stoich.shape[1]
        cdef unsigned num_timepoints = c_timepoints.shape[0]

        cdef double current_time = sim.get_initial_time()
        cdef double final_time = c_timepoints[num_timepoints-1]
        cdef double proposed_time = 0.0
        cdef double Lambda = 0.0
        cdef np.ndarray[np.double_t,ndim=2] c_results = np.zeros((num_timepoints,num_species))
        cdef np.ndarray[np.double_t,ndim=1] c_propensity = np.zeros(num_reactions)
        cdef np.ndarray[np.double_t,ndim=1] c_volume_trace = np.zeros(num_timepoints,)

        cdef unsigned current_index = 0
        cdef unsigned reaction_index = 4294967295
        cdef unsigned reaction_choice = 4294967295
        cdef unsigned species_index = 4294967295

        cdef double delta_t = sim.get_dt()
        cdef double next_vol_time = delta_t+current_time
        cdef double current_volume = v.get_volume()
        cdef unsigned cell_divided = 0

        cdef double next_queued_reaction_time = 1000000000.0
        cdef np.ndarray[np.double_t,ndim=1] c_delay_rxns = np.zeros(num_reactions)
        cdef double computed_delay = 0.0

        cdef unsigned short step_type = 0 # this variable is 0 if reaction, 1 if volume, 2 if delayed reaction

        # Do the SSA part now

        while current_index < num_timepoints:
            # Compute the propensity in place
            sim.apply_repeated_volume_rules(<double*> c_current_state.data, current_volume, current_time)
            sim.compute_stochastic_volume_propensities(<double*> (c_current_state.data), <double*> (c_propensity.data),
                                            current_volume, current_time)
            Lambda = cyrandom.array_sum(<double*> (c_propensity.data), num_reactions)

            # Either we are going to move to the next volume time, delay time, or reaction time.

            if Lambda == 0:
                proposed_time = final_time + 1
            else:
                proposed_time = current_time + cyrandom.exponential_rv(Lambda)
            next_queued_reaction_time = q.get_next_queue_time()

            if proposed_time < next_vol_time and proposed_time < next_queued_reaction_time:
                current_time = proposed_time
                step_type = 0
            elif next_vol_time < next_queued_reaction_time:
                current_time = next_vol_time
                next_vol_time += delta_t
                step_type = 1
            else:
                current_time = next_queued_reaction_time
                step_type = 2


            # Update the results array with the state for the time period that we just jumped through.
            while current_index < num_timepoints and c_timepoints[current_index] < current_time:
                for species_index in range(num_species):
                    c_results[current_index,species_index] = c_current_state[species_index]
                c_volume_trace[current_index] = current_volume
                current_index += 1

            # Handle the different cases here.

            # 1. If we had a regular reaction fire, add it and potentially add it to the queue.
            if step_type == 0:
                # select a reaction
                reaction_choice = cyrandom.sample_discrete(num_reactions, <double*> c_propensity.data, Lambda )

                # Compute the delay for the reaction.
                computed_delay = sim.compute_delay(<double *> c_current_state.data, reaction_choice)

                # Do the reaction's initial stoichiometry.
                for species_index in range(num_species):
                    c_current_state[species_index] += c_stoich[species_index,reaction_choice]

                # Do the delayed stoichiometry if needed.

                # Check if the delay is real, and then if needed add it to the queue.
                if computed_delay > 0.0:
                    q.add_reaction(current_time+computed_delay,reaction_choice,1.0)
                else:
                    for species_index in range(num_species):
                        c_current_state[species_index] += c_delay_stoich[species_index,reaction_choice]

            # 2. If we have a volume step occur instead.
            elif step_type == 1:
                # Update the volume
                current_volume += v.get_volume_step(<double*>(c_current_state.data),
                                                    <double*> sim.get_param_values(),
                                                    current_time, current_volume, delta_t)
                v.set_volume(current_volume)

                # IF the cell divided, just return and bounce from here!!!!
                if v.cell_divided(<double*>(c_current_state.data), <double*> sim.get_param_values(),
                                  current_time,current_volume,delta_t):
                    cell_divided = True
                    break
            # 3. If we have a delay reaction come on instead.
            else:
                # find out how much of each reaction happened from the delay queue
                q.get_next_reactions(<double*> (c_delay_rxns.data))

                # update the state
                for reaction_index in range(num_reactions):
                    for species_index in range(num_species):
                        c_current_state[species_index] += c_delay_rxns[reaction_index]*c_delay_stoich[species_index,reaction_index]

                # advance the queue in time.
                q.advance_time()

        # must account for cell_divided part by truncating up to the time where we were at when it divided.
        if cell_divided:
            c_timepoints = c_timepoints[:(current_index)]
            c_volume_trace = c_volume_trace[:(current_index)]
            c_results = c_results[:current_index,:]

        return DelayVolumeSSAResult(c_timepoints,c_results,c_volume_trace,q,cell_divided)


#A wrapper function to allow easy simulation of Models
def py_simulate_model(timepoints, Model = None, Interface = None, stochastic = False, 
                    delay = None, safe = False, volume = False, return_dataframe = True):
    #Check model and interface
    if Model is None and Interface is None:
        raise ValueError("py_simulate_model requires either a Model or CSimInterface to be passed in.")
    elif not Model is None and not Interface is None:
        raise ValueError("py_simulate_model requires either a Model OR a CSimInterface to be passed in. Not both.")
    elif Interface is None:
        if safe:
            Interface = SafeModelCSimInterface(Model)
        else:
            Interface = ModelCSimInterface(Model)
    elif not Interface is None and safe:
        logging.info("Cannot gaurantee that the interface passed in is safe. Simulating anyway.")

    #Create Volume (if necessary)
    if isinstance(volume, Volume):
        pass
    elif volume == False:
        v = None
    else:
        if volume == True:
            v = Volume()
            v.py_set_volume(1.0)
        else:
            try:
                if volume > 0:
                    v = Volume()
                    v.py_set_volume(volume)
            except TypeError:
                warnings.warn("Caught TypeError: invalid volume keyword. Setting volume to 1.")
                v = Volume()
                v.py_set_volume(1.0)
            
    #Create Simulator and Simulate    
    if delay:
        if not stochastic:
            warnings.warn("Delay Simulators only exist for stochastic simulations. Defaulting to Stochastic simulation")

        q = ArrayDelayQueue.setup_queue(Interface.py_get_num_reactions(),len(timepoints),timepoints[1]-timepoints[0])
        if v == None:
            Sim = DelaySSASimulator()
            result = Sim.py_delay_simulate(Interface, q, timepoints)
        else :
            Sim = DelayVolumeSimulator()
            result = Sim.py_delay_volume_simulate(Interface, q, v, timepoints)
    elif stochastic:
        if v == None:
            Sim = SSASimulator()
            result = Sim.py_simulate(Interface, timepoints)

        else:
            Sim = VolumeSSASimulator()
            result = Sim.py_volume_simulate(Interface, v, timepoints)
    else:
        if v != None:
            logging.info("uncessary volume parameter for deterministic simulation.")
        Sim = DeterministicSimulator()
        Interface.py_prep_deterministic_simulation()
        result = Sim.py_simulate(Interface, timepoints)

    if return_dataframe:
        return result.py_get_dataframe(Model = Model)
    else:
        return result


cdef list propagate_cell(ModelCSimInterface sim, VolumeCellState cell, double end_time,
                               VolumeSimulator vsim, VolumeSplitter vsplit):
    cells_to_simulate = [cell]
    cell.set_volume_object(cell.get_volume_object().copy())
    cells_to_return = []
    index = 0

    cdef VolumeCellState cs,d1,d2
    cdef np.ndarray timepoints
    cdef VolumeSSAResult r
    cdef np.ndarray daughters

    while index < len(cells_to_simulate):
        cs = cells_to_simulate[index]
        index += 1

        sim.set_initial_state(cs.get_state())
        sim.set_initial_time(cs.get_time())
        timepoints = np.linspace(cs.get_time(),end_time,int( (end_time-cs.get_time())/sim.get_dt() )+10)
        r = vsim.volume_simulate(sim, cs.get_volume_object(), timepoints)

        cs = r.get_final_cell_state()
        if cs.get_time() >= end_time - sim.get_dt() - 1E-8:
            cs.set_time(end_time)
            cells_to_return.append(cs)

        else:
            daughters = vsplit.partition(cs)
            d1 = <VolumeCellState>(daughters[0])
            d2 = <VolumeCellState>(daughters[1])

            d1.set_volume_object(cs.get_volume_object().copy())
            d1.get_volume_object().initialize(<double*> d1.get_state().data,sim.get_param_values(),
                                              d1.get_time(), d1.get_volume())
            d2.set_volume_object(cs.get_volume_object().copy())
            d2.get_volume_object().initialize(<double*> d2.get_state().data,sim.get_param_values(),
                                              d2.get_time(), d2.get_volume())
            cells_to_simulate.append(d1)
            cells_to_simulate.append(d2)

    return cells_to_return

def py_propagate_cell(ModelCSimInterface sim, VolumeCellState cell, double end_time,
                               VolumeSimulator vsim, VolumeSplitter vsplit):
    """
    Move a cell forward in time including division and partitioning
    :param sim: (ModelCSimInterface) the model to use to simulate
    :param cell: (VolumeCellState) the initial cell state to start from
    :param end_time: (double) the final time of the simulation
    :param vsim: (VolumeSimulator) simulator object that handles volume and division
    :param vsplit: (VolumeSplitter) splitter for cell division
    :return: (list) VolumeCellState corresponding to this cell (or its lineage) in the future. If
              the cell divided there will be muliple entries
    """
    return propagate_cell(sim,cell,end_time, vsim,vsplit)

cdef list propagate_cells(ModelCSimInterface sim, list cells, double end_time,
                          VolumeSimulator vsim, VolumeSplitter vsplit):
    cdef VolumeCellState cell
    cdef unsigned index
    cdef unsigned size = len(cells)
    cdef list cells_to_return = []

    for index in range(size):
        cell = cells[index]
        cells_to_return.extend(propagate_cell(sim,cell,end_time,vsim,vsplit))

    return cells_to_return


def py_propagate_cells(ModelCSimInterface sim, list cells, double end_time,
                          VolumeSimulator vsim, VolumeSplitter vsplit):
    """
    Move a population of cells forward in time.
    :param sim: (ModelCSimInterface) the model to use to simulate
    :param cells: list of (VolumeCellState) the initial cell states to start from
    :param end_time: (double) the final time of the simulation
    :param vsim: (VolumeSimulator) simulator object that handles volume and division
    :param vsplit: (VolumeSplitter) splitter for cell division
    :return: (list) VolumeCellState corresponding to this cells in the future. If
              the cells divided there will be muliple entries for each original cell
    """
    return propagate_cells(sim,cells,end_time,vsim,vsplit)


