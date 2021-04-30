from bioscrape.lineage import *
import numpy as np
import pytest
import copy
import test_utils
from random import random

SEED = 135334440

def test_isheap_check():
	heap_nums = [
		[22, 25, 46, 38, 34, 63, 67, 101, 59],
		[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
		[5, 10, 15, 11, 12, 16, 17],
		[5, 15, 10, 16, 17, 11, 12],
	]

	not_heap_nums = [
		[63, 34, 46, 67, 59, 63],
		[25, 22, 46, 38, 34, 63, 67, 101, 59],
		[22, 25, 46, 38, 101, 63, 67, 101, 34],
		[12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1],
		[5, 15, 10, 11, 12, 16, 17]
	]

	heaps = [[LineageVolumeCellState(t0 = i) for i in ns] for ns in heap_nums]
	not_heaps = [[LineageVolumeCellState(t0 = i) for i in ns] 
												for ns in not_heap_nums]

	for heap in heaps:
		q = CappedStateQueue(12)
		q.py_set_underlying_array(heap)
		assert q.is_still_heap(), f"CappedStateQueue incorrectly claimed that {q} is NOT a heap"

	for fake_heap in not_heaps:
		q = CappedStateQueue(12)
		q.py_set_underlying_array(fake_heap)
		assert not q.is_still_heap(), f"CappedStateQueue incorrectly claimed that {q} IS a heap"


def test_construction():
	test_utils.set_seed(SEED)
	n_iters = 100

	for i in range(n_iters):
		q_size = int(np.random.random()*100)
		q = CappedStateQueue(q_size)
		for n in range(q_size):
			q.py_push_event(LineageVolumeCellState(t0 = np.random.random()))
		assert q.is_still_heap(), f"Pushing onto CappedStateQueue produced non-heap {q}."

def test_capping():
	test_utils.set_seed(SEED)
	n_iters = 100
	n_cap_pushes = 101

	for i in range(n_iters):
		q_size = int(np.random.random()*100)
		q = CappedStateQueue(q_size)
		# Fill queue
		for n in range(q_size):
			q.py_push_event(LineageVolumeCellState(t0 = n * 1/q_size))
		# Push more stuff on the queue
		for n in range(n_cap_pushes):
			new_num = np.random.random()
			q.py_push_event(LineageVolumeCellState(t0 = new_num))
			assert len(q) == q_size, f"Pushing {new_num} onto a full CappedStateQueue changed its size ({q})."
			assert q.is_still_heap(), f"Pushing {new_num} onto a full CappedStateQueue produced non-heap {q}."

def test_removal():
	test_utils.set_seed(SEED)
	n_iters = 100

	for i in range(n_iters):
		# Build a queue and fill it
		q_size = int(np.random.random()*100)
		nums = [np.random.random() for j in range(q_size)]
		sorted_nums = copy.copy(nums)
		sorted_nums.sort()
		q = CappedStateQueue(q_size)
		for n in nums:
			q.py_push_event(LineageVolumeCellState(t0 = n))
		assert q.is_still_heap(), f"During setup, {q} somehow became not a heap."

		# Check that the queue returns the right numbers in the right order
		# (and maintains heapness).
		for n in sorted_nums:
			e = q.py_pop_event()
			assert e.py_get_initial_time() == n, \
				f"CappedStateQueue {q} should have returned {n}, returned {e.py_get_initial_time()}."
			assert q.is_still_heap(), \
				f"Popping from CappedStateQueue produced non-heap {q}."