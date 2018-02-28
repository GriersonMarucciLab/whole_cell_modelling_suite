from abc import ABCMeta, abstractmethod
import random
from connections import Bc3
from batch_jobs import JobSubmission, ManageSubmission, SimpleManageSubmission
import operator
import numpy as np
#from multiprocessing import Pool
#import contextlib
from concurrent.futures import ProcessPoolExecutor as Pool

class MGA(metaclass=ABCMeta):
	"""This is an abstract class that all multi-generation algorithms inherit from."""
	def __init__(self, dict_of_cluster_instances, MGA_name, checkStop_func, runSimulations_func, getGenerationName_func, func_name_to_variables_dict):
		"""
				
				"""
		self.cluster_instances_dict = dict_of_cluster_instances
		self.generation_counter = None
		self.MGA_name = MGA_name
                self.checkStop_func
                self.runSimulations_func
                self.getGenerationName_func
                self.func_name_to_variables_dict

	# instance methods
	def run(self):
		if self.generation_counter == None:
			self.generation_counter = 0
		while self.checkStop_func() != True:
			self.run_sim_out = self.runSimulations()
			self.generation_counter += 1

	def checkStop(self, in_dict):
		if 'max_generation' in in_dict:
			if self.generation_counter <= in_dict['max_generation'][0]:
				output = False
			else:
				output = True
		else:
			raise ValueError('You have not input a valid type of stopping the multi-generation algorithm. Here in_dict[\'stop_type\'] =', in_dict['stop_type'])

		return output
 
	# abstract methods
	@abstractmethod
	def getGenerationName(self):
		pass

	@abstractmethod
	def getNewGeneration(self):
		pass

	@abstractmethod
	def runSimulations(self):
		pass
