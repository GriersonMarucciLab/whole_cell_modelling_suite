import time
import subprocess
import os
import datetime
#from multiprocessing import Pool
#import contextlib
from concurrent.futures import ProcessPoolExecutor as Pool
class JobSubmission():
	def __init__(self, submission_name, cluster_connection, ko_name_to_set_dict, simulation_output_path, errorfile_path, outfile_path, runfiles_path, reps_of_unique_ko, wholecell_model_master_dir):
		self.submission_name = submission_name
		self.submission_file_name = self.submission_name + '_sub_file.sh'
		self.temp_storage_path = '/space/oc13378/myprojects/github/uob/wc/mg/oc2/whole_cell_modelling_suite/tmp_storage'
		self.ko_sets_file_name = 'ko_sets.list'
		self.ko_set_names_file_name = 'ko_set_names.list'
		self.unique_job_name = self.createUniqueJobName(self.submission_name + '_')
		os.makedirs(self.temp_storage_path + '/' + self.unique_job_name)
		self.temp_storage_path = self.temp_storage_path + '/' + self.unique_job_name
		self.cluster_connection = cluster_connection
		self.ko_details_dict = ko_name_to_set_dict
		self.simulation_output_path = simulation_output_path
		self.errorfile_path = errorfile_path
		self.outfile_path = outfile_path
		self.runfiles_path = runfiles_path
		self.ko_set_tmp_storage = self.temp_storage_path + '/' + self.ko_sets_file_name
		self.ko_set_names_tmp_storage = self.temp_storage_path + '/' + self.ko_set_names_file_name
		self.submission_script_tmp_storage = self.temp_storage_path + '/' + self.submission_file_name
		self.order_of_keys_written_to_file = self.cluster_connection.convertKosAndNamesToFile(self.ko_details_dict, self.ko_set_tmp_storage, self.ko_set_names_tmp_storage)
		self.resource_usage_dict = self.cluster_connection.createStandardKoSubmissionScript(self.submission_script_tmp_storage, self.submission_name, len(ko_name_to_set_dict.values()), self.runfiles_path + '/' + self.ko_set_names_file_name, reps_of_unique_ko, wholecell_model_master_dir, self.simulation_output_path, self.runfiles_path + '/' + self.ko_sets_file_name, self.outfile_path, self.errorfile_path)
		self.cluster_job_number = None
		self.list_of_folders_to_make_on_cluster = [self.simulation_output_path, self.errorfile_path, self.outfile_path, self.runfiles_path]
		self.cluster_job_number = None
		self.time_of_submission = None

	def prepareForSubmission(self):
		# create folders on the cluster
		self.list_of_folders_to_make_on_cluster
		makedir_commands_list = ['mkdir -p ' + dir for dir in self.list_of_folders_to_make_on_cluster]
		makedir_output_dict = self.cluster_connection.checkSuccess(self.cluster_connection.sendCommand, makedir_commands_list)

		# transfer all neccessary files to the cluster
		sendFiles_output_dict = self.cluster_connection.checkSuccess(self.cluster_connection.rsyncFile, self.ko_set_tmp_storage + ' ' + self.ko_set_names_tmp_storage + ' ' + self.submission_script_tmp_storage, self.runfiles_path)

		return [makedir_output_dict, sendFiles_output_dict]

	def submitJobToCluster(self):
		submit_command = self.cluster_connection.submit_command + ' ' + self.runfiles_path + '/' + self.submission_file_name
		now = datetime.datetime.now()
		submit_job_ouput_dict = self.cluster_connection.checkSuccess(self.cluster_connection.sendCommand, [submit_command])
		self.cluster_job_number = self.cluster_connection.getJobIdFromSubStdOut(submit_job_ouput_dict['stdout'])
		self.time_of_submission = {'day': now.day, 'month': now.month, 'year': now.year}
		# tidy up the tmp storage area if the submission was successful
		if submit_job_ouput_dict['return_code'] == 0:
			tmp_return_code = subprocess.call("rm -r " + self.temp_storage_path, shell=True)
			if tmp_return_code != 0:
				print("WARNING!!!! Could not remove the temporary files from ", self.temp_storage_path, " please fix this problem ASAP since if it carries on repeating it will filll the entire computer up until it breaks!")

		return submit_job_ouput_dict

	def createUniqueJobName(self, prefix):
		valid_name = False
		while valid_name == False:
#			unique_name_end = str(print("%.9f" % time.time()))
			unique_name_end = str('{:.9f}'.format(time.time()))
			unique_name_end = unique_name_end.replace(".","")
			unique_path = self.temp_storage_path + '/' + prefix + unique_name_end
			if os.path.isdir(unique_path) == False:
				valid_name = True

		return prefix + unique_name_end

class ManageSubmission():
	def __init__(self, input_tuple, test_mode = False):
		(submission_instance, misc_data) = input_tuple
		self.submission = submission_instance
		self.misc_data = misc_data # thiss is added just in case there needs to be some kind of unforeseen data transfer between library classes (fro example to save connecting to static.db to often the MGA class just passes the gene_code_to_id_dict to ManageSubmission through misc_data)
		self.simulation_data_dict = {}
		if test_mode == True:
			print("WARNING: This is in TEST mode so no files will be transfered and no job will be submitted.")
			self.submission.time_of_submission = {}
			self.submission.time_of_submission['day'] = 8
			self.submission.time_of_submission['month'] = 1
			self.submission.time_of_submission['year'] = 2018
		else:
			outs =  self.submission.prepareForSubmission()
			list_of_return_codes = [output['return_code'] for output in outs]
			if sum(list_of_return_codes) > 0:
				raise ValueError('There has been a problem preparing some jobs for submission. The return codes from JobSubmission.prepareForSubmission() are: ', list_of_return_codes, '. The submission file name is ', self.submission.submission_file_name, ' and the temp_storage_path was ', self.submission.temp_storage_path)

			submission_out_dict = self.submission.submitJobToCluster()
			if submission_out_dict['return_code'] > 0:
				raise ValueError('There has been a problem submitting a job. The return code from JobSubmission.submitJobToCluster() is: ', submission_out_dict['return_code'], ". The submission file name is ", self.submission.submission_file_name, " and the temp_storage_path was ", self.submission.temp_storage_path)

		print("preparing ko dict!")
		self.data_dict = self.prepareDictForKoDbSubmission()
		# there is no output given here because if any return codes come back none xero it raises an error and no other output is needed. If any output is needed for debugging it is suggested that you do that inside the self.monitor submission function if possible.
		if test_mode == True:
			pass
		else:
			print("Starting monitor submission function!")
			self.monitorSubmission(self.submission)

		return 

	def prepareDictForKoDbSubmission(self):
		# initialise dict that KoDb class uses to submit data to the database
		data_dict = {'people': {'id': None, 'first_name': self.submission.cluster_connection.forename_of_user, 'last_name': self.submission.cluster_connection.surename_of_user, 'user_name': self.submission.cluster_connection.user_name}, 'batchDescription': {'id': None, 'name': self.submission.submission_name, 'simulation_initiator_id': None, 'description': 'This batch is automatically generated by Oliver Chalkleys whole-cell modelling suite and has the name: ' + self.submission.submission_name, 'simulation_day': self.submission.time_of_submission['day'], 'simulation_month': self.submission.time_of_submission['month'], 'simulation_year': self.submission.time_of_submission['year'], 'cluster_info': self.submission.cluster_connection.information_about_cluster}, 'koIndex': {'id': None, 'number_of_kos': None}, 'simulationDetails': {'id': None, 'ko_index_id': None, 'batch_description_id': None, 'average_growth_rate': None, 'time_when_pinchedDiameter_is_first_zero': None}}

		return data_dict

	def prepareSimulationDictForKoDbSubmission(self, directory_of_specific_ko_repetition):
#		codes not IDs! simulation_dict = {(1, 2, 3): [(13, 13446), (17, 14532)], (1, 7, 39, 301): [(321, 251637), (536127, 36128), (5632, 637)], (34, 432, 12, 19, 234): [(432, 654), (432, 432), (324, 234), (543, 675)]}
		ko_name_to_ko_set_dict = self.submission.ko_details_dict 
		# the ko name is used as the name of the directory (remember that the last dir name is the repetition number so we want the second last)
		dirs_as_list = directory_of_specific_ko_repetition.split("/")
		sim_name1 = dirs_as_list[-2]
		sim_name2 = dirs_as_list[-1]

		# get the pinch time and average growth rate of this simulation in a tuple
		cmds_to_get_sim_details = self.submission.cluster_connection.activate_venv_list
		cmds_to_get_sim_details = cmds_to_get_sim_details + ['cd ' + directory_of_specific_ko_repetition, 'python -c "import pandas as pd;basic_summary = pd.io.pickle.read_pickle(\'basic_summary_' + sim_name1 + '_' + sim_name2 + '.pickle\');growth_mean = basic_summary[\'metabolicReaction_growth\'].mean();pinch_data = basic_summary[\'geometry_pinchedDiameter\'];pinched_times = list(pinch_data[pinch_data == 0].index);first_pinch = (pinched_times[0] if pinched_times else 0);print(growth_mean);print(first_pinch)"']
		sim_detail_out = self.submission.cluster_connection.checkSuccess(self.submission.cluster_connection.sendCommand, cmds_to_get_sim_details)
		list_of_mgr_and_pinch = sim_detail_out['stdout'].split("\n")
		del list_of_mgr_and_pinch[-1]
		tuple_of_sim_details = (list_of_mgr_and_pinch[0], list_of_mgr_and_pinch[1])

		# because this is being done in parallel we cannot update the self.simulation_data_dict until the parallel part is complete so we create a temp variable to store it in and pass it back to the create pandas function and then back through to monitor submission function where we are finall out of the  parallel bit and we can undate self.simulation_data_dict.
		tmp_sim_data_dict = {}
		if ko_name_to_ko_set_dict[sim_name1] not in tmp_sim_data_dict:
			tmp_sim_data_dict[ko_name_to_ko_set_dict[sim_name1]] = []

		tmp_sim_data_dict[ko_name_to_ko_set_dict[sim_name1]].append(tuple_of_sim_details)
		print("tmp_sim_data_dict = ", tmp_sim_data_dict)

		return tmp_sim_data_dict

	def monitorSubmission(self, submission):
		print("In monitorSubmission!")
		all_outputs_dict = {}
		time.sleep(3600) # wait 3,600 seconds or 1 hour because definitely nothing will finish in the first hour
#		time.sleep(15) # wait 14,400 seconds or 4 hours because definitely nothing will finish in the first four hours
		#useful info about the job
		dir_to_ko_set_dict = submission.ko_details_dict
		job_resource_allocation_dict = submission.resource_usage_dict
		no_of_arrays = job_resource_allocation_dict['no_of_arrays']
		no_of_unique_kos_per_array_job = job_resource_allocation_dict['no_of_unique_kos_per_array_job']
		no_of_repetitions_of_each_ko = job_resource_allocation_dict['no_of_repetitions_of_each_ko']
		# get details about the queue
		queue_output_dict = submission.cluster_connection.checkSuccess(submission.cluster_connection.checkQueue, submission.cluster_job_number)
		stdout_tmp = queue_output_dict['stdout'].split("\n")
		del stdout_tmp[-1]
		jobs_still_running_list = [int(digit) for digit in stdout_tmp]
                # add a faux job just so it goes round one whole time just in case all the jobs are finished (which they shouldn't but should the unlikely happen I THINK it would be better to get inside the while loop)
		jobs_still_running_list.append(1)
		# need to make sure that we record the data that has been converted so that we don't do it twice
		job_done_dict = {no + 1: False for no in range(no_of_arrays)}
		# keep looping until all the jobs are finished
		missing_jobs = []
		while len(jobs_still_running_list) > 0:
			print("Job array numbers still running: ", jobs_still_running_list)
			# get details about the queue
			queue_output_dict = submission.cluster_connection.checkSuccess(submission.cluster_connection.checkQueue, submission.cluster_job_number)
			stdout_tmp = queue_output_dict['stdout'].split("\n")
			del stdout_tmp[-1]
			jobs_still_running_list = [int(digit) for digit in stdout_tmp]
			# find all array jobs that have completed i.e. make a list of zeros the length of the amount of array jobs then change all zeros to ones that have a job running still. Finding all the zeros in thelist will now give you the array job numbers i.e. if total_list_of_jobs[idx] == 0 then job array idx + 1 has finished.
			total_list_of_jobs = [0]*no_of_arrays
			for idx in jobs_still_running_list:
				total_list_of_jobs[idx-1] = 1

			zero_idxs = [i for i, e in enumerate(total_list_of_jobs) if e == 0]
			finished_array_job_nos = [idx + 1 for idx in zero_idxs]
			print("finished_array_job_nos = ", finished_array_job_nos)
			# remove the jobs that have already had their data converted
			job_array_nos_to_convert = [job for job in finished_array_job_nos if job_done_dict[job] == False]
			print("job_array_nos_to_convert = ", job_array_nos_to_convert)
			# get the dirs corresponding to the arrays
			all_ordered_unique_gene_dirs = submission.order_of_keys_written_to_file
			array_no_to_dirs_to_convert_dict = {}
			for array_job in range(no_of_arrays):
				array_no_to_dirs_to_convert_dict[array_job + 1] = [all_ordered_unique_gene_dirs[no_of_unique_kos_per_array_job*(array_job) + unique_ko_idx] for unique_ko_idx in range(no_of_unique_kos_per_array_job)]
				
			
			# convert data of finished jobs
			if len(job_array_nos_to_convert) > 0:
				# create tuple of dirs (THIS ASSUMES THAT THE MAT DESTINATION IS THE SAME AS THE PANDAS DESTINATION)
				initial_tuple_of_dirs = [submission.simulation_output_path + '/' + dir + '/' + str(sub_idx) for job_id in job_array_nos_to_convert for dir in array_no_to_dirs_to_convert_dict[job_id] for sub_idx in submission.resource_usage_dict['list_of_rep_dir_names']]

				# because job arrays go missing and and simulation crash etc we remove all the dirs that don't have a summary.mat
				print("initial_tuple_of_dirs = ", initial_tuple_of_dirs)
				print("missing_jobs = ", missing_jobs)
				initial_tuple_of_dirs = initial_tuple_of_dirs + missing_jobs.copy()
				print("initial_tuple_of_dirs (plus missing_jobs) = ", initial_tuple_of_dirs)
				missing_jobs = []
				tuple_of_dirs = []
				for dir in initial_tuple_of_dirs:
					current_path = dir + "/summary.mat"
					cmds_to_check_file_exists = ['if [ -f ' + current_path + ' ]; then echo "yes";else echo "no";fi']
					does_file_exist = self.submission.cluster_connection.checkSuccess(self.submission.cluster_connection.sendCommand, cmds_to_check_file_exists)
					if does_file_exist['stdout'].strip()  == 'yes':
						tuple_of_dirs = tuple_of_dirs + [dir]
					else:
						missing_jobs = missing_jobs + [dir]

				print("tuple_of_dirs = ", tuple_of_dirs)
				print("New missing_jobs = ", missing_jobs)
#				tuple_of_dirs = tuple(tuple_of_dirs)
				if len(tuple_of_dirs) > 0:
					with Pool() as pool:
						print("converting to pandas!!")
						list_of_tuples = list(pool.map(self.convertDataToPandas, tuple_of_dirs))

					[list_of_data_conversion_output_dicts, list_of_tmp_sim_data_dicts] = zip(*list_of_tuples)
					# update self.simulation_data_dict
					for idx in range(len(list_of_tmp_sim_data_dicts)):
						for ko in list_of_tmp_sim_data_dicts[idx].keys(): # this shouldn't be neccessary as there should only be one KO set per dict but we don't know what the code is so resort to looping through the one key
							if ko not in self.simulation_data_dict:
								print("adding ko to simulation dict!")
								self.simulation_data_dict[ko] = []

							self.simulation_data_dict[ko] = self.simulation_data_dict[ko] + list_of_tmp_sim_data_dicts[idx][ko]
							print("self.simulation_data_dict[", ko, "] = ", self.simulation_data_dict[ko])

	#				data_coversion_ouput_dict = self.convertDataToPandas(tuple_of_dirs)
	#				all_outputs_dict['data_coversion_ouput_dict'] = data_coversion_ouput_dict
					# add coverted dirs to job_done_dict and delete matfiles
	#				for idx in range(len(data_conversion_output_dict)):
	#					if data_coversion_output_dict[idx]['return_code'] != 0:
	#						raise ValueError('There was an error converting the simulation output to Pandas dataframes! output_dict = ', data_coversion_ouput_dict)
	#				else:
				# add to job_done_dict
				for arr_no in job_array_nos_to_convert:
					job_done_dict[arr_no] = True

			time.sleep(900) # wait 900 seconds or 15 minutes
#			time.sleep(15) # wait 900 seconds or 15 minutes

		# update the KO database
		# convert self.simulation_data_dict from gene codes to gene ids
		simulation_dict = {}
		print("self.simulation_data_dict = ", self.simulation_data_dict)
		for ko in self.simulation_data_dict.keys():
			tmp_tuple = tuple([self.misc_data[code] for code in ko])
			simulation_dict[tmp_tuple] = self.simulation_data_dict[ko]

		cmds_to_update_db = self.submission.cluster_connection.db_connection.activate_venv_list
		cmds_to_update_db = cmds_to_update_db + ['python -c "import sys;sys.path.insert(0, \'' + self.submission.cluster_connection.db_connection.path_to_flex1 + '/database/ko_db/library\');import ko_db;path_to_ko_data = \'' + self.submission.cluster_connection.db_connection.path_to_flex1 + '/database/ko_db/new_test/ko.db\';path_to_static_db = \'' + self.submission.cluster_connection.db_connection.path_to_flex1 + '/database/staticDB/static.db\';ko_db_inst = ko_db.KoDb(path_to_ko_data, path_to_static_db);data_dict = ' + str(self.data_dict) + ';simulation_dict = ' + str(simulation_dict) + ';print(ko_db_inst.addNewSimulations(data_dict, simulation_dict))"']
		update_db_out = self.submission.cluster_connection.db_connection.checkSuccess(self.submission.cluster_connection.sendCommand, cmds_to_update_db)
		all_outputs_dict['update_db_out'] = update_db_out
		if update_db_out['return_code'] != 0:
			raise ValueError('There was an error converting the simulation output to Pandas dataframes! update_db_out = ', update_db_out)

		return

	def convertDataToPandas(self, directory_s):
#		if type(directory_s) is tuple:
#			for dir in directory_s:
#				# create simdir and save_dir
#				list_or_dict_of_simdir_and_save_dir = {}
#				list_or_dict_of_simdir_and_save_dir['simdir'] = dir
#				list_or_dict_of_simdir_and_save_dir['save_dir'] = dir
#				# create list of commands
#				cmd_list = self.submission.cluster_connection.activate_venv_list
#				cmd_list = cmd_list + ['python -c "import sys;sys.path.insert(0, \'' + self.submission.cluster_connection.path_to_flex1 + '/database/functions_for_extracting_data_from_matFiles/functions_for_joshuas_matFiles\');import extract_matFile_data_v73 as extract;list_or_dict_of_simdir_and_save_dir = ' + str(list_or_dict_of_simdir_and_save_dir) + ';extract.saveBasicDetailsAsPickle(list_or_dict_of_simdir_and_save_dir);extract.saveMatureRnaCountsAsPickle(list_or_dict_of_simdir_and_save_dir);extract.saveMetabolicReactionFluxsAsPickle(list_or_dict_of_simdir_and_save_dir);extract.saveMatureProteinComCountsAsPickle(list_or_dict_of_simdir_and_save_dir);extract.saveMatureProteinMonCountsAsPickle(list_or_dict_of_simdir_and_save_dir)"']
#				# submit commands
#				convert_data_output_dict = self.submission.cluster_connection.checkSuccess(self.submission.cluster_connection.sendCommand, cmd_list)
#
#				# remove .mat files
#				# find mat files
#				rm_mats_cmd = self.submission.cluster_connection.activate_venv_list
#				rm_mats_cmd = rm_mats_cmd + ['cd ' + str(list_or_dict_of_simdir_and_save_dir['simdir']), 'python -c "import glob;import os;list_of_matfiles = glob.glob(\'*.mat\');print(list_of_matfiles);print([os.remove(mfile) for mfile in list_of_matfiles])"']
#				self.submission.cluster_connection.checkSuccess(self.submission.cluster_connection.sendCommand, rm_mats_cmd)
#
#				# update KO database
#				# get simulation info
#				self.prepareSimulationDictForKoDbSubmission(dir)
#
#		else:
		dir = directory_s

		# create simdir and save_dir
		list_or_dict_of_simdir_and_save_dir = {}
		list_or_dict_of_simdir_and_save_dir['simdir'] = dir
		list_or_dict_of_simdir_and_save_dir['save_dir'] = dir
		# create list of commands
		cmd_list = self.submission.cluster_connection.activate_venv_list
		cmd_list = cmd_list + ['python -c "import sys;sys.path.insert(0, \'' + self.submission.cluster_connection.path_to_database_dir + '/database/functions_for_extracting_data_from_matFiles/functions_for_joshuas_matFiles\');import extract_matFile_data_v73 as extract;list_or_dict_of_simdir_and_save_dir = ' + str(list_or_dict_of_simdir_and_save_dir) + ';extract.saveBasicDetailsAsPickle(list_or_dict_of_simdir_and_save_dir);extract.saveMatureRnaCountsAsPickle(list_or_dict_of_simdir_and_save_dir);extract.saveMetabolicReactionFluxsAsPickle(list_or_dict_of_simdir_and_save_dir);extract.saveMatureProteinComCountsAsPickle(list_or_dict_of_simdir_and_save_dir);extract.saveMatureProteinMonCountsAsPickle(list_or_dict_of_simdir_and_save_dir)"']
		# submit commands
		convert_data_output_dict = self.submission.cluster_connection.checkSuccess(self.submission.cluster_connection.sendCommand, cmd_list)

		# remove .mat files
		# find mat files
		rm_mats_cmd = self.submission.cluster_connection.activate_venv_list
		rm_mats_cmd = rm_mats_cmd + ['cd ' + str(list_or_dict_of_simdir_and_save_dir['simdir']), 'python -c "import glob;import os;list_of_matfiles = glob.glob(\'*.mat\');print(list_of_matfiles);print([os.remove(mfile) for mfile in list_of_matfiles])"']
		self.submission.cluster_connection.checkSuccess(self.submission.cluster_connection.sendCommand, rm_mats_cmd)

		# update KO database
		# get simulation info
		tmp_sim_data_dict = self.prepareSimulationDictForKoDbSubmission(dir)

		return (convert_data_output_dict, tmp_sim_data_dict,)

class SimpleManageSubmission():
	def __init__(self, input_tuple):
		(submission_instance, misc_data) = input_tuple
		self.submission = submission_instance
		self.misc_data = misc_data # thiss is added just in case there needs to be some kind of unforeseen data transfer between library classes (fro example to save connecting to static.db to often the MGA class just passes the gene_code_to_id_dict to ManageSubmission through misc_data)
		self.simulation_data_dict = {}
		if misc_data == 'TEST':
			print("WARNING: This is in TEST mode so no files will be transfered and no job will be submitted.")
		else:
			outs =  self.submission.prepareForSubmission()
			list_of_return_codes = [output['return_code'] for output in outs]
			if sum(list_of_return_codes) > 0:
				raise ValueError('There has been a problem preparing some jobs for submission. The return codes from JobSubmission.prepareForSubmission() are: ', list_of_return_codes, '. The submission file name is ', self.submission.submission_file_name, ' and the temp_storage_path was ', self.submission.temp_storage_path)

			submission_out_dict = self.submission.submitJobToCluster()
			if submission_out_dict['return_code'] > 0:
				raise ValueError('There has been a problem submitting a job. The return code from JobSubmission.submitJobToCluster() is: ', submission_out_dict['return_code'], ". The submission file name is ", self.submission.submission_file_name, " and the temp_storage_path was ", self.submission.temp_storage_path)

		# there is no output given here because if any return codes come back none xero it raises an error and no other output is needed. If any output is needed for debugging it is suggested that you do that inside the self.monitor submission function if possible.
		if misc_data == 'TEST':
			pass
		else:
			self.monitorSubmission(self.submission)

		return 

	def monitorSubmission(self, submission):
		all_outputs_dict = {}
		time.sleep(3600) # wait 3,600 seconds or 1 hour because definitely nothing will finish in the first hour
#		time.sleep(15) # wait 14,400 seconds or 4 hours because definitely nothing will finish in the first four hours
		#useful info about the job
		dir_to_ko_set_dict = submission.ko_details_dict
		job_resource_allocation_dict = submission.resource_usage_dict
		no_of_arrays = job_resource_allocation_dict['no_of_arrays']
		no_of_unique_kos_per_array_job = job_resource_allocation_dict['no_of_unique_kos_per_array_job']
		no_of_repetitions_of_each_ko = job_resource_allocation_dict['no_of_repetitions_of_each_ko']
		# get details about the queue
		queue_output_dict = submission.cluster_connection.checkSuccess(submission.cluster_connection.checkQueue, submission.cluster_job_number)
		stdout_tmp = queue_output_dict['stdout'].split("\n")
		del stdout_tmp[-1]
		jobs_still_running_list = [int(digit) for digit in stdout_tmp]
                # add a faux job just so it goes round one whole time just in case all the jobs are finished (which they shouldn't but should the unlikely happen I THINK it would be better to get inside the while loop)
		jobs_still_running_list.append(1)
		# need to make sure that we record the data that has been converted so that we don't do it twice
		job_done_dict = {no + 1: False for no in range(no_of_arrays)}
		# keep looping until all the jobs are finished
		missing_jobs = []
		while len(jobs_still_running_list) > 0:
			print("Job array numbers still running: ", jobs_still_running_list)
			# get details about the queue
			queue_output_dict = submission.cluster_connection.checkSuccess(submission.cluster_connection.checkQueue, submission.cluster_job_number)
			stdout_tmp = queue_output_dict['stdout'].split("\n")
			del stdout_tmp[-1]
			jobs_still_running_list = [int(digit) for digit in stdout_tmp]

			time.sleep(900) # wait 900 seconds or 15 minutes
#			time.sleep(15) # wait 900 seconds or 15 minutes

		return
