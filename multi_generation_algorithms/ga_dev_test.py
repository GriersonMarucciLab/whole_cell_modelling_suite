from genetic_algorithms import GeneticAlgorithm
from connections import Bc3

cluster_user_name = 'oc13378'
ssh_config_alias = 'bc3'
path_to_key = '/home/oli/.ssh/uob/uob-rsa'
forename_of_user = 'Oliver'
surname_of_user = 'Chalkley'
user_email = 'o.chalkley@bristol.ac.uk'

bc3 = Bc3(cluster_user_name, ssh_config_alias, path_to_key, forename_of_user, surname_of_user, user_email)
dict_of_cluster_instances = {'bc3': bc3}
MGA_name = 'test_ga_18_03_01'
max_no_of_fit_individuals = 50
reps_of_unique_ko = 1
ga = GeneticAlgorithm(dict_of_cluster_instances, MGA_name, max_no_of_fit_individuals, reps_of_unique_ko, generation_num_to_gen_size_dict = {0: 1, -1: 1}, relative2clusterBasePath_simulation_output_path = 'kos/genetic_algorithms')
ga.run()
