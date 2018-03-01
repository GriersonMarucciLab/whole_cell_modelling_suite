from base_class import WholeCellModelBase
import random

class GeneticAlgorithm(WholeCellModelBase):
    def __init__(self, dict_of_cluster_instances, MGA_name, max_no_of_fit_individuals, reps_of_unique_ko, generation_num_to_gen_size_dict = {0: 150, -1: 50}, relative2clusterBasePath_simulation_output_path = 'kos/genetic_algorithms'):
        WholeCellModelBase.__init__(self, dict_of_cluster_instances, MGA_name)
        self.gene_code_to_id_dict = self.getDictOfJr358Codes(dict_of_cluster_instances[list(dict_of_cluster_instances.keys())[0]])
        self.gene_id_to_code_dict = self.invertDictionary(self.gene_code_to_id_dict)
        self.max_no_of_fit_individuals = max_no_of_fit_individuals
        self.fittest_individuals = []
        self.contemporary_generation = None
        self.relative2clusterBasePath_simulation_output_path = relative2clusterBasePath_simulation_output_path + '/' + self.MGA_name
        self.reps_of_unique_ko = reps_of_unique_ko
        ids_list = list(self.gene_code_to_id_dict.values())
        ids_list.sort()
        self.wt_genome = [1 for i in range(len(ids_list))]
        self.genome_idx_to_id_dict = {idx: ids_list[idx] for idx in range(len(ids_list))}
        self.id_to_genome_idx_dict = self.invertDictionary(self.genome_idx_to_id_dict)
        self.generation_num_to_gen_size_dict = generation_num_to_gen_size_dict

    def checkStop(self):
        if self.generation_counter > 2:
            return True
        else:
            return False

    def getGenerationName(self):
        generation_name = 'gen' + str(self.generation_counter)

        return generation_name

    def getPopulationSize(self):
        important_idxs = list(self.generation_num_to_gen_size_dict.keys())
        if important_idxs.count(self.generation_counter) > 0:
            generation_size = self.generation_num_to_gen_size_dict[self.generation_counter]
        else:
            generation_size = self.generation_num_to_gen_size_dict[-1]

        return generation_size

    def getRandomKos(self):
        submission_name = self.getGenerationName()
        population_size = self.getPopulationSize()
        jr_codes_to_ids_dict = self.gene_code_to_id_dict.copy()
        idx_to_ids_dict = self.createIdxToIdDict(jr_codes_to_ids_dict)
        list_of_jr_ids = list(jr_codes_to_ids_dict.values())
#        print("jr_codes_to_ids_dict = ", jr_codes_to_ids_dict)
        list_of_jr_ids.sort()
        all_possible_idxs_iter = range(len(list_of_jr_ids))
        # create the inverse dictionary
        jr_ids_to_codes_dict = self.invertDictionary(jr_codes_to_ids_dict)
        ko_name_to_set_dict = {}
        ko_names = [float('NaN') for i in range(population_size)]
        for ko_no in range(population_size):
            length_of_combo = random.randint(2, 5)
            ko_names[ko_no] = 'ko' + str(ko_no + 1)
            ko_idxs = self.random_combination(all_possible_idxs_iter, length_of_combo)
            # convert indexs into gene IDs
            ko_ids = [idx_to_ids_dict[idx] for idx in ko_idxs]
            ko_ids.sort() # this should be taken care of but to be safe always try and keep gene IDs in ascending order to avoid unforseen problems
         
            # convert gene IDs into codes
            ko_codes = [jr_ids_to_codes_dict[id] for id in ko_ids]
            ko_name_to_set_dict[ko_names[ko_no]] = tuple(ko_codes)

        return ko_name_to_set_dict

    def mateTheFittest(self):
        # get the fittest (note this is already in order of largest KO length at the top and smallest at the bottom)
        fittest = self.fittest_individuals.copy()
        print("fittest = ", fittest)
        # randomly pick a ko such that larger KOs are more likely to be picked
        ko_length_of_fittest = [len(ko) for ko in fittest]
        list_of_probabilities = [ko_length_of_fittest[idx]/sum(ko_length_of_fittest) for idx in range(len(ko_length_of_fittest))]
        # create new generation
        pop_size = self.getPopulationSize()
        list_of_children = [float('NaN') for i in range(pop_size)]
        ko_set_names = [float('NaN') for i in range(pop_size)]
        for child_idx in range(pop_size):
            parent1_codes = np.random.choice(fittest, p=list_of_probabilities)
            parent2_codes = parent1_codes
            while parent2_codes == parent1_codes:
                parent2_codes = np.random.choice(fittest, p=list_of_probabilities)

            # convert parent ko codes to ids
            parent1_ids = [self.gene_code_to_id_dict[code] for code in parent1_codes]
            parent2_ids = [self.gene_code_to_id_dict[code] for code in parent2_codes]
            # convert parent KO sets to genomes
            parent1_genome = self.wt_genome.copy()
            parent2_genome = self.wt_genome.copy()
            for id in parent1_ids:
                parent1_genome[self.id_to_genome_idx_dict[id]] = 0

            for id in parent2_ids:
                parent2_genome[self.id_to_genome_idx_dict[id]] = 0

            # mate the genomes
            # create empty child genome
            tmp_child = []
            # make sure the child has at least 2 KOs
            while tmp_child.count(0) < 2:
                # pick a idx to split the geneomes by
                split_idx = random.randint(0,len(parent1_genome) - 1)
#                # OLD WAY
#                tmp_child = parent1_genome[:split_idx] + parent2_genome[split_idx:]
                # NEW WAY
                # randomly create the gene indexs to take from parent1
                parent1_idxs_to_inherit = random.sample(range(len(parent1_genome)), split_idx)
                # create tmp child genome from randomly selected choice of genes from both parents
                tmp_child = [parent1_genome[idx] if parent1_idxs_to_inherit.count(idx) > 0 else parent2_genome[idx] for idx in range(len(parent1_genome))]

            # mutate genes randomly 10% of the time
            if random.random() < 0.1:
#                # OLD WAY
#                # randomly pick index from JRs 359 set
#                idx = random.randint(0, len(tmp_child) - 1)
#                # flip the gene
#                tmp_child[idx] = (tmp_child[idx] + 1)%2

                # NEW WAY
                # pick the amount of gene mutations from a exponentially distributed random number with parameter 2
                # exp R.V. can produce zero, we don't want zeros
                number_of_gene_mutations = 0
                while number_of_gene_mutations == 0:
                    number_of_gene_mutations = int(np.around(np.random.exponential(2)))

                # flip number_of_gene_mutations amount of genes randomly
                # create list of indexs to flip
                gene_idxs_to_flip = random.sample(range(len(tmp_child)), number_of_gene_mutations)
                # flip genes
                for idx in gene_idxs_to_flip:
                    tmp_child[idx] = (tmp_child[idx] + 1)%2

            # convert back into idxs, ids then codes
            tmp_child = [i for i,x in enumerate(tmp_child) if x == 0]
            tmp_child = [self.genome_idx_to_id_dict[idx] for idx in tmp_child]
            tmp_child.sort()
            tmp_child = tuple([self.gene_id_to_code_dict[gene_id] for gene_id in tmp_child])

            # update children
            list_of_children[child_idx] = tmp_child
            # create ko set names
            ko_set_names[child_idx] = 'ko' + str(child_idx + 1)
            
        ko_name_to_set_dict = {ko_set_names[idx]: list_of_children[idx] for idx in range(len(list_of_children))}

        return ko_name_to_set_dict

    def getNewGeneration(self):
        try:
            len(self.fittest_individuals)
            has_length = True
        except:
            has_length = False

        if has_length == True:
            if self.generation_counter == 0:
                print("generation 0!")
                ko_name_to_set_dict = self.getRandomKos()
            elif len(self.fittest_individuals) == 0:
                print("No survivors!")
                ko_name_to_set_dict = self.getRandomKos()
            elif len(self.fittest_individuals) == 1:
                print("Only one survivor!")
                ko_name_to_set_dict = self.getRandomKos()
                list_of_ko_numbers = [int(ko_name[2:]) for ko_name in ko_name_to_set_dict.keys()]
                new_ko_name = 'ko' + str(max(list_of_ko_numbers) + 1)
                ko_name_to_set_dict[new_ko_name] = self.fittest_individuals[0]
            else:
                print("Normal mating!")
                ko_name_to_set_dict = self.mateTheFittest()
        elif has_length == False:
            print("No length!")
            ko_name_to_set_dict = self.getRandomKos()
        else:
            raise ValueError("self.fittest_individuals must either have length or not have length here self.fittest_individuals = ", self.fittest_individuals)

        return ko_name_to_set_dict

    def updateFittestPopulation(self, submission_instance, submission_management_instance):
        # create dictionary which will be used to store all viable sims from this generation
        dict_of_this_gen_viable_sims = {}
        # extract all the simulations that divided from this generation and record them in a tmp dict
        tmp_sim_results_dict = {ko: submission_management_instance.simulation_data_dict[ko] for ko in submission_management_instance.simulation_data_dict.keys() if sum([int(submission_management_instance.simulation_data_dict[ko][rep][1]) for rep in range(len(submission_management_instance.simulation_data_dict[ko]))]) != 0}
        # add all results to dict_of_this_gen_viable_sims
        for ko in tmp_sim_results_dict.keys():
            if ko not in dict_of_this_gen_viable_sims:
                dict_of_this_gen_viable_sims[ko] = []

            dict_of_this_gen_viable_sims[ko].append(tmp_sim_results_dict[ko])

        #combine with the current fittest population
        all_viable_kos = list(dict_of_this_gen_viable_sims.keys())
        all_viable_kos = all_viable_kos + self.fittest_individuals.copy()
        # get unique ko sets
        all_viable_kos = list(set(all_viable_kos))
        # sort all_viable_ids in order of length so that we can pick the fittest 100
        # create dict of ko set to length
        ko_codes_to_len_dict = {ko: len(ko) for ko in all_viable_kos}
        ko_codes_to_len_list_sorted = sorted(ko_codes_to_len_dict.items(), key=operator.itemgetter(1), reverse=True)
        # create a list of the fittest individuals
        fittest_individuals = [ko_set[0] for ko_set in ko_codes_to_len_list_sorted]
        if len(fittest_individuals) > self.max_no_of_fit_individuals:
            fittest_individuals = fittest_individuals[:(self.max_no_of_fit_individuals)]

        self.fittest_individuals = fittest_individuals.copy()

        return
