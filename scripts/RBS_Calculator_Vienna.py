#Given an mRNA sequence, this Python class predicts the dG_total and translation initiation rate.

#This file is part of the Ribosome Binding Site Calculator.

#The Ribosome Binding Site Calculator is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#The Ribosome Binding Site Calculator is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with Ribosome Binding Site Calculator.  If not, see <http://www.gnu.org/licenses/>.
#Copyright 2008-2009 is owned by the University of California Regents. All rights reserved.

from ViennaRNA import *
import re
import math
import csv
from imthinkingRBS import module_searcher, calc_dG_distort, calc_dG_standby

class CalcError(Exception):
    """Base class for exceptions in this module."""

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

class RBS_Calculator:

    #From experimental characterization of ALL Predicted RBSs as of 1/20/08
    RT_eff = 2.222
    logK =  7.824
    K = 2500.0 #this have likely changed between 1.0 and 2.0 (or not!!)

    #Global parameters -- constants
    infinity = 1e12 #For all practical purposes, here.
    RNA_model = "rna1999"
    start_codon_energies = {"ATG":-1.194, "AUG": -1.194, "GTG": -0.0748, "GUG": -0.0748, "TTG":-0.0435, "UUG": -0.0435, "CTG": -0.03406, "CUG":-0.03406} #hybridization to CAT
    auto_dangles = True
    dangles_default = "all"
    temp = 37.0
    optimal_spacing = 5 #aligned spacing

    dG_spacing_constant_push = [12.2, 2.5, 2.0, 3.0]
    dG_spacing_constant_pull = [0.048, 0.24, 0.0]
    cutoff = 35 #number of nt +- start codon considering for folding
    standby_site_length = 4 #Number of nt before SD sequence that must be unpaired for ribosome binding
    energy_cutoff = 3.0
    start_codons = ["ATG", "AUG", "GTG", "GUG","TTG","UUG"] #substituted U for T in actual calcs. Ignores CTG/CUG
    #rRNA = "acctcctta" #These are the last 9 nt (3' end) of the 16S rRNA in E. coli

    footprint = 1000 #Footprint of the 30S complex that prevents formation of secondary structures downstream of the start codon. Here, we assume that the entire post-start RNA sequence does not form secondary structures once the 30S complex has bound.

    def __init__(self, mRNA, start_range, rRNA, constraints, verbose = False):
        """Initializes the RBS Calculator class with the mRNA sequence and the range of start codon positions considered."""

        #NuPACK.__init__(self,sequences,self.RNA_model)

        exp = re.compile('[ATGCU]',re.IGNORECASE)
        if exp.match(mRNA) == None:
            raise ValueError("Invalid letters found in sequence ""%s"". Only ATGCU accepted." % mRNA)

        if start_range[0] < 0: start_range[0] = 0
        if start_range[1] > len(mRNA): start_range[1] = len(mRNA)

        #self.name = name
        self.mRNA_input = mRNA.upper()
        self.constraints = constraints
        self.rRNA = rRNA
        self.rRNA_len = len(self.rRNA)
        self.mRNA_len = len(self.mRNA_input)
        self.total_sequence_length = len(mRNA) + len(self.rRNA)
        self.dG_rRNA = self.calc_dG_rRNA()
        self.run = 0
        self.start_range = start_range
        self.verbose = verbose

    def find_min(self,input_list):
        """Finds the minimum of a list of numbers."""

        min_item = self.infinity
        min_index = 0

        for i, item in enumerate(input_list):
            if item < min_item:
                min_item = item
                min_index = i

        return (min_item,min_index)


    def find_start_codons(self,sequence):
        """Finds all start codons in an mRNA sequence. Creates a list."""

        self.start_position_list = []
        self.start_codon_list = []

        seq_len = len(sequence)
        end = min(self.start_range[1],seq_len-2)
        begin = min(self.start_range[0],end)

        for i in range(begin,end+1):
            codon = sequence[i:i+3]
            if codon.upper() in self.start_codons:
                self.start_position_list.append(i)
                self.start_codon_list.append(codon)
                yield (i,codon)
            else:
                pass

    def calc_aligned_spacing(self,mRNA,start_pos,bp_x,bp_y):
        """Calculates the aligned spacing between the 16S rRNA binding site and the start codon."""

        #rRNA is the concatenated at the end of the sequence in 5' to 3' direction
        #first: identify the farthest 3' nt in the rRNA that binds to the mRNA and return its mRNA base pairer

        Ok = False
        seq_len = len(mRNA) + self.rRNA_len
        for (rRNA_nt) in range(seq_len,seq_len - self.rRNA_len,-1):

            if rRNA_nt in bp_y:
                rRNA_pos = bp_y.index(rRNA_nt)
                if bp_x[rRNA_pos] < start_pos:
                    Ok = True
                    farthest_3_prime_rRNA = rRNA_nt - len(mRNA)

                    mRNA_nt = bp_x[rRNA_pos]
                    distance_to_start = start_pos - mRNA_nt + 1 #start_pos is counting starting from 0 (python)

                    break
                else:
                    break

        if Ok:
            aligned_spacing = distance_to_start - farthest_3_prime_rRNA
        else:
            aligned_spacing = self.infinity

        return aligned_spacing

    def calc_dG_spacing(self, aligned_spacing):
        """Calculates the dG_spacing according to the value of the aligned spacing. This relationship was determined through experiments."""

        if (aligned_spacing < self.optimal_spacing):
            ds = aligned_spacing - self.optimal_spacing

            dG_spacing_penalty = self.dG_spacing_constant_push[0] / (1.0 + math.exp(self.dG_spacing_constant_push[1]*(ds + self.dG_spacing_constant_push[2] )))**self.dG_spacing_constant_push[3]

        else:
            ds = aligned_spacing - self.optimal_spacing
            dG_spacing_penalty = self.dG_spacing_constant_pull[0] * ds * ds + self.dG_spacing_constant_pull[1] * ds + self.dG_spacing_constant_pull[2]

        return dG_spacing_penalty

    def calc_dG_mRNA_rRNA(self,start_pos):
        """Calculates the dG_mRNA_rRNA from the mRNA and rRNA sequence. Considers all feasible 16S rRNA binding sites and includes the effects of non-optimal spacing."""

        begin = 0 #max(0,start_pos-self.cutoff)
        mRNA_len = min(len(self.mRNA_input),start_pos+self.cutoff)
        start_pos_in_subsequence = start_pos #min(start_pos, self.cutoff)
        startpos_to_end_len = mRNA_len - start_pos_in_subsequence - begin

        #print(begin, mRNA_len, start_pos_in_subsequence, startpos_to_end_len)

        #1. identify a list of rRNA-binding sites. Binding sites are hybridizations between the mRNA and rRNA and can include mismatches, bulges, etc. Intra-molecular folding is also allowed within the mRNA. The subopt program is used to generate a list of optimal & suboptimal binding sites.
        #Constraints: the entire rRNA-binding site must be upstream of the start codon
        mRNA = self.mRNA_input[begin:start_pos]

        if self.constraints is None:
            constraints = None
        else:
            constraints = self.constraints[begin:start_pos] #added so constraints file will be the same as the sequence

        if len(mRNA) == 0:
            raise CalcError("Warning: There is a leaderless start codon, which is being ignored.")

        #print "After exception"

        fold = ViennaRNA([mRNA,self.rRNA],material = self.RNA_model)
        #print(fold)
        fold.subopt([1, 2], constraints, self.energy_cutoff,dangles = 'all', Temp = self.temp)

        if len(fold["subopt_basepairing_x"]) == 0:
            raise CalcError("Warning: The 16S rRNA has no predicted binding site. Start codon is considered as leaderless and ignored.")

        #2. Calculate dG_spacing for each 16S rRNA binding site

        #Calculate the aligned spacing for each binding site in the list
        aligned_spacing = []
        for (bp_x, bp_y) in zip(fold["subopt_basepairing_x"], fold["subopt_basepairing_y"]):
            aligned_spacing.append(self.calc_aligned_spacing(mRNA, start_pos_in_subsequence, bp_x,bp_y))

        dG_spacing_list = []
        dG_mRNA_rRNA = []
        dG_mRNA_rRNA_withspacing = []

        #Calculate dG_spacing using aligned spacing value. Add it to dG_mRNA_rRNA.
        for (counter) in range(len(fold["subopt_basepairing_x"])):

            dG_mRNA_rRNA.append(fold["subopt_energy"][counter])
            val = self.calc_dG_spacing(aligned_spacing[counter])
            dG_spacing_list.append(val)
            dG_mRNA_rRNA_withspacing.append(val + fold["subopt_energy"][counter])

        #3. Find 16S rRNA binding site that minimizes dG_spacing+dG_mRNA_rRNA.
        [dG_mRNA_rRNA_folding, index] = self.find_min(dG_mRNA_rRNA_withspacing)
        dG_spacing_final = dG_spacing_list[index]

        dG_mRNA_rRNA_nospacing = dG_mRNA_rRNA[index]

        #Check: Is the dG spacing large compared to the energy gap? If so, this means the list of suboptimal 16S rRNA binding sites generated by subopt is too short.
        if dG_spacing_final > self.energy_cutoff:
            if self.verbose: print("Warning: The spacing penalty is greater than the energy gap. dG (spacing) = ", dG_spacing_final)


        #4. Identify the 5' and 3' ends of the identified 16S rRNA binding site. Create a base pair list.

        most_5p_mRNA = self.infinity
        most_3p_mRNA = -self.infinity

        #Generate a list of rRNA-mRNA base paired nucleotides
        bp_x_target = []
        bp_y_target = []

        bp_x = fold["subopt_basepairing_x"][index]
        bp_y = fold["subopt_basepairing_y"][index]
        for (nt_x, nt_y) in zip(bp_x, bp_y):
            if nt_y > len(mRNA): #nt is rRNA
                most_5p_mRNA = min(most_5p_mRNA, bp_x[bp_y.index(nt_y)])
                most_3p_mRNA = max(most_3p_mRNA, bp_x[bp_y.index(nt_y)])
                bp_x_target.append(nt_x)
                bp_y_target.append(nt_y)

        #The rRNA-binding site is between the nucleotides at positions most_5p_mRNA and most_3p_mRNA
        #Now, fold the pre-sequence, rRNA-binding-sequence and post-sequence separately. Take their base pairings and combine them together. Calculate the total energy. For secondary structures, this splitting operation is allowed.

        #We postulate that not all of the post-sequence can form secondary structures. Once the 30S complex binds to the mRNA, it prevents the formation of secondary structures that are mutually exclusive with ribosome binding. We define self.footprint to be the length of the 30S complex footprint. Here, we assume that the entire mRNA sequence downstream of the 16S rRNA binding site can not form secondary structures.

        mRNA_pre = self.mRNA_input[begin:begin+most_5p_mRNA-1]
        post_window_end = mRNA_len + 1
        post_window_begin = min(start_pos + self.footprint,post_window_end) #Footprint
        post_window_end = mRNA_len + 1
        mRNA_post = self.mRNA_input[post_window_begin:post_window_end]

        mRNA_pre_len = len(mRNA_pre)
        mRNA_post_len = len(mRNA_post)
        mRNA_rRNA_binding_len = most_3p_mRNA - most_5p_mRNA + 1
        total_folded_len = mRNA_pre_len + mRNA_post_len + mRNA_rRNA_binding_len

        total_bp_x = []
        total_bp_y = []

        if self.constraints is None:
            pre_constraints = None
            post_constraints = None
        else:
            pre_constraints = self.constraints[begin:begin+most_5p_mRNA-1]
            post_constraints = self.constraints[post_window_begin:post_window_end]

        #Calculate pre-sequence folding
        if len(mRNA_pre) > 0:
            fold_pre = ViennaRNA([mRNA_pre], material = self.RNA_model)
            fold_pre.mfe([1], pre_constraints, Temp = self.temp, dangles = self.dangles, )
            bp_x_pre = fold_pre["mfe_basepairing_x"][0]
            bp_y_pre = fold_pre["mfe_basepairing_y"][0]

        else:
            bp_x_pre = []
            bp_y_pre = []

        #Add pre-sequence base pairings to total base pairings
        offset = 0 #Begins at 0
        for (nt_x, nt_y) in zip(bp_x_pre, bp_y_pre):
            total_bp_x.append(nt_x + offset)
            total_bp_y.append(nt_y + offset)

        #Add rRNA-binding site base pairings to total base pairings
        offset = 0 #Begins at zero
        if startpos_to_end_len < self.cutoff:
            rRNA_offset = startpos_to_end_len
        else:
            rRNA_offset = startpos_to_end_len

        for (nt_x, nt_y) in zip(bp_x_target, bp_y_target):
            total_bp_x.append(nt_x + offset)
            total_bp_y.append(nt_y + rRNA_offset)

        #Calculate post-sequence folding
        if len(mRNA_post) > 0:
            fold_post = ViennaRNA([mRNA_post], material = self.RNA_model)
            fold_post.mfe([1], post_constraints, Temp = self.temp, dangles = self.dangles)
            bp_x_post = fold_post["mfe_basepairing_x"][0]
            bp_y_post = fold_post["mfe_basepairing_y"][0]
        else:
            bp_x_post = []
            bp_y_post = []

        offset = post_window_begin - begin
        for (nt_x, nt_y) in zip(bp_x_post, bp_y_post):
            total_bp_x.append(nt_x + offset)
            total_bp_y.append(nt_y + offset)

        mRNA = self.mRNA_input[begin:mRNA_len]
        fold = ViennaRNA([mRNA, self.rRNA], material = self.RNA_model)

        total_energy = fold.energy([1, 2], total_bp_x, total_bp_y, Temp = self.temp, dangles = self.dangles)

        energy_nowindows = dG_mRNA_rRNA_nospacing
        total_energy_withspacing = total_energy + dG_spacing_final

        structure = fold
        structure["program"] = "subopt"
        structure["mRNA"] = mRNA
        structure["MinStructureID"] = 0
        structure["dG_mRNA_rRNA"] = total_energy
        structure["dG_mRNA_rRNA_withspacing"] = total_energy_withspacing
        structure["dG_spacing"] = dG_spacing_final
        structure["subopt_energy"] = [total_energy_withspacing]
        structure["subopt_basepairing_x"] = [total_bp_x]
        structure["subopt_basepairing_y"] = [total_bp_y]
        structure["subopt_composition"] = [1, 2]
        structure["bp_x"] = total_bp_x
        structure["bp_y"] = total_bp_y
        structure['bracket_string'] = fold["energy_bracket_string"]

        return (total_energy_withspacing, structure)

    def calc_dG_standby_site(self,structure_old, rRNA_binding = True):
        """Calculates the dG_standby given the structure of the mRNA:rRNA complex"""

        #To calculate the mfe structure while disallowing base pairing at the standby site, we split the folded mRNA sequence into three parts: (i) a pre-sequence (before the standby site) that can fold; (ii) the standby site, which can not fold; (iii) the 16S rRNA binding site and downstream sequence, which has been previously folded.

        import copy
        structure = copy.deepcopy(structure_old)
        mRNA = structure["mRNA"]
        bp_x = structure["bp_x"]
        bp_y = structure["bp_y"]
        bracket_string = structure['bracket_string']
        energy_before = structure["dG_mRNA_rRNA"] #without spacing effects

        #Identify the most 5p mRNA nt that is bound to rRNA
        for (nt_x, nt_y) in zip(bp_x, bp_y):
            if nt_x <= len(mRNA) and nt_y > len(mRNA): #nt_x is mRNA, nt_y is rRNA, they are bound.
                most_5p_mRNA = nt_x #starts counting from 0
                break

        #Extract the base pairings that are 3' of the most_5p_mRNA base pairing
        bp_x_3p = []
        bp_y_3p = []
        for (nt_x, nt_y) in zip(bp_x, bp_y):
            if nt_x >= most_5p_mRNA:
                bp_x_3p.append(nt_x)
                bp_y_3p.append(nt_y)

        #Create the mRNA subsequence
        mRNA_subsequence = mRNA[0:most_5p_mRNA]   #[0:max(0,most_5p_mRNA - self.standby_site_length - 1)]
        if self.constraints is None:
            constraint_subsequence = None
        else:
            constraint_subsequence = self.constraints[0:most_5p_mRNA]   #0:max(0,most_5p_mRNA - self.standby_site_length - 1)]
        #standby_site = mRNA[most_5p_mRNA - self.standby_site_length - 1:most_5p_mRNA]

        #Fold it and extract the base pairings
        if (len(mRNA_subsequence)) > 0:
            fold = ViennaRNA([mRNA_subsequence], material = self.RNA_model)
            fold.mfe([1], constraint_subsequence, Temp = self.temp, dangles = self.dangles)
            energy_after_5p = fold["mfe_energy"][0]
            bp_x_5p = fold["mfe_basepairing_x"][0]   #[0] added 12/13/07
            bp_y_5p = fold["mfe_basepairing_y"][0]
            bracket_string = fold['bracket_string']
        else:
             bp_x_5p = []
             bp_y_5p = []
             energy_after_5p = 0.0

        #new commands for dG standby from v2.0 (AL)
        if len(bp_x_5p) > 2:
            module_range, modules = module_searcher(bracket_string, bp_x_5p, bp_y_5p)
            print(modules)
            dG_standby = calc_dG_standby(modules,mRNA_subsequence)
        else:
            dG_standby = 0

        #Put the sets of base pairings together
        #bp_x_after = []
        #bp_y_after = []
        #for (nt_x, nt_y) in zip(bp_x_5p, bp_y_5p):
        #    bp_x_after.append(nt_x)
        #    bp_y_after.append(nt_y)

        #for (nt_x, nt_y) in zip(bp_x_3p, bp_y_3p):
        #    bp_x_after.append(nt_x)
        #    bp_y_after.append(nt_y)

        #Calculate its energy
        #fold = ViennaRNA([mRNA, self.rRNA], material = self.RNA_model)
        #energy_after = fold.energy([1, 2], bp_x_after, bp_y_after, dangles = self.dangles, Temp = self.temp)

        #dG_standby_site = energy_before - energy_after

        #if (dG_standby_site > 0.0): dG_standby_site = 0.0

        #index = structure["MinStructureID"]
        #structure["bp_x"] = bp_x_after
        #structure["bp_y"] = bp_y_after
        #structure["subopt_basepairing_x"][index] = bp_x_after
        #structure["subopt_basepairing_y"][index] = bp_y_after
        #structure["subopt_energy"][index] = energy_after
        #structure["dG_mRNA_rRNA_corrected"] = energy_after

        return (dG_standby)

    def calc_dG_mRNA(self,start_pos):
        """Calculates the dG_mRNA given the mRNA sequence."""

        mRNA = self.mRNA_input[0:start_pos + self.cutoff] #max(0,start_pos-self.cutoff):min(len(self.mRNA_input),start_pos+self.cutoff)]
        if self.constraints is None:
            constraints = None
        else:
            constraints = self.constraints[0:start_pos+ self.cutoff] #max(0,start_pos-self.cutoff):min(len(self.mRNA_input),start_pos+self.cutoff)]
        
        fold = ViennaRNA([mRNA],self.RNA_model)
        fold.mfe([1], constraints, Temp = self.temp, dangles = self.dangles)

        structure = fold
        structure["mRNA"] = mRNA
        structure["bp_x"] = fold["mfe_basepairing_x"][0]
        structure["bp_y"] = fold["mfe_basepairing_y"][0]
        structure["dG_mRNA"] = fold["mfe_energy"][0]
        structure["MinStructureID"] = 0

        dG_mRNA_folding = fold["mfe_energy"][0]
        #print('this is the dG mRNA ' + str(fold['mfe_energy'][0]))
        return (dG_mRNA_folding, structure)

    def calc_dG_rRNA(self):
        """Calculates the dG of folding for the last 9 nt of the 16S rRNA. Not used in the free energy model."""
        fold = ViennaRNA([self.rRNA],self.RNA_model)
        fold.mfe([1], constraints = None, Temp = self.temp, dangles = "all")
        dG_rRNA_folding = fold["mfe_energy"][0]
        return dG_rRNA_folding

    def calc_dG_SDopen(self, mRNA_structure, mRNA_rRNA_structure): #Alex note: I havent converted this fxn since it is not used in calc_dG
        """Calculate the dG required to unfold the nucleotides in the 16S rRNA binding site."""

        mRNA = mRNA_structure["mRNA"]
        program = mRNA_structure["program"]
        index = mRNA_structure["MinStructureID"]
        dG_mRNA = mRNA_structure[program + "_energy"][index]

        index = mRNA_rRNA_structure["MinStructureID"]
        bp_x_1 = mRNA_rRNA_structure["subopt_basepairing_x"][index][:]
        bp_y_1 = mRNA_rRNA_structure["subopt_basepairing_y"][index][:]

        most_5p_mRNA = self.infinity
        most_3p_mRNA = -self.infinity

        for (nt_x, nt_y) in zip(bp_x_1, bp_y_1):
            if nt_y > len(mRNA): #nt is rRNA
                most_5p_mRNA = min(most_5p_mRNA, bp_x_1[bp_y_1.index(nt_y)])
                most_3p_mRNA = max(most_3p_mRNA, bp_x_1[bp_y_1.index(nt_y)])

        pre_mRNA = mRNA[0:most_5p_mRNA]
        pre_mRNA_constraints = self.constraints[0:most_5p_mRNA]
        post_mRNA = mRNA[most_3p_mRNA+1:len(mRNA)+1]
        post_mRNA_constraints = self.constraints[most_3p_mRNA+1:len(mRNA)+1]

        pre_fold = ViennaRNA([pre_mRNA],material = self.RNA_model)
        pre_fold.mfe([1],dangles = self.dangles, Temp = self.temp)
        dG_pre = pre_fold["mfe_energy"][0]

        post_fold = NuPACK([post_mRNA],material = self.RNA_model)
        post_fold.mfe([1],dangles = self.dangles, Temp = self.temp)
        dG_post = post_fold["mfe_energy"][0]

        energy = dG_pre + dG_post

        ddG_mRNA = energy - dG_mRNA #positive if work is required to unfold SD sequence
        return ddG_mRNA

    def calc_kinetic_score(self, structure = None, mRNA_in = None, bp_x_in = None, bp_y_in = None):
        """Calculate a "kinetic score", a heuristic measure of the maximum time required for the mRNA secondary structure to form. This is related to the RNA polymer model by David et. al. This heuristic should not be used in any way to quantify the folding kinetics of an mRNA sequence because it completely ignores cooperative RNA folding mechanisms, such as zipping or strand displacement. Here, we use it to eliminate mRNA sequences that MAY fold slowly."""

        if not (structure is None):
            program = structure["program"]
            mRNA = structure["mRNA"]
            index = structure["MinStructureID"]
            bp_x = structure[program + "_basepairing_x"][index]
            bp_y = structure[program + "_basepairing_y"][index]

        if not (bp_x_in is None) and not (bp_y_in is None) and not (mRNA_in is None):
            mRNA = mRNA_in[:]
            bp_x = bp_x_in[:]
            bp_y = bp_y_in[:]

        largest_range_helix = 0
        for (nt_x, nt_y) in zip(bp_x, bp_y):
            if nt_x <= len(mRNA) and nt_y <= len(mRNA):
                val = nt_y - nt_x
                largest_range_helix = max(val, largest_range_helix)

        kinetic_score = float(largest_range_helix) / float(len(mRNA))
        if float(largest_range_helix) > 0:
            min_bp_prob = float(largest_range_helix)**(-1.44) #From David et. al.
        else:
            min_bp_prob = 1.0

        return (kinetic_score, min_bp_prob)

    def calc_most_5p_mRNA(self,structure_old):
        """Calculates the most 5' nucleotide in the 16S rRNA binding site."""

        import copy
        structure = copy.deepcopy(structure_old)
        mRNA = structure["mRNA"]
        bp_x = structure["bp_x"]
        bp_y = structure["bp_y"]

        #Identify the most 5p mRNA nt that is bound to rRNA
        for (nt_x, nt_y) in zip(bp_x, bp_y):
            if nt_x <= len(mRNA) and nt_y > len(mRNA): #nt_x is mRNA, nt_y is rRNA, they are bound.
                most_5p_mRNA = nt_x
                break

        return most_5p_mRNA

    def calc_longest_helix(self,structure):
        """Calculate the longest helical structure (longest contiguous list of base pairings) in the secondary structure"""

        mRNA = structure["mRNA"]
        bp_x = structure["bp_x"]
        bp_y = structure["bp_y"]

        longest_helix = 0
        helix_length = 1

        for (nt_x, nt_y) in zip(bp_x, bp_y):
            if (bp_x.count(nt_x+1) > 0 and bp_y.count(nt_y-1) > 0):
                helix_length += 1
            else:
                longest_helix = max(longest_helix, helix_length)
                helix_length = 1

        return longest_helix

    def calc_longest_loop_bulge(self,structure,output_start_end=False,InRBSOnly=False,RBS=None):
        """Calculate the longest helical loop and bulge structure (longest contiguous list of un-base paired nucleotides starting and ending with a helix (loop -> same helix, bulge -> different helix) in the secondary structure"""

        mRNA = structure["mRNA"]
        bp_x = structure["bp_x"]
        bp_y = structure["bp_y"]

        loop_length = 0
        begin_helix = 1

        bulge_loop_list = []
        helical_loop_list = []

        if output_start_end:
            bulge_loop_start_end = []
            helical_loop_start_end = []

        if InRBSOnly and RBS is not None:
            RBS_begin = mRNA.find(RBS)
            RBS_end = RBS_begin + len(RBS)
            nucleotide_range = list(range(RBS_begin,RBS_end+1))

        else:
            nucleotide_range = list(range(1,len(mRNA)+1))

        #Find loops. Find bulges.
        for n in nucleotide_range:
            if bp_x.count(n) == 0 and bp_y.count(n) == 0: #nth nucleotide is not base-paired.

                #Determine if nearest neighbor nucleotides are base-paired
                (x1,x2,y1,y2) = (bp_x.count(n-1), bp_x.count(n+1), bp_y.count(n-1), bp_y.count(n+1))

                #print "#", n, (x1,x2,y1,y2)

                #middle unpaired nt
                if (x1,x2,y1,y2) == (0,0,0,0):
                    loop_length += 1

                #single mismatch -- loop
                elif (x1,x2,y1,y2) == (1,0,0,1) or (x1,x2,y1,y2) == (0,1,1,0):
                    loop_length += 1
                    begin_helix = n - 1
                    end_helix = n + 1

                #single mismatch -- bulge
                elif (x1,x2,y1,y2) == (1,1,0,0) or (x1,x2,y1,y2) == (0,0,1,1):
                    loop_length += 1
                    begin_helix = n - 1
                    end_helix = n + 1

                #starting unpaired nt
                elif (x1,x2,y1,y2) == (1,0,0,0) or (x1,x2,y1,y2) == (0,0,1,0):
                    loop_length += 1
                    begin_helix = n - 1

                #ending unpaired nt
                elif (x1,x2,y1,y2) == (0,1,0,0) or (x1,x2,y1,y2) == (0,0,0,1):
                    loop_length += 1
                    end_helix = n + 1

                #1,0,1,0 is impossible w/o psuedoknots
                #0,1,0,1 is impossible w/o psuedoknots
                #Also, all binary combinations with 3 or 4 true are impossible (n-1 or n+1 can not be in both bp_x and bp_y)


            elif loop_length > 0:
                #Bulge or loop?
                #loop

                #print "begin = ", begin_helix
                #print "end = ", end_helix

                if ( bp_x.count(begin_helix) > 0 and bp_y.count(end_helix) > 0 and bp_x.index(begin_helix) == bp_y.index(end_helix) ):
                    helical_loop_list.append(loop_length)
                    loop_length = 0

                    if output_start_end:
                        #Also return the starting and ending positions of each loop/bulge
                        helical_loop_start_end.append((begin_helix,end_helix))

                else:

                    bp_end = 0
                    bp_begin = 0

                    if (bp_x.count(end_helix) > 0): bp_begin = bp_y[bp_x.index(end_helix)]
                    if (bp_y.count(end_helix) > 0): bp_end = bp_x[bp_y.index(end_helix)]

                    if (bp_x.count(begin_helix) > 0): bp_end = bp_y[bp_x.index(begin_helix)]
                    if (bp_y.count(begin_helix) > 0): bp_begin = bp_x[bp_y.index(begin_helix)]

                    if bp_end > bp_begin:
                        bulge_loop_list.append(loop_length)
                        loop_length = 0

                        if output_start_end:
                            #Also return the starting and ending positions of each loop/bulge
                            bulge_loop_start_end.append((begin_helix,end_helix))

                    else:
                        loop_length = 0
        if output_start_end:
            return (helical_loop_list, bulge_loop_list, helical_loop_start_end, bulge_loop_start_end)
        else:
            return (helical_loop_list, bulge_loop_list)

    #**********************************************************************************************************
    #New standby site module calling functions + new dG_standby calculation (Alex Lukasiewicz 07/2020)
    #Note: should be used with defined start codon position. 

    def module_searcher(self,bracket_string,bp_x,bp_y):
        modules = {}
        module_ranges = []
        y_anchor = 0
        modulecount = 0

        #search for start and stop sites for sequence with multiple modules
        for i in range(len(bp_x)-1):
            if bp_y[i+1] > bp_y[i]:
                if len(module_ranges) == 0:
                    modules['module' +  str(modulecount)] = {}
                    module_ranges.append([0,bp_x[i+1]])
                    modules['module'+str(modulecount)]['module_start'] = 0
                    modules['module'+str(modulecount)]['module_end'] = bp_x[i+1]
                    modules['module'+str(modulecount)]['module_structure'] = bracket_string[0:bp_x[i+1]-1]
                    modules['module'+str(modulecount)]['bp_x'] = bp_x[0:i+1]                
                    modules['module'+str(modulecount)]['bp_y'] = bp_y[0:i+1]
                    y_anchor = bp_y[0]
                    x_anchor = bp_x[i+1]
                    modulecount += 1
                    #print(modulecount)
                elif len(module_ranges) >= 1:
                    modules['module' +  str(modulecount)] = {}
                    module_ranges.append([y_anchor,bp_x[i+1]])
                    modules['module' +  str(modulecount)]['module_start'] = y_anchor
                    modules['module' +  str(modulecount)]['module_end'] = bp_x[i+1]
                    modules['module'+ str(modulecount)]['module_structure'] = bracket_string[y_anchor:bp_x[i+1]]
                    modules['module' +  str(modulecount)]['bp_x'] = bp_x[bp_x.index(x_anchor):i]                
                    modules['module' +  str(modulecount)]['bp_y'] = bp_y[bp_x.index(x_anchor):i] 
                    y_anchor=bp_x.index(x_anchor)
                    x_anchor = bp_x[i+1]
                    modulecount +=1
            if modulecount >= 1 and i == len(bp_x) - 2:
                modules['module' +  str(modulecount)] = {}
                module_ranges.append([y_anchor,len(bracket_string)])
                modules['module' +  str(modulecount)]['module_start'] = y_anchor
                modules['module' +  str(modulecount)]['module_end'] = len(bracket_string)
                modules['module'+str(modulecount)]['module_structure'] = bracket_string[y_anchor:len(bracket_string)]
                modules['module' +  str(modulecount)]['bp_x'] = bp_x[bp_x.index(x_anchor):i+1]                
                modules['module' +  str(modulecount)]['bp_y'] = bp_y[bp_x.index(x_anchor):i+1] 
                #y_anchor=bp_y[i]
                #x_anchor = bp_x[i+1]
                modulecount +=1
                #print(modules)
            else:
                continue 
            
        #for standby sites with only one module 
        if len(module_ranges) == 0:
            modules['module' +  str(modulecount)] = {}
            module_ranges.append([0,bp_x[i+1]])
            modules['module'+str(modulecount)]['module_start'] = 0
            modules['module'+str(modulecount)]['module_end'] = len(bracket_string)
            modules['module'+str(modulecount)]['module_structure'] = bracket_string[0:len(bracket_string)]
            modules['module'+str(modulecount)]['bp_x'] = bp_x                
            modules['module'+str(modulecount)]['bp_y'] = bp_y
            y_anchor = bp_y[0]
            x_anchor = bp_x[i+1]
            #print(modules)
        if len(module_ranges) > 1:
            modules = orderedDict(sorted(modules.items()))

        return module_ranges, modules


    def calc_dG_distort(self,structure,start,stop,bp_x,bp_y):
        #determine hairpin height
        distance = []
        if len(bp_x) > 0: 
            for i in range(len(bp_x)):
                distance.append(bp_y[i]-bp_x[i]-1)
            H = (bp_x.pop() - bp_x[0]+1) + (min(distance)/2)

            #calculate P and D
            P = bp_x[0] - (start + 1)
            D = stop - (bp_y[0] + 1)
            print('proximal: ' + str(P) + " hairpin: " + str(H) + " distal: " + str(D))

            #from Espah-borujeni et. al., 2014
            As = 15 + P + D - H 
            dG_distort = 0.038*(As**2) + -1.629*(As) + 17.359  
        else:
            dG_distort = 20 #apply penalty for dissolution of module 
            H = 0

        return dG_distort, H

    def calc_dG_standby(self,module_dict, sequence):  #for compatibility with standby site modules as detailed in RBS calc v2.0 

        RNA_model = "rna37"
        dangles = "none"
        temp = "37.0"
        RT = 0.616

        all_dG_distort = []
        all_dG_unfold = []
        all_dG_standby = []
        num_modules = len(module_dict)
        #print(num_modules)

        for module,features in module_dict.items():
            dG_distort = []
            dG_unfold = []
            dG_standby = []

            if len(module_dict) == 1:
                #seq_slice = sequence[features['module_start']:(features['module_end']+1)] 
                mod_bp_x = features['bp_x']
                mod_bp_y = features['bp_y']
                #print(mod_bp_y,mod_bp_x)

                #calculate dG distort and dG unfold for initial structure 
                dG_distort.append(calc_dG_distort(features['module_structure'],features['module_start'],features['module_end'],mod_bp_x,mod_bp_y))
                constraints = None
                fold = ViennaRNA([sequence], RNA_model)
                fold.mfe([1], constraints, Temp = temp, dangles = dangles)
                no_unfold = fold["mfe_energy"][0]
                dG_unfold.append(fold['mfe_energy'][0]-no_unfold)
                
                #begin modifying hairpin 
                constraints = '.' * (mod_bp_x[0]-1) + 'x' + '.' * (mod_bp_y[0] - (mod_bp_x[0]+1)) + 'x' + '.' * (len(sequence)- mod_bp_y[0])
                #print(constraints)
                fold = ViennaRNA([sequence], RNA_model)
                fold.mfe([1], constraints, Temp = temp, dangles = dangles)
                dG_unfold.append(fold["mfe_energy"][0]- no_unfold)
                dG_distort.append(calc_dG_distort(constraints, features['module_start'],features['module_end'], fold['mfe_basepairing_x'][0], fold['mfe_basepairing_y'][0]))
                #print(constraints, fold['mfe_basepairing_x'][0], fold['mfe_basepairing_y'][0])

                #iterate through possible dG distort + unfolding penalties
                for i in range(1,4):
                    constraints = list(constraints)
                    constraints[mod_bp_x[i]-1] = 'x'
                    constraints[mod_bp_y[i]-1] = 'x'
                    constraints = "".join(constraints)
                    fold = ViennaRNA([sequence], RNA_model)
                    fold.mfe([1], constraints, Temp = temp, dangles = dangles)
                    dG_unfold.append(fold["mfe_energy"][0] - no_unfold)

                    fold_bp_x = []
                    fold_bp_y = []

                    for i in range(len(fold['mfe_basepairing_x'][0])):
                        if features['module_start'] <= fold['mfe_basepairing_x'][0][i] and fold['mfe_basepairing_y'][0][i] < features['module_end']:
                            fold_bp_x.append(fold['mfe_basepairing_x'][0][i])
                            fold_bp_y.append(fold['mfe_basepairing_y'][0][i])

                    if len(fold['mfe_basepairing_x']) > 0:
                        dG_distort.append(calc_dG_distort(constraints,features['module_start'],features['module_end'], fold_bp_x, fold_bp_y))
                    else:
                        break
                    #print(constraints, fold_bp_x,fold_bp_y)

                for i in range(len(dG_distort)):
                    dG_standby.append(dG_unfold[i] + dG_distort[i][0])
                dG_standby_min = min(dG_standby)
            
                all_dG_standby.append(dG_standby_min)
                #print(dG_unfold,dG_distort,dG_standby)
            
            else: #calculate dG standby with dG sliding penalty 
                #seq_slice = sequence[features['module_start']:(features['module_end']+1)] 
                mod_bp_x = features['bp_x']
                mod_bp_y = features['bp_y']
                #print(mod_bp_y,mod_bp_x)

                #calculate dG distort and dG unfold for initial structure 
                dG_distort.append(calc_dG_distort(features['module_structure'],features['module_start'],features['module_end'],features['bp_x'],features['bp_y']))
                constraints = None
                fold = ViennaRNA([sequence], RNA_model)
                fold.mfe([1], constraints, Temp = temp, dangles = dangles)
                no_unfold = fold["mfe_energy"][0]
                dG_unfold.append(fold['mfe_energy'][0]-no_unfold)
                
                #begin modifying hairpin 
                constraints = '.' * (mod_bp_x[0]-1) + 'x' + '.' * (mod_bp_y[0] - (mod_bp_x[0]+1)) + 'x' + '.' * (len(sequence)- mod_bp_y[0])
                #print('this is the initial constraint string ' + str(constraints))
                fold = ViennaRNA([sequence], RNA_model)
                fold.mfe([1], constraints, Temp = temp, dangles = dangles)
                dG_unfold.append(fold["mfe_energy"][0]- no_unfold)

                fold_bp_x = []
                fold_bp_y = []

                for i in range(len(fold['mfe_basepairing_x'][0])):
                    if features['module_start'] <= fold['mfe_basepairing_x'][0][i] and fold['mfe_basepairing_y'][0][i] < features['module_end']:
                        fold_bp_x.append(fold['mfe_basepairing_x'][0][i])
                        fold_bp_y.append(fold['mfe_basepairing_y'][0][i])
                    else:
                        continue

                #print(fold_bp_x,fold_bp_y)
                dG_distort.append(calc_dG_distort(constraints, features['module_start'],features['module_end'],fold_bp_x, fold_bp_y))
                #print(constraints, fold['mfe_basepairing_x'][0], fold['mfe_basepairing_y'][0])

                #iterate through possible dG distort + unfolding penalties
                for i in range(1,4):
                    constraints = list(constraints)
                    #print(mod_bp_x[i],mod_bp_y[i])
                    constraints[mod_bp_x[i]-1] = 'x'
                    constraints[mod_bp_y[i]-1] = 'x'
                    constraints = "".join(constraints)
                    #print('the constraint string passed is ' + str(constraints))
                    fold = ViennaRNA([sequence], RNA_model)
                    fold.mfe([1], constraints, Temp = temp, dangles = dangles)
                    dG_unfold.append(fold["mfe_energy"][0] - no_unfold)
                    fold_bp_x = []
                    fold_bp_y = []

                    #print(fold['mfe_basepairing_x'][0])
                    #print(fold['mfe_basepairing_y'][0])


                    for i in range(len(fold['mfe_basepairing_x'][0])):
                        if features['module_start'] <= fold['mfe_basepairing_x'][0][i] and fold['mfe_basepairing_y'][0][i] < features['module_end']:
                            fold_bp_x.append(fold['mfe_basepairing_x'][0][i])
                            fold_bp_y.append(fold['mfe_basepairing_y'][0][i])
                        else:
                            continue 

                    #print(fold['mfe_basepairing_x'][0])
                    #print(fold['mfe_basepairing_y'][0])
                    #print(fold_bp_x)
                    #print(fold_bp_y)
                    if len(fold_bp_x) > 1:
                        dG_distort.append(calc_dG_distort(constraints,features['module_start'],features['module_end'], fold_bp_x, fold_bp_y))
                        #print('the dG distort and haripin height are ' + str(dG_distort))
                    else:
                        break
                all_dG_distort.append(dG_distort)
                all_dG_unfold.append(dG_unfold)
                #calculate minimum dG_standby for individual module
                for i in range(len(dG_distort)):
                    dG_standby.append(dG_unfold[i] + dG_distort[i][0])
                    #+ sum(all_dG_distort[][1])*0.02)
                
                dG_standby_min = min(dG_standby)
                
                all_dG_standby.append(dG_standby_min)
                #print(all_dG_standby)
        #after iterating through all modules, find minimum dG distort + unfold + sliding 
        find_min_dG_standby = []
        heights = []
        for i in range(len(all_dG_distort)):
            heights.append(all_dG_distort[i:][0][0][1])
        #print(heights)

        for i in range(len(all_dG_distort)-1):
            find_min_dG_standby.append(all_dG_standby[i] + sum(heights[i+1:])*0.02)
        find_min_dG_standby.append(all_dG_standby.pop()) # add last value with no downstream modules, lazily 
        dG_standby = min(find_min_dG_standby)
        return dG_standby

    def cpu_time(self):
        import resource
        return resource.getrusage(resource.RUSAGE_SELF)[0]

    def combsort(self, input_list):
        """Sorts the input_list in increasing order (minimum first) according to the comb sort algorithm from Wikipedia. Outputs the corresponding sorted index. """

        index = list(range(len(input_list)))

        #Implementation of comb sort (from Wikipedia)
        shrink_factor = 1.24733095

        #Init
        gap = len(input_list)
        swaps = False

        while not (gap <= 1 and not swaps):

            #update the gap value for a next comb
            if gap > 1:
                gap = int(gap / shrink_factor)
                if gap == 10 or gap == 9:
                    gap = 11

            i = 0
            swaps = False #see bubblesort for an explanation

            #a single "comb" over the input list
            while gap + i < len(input_list):
                val = input_list[i]
                if val > input_list[i+gap]:
                    #Swap values -- verify that lists are being altered correctly
                    input_list[i] = input_list[gap + i]
                    input_list[gap + i] = val

                    #Swap indices of index list accordingly
                    temp = index[i]
                    index[i] = index[gap + i]
                    index[gap + i] = temp

                    swaps = True

                i += 1

        return index

    def calc_dG(self):
        """Calculates each dG term in the free energy model and sums them together to create dG_total"""

        start = self.cpu_time()

        #Initialization of data structures
        self.start_pos_list = []
        self.dG_total_list = []
        self.dG_mRNA_list = []
        self.dG_mRNA_rRNA_list = []
        self.fold_x_list = []
        self.fold_y_list = []
        self.dG_start_energy_list = []
        self.dG_spacing_list = []
        self.mRNA_structure_list = []
        self.mRNA_rRNA_uncorrected_structure_list = []
        self.mRNA_rRNA_corrected_structure_list = []
        self.dG_standby_site_list = []
        self.dG_spacer_site_list = []
        self.kinetic_score_list = []
        self.min_bp_prob_list = []
        self.longest_helix_list = []
        self.three_state_indicator_list = []
        self.helical_loop_list_list = []
        self.bulge_loop_list_list = []

        self.dS1_list = []
        self.dS2_list = []
        self.most_5p_mRNA_list = []
        self.Expression_list = []

        for (start_pos, codon) in self.find_start_codons(self.mRNA_input):

            try:

                #print "Top of calc_dG here"

                #Set dangles based on length between 5' end of mRNA and start codon
                if self.auto_dangles:

                    if start_pos > self.cutoff:
                        self.dangles = "none"

                    else:
                        self.dangles = "all"

                else:
                    self.dangles = self.dangles_default
                    #print "Auto Dangles set to ", self.dangles

                #Start codon energy
                dG_start_codon = self.start_codon_energies[codon]

                #Energy of mRNA folding
                [dG_mRNA,mRNA_structure] = self.calc_dG_mRNA(start_pos)

                #Energy of mRNA:rRNA hybridization & folding
                [dG_mRNA_rRNA_withspacing,mRNA_rRNA_structure] = self.calc_dG_mRNA_rRNA(start_pos)

                dG_mRNA_rRNA_nospacing = mRNA_rRNA_structure["dG_mRNA_rRNA"]

                #Standby site correction:
                dG_standby_site = self.calc_dG_standby_site(mRNA_rRNA_structure, rRNA_binding = True)

                #Total energy is mRNA:rRNA + start - rRNA - mRNA - standby_site
                dG_total = dG_mRNA_rRNA_withspacing + dG_start_codon + dG_standby_site - dG_mRNA #edited for new standby site functions

                #Calculate 'kinetic score': directly related to probability of base pair formation
                (kinetic_score,min_bp_prob) = self.calc_kinetic_score(mRNA_structure)

                #Calculate dG to open SD sequence
                #ddG_SD_open = self.calc_dG_SDopen(mRNA_structure, mRNA_rRNA_structure)


                self.mRNA_structure_list.append(mRNA_structure)
                #self.mRNA_rRNA_uncorrected_structure_list.append(mRNA_rRNA_structure)
                #self.mRNA_rRNA_corrected_structure_list.append(corrected_structure)
                self.most_5p_mRNA_list.append(self.calc_most_5p_mRNA(mRNA_rRNA_structure))

                (helical_loop_list, bulge_loop_list) = self.calc_longest_loop_bulge(structure=mRNA_structure)
                self.longest_helix_list.append(self.calc_longest_helix(structure=mRNA_structure))

                self.dG_start_energy_list.append(dG_start_codon)
                self.dG_mRNA_list.append(dG_mRNA)
                self.dG_mRNA_rRNA_list.append(dG_mRNA_rRNA_nospacing)
                self.dG_spacing_list.append(mRNA_rRNA_structure["dG_spacing"])
                self.dG_standby_site_list.append(dG_standby_site)
                self.dG_total_list.append(dG_total)

                self.helical_loop_list_list.append(helical_loop_list)
                self.bulge_loop_list_list.append(bulge_loop_list)

                self.min_bp_prob_list.append(min_bp_prob)
                self.kinetic_score_list.append(kinetic_score)
                #self.three_state_indicator_list.append(ddG_SD_open)

                #Start positions
                self.start_pos_list.append(start_pos)

                #Expression levels
                self.Expression_list.append(self.calc_expression_level(dG_total))

                #For exporting the relevant structure to a PDF
                #index = mRNA_rRNA_structure["MinStructureID"]
                #mRNA_rRNA_structure.export_PDF(index, name = self.name + ": Before standby site", filename =  self.name + "_Before_Standby_rRNA.pdf", program = "subopt")

                #index = corrected_structure["MinStructureID"]
                #corrected_structure.export_PDF(index, name = self.name + ": After standby site", filename =  self.name + "_After_Standby_rRNA.pdf", program = "subopt")

            except CalcError as msg:
                print(msg)
                self.mRNA_structure_list.append([])
                #self.mRNA_rRNA_uncorrected_structure_list.append([])
                #self.mRNA_rRNA_corrected_structure_list.append([])

                self.most_5p_mRNA_list.append(self.infinity)
                self.longest_helix_list.append(self.infinity)

                self.dG_start_energy_list.append(self.infinity)
                self.dG_mRNA_list.append(self.infinity)
                self.dG_mRNA_rRNA_list.append(self.infinity)
                self.dG_spacing_list.append(self.infinity)
                self.dG_standby_site_list.append(self.infinity)
                self.dG_total_list.append(self.infinity)

                self.helical_loop_list_list.append(self.infinity)
                self.bulge_loop_list_list.append(self.infinity)

                self.kinetic_score_list.append(self.infinity)
                #self.three_state_indicator_list.append(self.infinity)

        self.run = 1

        end = self.cpu_time()
        self.run_time = end - start

    def calc_expression_level(self, dG):

        import math
        return RBS_Calculator.K * math.exp(-dG / RBS_Calculator.RT_eff)

    def print_dG(self,max_dG = 1e12,brief = 0, return_string = False, print_expression = False):
        '''Print out useful information about the mRNA sequence'''

        import math

        print_string = ""
        if self.run == 1:

            print("mRNA Sequence: ", self.name)
            if return_string: print_string = print_string + "mRNA Sequence: " + self.name + "\n"

            if len(self.start_position_list) == 0:
                print("No start codons found in input mRNA sequence")
                print("----------------------------------------------------------------------------------------")
                if return_string: print_string = print_string + "No start codons found in input mRNA sequence" + "\n" + "----------------------------------------------------------------------------------------" + "\n"
            elif len( [dG_total for dG_total in self.dG_total_list if dG_total < max_dG] ) == 0:
                print("No RBSs found with dG <", str(max_dG), ".")
                print("----------------------------------------------------------------------------------------" + "\n")
                if return_string: print_string = print_string + "No RBSs found with dG <" + str(max_dG) + "." + "\n" + "----------------------------------------------------------------------------------------" + "\n"

            else:

                if print_expression:
                    Headers = ("Start", "(pos)", "Expression level","Kinetic score")
                    format = "%4s %4s %15s %15s"
                else:
                    Headers = ("Start","(pos)", "dG total", "dG (rRNA:mRNA)", "dG (mRNA)", "dG (spacing)", "dG (standby)", "Kinetic Score")
                    format = "%4s %4s %12s %12s %12s %12s %12s %12s"

                print(format % Headers)

                if return_string: print_string = print_string + format % Headers + "\n"

                for (start,counter) in zip(self.start_position_list,list(range(len(self.start_position_list)))):

                    dG_total = self.dG_total_list[counter]
                    Expression = RBS_Calculator.K * math.exp(-dG_total / RBS_Calculator.RT_eff)

                    if (dG_total < max_dG):

                        if (dG_total > self.infinity):
                            dG_total = "Inf"
                            Expression = 0.0

                        dG_mRNA = self.dG_mRNA_list[counter]
                        if (dG_mRNA > self.infinity): dG_mRNA = "Inf"

                        dG_mRNA_rRNA = self.dG_mRNA_rRNA_list[counter]
                        if (dG_mRNA_rRNA > self.infinity): dG_mRNA_rRNA = "Inf"

                        dG_spacing = self.dG_spacing_list[counter]
                        if (dG_spacing > self.infinity): dG_spacing = "Inf"

                        dG_standby_site = self.dG_standby_site_list[counter]
                        if (dG_standby_site > self.infinity): dG_standby_site = "Inf"

                        start_codon = self.start_codon_list[counter]

                        dG_rRNA = self.dG_rRNA

                        kinetic_score = self.kinetic_score_list[counter]

                        if print_expression:
                            print(format % (start_codon, str(start), str(round(Expression,2)), str(round(kinetic_score,2))))

                        else:
                            print(format % (start_codon, str(start), str(dG_total), str(dG_mRNA_rRNA), str(dG_mRNA), str(dG_spacing), str(dG_standby_site), str(round(kinetic_score,2))))

                        if return_string: print_string = print_string + format % (start_codon, str(start), str(dG_total), str(dG_mRNA_rRNA), str(dG_mRNA), str(dG_spacing), str(dG_standby_site), str(round(kinetic_score,2))) + "\n"

                print("----------------------------------------------------------------------------------------")
                if return_string: print_string = print_string + "----------------------------------------------------------------------------------------" + "\n"

                #print "Computation Time: ", str(self.run_time), " seconds."

                if return_string: return print_string
        else:
            raise RuntimeError("The RBS Calculator has not been run yet. Call the 'calc_dG' method.")
    #def run(self):


    #def output(self):


    def save_data(self, handle, header = False):

        infinity = 1e6

        if self.run == 1:

            if (header):
                #Print header information --> all the parameters
                parameters = "RNA model: \t%s \nDangles: \t%s \nTemperature: \t%s \nOptimal Spacing \t%s \nSpacing constant (push) \t%s\nSpacing constant (pull) \t%s\nCutoff \t%s\nStandby Site Length \t%s\nEnergy Cutoff \t%s\nrRNA sequence \t%s\n" % (self.RNA_model, str(self.dangles), str(self.temp), str(self.optimal_spacing), str(self.dG_spacing_constant_push), str(self.dG_spacing_constant_pull), str(self.cutoff), str(self.standby_site_length), str(self.energy_cutoff), self.rRNA)

                handle.writelines(parameters)

            #Name, dG total, dG rRNA:mRNA, dG mRNA, dG spacing, dG standby, dG start codon, kinetic score, longest helix, longest loop, start position
            #format = "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n"
            format = "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n"

            for counter in range(len(self.start_position_list)):

                dG_total = self.dG_total_list[counter]
                if (dG_total >= infinity): dG_total = "Inf"

                dG_mRNA = self.dG_mRNA_list[counter]
                if (dG_mRNA >= infinity): dG_mRNA = "Inf"

                dG_mRNA_rRNA = self.dG_mRNA_rRNA_list[counter]
                if (dG_mRNA_rRNA >= infinity): dG_mRNA_rRNA = "Inf"

                dG_spacing = self.dG_spacing_list[counter]
                if (dG_spacing >= infinity): dG_spacing = "Inf"

                dG_standby_site = self.dG_standby_site_list[counter]
                if (dG_standby_site >= infinity): dG_standby_site = "Inf"

                kinetic_score = self.kinetic_score_list[counter]

                longest_helix = self.longest_helix_list[counter]
                if len(self.helical_loop_list_list[counter])>0:
                    longest_loop = max(self.helical_loop_list_list[counter])
                else:
                    longest_loop = 0

                most_5p_mRNA = self.most_5p_mRNA_list[counter]

                #three_state = self.three_state_indicator_list[counter]

                handle.writelines(format % (self.name.split(" ")[0],dG_total, dG_mRNA_rRNA,dG_mRNA,dG_spacing,dG_standby_site, self.start_codon_energies[self.start_codon_list[counter]],round(kinetic_score,2), str(longest_helix), str(self.start_position_list[counter]), str(most_5p_mRNA), str(three_state) ))

        else:
            raise RuntimeError("The RBS Calculator has not been run yet. Call the 'calc_dG' method.")
    

#----------------------------------------------------------------------------------------------------------
#End RBS_Calculator class
#----------------------------------------------------------------------------------------------------------
def calc_dG_from_file(handle, output, verbose = True, parameters = {}):
    from Bio import SeqIO
    from RBS_Calculator import RBS_Calculator

    records = SeqIO.parse(handle,"fasta")
    First = True
    export_PDF = True

    #for i in range(30):
    #    records.next()

    for record in records:

        mRNA = record.seq.tostring().upper()

        #Set any defaults
        start_range = [0, len(mRNA)]
        name = record.description.split(" ")[0]

        #Create instance of RBS Calculator
        test = RBS_Calculator(mRNA, start_range, name)

        #Examine kvars dictionary and pull out any options. Assign them to instanced class.
        for (key,value) in list(parameters.items()):

            if key == "cutoff":
                test.cutoff = value
            elif key == "start_range":
                test.start_range = value
            elif key == "rRNA":
                test.rRNA = value
            elif key == "energy_cutoff":
                test.energy_cutoff = value
            elif key == "standby_site_length":
                test.standby_site_length = value
            elif key == "dangles":
                test.dangles = value
            elif key == "export_PDF":
                export_PDF = value

        test.calc_dG()
        test.print_dG(test.infinity,print_expression=verbose)
        test.save_data(output, First)

        if First:
            First = False

        if export_PDF:
            num_structs = len(test.mRNA_rRNA_uncorrected_structure_list)
            for (structure,counter) in zip(test.mRNA_rRNA_uncorrected_structure_list,list(range(num_structs))):
                index = structure["MinStructureID"]
                structure.export_PDF(index, name, filename =  name + "_rRNA" + "_" + str(counter) + ".pdf", program = "subopt")

            num_structs = len(test.mRNA_structure_list)
            for (structure,counter) in zip(test.mRNA_structure_list,list(range(num_structs))):
                structure.export_PDF(0, name, filename = name + "_mRNA" + "_" + str(counter) + ".pdf")

    output.close()

def calc_dG_pre_post_RBS(pre_list,post_list,RBS_list,name_list,output,verbose = True, parameters = {}):

    from RBS_Calculator import RBS_Calculator

    First = True

    for (pre,post,RBS,name) in zip(pre_list,post_list,RBS_list,name_list):

        mRNA = pre + RBS + post

        start_range = [0, len(mRNA)]

        #Create instance of RBS Calculator
        test = RBS_Calculator(mRNA, start_range, name)

        #Examine kvars dictionary and pull out any options. Assign them to instanced class.
        for (key,value) in list(parameters.items()):

            if key == "cutoff":
                test.cutoff = value
            elif key == "start_range":
                test.start_range = value
            elif key == "rRNA":
                test.rRNA = value
            elif key == "energy_cutoff":
                test.energy_cutoff = value
            elif key == "standby_site_length":
                test.standby_sitRBSe_length = value
            elif key == "dangles":
                test.dangles = value
            elif key == "export_PDF":
                export_PDF = value

        test.calc_dG()
        if verbose: test.print_dG(test.infinity)

        test.save_data(output, First)
        if First:
            First = False

    output.close()

if __name__ == "__main__":

    #from Bio import SeqIO



    #filename = "/common/RBS_Calculator/DataSets/Amin_RBSs.txt"
    #filename = "/common/RBS_Calculator/DataSets/Forward_Predictions.txt"
    #filename = "/common/RBS_Calculator/DataSets/Context_Tests_RBSs_2nd.txt"
    #filename = "/common/RBS_Calculator/DataSets/DNA20_CDSs.fasta"
    #output_filename = "/common/RBS_Calculator/Output_Forward_Predictions_dangles_all_TTC.txt"
    #output_filename = "/common/RBS_Calculator/Output_Context_Tests_RBSs_test.txt"
    #filename = "/common/RBS_Calculator/DataSets/All_Tested_RBSs.txt"
    #filename = "E:\RBS_Calculator\DataSets\All_Tested_RBSs.txt"
    #filename = "/common/RBS_Calculator/DataSets/Testing_RBSs.txt"
    #filename = "/common/RBS_Calculator/DataSets/Karsten_nif_RBSs.txt"
    #filename = "/common/RBS_Calculator/Dan_RBS_silk_sequences.fasta"
    #output_filename = "/common/RBS_Calculator/Output_All_Tested_RBSs_newspacing.txt"
    #output_filename = "/common/RBS_Calculator/Output_DNA20_Designed_RBS_1.txt"
    #output_filename = "/common/RBS_Calculator/Output_Context_Tests_RBSs_3.txt"
    #output_filename = "/common/RBS_Calculator/Output_Karsten_nif_RBSs.txt"
    #output_filename = "/common/RBS_Calculator/Amin_RBSs_Output.txt"

    #filename = "/common/RBS_Calculator/DataSets/Crt_version_RBSs.txt"

    #filename = "/common/RBS_Calculator/DataSets/Output_Crt_natural_RBSs.txt"
    #filename = "/common/RBS_Calculator/DataSets/RBS_Spacing.txt"
    #filename = "/common/RBS_Calculator/DataSets/Ethan_RBSs_all.txt"
    #filename = "/common/RBS_Calculator/DataSets/Dan_sicP.fasta"
    #filename = "/common/RBS_Calculator/DataSets/Natural_prg_org_operon.str"


    #output_filename = "/common/RBS_Calculator/Dan_Output_silk_RBSs.txt"
    #handle = open(filename,"rU")
    #output = open(output_filename,"w")

    #verbose = False

    #pars = {"start_range":[18,71],"pre_window_length":0,"post_window_length":0,"export_PDF":False}
    #calc_dG_from_file(handle, output, verbose,pars)
    #handle.close()
    #output.close()

    #window_list = list(range(10,65))
    #for window in window_list:
    #    print("Cutoff = ", window)
    #    output_filename = "../data/RBS_Calculator/DataSets/cutoff/Output_All_Tested_RBSs_cutoff_" + str(window) + ".txt"

    #    handle = open(filename,"rU")
    #    output = open(output_filename,"w")

    #   verbose = True

    #   pars = {"start_range":[18,44],"export_PDF":False,"cutoff":window}
    #   calc_dG_from_file(handle, output, verbose,pars)
    #   handle.close()
    #   output.close()

    #pre = "ggggaattgtgagcggataacaattcccctctagaa"
    #RBS = "GAGGAAGTTAGTAAGGAGGTCAGGCGA"

    #pre_list = []
    #RBS_list = []
    #CDS_list = []
    #name_list = []
    ##Get protein coding sequence
    #records = SeqIO.parse(handle,"fasta")
    #for record in records:
        #CDS_list.append(record.seq.tostring().upper())
        #pre_list.append(pre)
        #RBS_list.append(RBS)
        #name_list.append(record.description)

    #calc_dG_pre_post_RBS(pre_list,CDS_list,RBS_list,name_list,output,verbose = True, parameters = pars)

    #with open('eb_summary.csv') as f: 
    #    entries = csv.reader(f,delimiter= ',')
        #next(entries)
        #print(entries)
    #    start_codons = {}
     #   sequences = {}
    #    for row in entries:
    #        try:
    #           start_codons[row[0]] = int(row[1])
    #           sequences[row[0]] = str(row[2])
    #        except(ValueError):
     #           break

    #return start_codons, sequences
    #standby_site_list = []
    #utr_list = []
    #dG_total_list = []
    #for key, value in start_codons.items():
    #    print('this is the start site ' + str(start_codons[key]))
    #    print('this is the sequence ' + str(sequences[key]))
    #    print('this is the UTR name ' + str(key))
    #    start_range = [start_codons[key] - 1, start_codons[key]]
    #    constraint_str = None
    #    name = key
    #    calcObj = RBS_Calculator(sequences[key],start_range,'ACCTCCTTA', constraint_str, name)
    #    calcObj.calc_dG()
        
    #    dG_total_list.append(calcObj.dG_total_list[0])
    #    standby_site_list.append(calcObj.dG_standby_site_list[0])
    #    utr_list.append(name)
    
    #    expr_list = []
    #for dG in dG_total_list:
    #    expr_list.append(calcObj.K * math.exp(-dG/calcObj.RT_eff))
    #for utr, standby in zip(utr_list, standby_site_list):
    #    print(utr, standby)    
    
    name = "test"
    #mRNA = "TTCTAGAGGGGGGATCTCCCCCCAAAAAATAAGAGGTACACATGACTAAAACTTTCAAAGGCTCAGTATTCCCACTGAG"
    mRNA = "TCTAGAGGCCGACGCAAGCCCATATCGGGGCTTCCGTCGGCCAAAAAAAAAAAAAAATAAGGAGGTAAAAAATG"
    constraint_str = None
    start_range = [0,len(mRNA)]
    

    calcObj = RBS_Calculator(mRNA,start_range,'ACCTCCTTA', constraint_str, name)
    calcObj.calc_dG()

    dG_total_list = calcObj.dG_total_list[:]
    start_pos_list = calcObj.start_pos_list[:]
    kinetic_score_list = calcObj.kinetic_score_list[:]
    standby_site_list = calcObj.dG_standby_site_list[:]
    mRNA_structure_list = calcObj.mRNA_structure_list[:]
    dG_start_list = calcObj.dG_start_energy_list[:]
    dG_mRNA_list = calcObj.dG_mRNA_list[:]
    dG_mRNA_rRNA = calcObj.dG_mRNA_rRNA_list[:]

    expr_list = []
    for dG in dG_total_list:
        expr_list.append(calcObj.K * math.exp(-dG/calcObj.RT_eff))
        print(dG)#math.exp(-dG/calcObj.RT_eff))

    #print(len(expr_list))
    for (expr,start_pos,start,mfe,mRNA_rRNA, dG, dstandby) in zip(expr_list,start_pos_list,dG_start_list,dG_mRNA_list,dG_mRNA_rRNA, dG_total_list, standby_site_list):
        print(expr,start_pos,start,mfe,mRNA_rRNA, dG, dstandby)

    ##test.post_window_length = 19
    #test.auto_dangles = False
    #test.default_dangles = "all"
    #test.pre_window_length = 1000
    #test.post_window_length = 20
    #test.cutoff = 50
    #test.calc_dG()
    #test.print_dG(test.infinity)


    #if export_PDF:
            #num_structs = len(test.mRNA_rRNA_corrected_structure_list)
            #for (structure,counter) in zip(test.mRNA_rRNA_corrected_structure_list,range(num_structs)):

                #index = structure["MinStructureID"]
                #structure.export_PDF(index, name, filename =  name + "_rRNA" + "_" + str(counter) + ".pdf", program = "energy")

            #for (structure,counter) in zip(test.mRNA_rRNA_uncorrected_structure_list,range(num_structs)):

                #index = structure["MinStructureID"]
                #structure.export_PDF(index, name, filename =  name + "_rRNA_uncorrected" + "_" + str(counter) + ".pdf", program = "energy")

            #num_structs = len(test.mRNA_structure_list)
            #for (structure,counter) in zip(test.mRNA_structure_list,range(num_structs)):
                #structure.export_PDF(0, name, filename = name + "_mRNA" + "_" + str(counter) + ".pdf")
