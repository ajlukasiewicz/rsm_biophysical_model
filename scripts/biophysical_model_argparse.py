import sys, os, math, xlwt, xlrd, random, traceback, pickle, time 

#sys.path.append('../../')
#from DNAc import *
import argparse
import itertools
#import RNA
import pandas as pd
from RBS_Calculator import *
#from RBS_Calculator_Vienna import *
from ViennaRNA import *

RNAEnergyModel = ViennaRNA

beta = 0.45

#Calculate the CsrA-RNA free energy interactions on arbitrary RNA sequences
#Identify CsrA binding sites

#Model variant A:  on line 50, include/exclude dG_dimerization 
#Model variant B:  lines 156-162: include/exclude dG penalty for having a CsrA motif outside a hairpin loop 


#ALL FREE ENERGIES IN UNITS OF RT / kT

#binding (AAGGA) = -2.63 + -2.2 + -3.14 + -3.14 + -1.65 in units of RT
CsrA_freeEnergyMatrix = [   {'A' : -2.63, 'T' : 0, 'C' : 0, 'G' : 0}, #1st nucleotide position
                            {'A' : -2.20, 'T' : 0, 'C' : 0, 'G' : 0}, #2nd nucleotide position
                            {'G' : -3.14, 'T' : 0, 'C' : 0, 'A' : 0}, #3rd nucleotide position
                            {'G' : -3.14, 'T' : 0, 'C' : 0, 'A' : 0}, #4th nucleotide position
                            {'A' : -1.65, 'T' : 0, 'C' : 0, 'G' : 0}, #5th nucleotide position
                        ]

CsrA_MDgen_freeEnergyMatrix = [ {'A' : -2.5, 'T' : 0, 'C' : -0.42, 'G' : -0.57}, #1st nucleotide position
                            {'A' : -2.34, 'T' : -0.64, 'C' : 0, 'G' : -0.75}, #2nd nucleotide position
                            {'G' : -3.31, 'T' : 0, 'C' : -0.13, 'A' : -1.36}, #3rd nucleotide position
                            {'G' : -2.83, 'T' : -1.26, 'C' : 0, 'A' : -0.54}, #4th nucleotide position
                            {'A' : 0.0, 'T' : -1.64, 'C' : -1.65, 'G' : -1.78}, #5th nucleotide position
                        ]

RsmA_MD_freeEnergyMatrix =  [ {'A' : -2.48, 'T' : -1.07, 'C' : -0.87, 'G' : 0}, #1st nucleotide position
                            {'A' : -3.46, 'T' : -2.26, 'C' : -1.15, 'G' : 0}, #2nd nucleotide position
                            {'G' : -5.63, 'T' : 0, 'C' : -0.31, 'A' : -1.29}, #3rd nucleotide position
                            {'G' : -5.81, 'T' : 0, 'C' : -0.95, 'A' : -1.96}, #4th nucleotide position
                            {'A' : -3.32, 'T' : 0, 'C' : -0.84, 'G' : -1.60}, #5th nucleotide position
                        ]

CsrA_siteFoldingConstraint = "xxxxx" #added additional nucleotide to see how this changes binding site predictions 

MAX_PRE_CUTOFF = 400
MAX_POST_CUTOFF = 400
    
def CsrA_dimerizationEnergyFunction(dist):

    #Calculate free energy of dimerization between two CsrA units,
    #depending on the distance between their binding sites
    #Alex note: what is the purpose of this fxn, if it only returns
    #a value of 0 and doesnt offer any wanrnings?
    
    MINIMUM_DISTANCE = 5        #closest possible distance between binding sites before clashing
    MININUM_FREE_ENERGY = -5.0  #minimum possible dG_dimerization
    STERIC_PENALTY = 20.0       #positive free energy of steric clashing between 2 CsrA units
    
    c0 = MININUM_FREE_ENERGY 
    c1 = 0.05   #tunable  => currently set so that dG == 0 @ dist = 50 nt
    c2 = 0.001  #tunable  => currently set so that dG == 0 @ dist = 50 nt
    
    if dist >= MINIMUM_DISTANCE:
        dG_dimerization = c2*(dist^2) + c1*dist + c0
    else:
        dG_dimerization = STERIC_PENALTY
        
    return dG_dimerization #removed remove dG_dimerization

def calculateSortedBindingSitesbyEnergy(sequence, freeEnergyMatrix):

    #Walk across sequence. 
    #Calculate binding free energy.
    #Calculate dG_dimerization with neighboring sites. 
    #Return list of CsrA binding sites (pairs) and free energies
    
    siteLen = len(freeEnergyMatrix)
    #print(siteLen)
    seqLen = len(sequence)
    #print(seqLen)

    singleBindingSiteList = []   #list of tuples, (beginning position, dG)
    
    for n in range(seqLen - siteLen):
        
        #Read subsequence & calculate energy
        dG_sum = 0.0
        for k in range(siteLen):
            
            letter = sequence[n + k].upper()
            try:
                dG_sum += freeEnergyMatrix[k][letter]
            except:
                print("Invalid letter %s at position %s in sequence %s" % (letter, n+k, sequence))
                pass
        #print("this is the dG sum: " + str(dG_sum))
        if dG_sum < 0: singleBindingSiteList.append(  (n, dG_sum) )
    
    numSites = len(singleBindingSiteList)
    #print(numSites)
    min_dG_single = min([dG for (pos, dG) in singleBindingSiteList])
    #print(min_dG_single, dG_sum)
    #print(singleBindingSiteList)
    
    doubleBindingSiteList = []  #list of tuples, (begin pos1, begin pos2, dG12, dG1, dG2, dG_dimerization)
    for combinations in itertools.combinations(range(numSites),2):
        #for k in range(numSites):
        #if not (n == k):
        site1 = singleBindingSiteList[combinations[0]]
        site2 = singleBindingSiteList[combinations[1]]
        #print(site1,site2)
        if site1[0] < site2[0]:
            (pos1, dG1) = site1
            (pos2, dG2) = site2
        else: 
            (pos2, dG2) = site1
            (pos1, dG1) = site2

        dist = pos2 - pos1 - siteLen #this is generating negative values when the two positions arent in descending order
        dG_dimerization = CsrA_dimerizationEnergyFunction(dist)

        dG12 = dG1 + dG2 + dG_dimerization 
        
        #Stop adding site pairs if the dG is not less than the best single site's dG.
        if dG12 < min_dG_single: 
            doubleBindingSiteList.append((pos1, pos2, dG12, dG1, dG2, dG_dimerization))
    #print(doubleBindingSiteList[0:10])    
    #Sort lists
    sortedSingleBindingSiteList = sorted(singleBindingSiteList, key = lambda x: x[1])  #ascending order
    #print(sortedSingleBindingSiteList)
    sortedDoubleBindingSiteList = sorted(doubleBindingSiteList, key = lambda x: x[2])  #ascending order
    #print(sortedDoubleBindingSiteList)
    
    return (sortedSingleBindingSiteList, sortedDoubleBindingSiteList)

def calc_dG_mRNA(sequence, constraints):
    
    RNA_model = "rna37"
    dangles = "none"
    temp = "37.0"
    RT = 0.616

    #if constraints is None:
    #    constraints = "." * len(sequence)

    #Abby added note: Vienna RNA functions here. then this calc_dG_mRNA function defined here is used later    
    fold = RNAEnergyModel([sequence], RNA_model)
    fold.mfe([1], constraints, Temp = temp, dangles = dangles)
    
    bp_x = fold["mfe_basepairing_x"][0]
    #print(bp_x)
    bp_y = fold["mfe_basepairing_y"][0]
    #print(bp_y)
    dG = fold["mfe_energy"][0] / RT  #now in units of RT
    structure = fold.convert_numbered_pairs_to_bracket([len(sequence)],bp_x,bp_y)
    #rint(structure)
    
    output = {'bp_x' : bp_x, 'bp_y' : bp_y, 'dG' : dG, 'structure' : structure}        
    return output

def calculateSingleBindingSite(inputList):
    
    seqLen = len(inputList[0]) #sequence
    #print("this is the site length: " + str(siteLen))
    sequence = inputList[0]
    siteLen = len(inputList[1]) #sitefoldingconstraint
    siteFoldingConstraint = inputList[1]
    (begin_pos, dG_protein) = inputList[2] #site, site energy 
    #print(begin_pos,dG_protein)
    
    cutoffs = [min(MAX_PRE_CUTOFF, begin_pos), min(MAX_POST_CUTOFF, seqLen - begin_pos - siteLen)]
    #print("cutoffs are : " + str(cutoffs))
    subSequence = sequence[begin_pos - cutoffs[0] : begin_pos + siteLen + cutoffs[1]]
    #print(subSequence)
    foldingConstraint = cutoffs[0] * "." + siteFoldingConstraint + cutoffs[1] * "."
    #print(foldingConstraint)

    #print "Calculating RNA structures in free & bound states for single CsrA binding sites"
    #these call upon vienna to determine the mfe of the subsequence
    RNA_ref = calc_dG_mRNA(subSequence, constraints = None)
    RNA_bound = calc_dG_mRNA(subSequence, constraints = foldingConstraint)
    
    RNA_ref['structure'] = (begin_pos - cutoffs[0]) * " " + RNA_ref['structure'] + (seqLen - begin_pos - siteLen - cutoffs[1]) * " "        #what does the * sign mean? put that many in a row?
    RNA_bound['structure'] = (begin_pos - cutoffs[0]) * " " + RNA_bound['structure'] + (seqLen - begin_pos - siteLen - cutoffs[1]) * " "
    RNA_bound['foldingConstraint'] = (begin_pos - cutoffs[0]) * "." + foldingConstraint + (seqLen - begin_pos - siteLen - cutoffs[1]) * "."
    
    #Check whether the site motif is located inside a hairpin loop or an internal loop
    #This will be true when base pairings exist from the left to right side of the motif
    bp_x = RNA_bound['bp_x']
    bp_y = RNA_bound['bp_y']
    insideLoop = False

    for (nt_x, nt_y) in zip(bp_x, bp_y):
        if nt_x < begin_pos and nt_y > (begin_pos + siteLen):
            insideLoop = True
            break
    
    if insideLoop:
        #when binding motif exists inside a loop, there is no binding penalty
        dG_insideLoop = 0.0 
    else:
        #when binding motif does not exist inside a loop, there is a binding penalty 
        #dG_insideLoop = 4.0 #kcal/mol #what is the binding penalty? Total guess at this point.
        dG_insideLoop = 0.0 #no binding penalty in this model variant 
    
    dG_ref = RNA_ref['dG']
    #print(dG_ref)
    dG_bound = RNA_bound['dG']
    #print(dG_bound)

    ddG_RNA =  dG_bound - dG_ref
    
    dG_total = dG_protein + ddG_RNA + dG_insideLoop
    
    return (begin_pos, dG_total, dG_protein, ddG_RNA, RNA_ref, RNA_bound, cutoffs)

def calculateDoubleBindingSite(inputList):
        
    seqLen = len(inputList[0]) 
    sequence = inputList[0]
    siteLen = len(inputList[1])
    siteFoldingConstraint = inputList[1]
    (pos1, pos2, dG12, dG1, dG2, dG_dimerization) = inputList[2]
    
    #calculate distance between binding regions
    dist = pos2 - pos1 - siteLen
    #print("this is the distance: " + str(dist))
    cutoffs = [min(MAX_PRE_CUTOFF, pos1), min(MAX_POST_CUTOFF, seqLen - pos2 - siteLen)]
    subSequence = sequence[pos1 - cutoffs[0] : pos2 + siteLen + cutoffs[1]]
    
    #pre_cutoff -- site1 -- spacer -- site 2 -- post_cutoff
    doubleBindingConstraint = cutoffs[0] * "." + siteFoldingConstraint + dist * "." + siteFoldingConstraint + cutoffs[1] * "."  #cheap workaround here 
    foldingConstraint = doubleBindingConstraint[0:len(subSequence)]
    #print(foldingConstraint)
    
    print("Calculating RNA structures in free & bound states for double CsrA binding sites")
    
    RNA_ref = calc_dG_mRNA(subSequence, constraints = None)
    RNA_bound = calc_dG_mRNA(subSequence, constraints = foldingConstraint)
    
    RNA_ref['structure'] = (pos1 - cutoffs[0]) * " " + RNA_ref['structure'] + (seqLen - pos2 - siteLen - cutoffs[1]) * " "
    RNA_bound['structure'] = (pos1 - cutoffs[0]) * " " + RNA_bound['structure'] + (seqLen - pos2 - siteLen - cutoffs[1]) * " "
    RNA_bound['foldingConstraint'] = (pos1 - cutoffs[0]) * "." + foldingConstraint + (seqLen - pos2 - siteLen - cutoffs[1]) * "."
    
    dG_ref = RNA_ref['dG']
    dG_bound = RNA_bound['dG']
    
    ddG_RNA =  dG_bound - dG_ref
    dG_total = dG12 + ddG_RNA
    
    return (pos1, pos2, dG_total, dG1, dG2, dG_dimerization, ddG_RNA, RNA_ref, RNA_bound, cutoffs)
        
def calculateFoldingFreeEnergies(use_MPI, sequence, siteFoldingConstraint, sortedSingleBindingSiteList, sortedDoubleBindingSiteList):

    #Iterate through single and double binding sites.
    #Calculate free energy penalty for RNA shape-change to make binding site accessible

    singleBindingSiteList_RNAFolding = []       #list of tuples, (begin pos1, dG_total, dG1, ddG_RNA)
    inputList = [[sequence, siteFoldingConstraint, site] for site in sortedSingleBindingSiteList]
    print("the size of the inputList is " + str(sys.getsizeof(inputList)))
    singleBindingSiteList_RNAFolding = map(calculateSingleBindingSite, inputList)
    
    doubleBindingSiteList_RNAFolding = []       #list of tuples, (begin pos1, begin pos2, dG_total, dG1, dG2, dG_dimerization, ddG_RNA)
    inputList = [[sequence, siteFoldingConstraint, sitepair] for sitepair in sortedDoubleBindingSiteList]
    #print(inputList[0:10])
    doubleBindingSiteList_RNAFolding = map(calculateDoubleBindingSite, inputList)
    #print(doubleBindingSiteList_RNAFolding[0:10])

    #Sort lists
    sortedSingleBindingSiteList_RNAFolding = sorted(singleBindingSiteList_RNAFolding, key = lambda x: x[1])  #ascending order
    sortedDoubleBindingSiteList_RNAFolding = sorted(doubleBindingSiteList_RNAFolding, key = lambda x: x[2])  #ascending order
    
    return (sortedSingleBindingSiteList_RNAFolding, sortedDoubleBindingSiteList_RNAFolding)

def inputSequencesFromExcel(fileName, sheetName):

    #Reads Excel file and imports sequence information
    #Format:  1st column are gene names.  2nd column are sequences. 
    #Terminates at first empty row.
    
    sequences = {}
    startCodons = {}
    
    wb = xlrd.open_workbook(fileName)
    sheets = wb.sheets()
    for sheet in sheets:
        if sheet.name == sheetName:
        
            for row in range(sheet.nrows):
                gene_name = sheet.cell(row, 0).value
                sequence = sheet.cell(row, 1).value
                start_codons = sheet.cell(row, 2).value
                
                if start_codons is str:
                    start_codon_list = [int(x) for x in start_codons.split(" ")]
                else:
                    start_codon_list = [int(start_codons)]
                
                sequences[gene_name] = sequence
                startCodons[gene_name] = start_codon_list
    return (sequences, startCodons)

def predictTranslationRatesSingleSites(inputList):
    sequence = inputList[0]
    start_codon_list = inputList[1]
    CsrA_begin_pos = inputList[2]
    dG_total = inputList[3]
    dG_protein = inputList[4]
    ddG_RNA = inputList[5]
    RNA_ref = inputList[6]
    RNA_bound = inputList[7]
    cutoffs = inputList[8]
    print(sequence, start_codon_list, CsrA_begin_pos, dG_total, dG_protein, ddG_RNA, RNA_ref, RNA_bound, cutoffs)
    #(sequence, start_codon_list, CsrA_begin_pos, dG_total, dG_protein, ddG_RNA, RNA_ref, RNA_bound, cutoffs) = inputs
    foldingConstraint = RNA_bound['foldingConstraint']
    siteLen = len(CsrA_siteFoldingConstraint)

    #print "single site foldingConstraint: ", foldingConstraint
    outputList = []
    for start in start_codon_list:
        calc = RBS_Calculator(sequence, start_range = [start-1, start], rRNA = 'ACCTCCTTA', constraints = foldingConstraint, verbose = False)
        calc.calc_dG()
        #outputData = calc.output()
        #RBS = calc.RBS_list[0]        
        start_pos = calc.start_pos_list[0]
        most_5p_SD_nt_pos = calc.most_5p_mRNA_list[0] #RBS.most_5p_paired_SD_mRNA
        #start_pos = RBS.start_position
        footprint_pos = start_pos + 15
        tir_OFF = calc.Expression_list[0]
        
        if CsrA_begin_pos >= most_5p_SD_nt_pos and (CsrA_begin_pos + siteLen) <= footprint_pos:
            CsrA_steric_repression = math.exp(-beta * dG_protein)
            tir_ON = tir_OFF / CsrA_steric_repression
        else:
            tir_ON = tir_OFF
            CsrA_steric_repression = 1.0
        
        print("Pos %s: TIR changed to %s (steric repression was %s)" % (str(start_pos), str(round(tir_ON,2)), str(round(CsrA_steric_repression,2))))
        
        outputList.append((start_pos, tir_ON, tir_OFF / tir_ON))
    return outputList
        
def predictTranslationRatesDoubleSites(inputList):
    sequence = inputList[0]
    start_codon_list = inputList[1]
    CsrA_site1_begin_pos = inputList[2]
    CsrA_site2_begin_pos = inputList[3]
    dG_total = inputList[4]
    dG1 = inputList[5]
    dG2 = inputList[6]
    dG_dimerization = inputList[7]
    ddG_RNA = inputList[8]
    RNA_ref = inputList[9]
    RNA_bound = inputList[10]
    cutoffs = inputList[11]   
    #(sequence, start_codon_list, CsrA_site1_begin_pos, CsrA_site2_begin_pos, dG_total, dG1, dG2, dG_dimerization, ddG_RNA, RNA_ref, RNA_bound, cutoffs) = inputs
    foldingConstraint = RNA_bound['foldingConstraint']
    siteLen = len(CsrA_siteFoldingConstraint)
    
    #print "double site foldingConstraint: ", foldingConstraint
        
    outputList = []
    for start in start_codon_list:
        calc = RBS_Calculator(sequence, start_range = [start-1, start], rRNA = 'ACCTCCTTA', constraints = foldingConstraint, verbose = False)
        calc.calc_dG()

        #outputData = calc.output()
        #RBS = calc.RBS_list[0]        
        start_pos = calc.start_pos_list[0]
        most_5p_SD_nt_pos = calc.most_5p_mRNA_list[0] #RBS.most_5p_paired_SD_mRNA
        #start_pos = RBS.start_position
        footprint_pos = start_pos + 15
        tir_OFF = calc.Expression_list[0]
        
        if CsrA_site1_begin_pos >= most_5p_SD_nt_pos and (CsrA_site2_begin_pos + siteLen) <= footprint_pos:
            CsrA_steric_repression = math.exp( -beta * (dG1 + dG2 + dG_dimerization) )
            tir_ON = tir_OFF / CsrA_steric_repression
        elif CsrA_site1_begin_pos >= most_5p_SD_nt_pos and (CsrA_site1_begin_pos + siteLen) <= footprint_pos:
            CsrA_steric_repression = math.exp( -beta * (dG1 + dG_dimerization) )
            tir_ON = tir_OFF / CsrA_steric_repression
        elif CsrA_site2_begin_pos >= most_5p_SD_nt_pos and (CsrA_site2_begin_pos + siteLen) <= footprint_pos:
            CsrA_steric_repression = math.exp( -beta * (dG2 + dG_dimerization) )
            tir_ON = tir_OFF / CsrA_steric_repression
        else:
            tir_ON = tir_OFF
            CsrA_steric_repression = 1.0
        
        print("Pos %s: TIR changed to %s (steric repression was %s)" % (str(start_pos), str(round(tir_ON,2)), str(round(CsrA_steric_repression,2))))
        outputList.append( (start_pos, tir_ON, tir_OFF / tir_ON) )
    return outputList

def predictTranslationRates(sequence, start_codon_list, singleSites, doubleSites):

    siteLen = len(CsrA_siteFoldingConstraint)
    
    print("Starting Translation Rate Calculations")
    
    #Calculate mRNA translation rates in the absence of CsrA binding
    TranslationRates_free = []
    for start in start_codon_list:
        calc = RBS_Calculator(sequence, start_range = [start-1, start], rRNA = 'ACCTCCTTA', constraints = None, verbose = False)
        calc.calc_dG()      
        RBS = calc.start_pos_list[0]
        #most_5p_SD_nt_pos = calc.aligned_most_5p_SDmost_5p_mRNA_list[0] #RBS.most_5p_paired_SD_mRNA
        #start_pos = RBS.start_position
        #footprint_pos = start_pos + 15
        tir = calc.Expression_list[0]
        TranslationRates_free.append( (tir, 1.0) )
        #print "TIR %s at pos %s" % (str(round(RBS.tir,2)), str(pos))
    
    inputList = [ [sequence, start_codon_list] + list(inputs) for inputs in singleSites][0:15]
    print("Number of Single Site translation rate predictions: ", len(inputList))
    if use_MPI:
        print("Using MPI")
        translationRatesSingleSites = pool.map(predictTranslationRatesSingleSites, inputList)
    else:
        translationRatesSingleSites = map(predictTranslationRatesSingleSites, inputList) #list of lists
        translationRatesSingleSites = list(translationRatesSingleSites)
    
    inputList = [ [sequence, start_codon_list] + list(inputs) for inputs in doubleSites][0:15]
    print("Number of Double Site translation rate predictions: ", len(inputList))
    if use_MPI:
        translationRatesDoubleSites = pool.map(predictTranslationRatesDoubleSites, inputList) 
    else:
        translationRatesDoubleSites = map(predictTranslationRatesDoubleSites, inputList) #list of lists
        translationRatesDoubleSites = list(translationRatesDoubleSites)


    return (sequence, TranslationRates_free, translationRatesSingleSites, translationRatesDoubleSites)

def exportCalculationsToExcel(filename, calculationsDict):

    #singleSites
    #(begin_pos, dG_total, dG_protein, ddG_RNA, RNA_ref, RNA_bound, cutoffs)
    
    #doubleSites
    #(pos1, pos2, dG_total, dG1, dG2, dG_dimerization, ddG_RNA, RNA_ref, RNA_bound, cutoffs)
    
    workbook = xlwt.Workbook()
    
    courier_new = xlwt.easyxf('font: name Courier New, height 200;')
    
    summary = workbook.add_sheet("Summary")
    
    summary.write(0,0, 'Gene Name')
    summary.write(0,1, 'Start Codon')
    summary.write(0,2, 'Basal Translation Rate, no CsrA (au)')
    summary.write(0,3, 'Translation Rate, single bound CsrA (au)')
    summary.write(0,4, 'Translation Rate, double bound CsrA (au)')        
    summary.write(0,5, 'Single CsrA dG_total (RT)')
    summary.write(0,6, 'Single CsrA dG_protein (RT)')
    summary.write(0,7, 'Single CsrA ddG_RNA (RT)')
    summary.write(0,8, 'Single CsrA Binding Sequence')
    summary.write(0,9, 'Single CsrA Binding Structure')
    summary.write(0,10, 'Single CsrA Binding site ')
    summary.write(0,11, 'Double CsrA dG_total (RT)')
    summary.write(0,12, 'Double CsrA dG_protein (1) (RT)')
    summary.write(0,13, 'Double CsrA dG_protein (2) (RT)')
    summary.write(0,14, 'Double CsrA dG_dimerization (RT)')
    summary.write(0,15, 'Double CsrA ddG_RNA (RT)')
    summary.write(0,16, 'Double CsrA Binding Sequence')
    summary.write(0,17, 'Double CsrA Binding Structure')
    summary.write(0,18, 'Double CsrA Site 1')
    summary.write(0,19, 'Double CsrA Site 2')
    
    rowCounter = 0
    
    for (gene, (sequence, singleSites, doubleSites) ) in calculationsDict['sites'].items():
        try:
            rowCounter += 2
            summary.write(rowCounter-1,0, gene)


            #These are all the RBS calculator specific summaries, add back when this portion is fixed.  
            (sequence, start_codon_list, TranslationRates_free, TranslationRatesSingleSites, TranslationRatesDoubleSites) = calculationsDict['translation'][gene]
         
            start = start_codon_list[0]
            (TIR, regulation) = TranslationRates_free[0]
            summary.write(rowCounter-1,1, start)
            summary.write(rowCounter-1,2, TIR)
        
            (start_pos, tir_ON, ratio) = TranslationRatesSingleSites[0][0]
            summary.write(rowCounter-1,3, tir_ON)
        
            (start_pos, tir_ON, ratio) = TranslationRatesDoubleSites[0][0]
            summary.write(rowCounter-1,4, tir_ON)
        
        
            #Add highest affinity site(s) to summary sheet
            if len(singleSites) > 0:
                (begin_pos, dG_total, dG_protein, ddG_RNA, RNA_ref, RNA_bound, cutoffs) = singleSites[0]
                structure = RNA_bound['structure']
        
                sequence1 = sequence[0:begin_pos] + "*" + sequence[begin_pos:]
                structure = structure[0:begin_pos] + "*" + structure[begin_pos:]
        
                summary.write(rowCounter-1,5, dG_total)
                summary.write(rowCounter-1,6, dG_protein)
                summary.write(rowCounter-1,7, ddG_RNA)
                summary.write(rowCounter-1,8, sequence1, courier_new)
                summary.write(rowCounter-1,9, structure, courier_new)
                summary.write(rowCounter-1,10, begin_pos,courier_new)
        
            if len(doubleSites) > 0:
                (pos1, pos2, dG_total, dG1, dG2, dG_dimerization, ddG_RNA, RNA_ref, RNA_bound, cutoffs) = doubleSites[0]
                structure = RNA_bound['structure']
            
                sequence2 = sequence[0:pos1] + "*" + sequence[pos1:pos2] + "*" + sequence[pos2:]
                structure = structure[0:pos1] + "*" + structure[pos1:pos2] + "*" + structure[pos2:]
            
                summary.write(rowCounter-1,11, dG_total)
                summary.write(rowCounter-1,12, dG1)
                summary.write(rowCounter-1,13, dG2)
                summary.write(rowCounter-1,14, dG_dimerization)
                summary.write(rowCounter-1,15, ddG_RNA)
                summary.write(rowCounter-1,16, sequence2, courier_new)
                summary.write(rowCounter-1,17, structure, courier_new)
                summary.write(rowCounter-1,18, pos1,courier_new)
                summary.write(rowCounter-1,19, pos2,courier_new)
        except:
            print( "ERROR! Skipping %s gene" % gene)
            print(traceback.format_exc())

    for (gene, (sequence, singleSites, doubleSites) ) in calculationsDict['sites'].items():
        try:
            sheet = workbook.add_sheet(gene)
        
            (sequence, start_codon_list, TranslationRates_free, TranslationRatesSingleSites, TranslationRatesDoubleSites) = calculationsDict['translation'][gene]
        
            sheet.write(0,0, 'List of CsrA Binding Configurations for %s mRNA' % gene)
            sheet.write(1,0, 'Sequence:   ')
            sheet.write(1,4, sequence, courier_new)
            sheet.write(1,10, sequence, courier_new)
        
            sheet.write(2,0, '#')
            sheet.write(2,1, 'Start Codon (nt pos)')
            sheet.write(2,2, 'Basal Translation Rate, no CsrA (au) ')
            sheet.write(2,3, 'Translation Rate, single bound CsrA (au)')
            sheet.write(2,4, 'Translation Rate, double bound CsrA (au)')    
            sheet.write(2,5, 'Single CsrA dG_total (RT)')
            sheet.write(2,6, 'Single CsrA dG_protein (RT)')
            sheet.write(2,7, 'Single CsrA ddG_RNA (RT)')
            sheet.write(2,8, 'Single CsrA Binding Structure')
            sheet.write(2,9, 'Single CsrA Binding Site')
            sheet.write(2,10, 'Double CsrA dG_total (RT)')
            sheet.write(2,11, 'Double CsrA dG_protein (1) (RT)')
            sheet.write(2,12, 'Double CsrA dG_protein (2) (RT)')
            sheet.write(2,13, 'Double CsrA dG_dimerization (RT)')
            sheet.write(2,14, 'Double CsrA ddG_RNA (RT)')
            sheet.write(2,15, 'Double CsrA Binding Structure')
            sheet.write(2,16, 'Double CsrA Binding Site 1')
            sheet.write(2,17, 'Double CsrA Binding Site 2')

        
        
            for (i, start) in enumerate(start_codon_list):
                (TIR_OFF, regulation) = TranslationRates_free[i]
        
                rowCounter = 2
                for predictions in TranslationRatesSingleSites:
                    rowCounter += 1
                    (start_pos, tir_ON, ratio) = predictions[i]
                    sheet.write(rowCounter,1, start)
                    sheet.write(rowCounter,2, TIR_OFF)
                    sheet.write(rowCounter,3, tir_ON)
                
                rowCounter = 2
                for predictions in TranslationRatesDoubleSites:
                    rowCounter += 1
                    (start_pos, tir_ON, ratio) = predictions[i]
                    sheet.write(rowCounter,4, tir_ON)  
        
            rowCounter = 2
            for (begin_pos, dG_total, dG_protein, ddG_RNA, RNA_ref, RNA_bound, cutoffs) in singleSites:
                structure = RNA_bound['structure']
                structure = structure[0:begin_pos] + "*" + structure[begin_pos:]
            
                rowCounter += 1
                sheet.write(rowCounter,5, dG_total)
                sheet.write(rowCounter,6, dG_protein)
                sheet.write(rowCounter,7, ddG_RNA)
                sheet.write(rowCounter,8, structure, courier_new)
                sheet.write(rowCounter,9, begin_pos, courier_new)

        
            rowCounter = 2
            for (pos1, pos2, dG_total, dG1, dG2, dG_dimerization, ddG_RNA, RNA_ref, RNA_bound, cutoffs) in doubleSites:
                structure = RNA_bound['structure']
                structure = structure[0:pos1] + "*" + structure[pos1:pos2] + "*" + structure[pos2:]
            
                rowCounter += 1
                sheet.write(rowCounter,10, dG_total)
                sheet.write(rowCounter,11, dG1)
                sheet.write(rowCounter,12, dG2)
                sheet.write(rowCounter,13, dG_dimerization)
                sheet.write(rowCounter,14, ddG_RNA)
                sheet.write(rowCounter,15, structure, courier_new)
                sheet.write(rowCounter,16, pos1, courier_new)
                sheet.write(rowCounter,17, pos2, courier_new)

        except:
            print("ERROR! Skipping %s gene" % gene)
            print(traceback.format_exc())
        
    workbook.save(filename)

def split_merge_binding_dict(gene,dictionary,binding_sites):
    if binding_sites == 1:
        site = []
        dG_total_ss = []
        dG_bs1_ss = []
        dG_bs2_ss = []
        ddG_mRNA = []
        struct_dict = []
        gene_name = []
        for i in range(len(dictionary)):
            site.append(dictionary[i][0])
            dG_total_ss.append(dictionary[i][1])
            dG_bs1_ss.append(dictionary[i][2])
            ddG_mRNA.append(dictionary[i][3])
            struct_dict.append(dictionary[i][4])
        struct_dict_df = pd.DataFrame(struct_dict)
        gene_df = pd.DataFrame()
        gene_name = []
        for n in range(len(site)):
            gene_name.append(str(gene))
        print(gene_name)
        gene_df['gene'] = gene_name
        print(gene)
        gene_df['site'] = site
        gene_df['dG_total_ss'] = dG_total_ss
        gene_df['dG_bs1_ss'] = dG_bs1_ss
        gene_df['ddG_mRNA'] = ddG_mRNA
        merged_gene = pd.concat([gene_df, struct_dict_df], axis = 1)
    if binding_sites == 2:
        site1 = []
        site2 = []
        dG_total_ds = []
        dG_bs1 = []
        dG_bs2 = []
        ddG_mRNA = []
        dG_coop = []
        struct_dict = []
        for i in range(len(dictionary)):
            site1.append(dictionary[i][0])
            site2.append(dictionary[i][1])
            dG_total_ds.append(dictionary[i][2])
            dG_bs1.append(dictionary[i][3])
            dG_bs2.append(dictionary[i][4])
            ddG_mRNA.append(dictionary[i][5])
            dG_coop.append(dictionary[i][6])
            struct_dict.append(dictionary[i][7])
        struct_dict_df = pd.DataFrame(struct_dict)
        gene_df = pd.DataFrame()
        gene_df['site1'] = site1
        gene_df['site2'] = site2
        gene_df['dG_total_ds'] = dG_total_ds
        gene_df['dG_bs1'] = dG_bs1
        gene_df['dG_bs2'] = dG_bs2
        gene_df['ddG_mRNA'] = ddG_mRNA
        gene_df['dG_coop'] = dG_coop
        merged_gene = pd.concat([gene_df, struct_dict_df], axis = 1)
    return merged_gene


if __name__ == "__main__":

    print(80*"*")
    print("    .----.   @   @ biophysical")
    print('   / .-"-.`.  \ /             ')
    print("  | | ''\ \ \_/ )       model ")
    print(" ,-\ `-.' /.'  /              ")
    print(" ---`----'----'        AL,2020")
    print(80*"*")

    parser = argparse.ArgumentParser(description='import csv with constraints for folding')
    parser.add_argument('-i',  
        action='store', 
        dest='i',
        required=True,
        type=str,
        help="input using '-i' *.xls file")

    parser.add_argument('-o',
        action = 'store',
        dest = 'o',
        required = True,
        type = str,
        help = "prefix for outfile summary")
    
    parser.add_argument('-pwm',
        action = 'store',
        dest = 'pwm',
        required = False,
        default = "CsrA",
        type = str,
        help = "Position weight matrix used in model. Options are CsrA, CsrA_MD, or RsmA_MD")

    parser.add_argument('-export',
        action = 'store',
        dest = 'exp',
        required = False,
        default = 'excel',
        type = str,
        help = "define export type, can export to excel or to csv")

    options = parser.parse_args()

    if options.pwm == "CsrA":
        energy_matrix = CsrA_freeEnergyMatrix
    if options.pwm == "CsrA_MD":
        energy_matrix = CsrA_MDgen_freeEnergyMatrix
    if options.pwm == "RsmA_MD":
        energy_matrix = RsmA_MD_freeEnergyMatrix

    #import excel file 

    filename = options.i 

    outputFilename = '../data/' + options.o + ".xls"
    print(filename, outputFilename)
    sheetName = 'Sheet1'
    
    doRandom = False
    
    if doRandom:
        sequenceDict = {}
        seqlen = 200
        for n in range(100):
            sequenceDict['random%s' % n] = "".join( [random.choice( ('A','G','C','T') ) for x in range(seqlen)] )
    else:
        (sequenceDict, startCodonDict) = inputSequencesFromExcel(filename, sheetName)
        print(sequenceDict, startCodonDict)
    
    #parallel processing functions
    #try:
    #    from mpi4py import MPI
    #    comm = MPI.COMM_WORLD
    #    print("Number of nodes is %s" % comm.Get_size())
    #    if comm.Get_size() > 1:
    #        use_MPI = True
    #        from MPI_pool import Pool
    #        pool = Pool(MPI.COMM_WORLD)
    #        pool.start()
    #    else:
    #        use_MPI = False
    #except:
    #    print("Could not start MPI Pool. Stopping now.")
    #    print(traceback.format_exc())
    #    use_MPI = False
    
    #if (use_MPI and pool.is_master()):
    calculationsDict = {}
    calculationsDict['sites'] = {}
    calculationsDict['translation'] = {}
    
    use_MPI = False

    geneCounter = 0
    print("Number of sequences: %s" % len(sequenceDict.keys() ))
    for (gene, sequence) in sequenceDict.items():
        try:
            geneCounter += 1
            #print("GENE #%s: %s" % geneCounter, gene)
            sequence.upper()
            (sortedSingleBindingSiteList, sortedDoubleBindingSiteList) = calculateSortedBindingSitesbyEnergy(sequence, energy_matrix)
            print("%s SINGLE BINDING SITES (CsrA only)" % len(sortedSingleBindingSiteList))
            print("TOP 5: ", sortedSingleBindingSiteList[0:5])
            print("%s DOUBLE BINDING SITES (CsrA only)" % len(sortedDoubleBindingSiteList))
            print("TOP 5: ", sortedDoubleBindingSiteList[0:5])
            (sortedSingleBindingSiteList_RNAFolding, sortedDoubleBindingSiteList_RNAFolding) = calculateFoldingFreeEnergies(use_MPI, sequence, CsrA_siteFoldingConstraint, sortedSingleBindingSiteList, sortedDoubleBindingSiteList) 
            print("TOP 5: SINGLE BINDING SITES (CsrA + RNA shape change)")
            print([(begin_pos, dG_total, ddG_RNA, RNA_bound['structure']) for (begin_pos, dG_total, dG_protein, ddG_RNA, RNA_ref, RNA_bound, cutoffs) in sortedSingleBindingSiteList_RNAFolding[0:5]])
            print("TOP 5: DOUBLE BINDING SITES (CsrA + RNA shape change)")
            print([(pos1, pos2, dG_total, dG1, dG2, dG_dimerization, ddG_RNA, RNA_bound['structure']) for (pos1, pos2, dG_total, dG1, dG2, dG_dimerization, ddG_RNA, RNA_ref, RNA_bound, cutoffs) in sortedDoubleBindingSiteList_RNAFolding[0:5]])
            print(80 * "*")
            calculationsDict['sites'][gene] = (sequence, sortedSingleBindingSiteList_RNAFolding, sortedDoubleBindingSiteList_RNAFolding)
            print(calculationsDict)
            print('single binding site list ' + str(sortedSingleBindingSiteList_RNAFolding))
            print('double binding site list ' + str(sortedDoubleBindingSiteList_RNAFolding))

            print("Starting RBS Calculator")
            start_codon_list = startCodonDict[gene]
            print(start_codon_list)
            (sequence, TranslationRates_free, TranslationRatesSingleSites, TranslationRatesDoubleSites) = predictTranslationRates(sequence, start_codon_list, sortedSingleBindingSiteList_RNAFolding, sortedDoubleBindingSiteList_RNAFolding)
            calculationsDict['translation'][gene] = (sequence, start_codon_list, TranslationRates_free, TranslationRatesSingleSites, TranslationRatesDoubleSites)
            print(calculationsDict)
            
        except:
            print("ERROR. Skipping %s gene." % gene)
            print (traceback.format_exc())

    handle = open(outputFilename + '.pkl','wb')
    pickle.dump(calculationsDict, handle)
    handle.close()

    #export options----------------------------------------------------------
    if options.exp == 'excel':
        exportCalculationsToExcel(outputFilename, calculationsDict)

    if options.exp == 'csv':
        binding_sites = calculationsDict['sites']
        translation_rates = calculationsDict['translation']
        all_gene_list = []

        #process binding site data 
        for gene in binding_sites:
            #process initial weird formatting in dictionary
            sequence = binding_sites[gene][0]
            single_site_dict = binding_sites[gene][1]
            double_site_dict = binding_sites[gene][2]
            #split dictionaries and create dataframe for single sites 
            test_df = split_merge_binding_dict(gene,single_site_dict,1)
            test_df2 = split_merge_binding_dict(gene,double_site_dict,2)
            merged_df = pd.concat([test_df,test_df2], axis = 1)
            all_gene_list.append(merged_df)
        all_gene_df = pd.concat(all_gene_list)
        all_gene_df.to_csv("../data/binding_sites_" + options.o + '.csv')

        #process translation rate data 
        translation_rate_list = []
        for gene in translation_rates:
            translation_df = pd.DataFrame()
            gene_name = []
            for n in range(len(translation_rates[gene][4])):
                gene_name.append(str(gene))
            translation_df['gene'] = gene_name
            translation_df['sequence'] = translation_rates[gene][0]
            translation_df['start_codon'] = translation_rates[gene][1][0]
            TIR_off = []
            TIR_on_ss = []
            TIR_on_ds = []
            for i in range(len(translation_rates[gene][3])):
                TIR_on_ss.append(translation_rates[gene][3][i][0][1])
                TIR_on_ds.append(translation_rates[gene][4][i][0][1])
            translation_df['TIR_off'] = translation_rates[gene][2][0][0]
            translation_df['TIR_on_ss'] = TIR_on_ss
            translation_df['TIR_on_ds'] = TIR_on_ds
            translation_rate_list.append(translation_df)
        all_translation_df = pd.concat(translation_rate_list)
        all_translation_df.to_csv("../data/translation_rates_" + options.o + '.csv')

    #if use_MPI and pool.is_master(): pool.stop()
