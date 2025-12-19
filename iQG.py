import os
import pandas as pd
import re
import itertools
import numpy as np
from numpy import log as ln
import math
from functools import reduce
from statistics import mean
from ast import literal_eval

'''
pyVCF requires manual installation through the python interpreter, 
uncomment the import vcf line below if original vcf files still need to be converted
'''
# import vcf


#########################################**** VCF Conversions ****#####################################################


def splitINFO(file): # file would have to be a dataframe
    csq = []
    for row in range(file.shape[0]):
        info1 = str(file["info"][row])
        start = info1.find("CSQ")
        end = info1.find("SOMATIC")
        if end > start:
            info = info1[start+6:end-3]
            info = info.replace("'", "")
            info = info.replace("[", "")
            info = info.replace("]", "")
            info = info.replace(" ", "")
            info = info.split(",")
        else:
            info = info1[start+6:-2]
            info = info.replace("'", "")
            info = info.replace("[", "")
            info = info.replace("]", "")
            info = info.replace(" ", "")
            info = info.split(",")

        csq.append(info)
    file.insert(file.shape[1], "CSQ", csq)
    return file



def parseCSQ (file, header): #file needs to be a dataframe, header needs to be based on the information within the vcf
    newColheader = header.split("|")

    # creates new column headers based on inputted list where the entries are empty list that can be appended to:
    for col in newColheader:
        file[col] = [list() for x in range(len(file.index))]


    for row in range(file.shape[0]):
        info1 = str(file["CSQ"][row]).replace("'", "").replace("[", "").replace("]", "").replace('"', "")
        info1 = info1.split(",")

        for ele in info1:
            count = 0
            ele = ele.split("|")

            while count < len(ele):
                file[newColheader[count]][row].append(ele[count])
                count += 1
    return file


def PatientID(vcffile):
    temp = vcffile.split(".")

    with open(vcffile, mode='r') as vcf:
        for line in vcf:

            if line.startswith('##INDIVIDUAL'):

                ind1 = line.find("ID=")
                ind2 = line.find("NAME=")

                if ind1 == -1:
                    print("line Problem with: " + temp[0])
                    return
                else:

                    ele = line[ind2+5:ind1-1] #Patient ID

    return ele

'''
metadataParsing Parameters: 
info > "yes" or "no" string entries only, user needs to identify whether or not the VCF Files contain the CSQ-INFO 

pLine > In the line where the patient ID is found, string entry needs to be what the line begins with (i.e. "##INDIVIDUAL", "##NORMAL", or "##SAMPLE")
pStart > Within the pLine, this string entry is to locate where the begining of the ID is in the line
pEnd > Within the pLine, this string entry is to locate where the end of the ID is in the line
'''


def metadataParsing(vcffile, info, pLine, pStart, pEnd):
    temp = vcffile.split(".")
    pID = ""
    header = ""

    with open(vcffile, mode='r') as vcf:
        for line in vcf:

            if line.startswith(pLine):

                ind1 = line.find(pStart)
                ind2 = line.find(pEnd)

                if ind1 == -1:
                    print("Patient ID error in: " + temp[0])
                    return
                else:
                    pID = line[ind1 + len(pStart):ind2 - 1]  # Patient ID
            if info == "yes" and line.startswith("##INFO=<ID=CSQ"):
                hStart = line.find("Format:")
                hEnd = line.find('">')
                if hStart == -1 or hEnd == -1:
                    print("CSQ header error in:" + temp[0])
                    return
                else:
                    header = line[hStart + len("Format:"):hEnd - 1]
        return pID, header

def addALT_Seq(df):
    df['alt_seq'] = np.where((df['ref_seq'] != df['var_seq1']) & (df['ref_seq'] == df['var_seq2']),
                         df['var_seq1'], df['var_seq2'])
    return df

def writecsv(inputfile, csqInfo):
    vcf_reader = vcf.Reader(open(inputfile,'r'))  # read vcf file
    #use nested list to store variant information, one for tumor, one for normal
    normal_var=[]
    tumor_var=[]
    var_score = [28]
    refct = []
    altct = []
    refct1 = []
    altct1 = []

    for record in vcf_reader:
        curr_var = []
        left = []
        right = []
        ref_seq = []
        alt_seq = []
        info = []
        # Because the AML and ALL have a reverse order for which comes first in the VCF, the sample[0/1] needed to be changed in order to make sure
            # the tumor and normal variants were saved to the correct csv. Can be compared to the original PrCa code
        tumorAD = record.samples[1]["AD"]
        normalAD = record.samples[0]["AD"]
        ########################################### tumor_var #############################################
        #######check for 2 ALT options:
        if len(tumorAD) == 3:
            sum = tumorAD[0] + tumorAD[1] + tumorAD[2]
            x = tumorAD[0]
            y = tumorAD[1]
            z = tumorAD[2]
            # if two out of the three are not greater than 0.05 throw away variant
            if ((x / sum) < 0.05 and (y / sum) < 0.05) or ((x / sum) < 0.05 and (z / sum) < 0.05) or (
                    (z / sum) < 0.05 and (y / sum) < 0.05):
                pass
            # means two are going to be good to decide which to use
            else:
                print(record)

        ########only 1 base for the ALT sequence
        else:
            sum = tumorAD[0] + tumorAD[1]
            refcount = tumorAD[0]
            altcount = tumorAD[1]

            # both the REF and ALT have to be greater
            if (refcount / sum) > 0.05 and (altcount / sum) > 0.05:
                curr_var.append(record.CHROM)
                left.append(record.POS)
                right.append(record.POS + len(record.ALT))
                ref_seq.append(record.REF)
                alt_seq.append(str((record.ALT[0])))
                refct.append(refcount)
                altct.append(altcount)
                info.append(str(record.INFO))
                lists = curr_var + left + right + ref_seq + alt_seq + var_score + info
                tumor_var.append(lists)


            # throw away variant if not
            else:
                pass

        ############################# normal_var ##############################3
        curr_var1 = []
        left1 = []
        right1 = []
        ref_seq1 = []
        alt_seq1 = []
        info = []
        if len(normalAD) == 3:
            sum = normalAD[0] + normalAD[1] + normalAD[2]
            x = normalAD[0]
            y = normalAD[1]
            z = normalAD[2]
            # if two out of the three are not greater than 0.05 throw away variant
            if ((x / sum) < 0.05 and (y / sum) < 0.05) or ((x / sum) < 0.05 and (z / sum) < 0.05) or (
                    (z / sum) < 0.05 and (y / sum) < 0.05):
                pass
            # means two are going to be good to decide which to use
            else:
                print(record)

        ########only 1 base for the ALT sequence
        else:
            sum = normalAD[0] + normalAD[1]
            refcount = normalAD[0]
            altcount = normalAD[1]

            # both the REF and ALT have to be greater
            if (refcount / sum) > 0.05 and (altcount / sum) > 0.05:
                curr_var1.append(record.CHROM)
                left1.append(record.POS)
                right1.append(record.POS + len(record.ALT))
                ref_seq1.append(record.REF)
                alt_seq1.append(str((record.ALT[0])))
                info.append(str(record.INFO))
                refct1.append(refcount)
                altct1.append(altcount)

                lists = curr_var1 + left1 + right1 + ref_seq1 + alt_seq1 + var_score + info
                normal_var.append(lists)


            # throw away variant if not
            else:
                pass


    #transfer list into dataframe, it would be easier for following manipulation
    normal=pd.DataFrame(normal_var,columns=['chrom','left','right','ref_seq','var_seq1','var_score', 'info'])
    tumor=pd.DataFrame(tumor_var,columns=['chrom','left','right','ref_seq','var_seq1','var_score', 'info'])

    normal_seq2 = []
    tumor_seq2 = []

    normal.insert(5, "count1", altct1)
    normal.insert(6, "count2", refct1)
    for i in range(normal.shape[0]):

        if normal['ref_seq'][i] > normal['var_seq1'][i]:
            normal_seq2.append(normal['ref_seq'][i])

        else:
            normal_seq2.append(normal['var_seq1'][i])
            normal['var_seq1'][i] = normal['ref_seq'][i]
            # switching count values because the ref order is now var_seq1
            temp = normal['count1'][i]
            normal['count1'][i] = normal['count2'][i]
            normal['count2'][i] = temp

    normal.insert(5, "var_seq2", normal_seq2)  # inserts column with information given

    tumor.insert(5, "count1", altct)
    tumor.insert(6, "count2", refct)
    for j in range(tumor.shape[0]):
        if tumor['ref_seq'][j] > tumor['var_seq1'][j]:
            tumor_seq2.append(tumor['ref_seq'][j])

        else:
            tumor_seq2.append(tumor['var_seq1'][j])
            tumor['var_seq1'][j] = tumor['ref_seq'][j]
            # switching count values because the ref order is now var_seq1
            temp = tumor['count1'][j]
            tumor['count1'][j] = tumor['count2'][j]
            tumor['count2'][j] = temp

    tumor.insert(5, "var_seq2", tumor_seq2)

    normal = addALT_Seq(normal)
    tumor = addALT_Seq(tumor)

    #inserts and sets the var_index using a list of number within the range of the rows
    tumor.insert(0, "var_index", list(range(tumor.shape[0])))
    normal.insert(0, "var_index", list(range(normal.shape[0])))



    ################################
    # Inserts patient ID and extracts CSQ header info, see comment above defined function for parameter descriptions
    # Alter parameters as needed, these specifically apply for TCGA ThCa and OCa
    pID, header = metadataParsing(inputfile, "yes", "##INDIVIDUAL=", "NAME=", "ID=")

    tumor.insert(tumor.shape[1], "Patient_ID", pID)
    normal.insert(normal.shape[1], "Patient_ID", pID)
    ################################



    ################################
    if csqInfo == "yes":
        tumor = splitINFO(tumor)
        tumor = parseCSQ(tumor, header)

        normal = splitINFO(normal)
        normal = parseCSQ(normal, header)

    # Remove the Info and CSQ Column, too long with them kept
        tumor = tumor.drop(columns=["info", "CSQ"])
        normal = normal.drop(columns=["info","CSQ"])

    else:
        tumor = tumor.drop(columns=["info"])
        normal = normal.drop(columns=["info"])
    ##################################



    ############### export the CSV file ###############
    # tumor.to_csv(id + "_tumor.csv")
    # normal.to_csv(id + "_normal.csv")
    return tumor, normal



def multi_vcf2csv(vcfPath, csqInfo, outPath, cancerType):
    os.chdir(vcfPath)
    files = os.listdir()
    Merge = pd.DataFrame()
    Tumor = pd.DataFrame()
    Normal = pd.DataFrame()

    for file in files:
        T, N = writecsv(file, csqInfo)
        Tumor = pd.concat([Tumor, T], axis=0)
        Normal = pd.concat([Normal, N], axis=0)
        Merge = pd.concat([Merge, T], axis=0)
        Merge = pd.concat([Merge, N], axis=0)


    # Remove the chrM variants, those are to be ignored
    Tumor = Tumor[Tumor["chrom"] != "chrM"]
    Normal = Normal[Normal["chrom"] != "chrM"]
    Merge = Merge[Merge["chrom"] != "chrM"]

    ## remove duplicates based on patient id
    Normal = Normal.drop_duplicates(subset=["chrom", "left", "ref_seq", "alt_seq", "Patient_ID"])
    Tumor = Tumor.drop_duplicates(subset=["chrom", "left", "ref_seq", "alt_seq", "Patient_ID"])
    Merge = Merge.drop_duplicates(subset=["chrom", "left", "ref_seq", "alt_seq"])

####### Mutational Frequency: Patient's Dataframe #################
    ## subset the variant info and patients info only
    ptN = Normal[["chrom", "left", "ref_seq", "alt_seq", "Patient_ID"]]
    ptT = Tumor[["chrom", "left", "ref_seq", "alt_seq", "Patient_ID"]]
    ptN = ptN.reset_index(drop=True)
    ptT = ptT.reset_index(drop=True)

    ## group by variant info and
    ptN = ptN.groupby(["chrom", "left", "ref_seq", "alt_seq"], as_index=False).agg(list)
    ptT = ptT.groupby(["chrom", "left", "ref_seq", "alt_seq"], as_index=False).agg(list)

    #insert N# and T# counts for each
    ptN["N#"] = ptN["Patient_ID"].apply(lambda x: len(x))
    ptT["T#"] = ptT["Patient_ID"].apply(lambda x: len(x))

    ptN = ptN.rename(columns={"Patient_ID": "Patient_ID_Normal"})
    ptT = ptT.rename(columns={"Patient_ID": "Patient_ID_Tumor"})

    # merge normal and tumor patient dataframes
    pt = pd.merge(ptT, ptN, how="outer", on=["chrom", "left", "ref_seq", "alt_seq"])

    pt[["N#", "T#"]] = pt[["N#", "T#"]].fillna(0)

############# Merged Dataframe with CSQ columns #################
    ## subset the dataframe to only the required info columns, reset  index (parameter drop=True, does not create an extra column with previous index numbers)
    if csqInfo == "yes":
        Merge = Merge[["chrom", "left", "right", "ref_seq", "alt_seq", "STRAND", "SYMBOL", "Consequence", "BIOTYPE", "Gene", "Feature", "Amino_acids", "cDNA_position", "CDS_position", "Protein_position", "Codons", "SIFT", "PolyPhen", "RefSeq", "Existing_variation"]]

    else:
        Merge = Merge[["chrom", "left", "right", "ref_seq", "alt_seq"]]

    Merge = Merge.reset_index(drop=True)

###### Combine for final, patient and variant CSQ #############
    final = pd.merge(pt, Merge, how="outer", on=["chrom", "left", "ref_seq", "alt_seq"])


    ## Exporting the 3 csv files 1) patient info, 2) variant/transcript info, 3) both patient and variant/transcript:
    # os.chdir(outPath)
    # outPT = cancerType + "_Patients.csv"
    # outVar = cancerType + "_Variants.csv"
    # outFinal = cancerType + "_Final.csv"
    #
    # pt.to_csv(outPT, sep=",", index=False)
    # Merge.to_csv(outVar, sep=",", index=False)
    # final.to_csv(outFinal, sep=",", index=False)

    return final, pt, Merge



#########################################**** Effect Analyzers ****#####################################################
'''
Before continuing: 
    1) Ensure that the effect analyzers have been generated and the output files are saved
    2) Identify whether the CSQ information is provided, if not an additional file and function (ensembl BioMart) will be needed to
        obtain the transcript lengths
'''


def addFATHMM(df, outputFile, outPath):
    fathmm = pd.read_csv(outputFile, sep=",", low_memory=False)
    fathmm = fathmm.drop_duplicates()
    fathmm = fathmm.reset_index(drop=True)

    fathmm["Chromosome"] = "chr" + fathmm["Chromosome"]
    ## or
    # fathmm["# Chromosome"] = "chr" + fathmm["# Chromosome"]


    newdf = pd.merge(df, fathmm, how="left", left_on=["chrom", "left", "ref_seq", "alt_seq"],
                    right_on=["Chromosome", "Position", "Ref. Base", "Mutant Base"])

    newdf.replace("--", 0, inplace=True)

    newdf["FATHMM_Score"] = newdf["Coding Score"].astype(float) + newdf["Non-Coding Score"].astype(float)

    newdf = newdf.drop(columns=["Coding Score", "Non-Coding Score", "Warning"])

    ## Exporting to ensure they are combining properly, comment out if successful:
    # os.chdir(outPath)
    # output = "FATHMM_Check.csv"
    # newdf.to_csv(output, sep=",", index=False)
    return newdf

##### Run for ThCa and OvCa #####


def cleanEnsemblVEP(outputFile):
    vep = pd.read_csv(outputFile, sep="\t", low_memory=False)

    vep = vep.drop_duplicates().reset_index(drop=True)

    vep.replace("-", "", inplace=True)

    vep[['chrom', 'left', "SNV"]] = vep["#Uploaded_variation"].str.split('_', expand=True)
    vep[['ref_seq', 'alt_seq']] = vep["SNV"].str.split('/', expand=True)
    vep["chrom"] = "chr" + vep["chrom"]
    vep = vep[vep["alt_seq"] == vep["Allele"]]
    vep[['Feature1', 'TX_Num']] = vep["Feature"].str.split('.', expand=True)

    vep = vep.drop(columns=["Feature"])
    vep = vep.rename(columns={"Feature1": "Feature"})

    ## Exporting to ensure they are combining properly, comment out if successful:
    # os.chdir(outPath)
    # output = "VEP_Check.csv"
    # vep.to_csv(output, sep=",", index=False)
    return vep


##### Run for ThCa and OvCa #####


def cleanEnsemblBioMart(outputFile):
    bm = pd.read_csv(outputFile, sep=",", low_memory=False)

    bm = bm[["Transcript stable ID version", "CDS Length"]]
    bm = bm.drop_duplicates().reset_index(drop=True)

    ## Exporting to ensure they are combining properly, comment out if successful:
    # output = "VEP_BM_Check.csv"
    # bm.to_csv(output, sep=",", index=False)
    return bm


'''
If Consequence data is provided within original VCF and only CADD score is needed use addEnsembl1, otherwise use addEnsembl2
'''


def addEnsembl1(df, vep):
    vep = vep[["chrom", "left", "ref_seq", "alt_seq", "CADD_PHRED", "SYMBOL", "Feature"]]

    df = df[
        ["chrom", "left", "ref_seq", "alt_seq", "T#", "N#", "SYMBOL", "Feature", "BIOTYPE", "CDS_position",
         "Amino_acids", "SIFT", "PolyPhen", "FATHMM_Score"]]


    col = ["SYMBOL", "Feature", "BIOTYPE", "CDS_position", "Amino_acids", "SIFT", "PolyPhen"]

    for ele in col:
        df[ele] = df[ele].apply(literal_eval)

    print(type(df["SYMBOL"][3]))
    df = df.explode(col, ignore_index=True)
    print(type(df["SYMBOL"][3]))

    df[['CDS Position', 'CDS Length']] = df["CDS_position"].str.split('/', expand=True)

    df["left"] = df["left"].astype(int)
    vep["left"] = vep["left"].astype(int)
    vep = vep.rename(columns={"CADD_PHRED": "CADD_Score"})

    newdf = pd.merge(df, vep, how="left", on=["chrom", "left", "ref_seq", "alt_seq", "SYMBOL", "Feature"])



    ## Exporting to ensure they are combining properly, comment out if successful:
    # output = "Check_addEnsembl1.csv"
    # newdf.to_csv(output, sep=",", index=False)
    return newdf


def addEnsembl2(df, vep, bm):
    vep = vep[["chrom", "left", "ref_seq", "alt_seq", "SYMBOL", "Feature", "BIOTYPE", "CDS_position",
         "Amino_acids", "SIFT", "PolyPhen", "CADD_PHRED"]]
    vep = vep.drop_duplicates().reset_index()
    vep = vep.rename(columns={"CDS_position": "CDS Position"})
    vep = vep.rename(columns={"CADD_PHRED": "CADD_Score"})

    ensembl = pd.merge(vep, bm, how="left", left_on="Feature", right_on="Transcript stable ID version")

    df["left"] = df["left"].astype(int)
    ensembl["left"] = ensembl["left"].astype(int)

    newdf = pd.merge(df, ensembl, how="left", on=["chrom", "left", "ref_seq", "alt_seq"])

    ## Exporting to ensure they are combining properly, comment out if successful:
    # output = "Check_addEnsembl2.csv"
    # newdf.to_csv(output, sep=",", index=False)
    return newdf


def addCGI_CADD(df, outputFile):
    cgi = pd.read_csv(outputFile, low_memory=False)
    cgi = cgi.drop_duplicates().reset_index(drop=True)

    cgi["Variation_Temp"] = cgi["Variation ID"].str.split(":")
    cgi["Variation"] = cgi["Variation_Temp"].apply(lambda x: ":".join(str(i) for i in x[:-1]))
    cgi = cgi.drop(columns="Variation_Temp")

    cgi = cgi[cgi["Deleterious (CADD score)"] != ""]
    cgi = cgi[cgi["Deleterious (CADD score)"] != "None"]
    cgi = cgi[cgi["Protein change"] != "None"]
    cgi = cgi[cgi["Protein change"] != "intron_variant"]
    cgi = cgi.drop_duplicates().reset_index(drop=True)

    cgi = cgi.groupby(["Chromosome", "Gene", "Transcript", "Protein change", "Variation"], as_index=False).agg(list)

    cgi["CADD_Score"] = cgi["Deleterious (CADD score)"].apply(lambda x: mean([float(i) for i in x]))
    cgi = cgi.rename(columns={"Gene": "Gene_CADD"})
    cgi = cgi.drop(columns="Transcript")


    df = df[
        ["chrom", "left", "ref_seq", "alt_seq", "T#", "N#", "SYMBOL", "Feature", "BIOTYPE", "CDS_position", "Protein_position",
         "Amino_acids", "SIFT", "PolyPhen", "FATHMM_Score"]]

    df["Variation"] = df["chrom"] + ":" + df["left"].astype(str) + ":" + df["ref_seq"] + "/" + df["alt_seq"]

    col = ["SYMBOL", "Feature", "BIOTYPE", "CDS_position", "Protein_position", "Amino_acids", "SIFT", "PolyPhen"]
    for ele in col:
        df[ele] = df[ele].apply(literal_eval)
    df = df.explode(col, ignore_index=True)

    df[['CDS Position', 'CDS Length']] = df["CDS_position"].str.split('/', expand=True)

    df = df[df["Amino_acids"].apply(len) == 3]
    df = df.reset_index(drop=True)

    df["temp_Amino_acids"] = df["Amino_acids"].str.split("/")
    df["temp_Protein_position"] = df["Protein_position"].str.split("/")
    df["Protein change"] = df["temp_Amino_acids"].str[0] + df["temp_Protein_position"].str[0] + df["temp_Amino_acids"].str[1]

    newdf = pd.merge(df, cgi, how="left", left_on=["Variation", "chrom", "Protein change"],
                  right_on=["Variation", "Chromosome", "Protein change"])

    return newdf


def cleanSIFT_PolyPhen(df):
    # SIFT:
    df[['SIFT_Description', 'SIFT_Score']] = df["SIFT"].str.split('(', expand=True)
    df["SIFT_Score"] = df["SIFT_Score"].str.replace(")", "")

    # PolyPhen:
    df[['PolyPhen_Description', 'PolyPhen_Score']] = df["PolyPhen"].str.split('(', expand=True)
    df["PolyPhen_Score"] = df["PolyPhen_Score"].str.replace(")", "")

    # Convert Score columns into float:
    df["SIFT_Score"] = df["SIFT_Score"].astype(float)
    df["PolyPhen_Score"] = df["PolyPhen_Score"].astype(float)

    # Inverse SIFT:
    df["Inv_SIFT_Score"] = 1.0 - df["SIFT_Score"]

    ## Exporting to ensure they are combining properly, comment out if successful:
    # os.chdir(outPath)
    # output = "Check_parseSIFT_PolyPhen.csv"
    # df.to_csv(output, sep=",", index=False)
    return df


def minMaxNorm(df, col):
    df[col] = df[col].astype(float)

    newCol = col + "_Normalized"
    ## Min Max normalization
    df[newCol] = (df[col] - df[col].min()) / (df[col].max() - df[col].min())

    ## Exporting to ensure they are combining properly, comment out if successful:
    # os.chdir(outPath)
    # output = "Check_Normalization.csv"
    # df.to_csv(output, sep=",", index=False)
    return df


'''
Parameters determine which variants to use for Q(G) scoring
Biotype (tx): "protein_coding", transcript type
AminoAcid (aa): 3, length of column 
Score (s): ScoreDic = {"FATHMM_Score":0.5, "SIFT_Score":.05, "PolyPhen_Score":.447, "CADD_Score":10.0}, dictionary of cutoffs
'''

def filterSNV(df, aa, tx, s, outPath):

    # Non-Synonymous mutations have length 3, Synonymous have 1
    df = df[df["Amino_acids"].str.len() == aa]

    # Only protein coding variants
    df = df[df["BIOTYPE"] == tx]

    # Keep rows with given score cutoffs
    newdf = df[((df["FATHMM_Score"].astype(float) >= s["FATHMM_Score"]) | (df["SIFT_Score"].astype(float) <= s["SIFT_Score"]) | (df["PolyPhen_Score"].astype(float) >= s["PolyPhen_Score"]) | (df["CADD_Score"].astype(float) >= s["CADD_Score"]))]

    ## Exporting to ensure they are combining properly, comment out if successful:
    # os.chdir(outPath)
    # output = "Check_Filters.csv"
    # newdf.to_csv(output, sep=",", index=False)
    return newdf





def qGene(csv, score, s, output, path):
    os.chdir(path)
    # Take away the unknown values from the Column that is being scored, depending if
    if score != "Average":
        scoreCol = score + "_Normalized"
    else:
        scoreCol = score

    csv = csv[csv[scoreCol] != "Unknown"]
    csv = csv[csv[scoreCol] != "None"]
    csv = csv[csv[scoreCol] != ""]
    csv = csv[csv[scoreCol] != "-"]
    csv = csv[csv[score] >= s[score]]
    csv = csv.dropna(subset=[scoreCol])
    csv = csv.reset_index(drop=True)


    print("Step 1")
    # Will need to remove any Transcripts that do not have a given length
    csv = csv[csv["CDS Length"] != ""]
    csv = csv[csv["CDS Length"] != "-"]
    csv = csv.dropna(subset=["CDS Length"])

        ## Whenever values are removed from dataframe, reset index:
    csv = csv.reset_index(drop=True)


    # First tier, unique genes with the transcript list, include print statement to see how many within this list:
    csv1 = csv[["SYMBOL", "Feature"]]
    csv1 = csv1.drop_duplicates().reset_index(drop=True)

    # Intermediate step, want to know how many unique variants per gene
    csv11 = csv[["SYMBOL", "chrom", "left", "ref_seq", "alt_seq"]]
    csv11 = csv11.drop_duplicates().reset_index(drop=True)
    csv11 = csv11.groupby("SYMBOL", as_index=False).agg(list)
    csv11["#Variants per Gene"] = csv11["chrom"].str.len()
    csv11 = csv11.drop(columns=["chrom", "left", "ref_seq", "alt_seq"])


    # Second tier,
    csv0 = csv[["chrom", "left", "ref_seq", "alt_seq", "Feature", "CDS Length", "T#", "N#", scoreCol]]
    csv0 = csv0.drop_duplicates(subset=["chrom", "left", "ref_seq", "alt_seq", "Feature"]).reset_index(drop=True)


    step2 = "Transcripts_" + output
    step3 = "Gene_" + output
    step4 = "PreQ_" + output

    csv2 = csv0.groupby("Feature", as_index=False).agg(list)
    csv2["#Variants per Transcript"] = csv2["chrom"].str.len()

    # print statement to see how many unique transcripts, need to be seen after grouping
    uniqueTx = str(csv2.shape[0])
    print("With the given parameters, this dataset contains " + uniqueTx + " unique transcripts")


    csv2.to_csv(step2, sep=",", index = False)

    csv3 = pd.merge(csv1, csv2, how="left", on="Feature")
    csv3.to_csv(step3, sep=",", index=False)

    csv4 = csv3.groupby("SYMBOL", as_index=False).agg(list)
    csv4 = pd.merge(csv4, csv11, how="left", on="SYMBOL")
    csv4["#Transcripts per Gene"] = csv4["Feature"].str.len()


    # print statement to see how many unique genes, need to be seen after grouping
    uniqueGenes = str(csv4.shape[0])
    print("With the given parameters, this dataset contains " + uniqueGenes + " unique genes")

    csv4.to_csv(step4, sep=",", index=False)


    print("Step 6")
    score = []
    for row in range(csv4.shape[0]):

        # number of unique transcripts per gene
        Ng = len(csv4["Feature"][row])

        #list of list that need to be iterated thorugh in next loop
        row_score = csv4[scoreCol][row]
        row_tumor = csv4["T#"][row]
        row_normal = csv4["N#"][row]
        row_length = csv4["CDS Length"][row]
        count = 1
        sum1 = 0
        for (i, j, k, g) in zip(row_score, row_tumor, row_normal, row_length):
            sum2 = 0
            for (s, t, n, l) in zip(i, j, k , g):
                count+=1
                # temp = 0
                # l = l.split("/")
                l1 = int(l)

                ##### Without the mutational data, t-n will always be equal to 1 ######
                # temp = s * (1)

                ###### With the mutational data, t-n is determined by mutational data ######
                temp = s * (t - n)


                temp = temp / math.log(l1)
                sum2 += temp
            sum1 += sum2

        score1 = sum1 / Ng

        score.append(score1)
    outCol = "Q(Gene)"
    csv4.insert(csv4.shape[1], outCol, score)
    csv4 = csv4[["SYMBOL", outCol, "#Variants per Transcript", "#Variants per Gene", "#Transcripts per Gene"]]
    # Export Q_Gene score dataframe, function will return new dataframe to merge (how="outer") all the results into a single dataframe

    csv4.to_csv(output, sep=",", index = False)
    return csv4


def rankAve2(csvList, knownGenes, cancer, s, outPath):
    os.chdir(outPath)
    outfile = s + "_Rank_Average_THCA.txt"
    file1 = open(outfile, "w+")
    temp = cancer + ":"
    file1.write(temp)
    file1.write("\n")
    iteration = 0

    for csv in csvList:
        sum = 0
        csv = csv[["SYMBOL", "Q(Gene)"]]
        csv["Q(Gene)"] = csv["Q(Gene)"].astype(float)

        csv = csv.sort_values("Q(Gene)", ascending=False)
        csv = csv.reset_index(drop=True)
        csv = csv.reset_index().rename(columns={"index": "Rank"})
        l = len(csv)
        length = len(csv)
        den = len(knownGenes)
        count = 0
        nCSV = pd.merge(csv, knownGenes, how="right")

        for row in range(nCSV.shape[0]):
            if pd.isnull(nCSV.at[row, "Rank"]) is True:
                nCSV["Rank"][row] = length
                length += 1
                count += 1
        nCSV["Rank"] = nCSV["Rank"] + 1
        # nCSV["ln(Rank)"] = nCSV["Rank"].apply(lambda x: math.log(x))
        sum = nCSV["Rank"].sum()
        # sum = sum * count
        average = sum/float(den)
        Lave = average/float(l)
        sar = Lave - (1.0/float(l))
        # file1.write(extLst[iteration])
        iteration += 1
        file1.write("\n")
        file1.write("Average Rank:\t")
        file1.write(str(sar))
        file1.write("\n")
        file1.write("Number of Known Genes without Q(Gene) Score:\t")
        file1.write(str(count))
        file1.write("\n")
        file1.write("\n")
    file1.write("\n")
    file1.write("\n")
    file1.close()
    return



'''
Complete code running the functions for the iQG analysis on both OvCa and ThCa, 
Put all files in the working directory where the python file is saved, 
insert the working directory pathway below
example of pathway format: "/Users/abataycan/"
'''
pathway = "(Insert Working Directory Pathway)"


## Thyroid
# thcaFinal, thcaPT, thcaMerge = multi_vcf2csv("/Users/abataycan/VCF23/raw/thyroid/Thyroid_VCF/", "yes", "/Users/abataycan/PycharmProjects/Manuscript/Takes/Take4/Thyroid", "ThCa")
# thcaFinal.to_csv("ThCa_variants.csv", sep=",", index=False)

thcaFinal = pd.read_csv("ThCa_variants.csv", sep=",")


thca1 = addFATHMM(thcaFinal, "ThCa_FATHMM.csv", pathway)
thcaVEP = cleanEnsemblVEP("ThCa_CADD.txt")
thca2 = addEnsembl1(thca1, thcaVEP)
thca2 = cleanSIFT_PolyPhen(thca2)

ScoreDic = {"FATHMM_Score":0.5, "SIFT_Score":.05, "Inv_SIFT_Score":.95, "PolyPhen_Score":.447, "CADD_Score":10.0, "Average":0.0}
thca3 = filterSNV(thca2, 3, "protein_coding", ScoreDic, pathway)

col = ["FATHMM_Score", "Inv_SIFT_Score", "PolyPhen_Score", "CADD_Score"]
for ele in col:
    thca3 = minMaxNorm(thca3, ele)

qCols = ["FATHMM_Score_Normalized", "Inv_SIFT_Score_Normalized", "PolyPhen_Score_Normalized", "CADD_Score_Normalized"]
thca3['Average'] = thca3[qCols].mean(skipna=True, axis=1)

sDic = {"FATHMM_Score":0.5, "Inv_SIFT_Score":.95, "PolyPhen_Score":.447, "CADD_Score":10.0, "Average":0.0}

thca4 = qGene(thca3, "FATHMM_Score", ScoreDic, "Q(Gene)_FATHMM.csv", pathway)
thca5 = qGene(thca3, "Inv_SIFT_Score", ScoreDic, "Q(Gene)_SIFT.csv", pathway)
thca6 = qGene(thca3, "PolyPhen_Score", ScoreDic, "Q(Gene)_PolyPhen.csv", pathway)
thca7 = qGene(thca3, "CADD_Score", ScoreDic, "Q(Gene)_CADD.csv", pathway)
thca8 = qGene(thca3, "Average", ScoreDic, "iQ(Gene).csv", pathway)



KGthca = pd.read_csv("ThCa_KnownGenes.csv", sep=",")
csvTHCA1 = [thca4]
csvTHCA2 = [thca5]
csvTHCA3 = [thca6]
csvTHCA4 = [thca7]
csvTHCA5 = [thca8]
#
rankAve2(csvTHCA1, KGthca, "THCA", "FATHMM", pathway)

rankAve2(csvTHCA2, KGthca, "THCA", "SIFT", pathway)
rankAve2(csvTHCA3, KGthca, "THCA", "PolyPhen", pathway)

rankAve2(csvTHCA4, KGthca, "THCA", "CADD", pathway)
rankAve2(csvTHCA5, KGthca, "THCA", "Average", pathway)




## Ovarian
# ocaFinal, ocaPT, ocaMerge = multi_vcf2csv("/Users/abataycan/VCF23/raw/ovary/Ovary_VCF/", "yes", "/Users/abataycan/PycharmProjects/Manuscript/Takes/Take6/Ovarian", "OCa")
# ocaFinal.to_csv("/Users/abataycan/PycharmProjects/Manuscript/Takes/OCA_Khodeza.csv", sep=",", index=False)

ocaFinal = pd.read_csv("OvCa_variants.csv", sep=",")
oca1 = addFATHMM(ocaFinal, "OvCa_FATHMM.csv", pathway)
oca2 = addCGI_CADD(oca1, "OvCa_CADD.csv")
oca2 = cleanSIFT_PolyPhen(oca2)

ScoreDic = {"FATHMM_Score":0.5, "SIFT_Score":.05, "Inv_SIFT_Score":.95, "PolyPhen_Score":.447, "CADD_Score":10.0, "Average":0.0}
oca3 = filterSNV(oca2, 3, "protein_coding", ScoreDic, pathway)

col = ["FATHMM_Score", "Inv_SIFT_Score", "PolyPhen_Score", "CADD_Score"]
for ele in col:
    oca3 = minMaxNorm(oca3, ele)

qCols = ["FATHMM_Score_Normalized", "Inv_SIFT_Score_Normalized", "PolyPhen_Score_Normalized", "CADD_Score_Normalized"]
oca3['Average'] = oca3[qCols].mean(skipna=True, axis=1)

oca4 = qGene(oca3, "FATHMM_Score", ScoreDic, "Q(Gene)_FATHMM.csv", pathway)
oca5 = qGene(oca3, "Inv_SIFT_Score", ScoreDic, "Q(Gene)_SIFT.csv",pathway)
oca6 = qGene(oca3, "PolyPhen_Score", ScoreDic, "Q(Gene)_PolyPhen.csv", pathway)
oca7 = qGene(oca3, "CADD_Score", ScoreDic, "Q(Gene)_CADD.csv", pathway)
oca8 = qGene(oca3, "Average", ScoreDic, "iQ(Gene).csv", pathway)



KGoca = pd.read_csv("OvCa_KnownGenes.csv", sep=",")
csvOCA1 = [oca4]
csvOCA2 = [oca5]
csvOCA3 = [oca6]
csvOCA4 = [oca7]
csvOCA5 = [oca8]
#
rankAve2(csvOCA1, KGoca, "OCA", "FATHMM", pathway)

rankAve2(csvOCA2, KGoca, "OCA", "SIFT", pathway)
rankAve2(csvOCA3, KGoca, "OCA", "PolyPhen", pathway)

rankAve2(csvOCA4, KGoca, "OCA", "CADD", pathway)
rankAve2(csvOCA5, KGoca, "OCA", "Average", pathway)





