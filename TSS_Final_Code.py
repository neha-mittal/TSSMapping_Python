# import pyranges as pr
# import pandas as pd
# import pybedtools
# from fuc import pybed
 
 # *********To Generate TSS Sites from the known gene ref file********** #

# gr = pr.read_gtf("/Users/neha/Desktop/Gene_Mapping/Reference_genomes/hg38.knownGene.gtf") #read the known gene ref file using pyranges

# print(gr.columns)

# gr = gr[["Chromosome",'Feature', "Start", 'End', 'Strand','gene_id']] # Filter the columns

# print(gr)

# gr = gr.features.tss() #to find the tss site using pyranges


# s = (gr.Chromosome == "chr21") #Filter Chr21

# gr[s].df.to_csv('tss_sites_knowngene.csv', header=True, index=False) # Saving tss sites to a csv file


#  # *********To convert csv file to bed file  ********** #

# df = pd.read_csv('tss_sites_knowngene.csv', usecols = ['Chromosome','Start', 'End','Strand',"gene_id"])
# df.columns = ['Chromosome', 'Start', 'End', 'Strand', "gene_id"]

# bf = pybed.BedFrame.from_frame(meta=[], data=df)
# bf.to_file('tss_sites_knowngene.bed')


#  # *********To extend TSS sites 5kb upstream & 5kb downstream ********** #


# tss_sites = pybedtools.BedTool('tss_sites_knowngene.bed') #get the tss sites bed file

# tss_extend_5kb= tss_sites.slop(b=5000, genome='hg38', output='tss_extend-5kb.bed') #extend 5kb upstream and downstream

#  # *********To find intersection of sample with 5kb upstream & 5kb downstream region of tss  ********** #

# bt_sample = pybedtools.BedTool("BT_sample.bed")   #read sample file


# #Fag - wo :- Write the original A and B entries plus the number of base pairs of overlap between the 
# #two features. Only A features with overlap are reported
# sample_intersect_5kb = tss_extend_5kb.intersect(bt_sample, wo=True, output="sample_intersect_5kb_wo.bed")

## **********Calculating the bin from TSS site Upstream and Downstream ## *********
# tss_Minus_1000 = tss_sites.slop(l=1000, r=0, s=True, genome='hg38', output="TSSMinus1000.bed")
# tss_Minus_2000 = tss_Minus_1000.slop(l=1000, r=-1000, s=True, genome='hg38', output="TSSMinus2000.bed")
# tss_Minus_3000 = tss_Minus_2000.slop(l=1000, r=-1000, s=True, genome='hg38', output="TSSMinus3000.bed")
# tss_Minus_4000 = tss_Minus_3000.slop(l=1000, r=-1000, s=True, genome='hg38', output="TSSMinus4000.bed")
# tss_Minus_5000 = tss_Minus_4000.slop(l=1000, r=-1000, s=True, genome='hg38', output="TSSMinus5000.bed")

# tss_plus_1000 = tss_sites.slop(l=0, r=1000, s=True, genome='hg38', output="TSSPlus1000.bed")
# tss_plus_2000 = tss_plus_1000.slop(l=-1000, r=1000, s=True, genome='hg38', output="TSSPlus2000.bed")
# tss_plus_3000 = tss_plus_2000.slop(l=-1000, r=1000, s=True, genome='hg38', output="TSSPlus3000.bed")
# tss_plus_4000 = tss_plus_3000.slop(l=-1000, r=1000, s=True, genome='hg38', output="TSSPlus4000.bed")
# tss_plus_5000 = tss_plus_4000.slop(l=-1000, r=1000, s=True, genome='hg38', output="TSSPlus5000.bed")

## **********Calculating the coverage count of reads from TSS site and Upstream/Downstream region in each bin ## *********

# bt_sample = pybedtools.BedTool("BT_sample.bed")
# bt_sample1 = pr.read_bed("BT_sample.bed")

# TSSMinus1000_count = tss_Minus_1000.coverage(bt_sample, counts=True, output="TSSMinus1000_coverage.bed")
# TSSMinus2000_count = tss_Minus_2000.coverage(bt_sample, counts=True, output="TSSMinus2000_coverage.bed")
# TSSMinus3000_count = tss_Minus_3000.coverage(bt_sample, counts=True, output="TSSMinus3000_coverage.bed")
# TSSMinus4000_count = tss_Minus_4000.coverage(bt_sample, counts=True, output="TSSMinus4000_coverage.bed")
# TSSMinus5000_count = tss_Minus_5000.coverage(bt_sample, counts=True, output="TSSMinus5000_coverage.bed")

# TSSPlus1000_count = tss_plus_1000.coverage(bt_sample, counts=True, output="TSSPlus1000_coverage.bed")
# TSSplus2000_count = tss_plus_2000.coverage(bt_sample, counts=True, output="TSSPlus2000_coverage.bed")
# TSSPlus3000_count = tss_plus_3000.coverage(bt_sample, counts=True, output="TSSPlus3000_coverage.bed")
# TSSPlus4000_count = tss_plus_4000.coverage(bt_sample, counts=True, output="TSSPlus4000_coverage.bed")
# TSSPlus5000_count = tss_plus_5000.coverage(bt_sample, counts=True, output="TSSPlus5000_coverage.bed")

## **********Calculating the coverage fraction of reads from TSS site and Upstream/Downstream region in each bin ## *********

# bt_sample = pybedtools.BedTool("BT_sample.bed")
# bt_sample1 = pr.read_bed("BT_sample.bed")

# TSSMinus1000_count = tss_Minus_1000.coverage(bt_sample, output="TSSMinus1000_coverage_f.bed")
# TSSMinus2000_count = tss_Minus_2000.coverage(bt_sample, output="TSSMinus2000_coverage_f.bed")
# TSSMinus3000_count = tss_Minus_3000.coverage(bt_sample, output="TSSMinus3000_coverage_f.bed")
# TSSMinus4000_count = tss_Minus_4000.coverage(bt_sample, output="TSSMinus4000_coverage_f.bed")
# TSSMinus5000_count = tss_Minus_5000.coverage(bt_sample, output="TSSMinus5000_coverage_f.bed")

# TSSPlus1000_count = tss_plus_1000.coverage(bt_sample, output="TSSPlus1000_coverage_f.bed")
# TSSplus2000_count = tss_plus_2000.coverage(bt_sample, output="TSSPlus2000_coverage_f.bed")
# TSSPlus3000_count = tss_plus_3000.coverage(bt_sample, output="TSSPlus3000_coverage_f.bed")
# TSSPlus4000_count = tss_plus_4000.coverage(bt_sample, output="TSSPlus4000_coverage_f.bed")
# TSSPlus5000_count = tss_plus_5000.coverage(bt_sample, output="TSSPlus5000_coverage_f.bed")

### ******* Coverage Total Sum Upstream ***********

# gr_1 = pr.read_bed("TSSMinus1000_coverage.bed")
# gr_1.df.to_csv("TSSMinus1000_coverage.csv")
# # print(gr_1)

# coverage_tssminus1000 = pd.read_csv("TSSMinus1000_coverage.csv", usecols = [6])
# # # print(coverage_tssminus1000)
# list_tssminus1000 = coverage_tssminus1000["Strand"].to_list()
# # print(list_tssminus1000)
# list_tssminus1000 = (sum(list_tssminus1000))
# print(list_tssminus1000)

# gr_2 = pr.read_bed("TSSMinus2000_coverage.bed")
# gr_2.df.to_csv("TSSMinus2000_coverage.csv")
# # print(gr_2)

# coverage_tssminus2000 = pd.read_csv("TSSMinus2000_coverage.csv", usecols = [6])
# # print(coverage_tssminus2000)
# list_tssminus2000 = coverage_tssminus2000["Strand"].to_list()
# # print(list_tssminus2000)
# list_tssminus2000 = (sum(list_tssminus2000))
# print(list_tssminus2000)

# gr_3 = pr.read_bed("TSSMinus3000_coverage.bed")
# gr_3.df.to_csv("TSSMinus3000_coverage.csv")
# # # print(gr_3)

# coverage_tssminus3000 = pd.read_csv("TSSMinus3000_coverage.csv", usecols = [6])
# # print(coverage_tssminus3000)
# list_tssminus3000 = coverage_tssminus3000["Strand"].to_list()
# # print(list_tssminus3000)
# list_tssminus3000 = (sum(list_tssminus3000))
# print(list_tssminus3000)

# gr_4 = pr.read_bed("TSSMinus4000_coverage.bed")
# gr_4.df.to_csv("TSSMinus4000_coverage.csv")
# # print(gr_4)

# coverage_tssminus4000 = pd.read_csv("TSSMinus4000_coverage.csv", usecols = [6])
# # print(coverage_tssminus4000)
# list_tssminus4000 = coverage_tssminus4000["Strand"].to_list()
# # print(list_tssminus4000)
# list_tssminus4000 = (sum(list_tssminus4000))
# print(list_tssminus4000)

# gr_5 = pr.read_bed("TSSMinus5000_coverage.bed")
# gr_5.df.to_csv("TSSMinus5000_coverage.csv")
# # print(gr_5)

# coverage_tssminus5000 = pd.read_csv("TSSMinus5000_coverage.csv", usecols = [6])
# # print(coverage_tssminus5000)
# list_tssminus5000 = coverage_tssminus5000["Strand"].to_list()
# # print(list_tssminus5000)
# list_tssminus5000 = (sum(list_tssminus5000))
# print(list_tssminus5000)

# Upstream_List = [86302, 55212, 33950, 45580, 48944]
# print(Upstream_List)


### ******* Coverage Total Sum Downstream ***********

# gr_1_d = pr.read_bed("TSSPlus1000_coverage.bed")
# gr_1_d.df.to_csv("TSSPlus1000_coverage.csv")
# print(gr_1_d)

# coverage_tssplus1000 = pd.read_csv("TSSPlus1000_coverage.csv", usecols = [6])
# print(coverage_tssplus1000)
# list_tssplus1000 = coverage_tssplus1000["Strand"].to_list()
# print(list_tssplus1000)
# list_tssplus1000 = (sum(list_tssplus1000))
# print(list_tssplus1000)

# gr_2_d = pr.read_bed("TSSPlus2000_coverage.bed")
# gr_2_d.df.to_csv("TSSPlus2000_coverage.csv")
# print(gr_2_d)

# coverage_tssplus2000 = pd.read_csv("TSSPlus2000_coverage.csv", usecols = [6])
# print(coverage_tssplus2000)
# list_tssplus2000 = coverage_tssplus2000["Strand"].to_list()
# print(list_tssplus2000)
# list_tssplus2000 = (sum(list_tssplus2000))
# print(list_tssplus2000)

# gr_3_d = pr.read_bed("TSSPlus3000_coverage.bed")
# gr_3_d.df.to_csv("TSSPlus3000_coverage.csv")
# print(gr_3_d)

# coverage_tssplus3000 = pd.read_csv("TSSPlus3000_coverage.csv", usecols = [6])
# # print(coverage_tssplus3000)
# list_tssplus3000 = coverage_tssplus3000["Strand"].to_list()
# print(list_tssplus3000)
# list_tssplus3000 = (sum(list_tssplus3000))
# print(list_tssplus3000)

# gr_4_d = pr.read_bed("TSSPlus4000_coverage.bed")
# gr_4_d.df.to_csv("TSSPlus4000_coverage.csv")
# print(gr_4_d)

# coverage_tssplus4000 = pd.read_csv("TSSPlus4000_coverage.csv", usecols = [6])
# # print(coverage_tssplus4000)
# list_tssplus4000 = coverage_tssplus4000["Strand"].to_list()
# print(list_tssplus4000)
# list_tssplus4000 = (sum(list_tssplus4000))
# print(list_tssplus4000)

# gr_5_d = pr.read_bed("TSSPlus5000_coverage.bed")
# gr_5_d.df.to_csv("TSSPlus5000_coverage.csv")
# print(gr_5_d)

# coverage_tssplus5000 = pd.read_csv("TSSPlus5000_coverage.csv", usecols = [6])
# # print(coverage_tssplus5000)
# list_tssplus5000 = coverage_tssplus5000["Strand"].to_list()
# print(list_tssplus5000)
# list_tssplus5000 = (sum(list_tssplus5000))
# print(list_tssplus5000)

# Downstream_List = [35050, 38668, 51121, 46646, 50926]
# print(Downstream_List)

# Final_List_Count = Upstream_List + Downstream_List
# print ("Concatenated list: " + str(Final_List_Count)) 
#####################*************************************************##################################
#### New Pipeline for TSS_Sites #######
import pyranges as pr
import pandas as pd
import pybedtools
from fuc import pybed
import genomepy
import numpy

### Identification of the TSS Sites ###

# genome = pr.read_gtf("/Users/neha/Desktop/TSS_Gene_Mapping/Reference_genomes/hg38.refGene.gtf")
# print(genome)
# print(genome.columns) #printing all the columns name from genome variable;

# genome1 = (genome.Chromosome == "chr21") & (genome.Feature == "transcript") #Filtered all the transcript rows on chromosome 21 from the reference file but still having all the data columns
# print(genome[genome1])

# genome1 = (genome.Feature == "transcript")
# # print(genome[genome1])
# genome2 = genome[genome1].gene_id.drop_duplicates()
# print(genome[genome2])

# genome1 = genome[genome1].gene_id.drop_duplicates() #removing the genes duplicates using the coulmn gene_name and storing the gene_name in a list
# # print(genome1)

# genome1 = genome[genome1].drop_duplicate_positions(strand = True, keep = 'first') 
# print(genome1)

# genome2 = genome[genome1].gene_id.drop_duplicates()
# print(genome2)

# genome1 = genome1.features.tss() #to find the tss site using pyranges
# print(genome1)

# genome1.df.to_csv('TSS_Sites_Chr21.csv', header=True, index=True) #storing the TSS sites in to the csv format

#  # *********To convert TSS csv file to bed file  ********** # #

# TSS_sites = pd.read_csv('TSS_Sites_Chr21.csv', usecols = ['Chromosome','Start', 'End','Strand',"gene_id"])
# TSS_sites_chr21 = TSS_sites.drop_duplicates(subset='gene_id', keep="first", ignore_index=True) #filtered the gene ID
# TSS_sites_chr21.to_csv('Tss_pyranges.csv', header=True, index=True)

# bf = pybed.BedFrame.from_frame(meta=[], data=TSS_sites_chr21) #data frame converted to bed file
# bf.to_file('TSS_Sites_chr21_Final.bed')


# df = pd.read_csv('tss_sites_knowngene.csv', usecols = ['Chromosome','Start', 'End','Strand',"gene_id"])
# df.columns = ['Chromosome', 'Start', 'End', 'Strand', "gene_id"]

# bf = pybed.BedFrame.from_frame(meta=[], data=df)
# bf.to_file('tss_sites_knowngene.bed')


# print(TSS_sites_chr21)
# TSS_sites = TSS_sites.drop_duplicates(subset='gene_id', keep="first", ignore_index=False)
# TSS_sites_plus_strand = print(TSS_sites.loc[0:201, :])
# TSS_sites_minus_strand = print(TSS_sites.loc[202:, :])
# TSS_sites.columns = ['Chromosome', 'Start', 'End', 'Strand', "gene_id"]
# print(TSS_sites.columns)

# TSS_Sites_Chr21 = pybed.BedFrame.from_frame(meta=[], data=TSS_sites)
# TSS_Sites_Chr21.to_file('Tss_Sites_Chr21.bed')

### Extension of the TSS Sites in to upstream and downstream regions among both strands ###

# tss_sites_chr21 = pybedtools.BedTool("TSS_Sites_chr21_Final.bed")   #read the converted bed file and saving in a separate variable;


# Tss_up_down_5k = tss_sites_chr21.slop(b=5000, genome= 'hg38', output="TSS_up_down_5k.bed")# upstream and downstream region identification
# # print(Tss_up_down_5k)

#Sample File#

Read1 = pd.read_csv('/Users/neha/Desktop/TSS_Gene_Mapping/Samples_Data/streck_cfDNA_205_S3__chr21__alignment_summary.csv') #Reading the sample file;

Read1.columns = ['Chromosome', 'Start', 'End', 'Mapping Quality', 'Fragment Length'] #Adding the name in each columns;
# print(Read1)
Read1 = pr.PyRanges(Read1) #Converting the Start, and End Columns in to the range using PyRanges function;
# print(Read1)
Read1 = Read1.drop_duplicate_positions(strand = None, keep = 'first') #Removing the duplicated range of start and end position;
# print(Read1)

# #Just taking the first three columns
# Read2 = pd.read_csv('/Users/neha/Desktop/TSS_Gene_Mapping/Samples_Data/streck_cfDNA_205_S3__chr21__alignment_summary.csv', usecols=[0,1,2]) #Reading the sample file;

# Read2.columns = ['Chromosome', 'Start', 'End'] #Adding the name in each columns;
# print(Read2)
# Read2 = pr.PyRanges(Read2) #Converting the Start, and End Columns in to the range using PyRanges function;
# print(Read2)
# Read2 = Read2.drop_duplicate_positions(strand = None, keep = 'first') #Removing the duplicated range of start and end position;
# print(Read2)

TSS1 = pd.read_csv('/Users/neha/Desktop/TSS_Gene_Mapping/Tss_pyranges.csv', usecols = [1,2,3,4]) #Reading the sample file;
# print (TSS1)
TSS1 = pr.PyRanges(TSS1) #Converting the Start, and End Columns in to the range using PyRanges function;
# print(TSS1)
TSS2 = TSS1.extend({"3": 5000, "5": 5000}) #Calculating the upstream and downstream 5kb region;
# TSS2 = pr.Pyranges(TSS2)
TSS2_1 = TSS2[TSS2.Start < 17912340]
# print(TSS2_1)
# TSS2.boundaries("gene_id")

#**** Count of reads overlapping within the TSS start, end range******** #
# Overlap_count = TSS2.count_overlaps(Read1, strandedness = None, overlap_col="Count")
# print(Overlap_count)

#****** Calculating the run length encodings using read data file using RLE function ********#

coverage_read1 = (Read1.to_rle()) #Calculating the coverage using read length encoding function of pyranges
# print(coverage_read1["chr21", "+"])


# coverage_read2 = (Read2.to_rle())
# print(coverage_read2)
# site = TSS2.Start.values[0],TSS2.End.values[0]
# site = TSS2_1
# print(site)
# print(coverage_read1[TSS_s_1])
# r1 = coverage_read1[["chr21"]][Start(site):End(site)]
# print(r1)






