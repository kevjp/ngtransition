import pandas as pd
import re, os
import argparse
import subprocess
from pandas.lib import astype_str

'''
Created on 4 Oct 2016

@author: kevin
'''

'''
HGVS class take a VLookup csv file and extract SNP regions. These SNP regions are fed in to Mutalyser batch upload to generate variant in HGVS format.
'''
class Hgvsconvert (object):
    
    
    def __init__(self):
        self.input = "/home/kevin/Documents/HGVS/VLookup_table_duplicates_removed.csv" # input file relates variant to Pathology class
        self.keyfile = "/home/kevin/Documents/HGVS/VLookup_tableNMaccession.csv" # key file relates gene symbol to a specific NM accession 
        self.output = []
        self.hgvs = []
        self.clas = []
        self.convertcount = 0 # counter for the number of variants sucessfully converted
        self.notconvertcount = 0
        self.singleexongenevariants = 0
        self.LOClist = []
        

    def gene2nm (self):
        inputfile = pd.read_table(self.input, header= 0) # takes the input file
        accessionfile = pd.read_table(self.keyfile, header= 0) # takes the key file
        # List of variants
        
        varlist = inputfile[inputfile.columns[0]].tolist() # generate a list of variants from vlookup file
        genelist = accessionfile[accessionfile.columns[0]].tolist() # generate a list of genes from key file
        # Open log file for run analysis
        log = open('/home/kevin/Documents/HGVS/Vlookupconversion_Log.txt','w') # Generate an activity log file 
        # iterate over key file
        for index, gene, accession, keypanel in accessionfile[accessionfile.columns[0:3]].itertuples():
            
            # Discount the following single exon genes 
            singleexongenes = ["FKRP", "GTDC2", "DOLK"]
            
            if gene in singleexongenes:
                g = None
                log.write("Single exon gene:" + gene + "\n")
                pass
            
            elif gene in genelist:
                g = gene
                # strip white space from gene input
                gene = gene.strip()
                # add backslash to protect bracketed values in gene symbols and remove white space this is required for downstream re.sub function
                gene = gene.replace("(", "\(").replace(")", "\)").replace(" ", "")
                # pattern from key file used to find corresponding gene symbol in vlookup file
                pattern = gene + "c."
                # replacement string for pattern    
                replace = accession + "(" + gene + ")" + ":c."
                # counter for the number of variant entries that were not recognised by the key file that are in the input file
                self.notconvertcount = 0
                # iterate over input file 
                for index , variant, clas, panel in inputfile[inputfile.columns[0:3]].itertuples():
                    
                    v = variant # v is used at end of for loop to remove the variant from the list once it has been sucessfully matched
                    # remove all white space from variant entry
                    variant = variant.replace(" ", "")
                    
                    # Perform match only for variants and genes from the same panel i.e. don't perform matching indiscrimiately
                    if panel == keypanel:                                  
                        # variant not present in panel then pass on to next item in input file (re.sub returns item not matched in this case variant is the item compared against pattern for a match and replace is the item to replace with if matched)
                        if re.sub(pattern, replace, variant) == variant:
                            self.notconvertcount += 1
                            pass
        
                        else:
                            # append matched values to hgvs list for panel and do the same for the classifcation column
                            self.hgvs.append(re.sub(pattern, replace, variant))
                            self.clas.append(clas)
                            self.convertcount += 1
                            self.notconvertcount = 0
                        
                            #Remove variant item from variant list to check if any variants have not been matched at the end of the process varlist contains all variants that were not matched
                            varlist.remove(v)
                            
                        
            if g is not None:
                # Remove gene from genelist if gene has been successfully processed
                genelist.remove(g)
                        
        #Write variants not processed by above pattern search to log file
        for varnotincluded in varlist:
            log.write("variant not included:" + varnotincluded + "\n")
            
            # Count the number of singleexon variants 
            if any(singleexons in varnotincluded for singleexons in ("FKRP", "GTDC2", "DOLK")) == True:
                self.singleexongenevariants += 1
                       
            
            
        # LOC gene symbols in vlookup but not in key        
        for index , variant, clas in inputfile[inputfile.columns[0:2]].itertuples():
            if variant.startswith("LOC") == True:
                self.LOClist.append(variant)
            
        log.write("The number of variants starting with LOC (no gene symbol assignment) equal:" + str(len(self.LOClist)) + "\n")
                 
        
        # write outputs to log file
        log.write("The number of variants not parsed equal:" + str(len(varlist)) + "\n")
        log.write("Of which the number of single exon variants not included equal:" + str(self.singleexongenevariants) + "\n")
        log.write("This means the number of other variants not included which are non single exon genes or a LOC gene equal:" + str(len(varlist) - self.singleexongenevariants - len(self.LOClist)) + "\n")
                    
    # Write output of gene2nm to output file
    def write (self):
        hgvs = pd.Series(self.hgvs)
        hgvs = hgvs.values
        clas = pd.Series(self.clas)
        clas = clas.values
        df = pd.DataFrame(zip(hgvs, clas),  columns = ["HGVS", "Classifiaction"])
        df.to_csv(path_or_buf="/home/kevin/Documents/HGVS/HGVS_VLookup.csv", sep='\t')
            
                    
""" File written by Hgvsconvert is put through Mutalyser. Parsebed class then takes input file generated by Mutalyser. 
    Mutalyser file contains a col of NM accession variants and a col with the corresponding of NC accession variants as input.
    Extracts SNP variants for this list and converts them over to a vcf like format
"""    
class Parsebed(object):
    
    def __init__(self):
        self.input=''
        self.output=''
        self.mutalyserinput = pd.read_table("/home/kevin/Documents/HGVS/HGVS_VLookup.csv", header= 0)
        
    # Retrieve Class of corresponding HGVS formatted NM accession. This method is called from within gencoords method  
    def matchupmutalyser(self, inputvar):
        
        # Extract Classification of inputted variant from HGVS_VLookup.csv 
        retrieveclass = self.mutalyserinput.loc[self.mutalyserinput['HGVS'] == inputvar, 'Classification'].to_string(index=False).encode("utf-8")
        # Any where there is no match assign the variant as having NoRecord for the classification
        retrieveclass = re.sub('Series\(\[\]\, \)', 'NoRecord', retrieveclass)
        # list of strings to exclude from class column
        exclusionlist = ["Unknown Significance", "Unknown significance", "Unlikely to be pathogenic", "Probably not pathogenic", "Likely to be pathogenic", "Clearly pathogenic", "Clearly Pathogenic", "tract", "Likely to be Pathogenic", "Clearly p...", "Clearly not pathogenic", "Unknown..."]
        exclusions = "|".join(exclusionlist)
        # Generate a list of classes removing those where the classification relates to any of the terms in the exclusionlist
        retrieveclass = re.sub(exclusions, "", retrieveclass)
        # Stanadardise Not real assignment to the term "NotReal"
        notreal = ["Not real", "not real", "Never Real", "Never real"]
        notreallist = "|".join(notreal)
        retrieveclass = re.sub(notreallist, "NotReal", retrieveclass)
        # Stanadardise Not classified assignment to the term "NotClassified"
        retrieveclass = re.sub("Not classified", "NotClassified", retrieveclass)
        retrieveclass = re.sub(r"\\|/", ",", retrieveclass)
        # Stanadardise not in ex26 assignment to the term "NotinEx26"
        retrieveclass = re.sub("not in ex26", "NotinEx26", retrieveclass)
        # Stanadardise class assignment to the term "Class"
        retrieveclass = re.sub("class", "Class", retrieveclass)
        # Standardise not done to "NotDone"
        retrieveclass = re.sub("not done", "NotDone", retrieveclass)
        retrieveclass = re.sub(" ", "", retrieveclass)
        return retrieveclass
        
    
        
                    
   
    def gencoords(self):
        # Input table of mutalyser outputs
        self.inputfrommutalyser = pd.read_table(self.input, header= 0)
        
        # Take the chromosomal variant column from the mutalyser file (file used was /home/kevin/Documents/HGVS/Recognised_by_mutalyser.csv)
        genomichgvs = self.inputfrommutalyser["ChromosomalVariant"]
        
        # Initialise vcf columns chr, coord, ref and alt
        chr = []
        coord = []
        ref = []
        alt = []
        varlist = [] # list of variants to be passed to matchupmutalyser function
        info = []
        # Iterate over the list of chromosomal coords defined in chromosomal variant column from Mutalyser output
        for hgvs in genomichgvs:
            #  If the coord contains a > character indicative of it being a SNP
            if ">" in hgvs:
                splitval = re.split('\>|\.',hgvs) # split on > or dot (.)
                chr.append(int(re.split('00',splitval[0])[-1].lstrip("0"))) # Take the first element of the splitval output eg 'NC_000023' and split it on 00 and append the last element which is the chromosome number to the list chr (strip off leading 0 from single digit values)
                #Genomic coordinate
                coord.append(re.split('[A-Za-z]',splitval[2])[0]) # Take the second element of splitval output eg. '149783173G' and split it on alphabet character and take the first element resulting in just the genomic coordinate being returned and append this to the list coord 
                ref.append(hgvs[-3]) # Take the third last element from the genomichgvs row output which corresponds to the Ref base
                alt.append(splitval[-1]) # Take the last element from the splitval output which corresponds to the Alt base
                # Extract the value from the InputVariant column of the Mutalyser input file where the variant matches the entry in the ChromosomalVariant column
                inputvar = self.inputfrommutalyser.loc[self.inputfrommutalyser["ChromosomalVariant"] == hgvs, 'InputVariant'].to_string(index=False).encode("utf-8")
                # Input the Input variant into the matchupmutalyser method to retrieve the Classification for that variant from the HGVS output file from gene2nm method
                retrieveclass = self.matchupmutalyser(inputvar)
                info.append("PreviousClassification=" + retrieveclass) # Generate a list of Classifications
             
        id = ["."] * len(coord) # Fill id column to be inputed into vcf as .
        qual = ["."] * len(coord) # Fill qual column to be inputed into vcf as .
        filter = ["."] * len(coord) # Fill filter column to be inputed into vcf as .
        #format = ["."] * len(coord) # Fill format column to be inputed into vcf as .
        print len(coord)
        
        # Convert all list generated to be zipped into the vcf into Series objects
        chr = pd.Series(chr)
        coord = pd.Series(coord)
        id  = pd.Series(id)
        ref = pd.Series(ref)
        alt = pd.Series(alt)
        qual = pd.Series(qual)
        filter = pd.Series(filter)
        info = pd.Series(info)
        #format = pd.Series(format)
        
        
        
        chr = chr.values
        chr = ["X" if c == 23 else "Y" if c == 24 else c for c in chr]
        coord = coord.values
        id = id.values
        print len(id)
        ref = ref.values
        alt = alt.values
        qual = qual.values
        filter = filter.values
        info = info.values
        #format = format.values
        
        df2 = pd.DataFrame(zip(coord, id, ref, alt, qual, filter, info),  columns = ["POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"], index = chr)
        df2.index.name = "#CHROM"
        df2.to_csv(path_or_buf=self.output, sep='\t')
        # Generate header for VCF file
        pd.set_option('display.max_colwidth', -1)
        dfheader = pd.read_table("/home/kevin/Documents/HGVS/vcfheader.csv", header=0)
        commentheader = str(dfheader.to_string(index=False))
        formattedcomment = commentheader.replace('\n', '\\n').replace(' ', '') # replace newlines with escaped newlines which sed can then read as newlines
        temp=os.path.splitext(self.output)[0] + "_temp.txt" # create temporary file to write to
        out = open(temp, "w")
        subprocess.call("(head -n 1 %s && tail -n +2 %s | sort -k1,1V -k2,2n)" % (self.output, self.output), shell=True, executable="/bin/bash", stdout=out) # sort vcf according to genomic location
        subprocess.call("mv %s %s" % (temp, self.output), shell=True)
        #subprocess.call(["/bin/bash", "-c", "(head", "-n", "1", "%s", "&&", "tail", "-n", "+2", "%s", "|", "sort", "-k1,1V", "-k2,2n)", ">", "%s"] % (self.output, self.output, self.output))
        subprocess.call("sed -i '1s/^/%s\\n/' %s" % (formattedcomment, self.output), shell=True) # add comment line to top of vcf
        
        
        
        subprocess.call("bgzip -f %s; tabix -f %s" % (self.output, self.output + ".gz"), shell=True)
        
        #self.matchupmutalyser()
        
        
    

def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, help="supply path to vlookup input file for script")
    parser.add_argument("--output", type=str, help="supply path to output file for script")
    args = parser.parse_args()
    if args.input:
        parsebed = Parsebed() # Instantiate Parsebed class
        parsebed.input = args.input # Define input file for Parsebed class
        parsebed.output = args.output # Define output file for Parsebed class
        parsebed.gencoords() # Call gencoords method
     
if __name__ == '__main__':
    
    #h = Hgvsconvert() # instantiate Hgvsconvert class
    #print h
    #h.gene2nm() # run gene2nm function to convert the variants in the vlookup into HGVS format for analysis using Mutalyser
    #h.write()
    arguments()
    