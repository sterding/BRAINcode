import argparse;
#import subprocess;
import pysam;
import vcf;
import calendar;
import time;
import random;
import re;

## python script to count reads per allele in a specific sample by parsing the samtools mpileup output
## modified based on https://github.com/secastel/allelecounter
## TODO: allow to work for multiple samples
## another option is to use ASEReadCounter in GATK

def main():
	#Arguments passed 
	parser = argparse.ArgumentParser()
	parser.add_argument("--vcf", help="Genotype VCF from SNPs.")
	parser.add_argument("--mpileup", help="mpileup Output.")
	parser.add_argument("--sample", help="Sample name")
	parser.add_argument("--samplelist", help="Samples list")
	parser.add_argument("--min_cov", type=int, help="Minimum total coverage to report read counts")
	parser.add_argument("--min_baseq", type=int, help="Minimum base quality to count reads")
	parser.add_argument("--min_mapq", type=int, help="Minimum map quality to count reads")
	parser.add_argument("--o", help="Output file")
	args = parser.parse_args()
	
	# deal with STAR alignments having mapq of 255
	# the highest mapq that can be encoded by samtools mpileup is 93
	# so if min_mapq > 93 set it to 93
	if args.min_mapq > 93: args.min_mapq = 93;
	
	if args.min_cov == "": args.min_cov = 5;
	if args.min_baseq == "": args.min_baseq = 0;
	if args.min_mapq == "": args.min_mapq = 0;
	
	start_timestamp = calendar.timegm(time.gmtime());
	
	# fine the rank of sample in the sample list
	nLine = 1
	inputFile = open(args.samplelist,'r')
	corpusLines = inputFile.readlines()
	for line_i, line in enumerate(corpusLines, 1): 
		if re.findall( args.sample, line ):
			nLine = line_i
			SAMPLE = re.sub(r'.*run_output/(.*'+args.sample+').*',r'\1',line).rstrip()
	inputFile.close()
	SAMPLE1 = SAMPLE
	if(SAMPLE == 'HC_UWA479'): SAMPLE1 = 'HC_UW479'
	if(SAMPLE == 'HC_NZ-H83'): SAMPLE1 = 'HC_NZ-H183'
	#print nLine, SAMPLE, SAMPLE1;
	
	#print("Loading the VCF ...");
	vcf_reader = vcf.Reader(open(args.vcf,'r'));
	
	#out_stream = open(args.o, "w");
	#out_stream.write("contig	position	variantID	refAllele	altAllele	refCount	altCount	totalCount	lowMAPQDepth	lowBaseQDepth	rawDepth	otherBases\n");
	#print "contig position variantID refAllele altAllele refCount altCount totalCount lowMAPQDepth lowBaseQDepth rawDepth otherBases";
	
	# 2 process the mpileup result
	#print("Processing pileup...");
	# B go through the pileup result line by line
	for line in open(args.mpileup, 'r'):
		cols = line.replace("\n","").split("\t");
		#chr	pos	REF	count	reads
		chr = cols[0];
		chr = chr.replace("chr","")
		pos = int(cols[1]);
		ref = "";
		alt = "";
		rsid = ".";
		# first retrieve the VCF record for this SNP
		#print chr, str(pos)
		records = vcf_reader.fetch(chr,pos-1,pos);
		
		for record in records:
			#print record, pos, record.genotype(args.sample)['GT']
			if int(record.POS) == pos:
				# only want bi-allelic SNPS
				if len(record.ALT) == 1 and len(record.REF) == 1 and len(record.ALT[0]) == 1:
					#only what heterozygous sites
					genotype = record.genotype(SAMPLE1)['GT']
					if genotype.count("1") == 1 and genotype.count("0") == 1:
						ref = str(record.REF);
						alt = str(record.ALT[0]);
						if record.ID != None:
							rsid = str(record.ID);
		#print rsid;			
		# if the site is biallelic and we have ref / alt genotype data go ahead and do the read counts
		if ref != "" and alt != "":
			reads = list(cols[4*nLine]);
			reads = [x for x in reads if (x != "$" and x != "^" and x != "~" and x != "\"" and x != "!")]
			block = 0;
			out_reads = [];
			for read in reads:
				if block == -1:
					block = int(read);
				elif read == "+" or read == "-":
					block = -1;
				elif block > 0:
					block -= 1;
				elif block == 0:
					out_reads.append(read);
			reads = out_reads;
			baseqs = list(cols[4*nLine+1]);
			mapqs = list(cols[4*nLine+2]);
			
			low_baseq = 0;
			low_mapq = 0;
			ref_count = 0;
			alt_count = 0;
			other_count = 0;
			raw_depth = len([x for x in reads if (x != ">" and x != "<")]);
			
			for read,baseq,mapq in zip(reads,baseqs,mapqs):
				if read != "<" and read != ">":
					basequal = ord(baseq)-33;
					mapqual = ord(mapq)-33;
					if basequal >= args.min_baseq and mapqual >= args.min_mapq:
						# count this base
						if read == "." or read == ",":
							ref_count += 1;
						elif read == alt.upper() or read == alt.lower():
							alt_count += 1;
						else:
							other_count += 1;
					if basequal < args.min_baseq:
						low_baseq += 1;
					if mapqual < args.min_mapq:
						low_mapq += 1;
			totalCount = ref_count+alt_count;
			
			if totalCount > args.min_cov:
				#out_stream.write("\t".join([cols[0],cols[1],rsid,ref,alt,str(ref_count),str(alt_count),str(totalCount),str(low_mapq),str(low_baseq),str(raw_depth),str(other_count),"\n"]));
				print SAMPLE,cols[0],cols[1],rsid,ref,alt,str(ref_count),str(alt_count),str(totalCount),str(low_mapq),str(low_baseq),str(raw_depth),str(other_count);
	
	#out_stream.close();
	
	stop_timestamp = calendar.timegm(time.gmtime());
	#print("Total time to complete: %d seconds"%(stop_timestamp-start_timestamp));
	
if __name__ == "__main__":
	main();
