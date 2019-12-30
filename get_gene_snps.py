#!/usr/bin/python

import re, os, time, sys

#path = sys.argv[1]
multi_vcf = sys.argv[1]
fasta_file = sys.argv[2]

exceptions = open("exceptions.txt", 'w')

path = ''
prefix = ''
for file in sys.argv:
	if re.search("\.fasta", file):
		path_file_l = re.split('/', file)
		file = path_file_l.pop(-1)
		path = '/'.join(path_file_l)
		prefix = re.split("\.", file)[0]

'''
"Prokka is loaded by the shell script so you need to either run this \
from the shell script or load Prokka prior to running this python script. \
The locations of the variants are relative to the gene not the genome. \
To get the genomic position of '+' oriented genes, simply add the SNP position \
to the beginning of the gene. If the gene is '-' oriented, subtract from the \
end of the gene and if you are checking the results for a revcomped sequence \
remember to switch 'G' and 'C' as well as 'A' and 'T'.
'''

print("Script written by John Miller incorporating useful suggestions by James Pettengill. \n")

def print_time():
	print (time.asctime( time.localtime(time.time())))	

def test_defline_len(fasta_file):
	'''
	Finds length of first defline from fasta file containing contigs.
	'''
	fa = open(fasta_file, 'r')
	for line in fa:
		if re.search('^>', line):
			length = len(line)
			break
	return length
			
def shorten_deflines(fasta_file):
	'''
	Shortens deflines.
	'''
	fa = open(fasta_file, 'r')
	new_fasta_file = re.sub('.fasta', '_short_deflines.fasta', fasta_file)
	new_fa = open(new_fasta_file, 'w')
	def_count = 0
	defline2short_defline = {}
	for line in fa:
		if re.search('^>', line):
			line = line.rstrip("\n")
			def_count += 1
			new_defline = line[:15]+'_'+str(def_count)
			id = re.sub('>', '', line)
			new_id = id[:14]+'_'+str(def_count)
			defline2short_defline[id] = new_id 
			new_fa.write(new_defline+"\n")
		else:
			new_fa.write(line)
	return defline2short_defline, new_fasta_file

def shorten_vcf_seq_ids(vcf_file, defline2short_defline):
	'''
	Replaces ids with those corresponding to shortened deflines from fasta file.
	'''
	new_vcf_file = re.sub('.vcf', '_short_ids.vcf', vcf_file)
	vcf = open(vcf_file, 'r')
	new_vcf = open(new_vcf_file, 'w')
	for line in vcf:
		if not re.search("^#", line):
			linel = re.split("\t", line)
			linel[0] = defline2short_defline[linel[0]]
			new_line = "\t".join(linel)
			new_vcf.write(new_line)
		else:
			new_vcf.write(line)
	return new_vcf_file

def gff2gene_pos(gff_file):
	'''
	Takes a gff3 file from Prokka file (Prokka does not add the 3 to the name).
	Generates dict mapping locus tags to gene [beginnings, endings].
	Also generates a dict mapping locus_tags to gene/product names.
	'''
	gff = open(gff_file, 'r')
	id2start_end = {}
	id2description = {}
	id2orientation = {}
	locus2dna_seg = {}
	for line in gff:
#		print line
		line = line.rstrip("\n")
		if not re.search("^#", line):
			linel = re.split("\t", line)
			begin_end = [linel[3], linel[4]]
			attributes = re.split(";", linel[8])
			id = ''
			product = ''
			if re.search('^ID=', attributes[0]):
				id = re.sub('ID=', '', attributes[0])
				product = re.sub('product=', '', attributes[-1])
#			print begin_end, id, product
			id2start_end[id] = begin_end
			id2description[id] = product
			id2orientation[id] = linel[6]
			locus2dna_seg[id] = linel[0]
		if re.search("FASTA", line):
			break
	gff.close()
	return(id2start_end, id2description, id2orientation, locus2dna_seg)

def parse_multi_vcf(multi_multi_vcf):
	'''
	Takes multiVCF.
	Produces dict snp_info keys=snp pos on chromosome/contig and values=[various info] and 
	a dict with key=col_label value=index to allow accessing information about a particular isolate.
	'''
	vcf = open(multi_multi_vcf, 'r')
	snp_info = {}
	label2index = {}
	isolateList = []
	for line in vcf:
		line = line.rstrip("\n")
		if re.search("^#CHROM", line):
			labels = re.split("\t", line)
			count = 0
			for col_lab in labels:
				col_lab = re.sub('a_', '', col_lab)
				label2index[col_lab] = count
				isolateList.append(col_lab)
				count += 1
		elif not re.search("^#", line):
			linel = re.split("\t", line)
			chr_pos = linel[0]+'_'+linel[1]
#			print chr_pos
			snp_info[chr_pos] = linel
	vcf.close()
	return snp_info, label2index, isolateList[9:]

def get_locus2snp_pos(locus_tag2pos, locus_tag2dna_seq_id, snp_dict, col_label2index):#, cfsan2bcw):	# isolate, 
	'''
	Associates locus_tags from PROKKA with SNPs from SAMTOOLs/VarScan.
	Returns a dict with locus tags as keys and a list of variant positions for each value.
	'''
	locus2pos = locus_tag2pos.keys()
	locus2snp_pos = {}
	for locus in locus2pos:						# just a list of keys (locus_tags) from gene2pos.
		chr = locus_tag2dna_seq_id[locus]		# find chromosome or contig on which locus resides.	####
		gene_start = int(locus_tag2pos[locus][0])
		gene_end = int(locus_tag2pos[locus][1])
		snp_pos = []
		for pos in range(gene_start, gene_end):	# try each position in each gene once.
			seq_id = chr+'_'+str(pos)
###			print seq_id
			'''
			if snp_dict[seq_id]:
				new_seq_id = snp_dict[seq_id][0]+'_'+snp_dict[seq_id][1]
				snp_pos.append(new_seq_id) #this is only the number
			'''
			try:
				new_seq_id = snp_dict[seq_id][0]+'_'+snp_dict[seq_id][1]
				snp_pos.append(new_seq_id) #this is only the number
			except:
#				print "no snp at position", pos
				pass
#		print len(snp_pos)
		locus2snp_pos[locus] = snp_pos
	return locus2snp_pos


def get_isolate_specific_snps(locus2snp_pos, snp_dict, col_label2index, isolate):
	'''
	creates dict with locus_tags as keys and list of snp_positions as values.
	snp_pos is a list.
	'''
	isolate_specific_snp_pos2nt = {}
	isolate_specific_locus2snp_pos = {}
	for locus in sorted(locus2snp_pos.keys()):	# loop thru locus tags
#		print locus
		snp_positions = []
#		print locus2snp_pos[locus]
		for snp_pos in locus2snp_pos[locus]:	# loop thru variant positions within one locus. 
			try:
				if re.search("\d", snp_dict[snp_pos][col_label2index[isolate]][0]):
					var_index = int(snp_dict[snp_pos][col_label2index[isolate]][0])-1
					if var_index >= 0:									# test whether isolate of interest has a SNP at this position.
						variants = re.split(',', snp_dict[snp_pos][4])	# the variant column in case there is more than one variant.
						snp_positions.append(snp_pos)
						isolate_specific_snp_pos2nt[snp_pos] = variants[var_index]
#						print 'monkey looking for bugs.'
			except:
#				print "cannot find ", snp_dict[snp_pos][col_label2index[isolate]][0]
				exceptions.write("cannot find snp at "+snp_pos+" for isolate "+isolate+'.')
		isolate_specific_locus2snp_pos[locus] = snp_positions
	return isolate_specific_snp_pos2nt, isolate_specific_locus2snp_pos

def orfs2dict(fasta):
	'''
	Puts contents of one fasta file into a dictionary.
	The keys are the deflines. The values are the sequences.
	'''
	fa = open(fasta, 'r')
	seq = ''
	defline = ''
	def2seq = {}
	for line in fa:
		line = line.rstrip("\n")
		if re.search("^>", line):
			def2seq[defline] = seq	# get seq prior to changing to new defline.
			defline = re.sub('>', '', line)
			defline = re.split("\s", defline)[0]
			seq = ''
		else:
			seq = seq + line.upper()
	def2seq[defline] = seq		# get the last one
	return def2seq

def find_alternate_seqs(locusTag2snpPos, locusTag2genePos, locusTags2orfs, locusTag2orientations, snpPos2ntState):
	'''
	Returns a dict with locus tags as keys and a list of two sequences the old one and the new one for each value.
	Also returns a dict with SNP positions (value) for each locus (key).
	'''
	locus2seqs = {}									# key is snp_pos. value is [old_seq, new_seq]
	locusTag2indelPos = {}
	for loc in locusTag2snpPos:						# loop thru locus tags
		start = int(locusTag2genePos[loc][0])-1
		end = int(locusTag2genePos[loc][1])
		gene_len = len(locusTags2orfs[loc])
		pos_syn = []
		positions = []
		indel_positions = []
		new_gene_seq = list(locusTags2orfs[loc])
		for snp_pos in locusTag2snpPos[loc]:				# loop thru SNPs associated with each locus tag
			pos = re.split('_', snp_pos)[-1]
			snp_pos_i = int(pos)
			if start <= snp_pos_i and snp_pos_i <= end:		# make sure variant is within gene of interest.
				snp_in_gene = 0								# the snp position relative to the gene start and end.
				if locusTag2orientations[loc] == '-':		# REVCOMPED SEQUENCES
					snp_in_gene = end - snp_pos_i
					# need to switch to complement nucleotide
					if snpPos2ntState[snp_pos] == 'G':
						complement_nt = re.sub('G', 'C', snpPos2ntState[snp_pos])
					elif snpPos2ntState[snp_pos] == 'C':
						complement_nt = re.sub('C', 'G', snpPos2ntState[snp_pos])
					elif snpPos2ntState[snp_pos] == 'A':
						complement_nt = re.sub('A', 'T', snpPos2ntState[snp_pos])
					elif snpPos2ntState[snp_pos] == 'T':
						complement_nt = re.sub('T', 'A', snpPos2ntState[snp_pos])
					elif snpPos2ntState[snp_pos] == '-':
						indel_positions.append([snp_pos, snp_in_gene])
					new_gene_seq[snp_in_gene] = complement_nt		# changing one element of a list.
				elif locusTag2orientations[loc] == '+':				# FORWARD SEQUENCES
					snp_in_gene = snp_pos_i - start - 1
					new_gene_seq[snp_in_gene] = snpPos2ntState[snp_pos]		# changing one nucleotide.
					if snpPos2ntState[snp_pos] == '-':
						indel_positions.append([snp_pos, snp_in_gene])
		new_gene_seq = ''.join(new_gene_seq)
		if len(indel_positions) > 0:
			locusTag2indelPos[loc] = indel_positions
		len1 = len(new_gene_seq)
		new_gene_seq = re.sub('-', '', new_gene_seq)		# is this in the correct place
		len2 = len(new_gene_seq)
		locus2seqs[loc] = [locusTags2orfs[loc], new_gene_seq]	# old_seq, new_seq
	return locus2seqs, locusTag2indelPos

def test_dN_dS(locusTag2seqs):
	'''
	Input locusTag2seqs is a dict with locusTags as keys and list [old_gene, new_gene] as values.
	Returns a list of lists in which element one is the SNP position and element two is information
	about the codons and amino acids
	'''
	trans_tab = {'TTT':'Phe', 'TCT':'Ser', 'TAT':'Tyr', 'TGT':'Cys', 'TTC':'Phe', 'TCC':'Ser', 'TAC':'Tyr', 'TGC':'Cys', 'TTA':'Leu', 'TCA':'Ser', 'TAA':'Ter', 'TGA':'Ter', 'TTG':'Leu', 'TCG':'Ser', 'TAG':'Ter', 'TGG':'Trp', 'CTT':'Leu', 'CCT':'Pro', 'CAT':'His', 'CGT':'Arg', 'CTC':'Leu', 'CCC':'Pro', 'CAC':'His', 'CGC':'Arg', 'CTA':'Leu', 'CCA':'Pro', 'CAA':'Gln', 'CGA':'Arg', 'CTG':'Leu', 'CCG':'Pro', 'CAG':'Gln', 'CGG':'Arg', 'ATT':'Ile', 'ACT':'Thr', 'AAT':'Asn', 'AGT':'Ser', 'ATC':'Ile', 'ACC':'Thr', 'AAC':'Asn', 'AGC':'Ser', 'ATA':'Ile', 'ACA':'Thr', 'AAA':'Lys', 'AGA':'Arg', 'ATG':'Met', 'ACG':'Thr', 'AAG':'Lys', 'AGG':'Arg', 'GTT':'Val', 'GCT':'Ala', 'GAT':'Asp', 'GGT':'Gly', 'GTC':'Val', 'GCC':'Ala', 'GAC':'Asp', 'GGC':'Gly', 'GTA':'Val', 'GCA':'Ala', 'GAA':'Glu', 'GGA':'Gly', 'GTG':'Val', 'GCG':'Ala', 'GAG':'Glu', 'GGG':'Gly', }
	locus_snp_info = []
	for loc in sorted(locusTag2seqs):	# loop thru locus tags.
#		print loc
		seq_pair = locusTag2seqs[loc]
		seq_len = 0
		# in the event of different length sequences find the shortest.
		if len(seq_pair[0]) >= len(seq_pair[1]):
			seq_len = len(seq_pair[1])
		elif len(seq_pair[0]) < len(seq_pair[1]):
			seq_len = len(seq_pair[0])
		# find the longest seq that is a multiple of 3 so the rest of the block will run.
		remainder = seq_len % 3
		if remainder != 0:
			seq_len = seq_len - remainder
		codon_start = 0
		for i in range(3, seq_len+3, 3):	# loop thru codons
			if re.search("N", seq_pair[1][codon_start:i]):
				# put together all possible codons and amino acids resulting from N
				A = re.sub('N', 'A', seq_pair[1][codon_start:i])
				C = re.sub('N', 'C', seq_pair[1][codon_start:i])
				G = re.sub('N', 'G', seq_pair[1][codon_start:i])
				T = re.sub('N', 'T', seq_pair[1][codon_start:i])
				refcodon = ':'.join([seq_pair[0][codon_start:i], trans_tab[seq_pair[0][codon_start:i]]]) 
				Acodon = ':'.join([A, trans_tab[A]]) 
				Ccodon = ':'.join([C, trans_tab[C]]) 
				Gcodon = ':'.join([G, trans_tab[G]])
				Tcodon = ':'.join([T, trans_tab[T]])
				line = "\t".join([str(pos), refcodon, ','.join([Acodon, Ccodon, Gcodon, Tcodon])])
				locus_snp_info.append([loc, line])
			elif seq_pair[0][codon_start:i] != seq_pair[1][codon_start:i]:
				if trans_tab[seq_pair[0][codon_start:i]] != trans_tab[seq_pair[1][codon_start:i]]:	# identify non-synonymous variants
					pos = ''
					for j in range(0, 3):	# find exact position of non-synonymous snp
						if seq_pair[0][codon_start+j] != seq_pair[1][codon_start+j]:
							pos = codon_start+j
							refcodon = ':'.join([seq_pair[0][codon_start:i], trans_tab[seq_pair[0][codon_start:i]]])
							varcodon = ':'.join([seq_pair[1][codon_start:i], trans_tab[seq_pair[1][codon_start:i]]])
							line = "\t".join([str(pos), refcodon, "\t", varcodon])
							locus_snp_info.append([loc, line])
			codon_start = i
	return locus_snp_info

def print_one_isolate(locus2snp_info, locus_tag2dna_seq_id, locus_tag2gene_pos, locus_tag2orientation, locus_tag2name, isolate, path):
	'''
	Prints the results for each isolate.
	'''
	results = path+'/gene_snp_results'
	path_file = results+'/'+isolate+'_gene_snps.txt'
	out = open(path_file, 'w')
	header = "\t".join(['locus_tag', 'contig/chromosome', 'gene_start', 'gene_end', 'snp_pos_in_gene', 'ref_codon:AA', 'variant_codon:AA', 'orientation', 'description'])
	out.write(header+"\n")
	for l in locus2snp_info:
		to_print = "\t".join([l[0], locus_tag2dna_seq_id[l[0]], locus_tag2gene_pos[l[0]][0], locus_tag2gene_pos[l[0]][1], l[1], locus_tag2orientation[l[0]], locus_tag2name[l[0]]])
		out.write(to_print+"\n")
#		print l[0], l[1]

def consolidate_codons(locus2snp_info, locus2snp_info_all_isolates):
	'''
	takes codon and amino acids variants from different strains and puts them into one list
	'''
	for locus in locus2snp_info:
		snp_info = re.split("\t+", locus[1])
		pos = snp_info.pop(0)
		locus_pos = locus[0]+'\t'+pos
		if locus_pos not in locus2snp_info_all_isolates.keys():
			locus2snp_info_all_isolates[locus_pos] = snp_info
		else:
			# don't need to test for snp_info[0] because if the key is present so is the reference state.
			if snp_info[1] not in locus2snp_info_all_isolates[locus_pos]:
				locus2snp_info_all_isolates[locus_pos].append(snp_info[1])	#adding additional codon:aa
	return locus2snp_info_all_isolates
	
def print_summary(locus2snp_info_all_isolates, locus_tag2dna_seq_id, locus_tag2gene_pos, locus_tag2orientation, locus_tag2name, path):
	'''
	Prints summary of variant information for all isolates without regard to isolate.
	'''
	out = open(path+'/summary.txt', 'w')
	header = "\t".join(['locus_tag', 'contig/chromosome', 'gene_start', 'gene_end', 'snp_pos_in_gene', 'codons:amino_acids', 'orientation', 'description'])
	out.write(header+"\n")
	for locus_pos in sorted(locus2snp_info_all_isolates.keys()):
		locus_posl = re.split("\t", locus_pos)
		locus = locus_posl[0]
		pos = locus_posl[1]
		print locus, pos
		to_print = "\t".join([locus, locus_tag2dna_seq_id[locus], locus_tag2gene_pos[locus][0], locus_tag2gene_pos[locus][1], pos, ','.join(locus2snp_info_all_isolates[locus_pos]), locus_tag2orientation[locus], locus_tag2name[locus]])
		out.write(to_print+"\n")

# ONLY DO ONCE.
print_time()
length = test_defline_len(fasta_file)
if length > 20:	# deciding whether to shorted deflines for PROKKA.
	defline2short, fasta_file = shorten_deflines(fasta_file)
	multi_vcf = shorten_vcf_seq_ids(multi_vcf, defline2short)

prokka_out = path+'/prokka_output'
if not os.path.exists(prokka_out):
	os.system('mkdir '+prokka_out)
	os.system('cp '+fasta_file+' '+prokka_out)

ffn_file = prokka_out+"/"+prefix+".ffn"
gff_file = prokka_out+"/"+prefix+".gff"

os.system('prokka --outdir '+prokka_out+' --force -norrna -notrna --prefix '+prefix+' '+fasta_file)
snp_dict, col_label2index, isolate_list = parse_multi_vcf(multi_vcf)
locus_tag2gene_pos, locus_tag2name, locus_tag2orientation, locus_tag2dna_seq_id = gff2gene_pos(gff_file)
locus2orfs = orfs2dict(ffn_file)
locus2snp_pos = get_locus2snp_pos(locus_tag2gene_pos, locus_tag2dna_seq_id, snp_dict, col_label2index)
print_time()

results = path+'/gene_snp_results'
if not os.path.exists(results):
	os.system('mkdir '+results)

# DO FOR EVERY STRAIN.
locus2snp_info_all_isolates = {}
for isolate in isolate_list:
	print isolate
	isolate_specific_snp_pos2nt, isolate_specific_locus2snp_pos = get_isolate_specific_snps(locus2snp_pos, snp_dict, col_label2index, isolate)
	locus2seqs, locus2indel_pos = find_alternate_seqs(isolate_specific_locus2snp_pos, locus_tag2gene_pos, locus2orfs, locus_tag2orientation, isolate_specific_snp_pos2nt)
	locus2snp_info = test_dN_dS(locus2seqs)
	locus2snp_info_all_isolates = consolidate_codons(locus2snp_info, locus2snp_info_all_isolates)
	print_one_isolate(locus2snp_info, locus_tag2dna_seq_id, locus_tag2gene_pos, locus_tag2orientation, locus_tag2name, isolate, path)
	print print_time()

print_summary(locus2snp_info_all_isolates, locus_tag2dna_seq_id, locus_tag2gene_pos, locus_tag2orientation, locus_tag2name, path)

print_time()