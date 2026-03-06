#! /bin/env python3

# ----------- IMPORTS ---------------------------- ||
import sys
import re
import argparse
import pkgutil
from pathlib import Path
from .tree import Tree

# ----------- GLOBALS ---------------------------- ||
__author__ = "Brandon Pickett"
__copyright_owner__ = "Brandon Pickett"
__copyright_year__ = "2019"
__version__ = str(pkgutil.get_data(__package__, "VERSION").decode(encoding="UTF-8")).rstrip('\n')

# ----------- CLASSES ---------------------------- ||
class CalcScoreException(Exception):
	pass

# ---------- FUNCTIONS --------------------------- ||
def handleArgs():
	"""
	Defines possible arguments user can include in a call to calcScore to influence input or output, or see other details about this project/file

	Args:
		none

	Returns: 
		parsed arguments object
	"""
	# define the main argument parser
	parser = argparse.ArgumentParser(prog=Path(sys.argv[0]).name, add_help=False, allow_abbrev=True, 
									formatter_class=argparse.RawTextHelpFormatter, 
									description="Calculate the taxa resiliency for a provided tree. Please see the README file\n" 
									"distributed with this projects repository for further explanation.\n")

	# define argument groups
	input_group = parser.add_argument_group("Input Options", "These options affect how input files are found and/or interpreted.")
	output_group = parser.add_argument_group("Output Options", "These options affect which and how output files are generated. Specify an output" 
																" filename as \"\" to skip outputting it.")
	misc_group = parser.add_argument_group("Misc. Options", "")

	# 	define input group options
	#		main tree
	input_group.add_argument("-m", "--main-tree", dest="main_tree", metavar="tree.nwk", action="store", type=str, required=False,  default="data/mainTree/tree.nwk", 
						help="The main tree in Newick format for which you wish to determine the resiliency\n" 
						"against the removal of taxa." 
						" [data/mainTree/tree.nwk]\n \n")
	#		jackknife trees
	input_group.add_argument("-t", "--jackknife-tree", dest="jack_tree_dir", metavar="path/to/jackknife/tree/", action="store", type=str, required=False, default="data/jackknife/tree", 
						help="The directory in which the jackknife tree data exists. The directory should\n"
						"contain one subdirectory per taxon in the main tree. Each directory name should\n"
						"match the taxon name from the main tree file (case-sensitive!). In turn, each\n"
						"subdirectory should contain tree files in Newick format. Unless specified  by\n"
						"other options, the file names should be tree-${num}.nwk, where ${num} ranges from\n"
						"1-${rep}, where ${rep} is the number of jackknife replicates. In thoery, your\n"
						"files should fit in a nice range from 1-${rep}, but that is not strictly\n"
						"necessary. Using this option, the program will search for all files in the\n"
						"directory specified by -jt/--jackknife-tree that match the reg. exp.\n"
						"\"tree-\d+\.nwk\"." 
						" [data/jackknife/tree]\n \n")
	input_group.add_argument("-e", "--tree-ext", dest="jack_tree_fn_ext", metavar=".ext", action="store", type=str, required=False, default="nwk", 
						help="When using -jt, the assumed filename extension for the tree files is\n"
						"\"nwk\". If your files are named differently you may specify a different file\n"
						"extension here. This is useful if you used IQ-TREE to generate your trees as the\n"
						"output Newick trees are in files with the extension \".treefile\". Note that this\n"
						"changes only the end of the expected filename pattern. The beginning of the\n"
						"filename must still match this reg. exp: \"tree-\d+\". For more control over the\n"
						"filenames, you will need to use other options." 
						" [nwk]\n \n")
	input_group.add_argument("-f", "--jackknife-tree-fofn", dest="jack_tree_fofn", metavar="path/to/trees.tsv", action="store", type=str, required=False, default=None, 
						help="Instead of using -t and -e to provide the jackknifed trees, you may use this\n"
						"option to provide a list of filenames. The file must be in tab-separated value\n"
						"(tsv) format. Exactly two columns must exist on every line. No header/title line\n"
						"should exist. The first column should contain the taxon name (case-sensitive!).\n"
						"The second column should contain the path (absolute -or- relative to the\n"
						"directory from which the program is being run) to a jackknifed tree file in\n"
						"Newick format. This option overrides -t and -e for the jackknife tree files.\n"
						"The following is an example of how the file should look:\n" 
						"\tCignobilis\tdata/jackknife/tree/Cignobilis/tree-1.treefile\n" 
						"\t...\n" 
						"\tCignobilis\tdata/jackknife/tree/Cignobilis/tree-50.treefile\n" 
						"\tRmuscosa\tdata/jackknife/tree/Rmuscosa/tree-1.treefile\n" 
						"\t...\n" 
						"\tRmuscosa\tdata/jackknife/tree/Rmuscosa/tree-50.treefile\n" 
						"\tAglossodonta\tdata/jackknife/tree/Aglossodonta/tree-1.treefile\n" 
						"\t...\n" 
						"\tAglossodonta\tdata/jackknife/tree/Aglossodonta/tree-50.treefile\n" 
						"\tOclarkiistomias\tdata/jackknife/tree/Oclarkiistomias/tree-1.treefile\n" 
						"\t...\n" 
						"\tOclarkiistomias\tdata/jackknife/tree/Oclarkiistomias/tree-50.treefile\n" 
						"\t...\n \n")

	# 	define output group options
	#		newick format
	output_group.add_argument("-n", "--output-nwk", dest="output_nwk", metavar="out.nwk", action="store", type=str, required=False, default="out.nwk", 
						help="The output tree in Newick format with the taxon resiliency score in a comment\n"
						"after the branch length. To replace the branch length with the taxon resiliency\n"
						"score, specify the flag -b/--replace-branch-len. To replace the label of the\n" 
						"internal nodes with the taxa resiliency score, specify the flag\n" 
						"-s/--replace-label. Use of either -b or -s will prevent comments from being\n" 
						"output." 
						" [out.nwk]\n \n")
	output_group.add_argument("-b", "--replace-branch-len", dest="replace_branch_len", action="store_true", required=False, 
						help="Specify this flag to replace the branch length with the taxon resiliency score in\n" 
						"the output Newick tree. It probably does not make sense to use this in\n" 
						"conjunction with \"-s\", though it is functionally possible.\n \n")
	output_group.add_argument("-s", "--replace-label", dest="replace_internal_labels", action="store_true", required=False, 
						help="Specify this flag to replace the internal node's labels with the taxon resiliency\n" 
						"score in the output Newick tree. If \"-s\" seems inobvious for this functionality,\n" 
						"think \"--replace-sobriquet\". It probably does not make sense to use this in\n" 
						"conjunction with \"-b\", though it is functionally possible.\n \n")
	#		json format
	output_group.add_argument("-j", "--output-json", dest="output_json", metavar="out.json", action="store", type=str, required=False, default="out.json", 
						help="The output tree in json format. Each node will have a label (which may be empty),\n" 
						"branch length (which may be empty), and metadata with name/value pairs. The taxon\n" 
						"resiliency score will be assigned with the name \"taxon-resiliency\"." 
						" [out.json]\n \n")
	output_group.add_argument("-p", "--output-json-pretty", dest="output_json_pretty", metavar="out_pretty.json", action="store", type=str, required=False, default="out_pretty.json", 
						help="The output tree in json format. Each node will have a label (which may be empty),\n" 
						"branch length (which may be empty), and metadata with name/value pairs. The taxa\n" 
						"resiliency score will be assigned with the name \"taxa-resiliency\".\n" 
						"[out_pretty.json]\n \n")
	#		mermaid format
	output_group.add_argument("-w", "--output-mmd", dest="output_mmd", metavar="out.mmd", action="store", type=str, required=False, default="", 
						help="The output tree in mermaid format. Each node will have a label (which may be\n"
						"empty). The internal labels may be replaced with the taxa resiliency score with\n" 
						"the \"-s\" option, similar to the Newick output. This file is not generated by\n" 
						"default.\n \n") 

	# 	define misc. group options
	misc_group.add_argument("-c", "--cite", dest="display_citation", action="store_true", required=False,
							help="Describe how to cite this program.\n \n")
	misc_group.add_argument("-h", "--help", action="help", help="Show this help message and exit.\n \n")
	misc_group.add_argument("-l", "--license", dest="display_license", action="store_true", required=False,
							help="Display the license under which this software is released (MIT).\n \n")
	misc_group.add_argument("-v", "--version", dest="display_version", action="store_true", required=False,
							help="Display the program version and exit.\n \n")
	
	# parse the arguments
	args = parser.parse_args()

	# process arguments, as needed
	# 	display some information, quit immediately
	if args.display_citation or args.display_license or args.display_version:
		if args.display_citation:
			print("Citation:\n" 
				"Please include a link to this repository and the following citation:\n" 
				"\n" 
				"Powell GS, Pickett BD, Ridge PG, Martin GJ, Whiting MF, Bybee SM. Paper title. _journal_. Year.\n" 
				"\tVolume(Issue):pages.\n", 
				file=sys.stdout)
		if args.display_license:
			print("License:\n" 
				f"Copyright (c) {__copyright_year__} {__copyright_owner__}\n" 
				"\n" 
				"=================================================================\n" 
				"Permission is hereby granted, free of charge, to any person\n" 
				"obtaining a copy of this software and associated documentation\n" 
				"files (the \"Software\"), to deal in the Software without\n" 
				"restriction, including without limitation the rights to use,\n" 
				"copy, modify, merge, publish, distribute, sublicense, and/or sell\n" 
				"copies of the Software, and to permit persons to whom the\n" 
				"Software is furnished to do so, subject to the following\n" 
				"conditions:\n" 
				"\n" 
				"The above copyright notice and this permission notice shall be\n" 
				"included in all copies or substantial portions of the Software.\n" 
				"\n" 
				"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND,\n" 
				"EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES\n" 
				"OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND\n" 
				"NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT\n" 
				"HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,\n" 
				"WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING\n" 
				"FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR\n" 
				"THE USE OR OTHER DEALINGS IN THE SOFTWARE.\n" 
				"=================================================================\n", 
				file=sys.stdout)
		if args.display_version:
			print(f"Version: {__version__}\n", file=sys.stdout)
		sys.exit(0)
	else: # don't display some info, don't quit (immediately)
		# sanity check on input paths
		if args.jack_tree_fofn is not None: # fofn is provided
			p = Path(args.jack_tree_fofn)
			if not (p.exists() and p.is_file()): # exists and is file?
				raise CalcScoreException(f"ERROR: You provided -f \"{args.jack_tree_fofn}\", but it either did not exist or was not a regular file.")
		else: # fofn not specificed
			p = Path(args.jack_tree_dir)
			if not (p.exists() and p.is_dir()): # exists and is dir?
				raise CalcScoreException(f"ERROR: Problem with argument used for -t, \"{args.jack_tree_dir}\" either did not exist or was not a directory.")

		# sanity check on output files
		for ofn in (args.output_nwk, args.output_json, args.output_json_pretty):
			if ofn != '':
				p = Path(ofn)
				p = p.resolve()
				if p.exists():
					if not p.is_file():
						raise CalcScoreException(f"ERROR: Output file \"{ofn}\" exists and is not a regular file. It cannot be\noverwritten since it is not a regular file.")
				else:
					d = p.parent
					if d.exists():
						if not d.is_dir():
							raise CalcScoreException(f"ERROR: Output file \"{ofn}\" cannot be created because its theoretical parent\ndirectory already exists as a non-directory file.")
					else:
						d.mkdir(parents=True)


	# return the parsed arguments object
	return args

def createTreeFromNewickFile(filename, treename):
	"""
	Creates a Tree object from a Newick file

	Args:
		filename (str): location/name of the newick file
		treename (str): name for the Tree object, will appear when user prints the tree

	Returns:
		Tree: object that conceptually represents the data in the Newick file
	"""
	nwk = ''
	with open(filename, 'r') as ifd:
		for line in ifd:
			nwk += line.rstrip('\n')
	return Tree(newick=nwk, name=treename)

def getJackknifedTreesFileNames(tree_dir, tree_ext, trees_fofn):
	"""
	Extracts the names of every jackknifed tree file, either from a file listing taxa names with file names (trees_fofn) or a folder that contains all the tree files (tree_dir)
	
	Args:
		tree_dir (str): name of directory containing jackknifed tree files, may be None if user chose a different way of identifying file names
		tree_ext (str): the file extension at the end of each jackknifed file's name
		trees_fofn (str): file location for a file containing a list of taxa names and names of the jackknifed files, may be None

	Returns:
		dict: dictionary with a key for each taxon name pointing to a file name for the jackknifed file containing all taxa except that one
	"""
	taxa_x_fns = {}
	if trees_fofn is not None: # user specified the fofn
		with open(trees_fofn, 'r') as ifd:
			for line in ifd:
				fields = line.rstrip('\n').split('\t')
				taxon = fields[0]
				fn = fields[1]
				if taxon not in taxa_x_fns:
					taxa_x_fns[taxon] = []
				taxa_x_fns[taxon].append(fn)

	else: # user did not specify the fofn
		match_pattern = rf"^tree-\d+\.{re.escape(tree_ext)}$"
		# we assume tree_dir exists and is a directory (handled during handleArgs)
		d = Path(tree_dir)
		for sd in d.iterdir(): # search for sub directories (one level, assume one dir per taxa)
			if sd.is_dir(): # look only at dirs
				taxon = sd.name
				for f in sd.iterdir(): # search for files in the dir
					if f.is_file() and re.match(match_pattern, f.name) is not None: # look only at files. they must match f"tree-\d+.{tree_ext}"
						if taxon not in taxa_x_fns:
							taxa_x_fns[taxon] = []
						taxa_x_fns[taxon].append(str(f.resolve()))
	
	return taxa_x_fns

def generateReplicatesHistogram(reps):
	"""
	Shows user how many replicates are in each taxon. Every taxon should have the same number of replicates, so this is only shown to the user when that assumption is broken. 

	Args:
		reps (list[int]): list containing the number of reps in each taxon

	Returns:
		str: human-readable representation of the number of reps in each taxon (not labelled by taxon)

	"""
	# do the counting
	counts = {}
	for rep in reps:
		if rep not in counts:
			counts[rep] = 0
		counts[rep] += 1
	
	# do some formatting
	counts_max_width = len(str(max(list(counts.keys()))))
	left_pad_char = ' '
	sep = "|"
	tick_mark = '='
	line_fmt = f"{{c:{left_pad_char}{counts_max_width}}}{sep}{{t}} ({{f}})\n".format # c: count, t: ticks, f: freq

	# create the visual repr of the hist
	output = ""
	for count in sorted(list(counts.keys())):
		freq = counts[count]
		ticks = tick_mark * freq
		output += line_fmt(c=count, t=ticks, f=freq)
	
	return output

def validateAndResolveJackknifedTrees(taxa_x_fns, taxa):
	"""
	Makes sure that the dictionary of jackknifed trees matches the original tree, the taxa in the dictionary match each other, and all the required files exist

	Args:
		taxa_x_fns (dict[str,list[str]]): keys are taxa names, values are a list of jackknifed file names
		taxa (list[str]): names of taxa present in the original tree

	Returns:
		none, but errors if any requirements are not met
	"""
	# confirms that all taxa from the original tree exist in taxa_x_fns, errors otherwise
	taxa_set = frozenset(list(taxa_x_fns.keys()))
	for taxon in taxa:
		if taxon not in taxa_set:
			raise CalcScoreException("ERROR: One or more taxa from the original/main tree were not present in the\njackknifed trees dir.")

	# confirms that all taxa in taxa_x_fns exist in the original tree, errors otherwise
	taxa_set = frozenset(taxa)
	for taxon in taxa_x_fns.keys():
		if taxon not in taxa_set:
			raise CalcScoreException("ERROR: One or more taxa from the jackknifed trees dir were not present in the\noriginal/main tree.")
	
	# confirms that all taxa in taxa_x_fns have the same number of replicates
	reps = [len(taxa_x_fns[taxon]) for taxon in taxa_set]
	first_len = reps[0]
	if not all(rep == first_len for rep in reps):
		hist = generateReplicatesHistogram(reps)
		raise CalcScoreException(f"ERROR: All taxa should have the same number of replicates. Here is the replicate\nhistogram:\n{hist}\n)")
	
	# confirms that all the required files exist
	for taxon in taxa_x_fns.keys():
		fns = taxa_x_fns[taxon]
		for i,fn in enumerate(fns):
			try:
				p = Path(fn)
				p = p.resolve()
				if not (p.exists() and p.is_file()):
					raise CalcScoreException(f"ERROR: \"{fn}\" does not exist or is not a regular file (taxon: {taxon}).")
				taxa_x_fns[taxon][i] = str(p)
			except:
				raise CalcScoreException(f"ERROR: while testing if \"{fn}\" was valid for {taxon}, failed to create or\nresolve Path object.")

def sortJackknifedTrees(taxa_x_fns):
	"""
	Sorts the list of jackknifed file names associated with each taxon based on the number on the filename

	Args:
		taxa_x_fns (dict[str,list[str]]): keys are taxa names, values are a list of jackknifed file names

	Returns:
		(dict[str,list[str]]): sorted version of taxa_x_fns
	"""
	for taxon in taxa_x_fns.keys():
		taxa_x_fns[taxon].sort(key=lambda x: int(re.sub(r"^.*(\d+).*$", r"\1", Path(x).stem)))

def buildJackknifedTreesFromFiles(taxa_x_fns):
	"""
	Creates trees from file names

	Args:
		taxa_x_fns (dict[str,list[str]]): keys are taxa names, values are a list of jackknifed file names

	Returns:
		(dict[str, list[Tree]]): keys are taxa names, values are a list of jackknifed trees
	"""
	taxa_x_trees = {}
	for taxon in taxa_x_fns.keys():
		if taxon not in taxa_x_trees:
			taxa_x_trees[taxon] = []
		fns = taxa_x_fns[taxon]
		for i,fn in enumerate(fns):
			try:
				taxa_x_trees[taxon].append(createTreeFromNewickFile(fn, f"{taxon}-{i}"))
			except:
				raise CalcScoreException(f"ERROR: failed to create Tree object from newick tree file \"{fn}\" (taxon: {taxon})")

	return taxa_x_trees

# ------------- MAIN ----------------------------- ||
def main():
	# handle the arguments
	args = handleArgs()

	# read in the main tree
	mt = createTreeFromNewickFile(args.main_tree, "main")

	# get a list of taxa
	taxa = sorted(mt.getLeafLabels())

	# obtain list of jackknifed tree files mapped to taxa names
	taxa_x_fns = getJackknifedTreesFileNames(args.jack_tree_dir, args.jack_tree_fn_ext, args.jack_tree_fofn)
	
	# validate and resolve jackknifed trees (paths, not tree objects)
	validateAndResolveJackknifedTrees(taxa_x_fns, taxa) # side-effect (arg1), no change (arg2), no return

	# sort jackknifed trees (individually sort each path list) 
	# HOW DO YOU SAY THIS MORE PROFESSIONALLY (arguably not necessary, but it feels nice)
	sortJackknifedTrees(taxa_x_fns) # side-effect, no return

	# build jackknifed trees from file
	taxa_x_trees = buildJackknifedTreesFromFiles(taxa_x_fns)

	# compare
	mt.scoreResiliency(taxa_x_trees) # changes mt, but not taxa_x_trees

	# generate output
	#	json
	#		ugly
	if args.output_json:
		with open(args.output_json, 'w') as ofd:
			ofd.write(mt.getJson())

	#		pretty
	if args.output_json_pretty:
		with open(args.output_json_pretty, 'w') as ofd:
			ofd.write(mt.getPrettyJson())

	# 	mmd
	if args.output_mmd:
		with open(args.output_mmd, 'w') as ofd:
			ofd.write(mt.getMermaid(replace_internal=args.replace_internal_labels))

	#	nwk
	if args.output_nwk:
		with open(args.output_nwk, 'w') as ofd:
			if args.replace_branch_len or args.replace_internal_labels:
				if args.replace_branch_len:
					mt.replaceBranchLenWithOtherValue("taxa-resiliency")
				if args.replace_internal_labels:
					mt.replaceInternalLabelsWithOtherValue("taxa-resiliency")
				ofd.write(mt.getNewick())
			else:
				ofd.write(mt.getNewickWithCommentedMetadata())
	
if __name__ == "__main__":
	main()
