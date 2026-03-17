#! /bin/env python3

__author__ = "Brandon Pickett"

# ----------- IMPORTS ---------------------------- ||
import sys
from .node import Node,MalformedNewickTree

# ---------- FUNCTIONS --------------------------- ||

# ----------- CLASSES ---------------------------- ||
class Tree:
	"""
	A class representing a phylogenetic tree.

	This class wraps a single root Node and provides methods for
	accessing, analyzing, and exporting the tree.

	Attributes:
		root (Node): The root node of the tree.
		name (str): Optional name for the tree,
	"""

	# constructor(s)
	def __init__(self, newick="", name=""):
		# "normal" "public" member fields
		self.root = Node()
		self.name = name

		self.__initializeNodes__(newick)

	# "normal" "public" member functions
	def getLeafLabels(self):
		return self.root.getLeafLabels()
	
	def getEachSubTreeLeafLabelSets(self):
		return self.root.getEachSubTreeLeafLabelSets()

	def getEachSubTreeLeafLabelSetStrs(self):
		return self.root.getEachSubTreeLeafLabelSetStrs()

	def containsSubtreeBasedOnSetOfLeafLabels(self, node):
		return self.root.containsSubtreeBasedOnSetOfLeafLabels(node)

	def containsSubtreeBasedOnPreFetchedSetOfLeafLabels(self, leaf_labels): # leaf_labels must be sorted
			return self.root.containsSubtreeBasedOnPreFetchedSetOfLeafLabels(leaf_labels)

	def generateNodesViaDepthFirstTraversal(self):
		yield from self.root.generateNodesViaDepthFirstTraversal()

	def scoreResiliency(self, taxa_x_trees):
		for node in self.generateNodesViaDepthFirstTraversal():
			node.scoreResiliency(taxa_x_trees)
		self.root.scoreResiliency(taxa_x_trees, meaningful=False) # force root score to 0
	
	def replaceBranchLenWithOtherValue(self, meta_key):
		for node in self.generateNodesViaDepthFirstTraversal():
			node.replaceBranchLenWithOtherValue(meta_key)

	def replaceInternalLabelsWithOtherValue(self, meta_key):
		for node in self.generateNodesViaDepthFirstTraversal():
			node.replaceInternalLabelWithOtherValue(meta_key)

	def getNewick(self):
		return self.root.getNewick() + ";\n"
	
	def getNewickWithCommentedMetadata(self):
		return self.root.getNewickWithCommentedMetadata() + ";\n"

	def getJson(self):
		return '{"name":"' + self.name + '","root":' + self.root.getJson() + '}'
	
	def getPrettyJson(self):
		#return json.loads(self.getJson())
		return '{\n\t"name": "' + self.name + '",\n\t"root":\n' + self.root.getPrettyJson(indent=2) + '\n}\n'
	
	def getAscii(self, prefix="", children_prefix=""):
		return self.root.getAscii(prefix=prefix, children_prefix=children_prefix)
	
	def getMermaid(self, replace_internal=False):
		#return "graph LR:\n" + self.root.getMermaid()
		mmd = "graph LR\n"
		leaf_ids = []
		long_ids = {}
		gfa = ''
		i = 0
		for node in self.generateNodesViaDepthFirstTraversal():
			long_ids[str(id(node))] = str(i)
			node_label = node.label if node.label else " "
			if node.isLeaf():
				mmd += '\t' + str(i) + "[" + node_label + "]\n"
				leaf_ids.append(i)
			else:
				if replace_internal and "taxa-resiliency" in node.metadata:
					node_label = str(node.metadata["taxa-resiliency"])
				mmd += '\t' + str(i) + "((" + node_label + "))\n"
				for child in node.children:
					gfa += '\t' + str(i) + " --- " + long_ids[str(id(child))] + '\n'
			i += 1
		mmd += gfa
		mmd += "\tclassDef nodes fill:#eee,stroke:#fff,stroke-width:0px,color:black;\n"
		mmd += "\tclassDef leaf-nodes fill:#fff;\n"
		mmd += "\tclass " + ','.join(map(str, range(0, i))) + " nodes;\n"
		mmd += "\tclass " + ','.join(map(str, leaf_ids)) + " leaf-nodes;\n"
		return mmd

	# "private" member functions
	def __initializeNodes__(self, newick):
		newick = self.__removeNewickComments__(newick).rstrip()
		index = self.root.initializeNode(newick)

		if index < len(newick):
			if newick[index] == ';':
				index += 1
				newick = newick[index:]
				if newick and not newick.isspace():
					raise MalformedNewickTree(f"Reached end of tree and found semi-colon, but found 1 or more non-space characters after semi-colon.")
			else:
				raise MalformedNewickTree(f"Reached end of tree and expected a semi-colon, but found {newick[index]} instead.")
		else:
			raise MalformedNewickTree("Reached end of tree without encountering a semi-colon")

	def __removeNewickComments__(self, newick):
		keep = ""
		i = 0
		while i < len(newick):
			if newick[i] == "[":
				found_r_bracket = False
				while i < len(newick):
					if newick[i] == ']':
						found_r_bracket = True
						break
					i += 1

				if not found_r_bracket:
					raise MalformedNewickTree("comment had no ending ']'")

				i += 1 # go to char after ']' (or if at len already, it won't matter that its past len)

			elif newick[i] == '"' or newick[i] == "'":
				# this section allows for '[' and ']' inside quoted labels
				found_end_quote = False
				q = newick[i]
				keep += q
				i += 1 # go to char after opening quote
				while i < len(newick):
					keep += newick[i]
					if newick[i] == q:
						found_end_quote = True
						break
					i += 1

				if not found_end_quote:
					raise MalformedNewickTree(f"quoted label had no ending quote ({q})")

				i += 1 # go to char after end quote (or if at len already, it won't matter that its past len)
			else:
				keep += newick[i]
				i += 1
		return keep

	# comparison operators
	#def __lt__(self, other):
	#	return self.root.__lt__(other.root)
	
	def isEqualBasedOnSetOfLeafLabels(self, other):
		return self.root.isEqualBasedOnSetOfLeafLabels(other.root)

	# make str(some_node) meaningful
	def __str__(self):
		return f'{{ name: "{self.name}", root: {str(self.root)} }}'
	
	# make print(some_node) meaningful
	def __repr__(self):
		return "Tree: " + self.__str__()
	

# ------------- MAIN ----------------------------- ||
if __name__ == "__main__":
	"""
	This module is intended to be imported as a library for phylogenetic tree operations,
	and running it directly is not currently supported.

	Attempting to run will procuce an error message and exit.
	"""
	sys.stderr.write(
		"ERROR: This is a module, it is meant to be imported, not run directly!\n"
		)
	sys.exit(1)
