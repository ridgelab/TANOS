#! /bin/env python3

__author__ = "Brandon Pickett"

# ----------- IMPORTS ---------------------------- ||
import sys
from .node import Node,MalformedNewickTree

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

	# ---------- CONSTRUCTOR --------------------- ||
	def __init__(self, newick="", name=""):
		"""
		Initialize a Tree from a Newick string.

		Args:
			newick (str, optional): Newick-formatted string representing the tree.
			name (str, optional): Name of the tree.
		"""
		self.root = Node()
		self.name = name
		self.__initializeNodes__(newick)

	# ---------- PUBLIC MEMBER FUNCTIONS --------- ||
	def getLeafLabels(self):
		"""
		Get all leaf labels in the tree.

		Returns:
			list[str]: List of leaf labels.
		"""
		return self.root.getLeafLabels()
	
	def getEachSubTreeLeafLabelSets(self):
		"""
		Get sets of leaf labels for each subtree.

		Returns:
			list[list[str]]: Nested lists of leaf labels for each subtree.
		"""
		return self.root.getEachSubTreeLeafLabelSets()

	def getEachSubTreeLeafLabelSetStrs(self):
		"""
		Get string representations of leaf label sets for each subtree.

		Returns:
			list[str]: List of concatonated leaf labels for each subtree.
		"""
		return self.root.getEachSubTreeLeafLabelSetStrs()

	def containsSubtreeBasedOnSetOfLeafLabels(self, node):
		"""
		Check whether the tree contains a subtree matching the leaf labels of another tree.

		Args:
			node (Node): Node whose leaf labels define the subtree of interest.
		
		Returns:
			bool: True if a matching subtree exits, False otherwise. 
		"""
		return self.root.containsSubtreeBasedOnSetOfLeafLabels(node)

	def containsSubtreeBasedOnPreFetchedSetOfLeafLabels(self, leaf_labels):
			"""
			Check whether the tree contains a subtree matching a pre-fetched set of leaf labels.

			Args:
				leaf_labels (list[str]): Sorted list of leaf labels to search for.

			Returns:
				bool: True if a matching subtree exists, False otherwise.
			"""
			return self.root.containsSubtreeBasedOnPreFetchedSetOfLeafLabels(leaf_labels)

	def generateNodesViaDepthFirstTraversal(self):
		"""
		Generate all nodes in the tree via depth-first traversal.

		Yields:
			Node: Each node in depth-first order.
		"""
		yield from self.root.generateNodesViaDepthFirstTraversal()

	def scoreResiliency(self, taxa_x_trees):
		"""
		Compute and store taxa-resiliency scores for all nodes in tree.

		Args:
			taxa_x_trees (dict[str, list[Node]]): Mapping from taxa to lists of trees.

		Sets self.metadata["taxa-resiliency"] in metadata of each node.
		"""
		for node in self.generateNodesViaDepthFirstTraversal():
			node.scoreResiliency(taxa_x_trees)
		self.root.scoreResiliency(taxa_x_trees, meaningful=False) # force root score to 0
	
	def replaceBranchLenWithOtherValue(self, meta_key):
		"""
		Replace branch lengths in all nodes with another metadata value.

		Args:
			meta_key (str): Metadata key whose value replaces branch_length.
		"""
		for node in self.generateNodesViaDepthFirstTraversal():
			node.replaceBranchLenWithOtherValue(meta_key)

	def replaceInternalLabelsWithOtherValue(self, meta_key):
		"""
		Replace labels of internal nodes with a metadata value.

		Args:
			meta_key (str): Metadata key to use for replacing internal node labels.
		"""
		for node in self.generateNodesViaDepthFirstTraversal():
			node.replaceInternalLabelWithOtherValue(meta_key)

	def getNewick(self):
		"""
		Export the tree in Newick format.

		Returns:
			str: Newick string representing the tree (end with';').
		"""
		return self.root.getNewick() + ";\n"
	
	def getNewickWithCommentedMetadata(self):
		"""Export the tree in Newick format including metadata comments.

		Returns:
			str: Newick string with metadata annotations (ends with ';').
		"""
		return self.root.getNewickWithCommentedMetadata() + ";\n"

	def getJson(self):
		"""
		Export the tree in compact JSON format.

		Returns:
			str: JSON string representing the tree.
		"""
		return '{"name":"' + self.name + '","root":' + self.root.getJson() + '}'
	
	def getPrettyJson(self):
		"""
		Export the tree in a human-readable, indented JSON format.

		Returns:
			str: Pretty JSON string.
		"""
		return '{\n\t"name": "' + self.name + '",\n\t"root":\n' + self.root.getPrettyJson(indent=2) + '\n}\n'
	
	def getAscii(self, prefix="", children_prefix=""):
		"""
		Generate an ASCII representation of the tree.

		Args:
			prefix (str, optional): Prefix for the current node line.
			children_prefix (str, optional): Prefix for child lines.

		Returns:
			str: ASCII art of the tree.
		"""
		return self.root.getAscii(prefix=prefix, children_prefix=children_prefix)
	
	def getMermaid(self, replace_internal=False):
		"""
		Generate a Mermaid diagram representation of the tree.
		
		Args:
			replace_internal (bool, optional): Whether to replace internal node labels
				with their taxa-resiliency values.

		Returns:
			str: Mermaid graph syntax representing the tree.
		"""
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

	# ---------- PRIVATE MEMBER FUNCTIONS ------------------ ||
	def __initializeNodes__(self, newick):
		"""
		Initialize the tree nodes from a Newick string, remove comments, validate.

		Args:
			newick (str): Newick-formatted string.
		"""
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
		"""
		Remove bracketed comments and preserve quoted labels in a Newick string.

		Args:
			newick (str): Raw Newick string.

		Returns:
			str: Cleaned Newick string without comments.
		"""
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

				i += 1 # go to char after ']'

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

				i += 1 # go to char after end quote
			else:
				keep += newick[i]
				i += 1
		return keep

	def isEqualBasedOnSetOfLeafLabels(self, other):
		"""
		Check if two trees are equal based on their leaf labels.

		Args:
			other (Tree): Another tree to compare.
		
		Returns:
			bool: True if both trees have identical leaf labels sets.
		"""
		return self.root.isEqualBasedOnSetOfLeafLabels(other.root)

	def __str__(self):
		"""
		String representation of the tree.

		Returns:
			str: Name and root node.
		"""
		return f'{{ name: "{self.name}", root: {str(self.root)} }}'
	
	def __repr__(self):
		"""
		Official representation of the tree.

		Returns:
			str: Debug string.
		"""
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
