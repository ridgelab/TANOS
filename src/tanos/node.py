#! /bin/env python3

__author__ = "Brandon Pickett"

# ----------- IMPORTS ---------------------------- ||
import sys

# ---------- FUNCTIONS --------------------------- ||

# ----------- CLASSES ---------------------------- ||
class MalformedNewickTree(Exception):
	pass

class Node:
	
	# class level variables (define once, not for every instance of the class)
	# for example:
	# PI = 3.14

	# constructor(s)
	def __init__(self):
		# "normal" "public" member fields

		# 	children: immediate children Nodes
		self.children = []
		#	label: used to store a taxon name for leaf Nodes or clade name for internal nodes
		self.label = ""
		#	metadata: Here we enable the storage
		#	of multiple name value pairs.
		#	The values can be numbers or strings.
		#	To get the standard branch length, simply use
		#	self.metadata["branch_length"].
		self.metadata = {}
		#if kwargs is not None:
		#	for k,v in kwargs.iteritems():
		#		self.metadata[str(k)] = str(v)
	
	def initializeNode(self, newick, index=0):
		# 1- conceptual str.lstrip() beginning at position index 
		index = self.__consumeNewickWhitespace__(newick, index=index)

		# 2- process children (recursively if neeeded)
		while index < len(newick) and ( newick[index] == '(' or newick[index] == ',' ):
			index += 1 # get past the recursive signal (left paren or comma)
			self.children.append(Node())
			index = self.children[-1].initializeNode(newick, index=index)
			index = self.__consumeNewickWhitespace__(newick, index=index)

		# 3- process label (quote or unquoted)
		if index < len(newick):
			if newick[index] == '"' or newick[index] == "'":
				# quoted label present
				self.label, index = self.__getQuotedNewickLabel__(newick, index)
				index = self.__consumeNewickWhitespace__(newick, index=index)
			elif not newick[index] in "),:;":
				# unquoted label present
				self.label, index = self.__getUnquotedNewickLabel__(newick, index)
				index = self.__consumeNewickWhitespace__(newick, index=index)
		else:
			raise MalformedNewickTree("Reached end of tree (while about to process a potential label) without encountering a semi-colon")

		# 4- process branch length
		if index < len(newick):
			index = self.__getFromNewickAndPossiblySetBranchLength__(newick, index)

		else:
			raise MalformedNewickTree("Reached end of tree (while about to process a potential branch length) without encountering a semi-colon")

		# 5- process end of node (possibly of entire tree)
		if index < len(newick):
			if newick[index] == ')':
				index += 1 # get past the )
				return index # position of char after )

			elif newick[index] == ',' or newick[index] == ';':
				return index # position of comma or semi-colon

			else:
				raise MalformedNewickTree("Reached what should have been the end of a node (and possibly the entire tree), but found a character other than a right paren, comma, or semi-colon")
					
		else:
			raise MalformedNewickTree("Reached end of tree (while expecting either a right paren, comma, or semi-colon) without encountering a semi-colon")

	def __consumeNewickWhitespace__(self, newick, index=0):
		while index < len(newick) and newick[index].isspace():
			index += 1
		return index

	def __getQuotedNewickLabel__(self, newick, index=0):
		if ( len(newick) - index ) > 1:
			if newick[index] == '"' or newick[index] == "'":
				starting_index = index # deep copy
				q = newick[index]
				index += 1 # get past opening quote
				while index < len(newick):
					if newick[index] == q:
						index + 1 # get past closing quote
						return newick[starting_index:index], index
					index += 1

				raise MalformedNewickTree("quoted label had no closing quote")
			else:
				raise MalformedNewickTree("attempted to extract the end position of a quoted label from a newick string, but the starting position provided was not a quote")
		else:
			raise MalformedNewickTree("attempted to extract the end position of a quoted label, but the newick string had <2 characters remaining")

	def __getUnquotedNewickLabel__(self, newick, index=0):
		starting_index = index # deep copy
		index += 1 # get past first label character (which is already validated)
		while index < len(newick) and not ( newick[index].isspace() or newick[index] in "),:;" ):
			index += 1
		return newick[starting_index:index], index
	
	def __getFromNewickAndPossiblySetBranchLength__(self, newick, index=0):
		if newick[index] == ':':
			# branch length present
			index += 1 # get past the ':'
			index = self.__consumeNewickWhitespace__(newick, index=index)
			if index < len(newick):
				branch_length_str = ''
				while index < len(newick) and not ( newick[index].isspace() or newick[index] in '),;' ):
					branch_length_str += newick[index]
					index += 1
				# index is either past end (we didn't find a semi-colon),
				# or it is a right paren, comma, semi-colon, or space
				if index < len(newick):
					try:
						branch_length = float(branch_length_str) if '.' in branch_length_str else int(branch_length_str)
						self.metadata["branch_length"] = branch_length
					except ValueError:
						raise MalformedNewickTree(f"Found branch length ({branch_length_str}), but it was not an int or a float.")
				else:
					raise MalformedNewickTree(f"Found branch length ({branch_length_str}) after ':', but reached end of tree without encountering a semi-colon")
			else:
				raise MalformedNewickTree("Expected branch length after ':', but found nothing.")

			# index is either a right paren, comma, semi-colon, or space
			# now let's skip past spaces to ensure it is one of the other three
			index = self.__consumeNewickWhitespace__(newick, index=index)

		return index

	def isEqualBasedOnSetOfLeafLabels(self, other):
		# it should already not have duplicates
		return self.isEqualBasedOnPreFetchedSetOfLeafLabels(sorted(other.getLeafLabels()))

	def isEqualBasedOnPreFetchedSetOfLeafLabels(self, leaf_labels): # leaf_labels must be sorted
		return sorted(self.getLeafLabels()) == leaf_labels
	
	# "normal" "public" member functions
	def isLeaf(self):
		return not self.hasChildren()
	
	def hasChildren(self): # true if 1+ children
		return bool(len(self.children))
	
	def hasGrandChildren(self): # true if 1+ grandchildren
		for child in self.children:
			if child.hasChildren():
				return True
		return False
	
	def getLeafLabels(self): # returns list leaf labels, e.g., [ "A", "B", "C", ... ]
		leaves = []
		if self.isLeaf():
			leaves.append(self.label)
		else:
			for child in self.children:
				leaves.extend(child.getLeafLabels())
		#assert sorted(list(frozen_set(leaves))) == sorted(leaves) # assuming the tree does _not_ have leaves w/ identical labels
		return leaves
	
	def getEachSubTreeLeafLabelSets(self): # returns list of lists of leaf labels for each subtree, e.g., [ ["A"], ["B"], ["A", "B"], ["C"], ["A", "B", "C"], ... ]
		leaves = []
		if self.isLeaf():
			leaves.append([self.label])
		else:
			for child in self.children:
				leaves.extend(child.getEachSubTreeLeafLabelSets())
			leaves.append(self.getLeafLabels())
		return leaves

	def getEachSubTreeLeafLabelSetStrs(self): # returns list of leaf labels for each subtree, e.g., [ "A", "B", "AB", "C", "ABC", ... ]
		leaves = []
		if self.isLeaf():
			leaves.append(self.label)
		else:
			for child in self.children:
				leaves.extend(child.getEachSubTreeLeafLabelSetStrs())
			leaves.append(''.join(self.getLeafLabels()))
		return leaves
	
	def containsSubtreeBasedOnSetOfLeafLabels(self, node):
		subtree_of_interest = sorted(node.getLeafLabels())
		return self.containsSubtreeBasedOnPreFetchedSetOfLeafLabels(subtree_of_interest)
		
	def containsSubtreeBasedOnPreFetchedSetOfLeafLabels(self, leaf_labels): # leaf_labels must be sorted
		for child in self.children:
			if child.containsSubtreeBasedOnPreFetchedSetOfLeafLabels(leaf_labels):
				return True
		return self.isEqualBasedOnPreFetchedSetOfLeafLabels(leaf_labels)
	
	def generateNodesViaDepthFirstTraversal(self):
		for child in self.children:
			yield from child.generateNodesViaDepthFirstTraversal()
		yield self
	
	def scoreResiliency(self, taxa_x_trees, meaningful=True):
		score = 0
		if meaningful: # root has no meaningful resiliency score
			if self.hasGrandChildren():
				taxa = sorted(self.getLeafLabels())
				total_possible = 0
				count = 0
				for i in range(0, len(taxa), 1):
					excluded_taxon = taxa[i]
					included_taxa = taxa[:i] + taxa[i+1:]
					total_possible += len(taxa_x_trees[excluded_taxon])
					for tree in taxa_x_trees[excluded_taxon]:
						if tree.containsSubtreeBasedOnPreFetchedSetOfLeafLabels(included_taxa):
							count += 1
				score = float(count) / total_possible
			else:
				score = 1 # nodes that have _no_ grandchildren have no meaningful resiliency score
			if score == 1 or count == 0:
				score = int(score)
		self.metadata["taxa-resiliency"] = score

	def replaceBranchLenWithOtherValue(self, meta_key):
		if meta_key in self.metadata:
			self.metadata["branch_length"] = self.metadata[meta_key]
	
	def replaceInternalLabelWithOtherValue(self, meta_key):
		if self.hasChildren():
			if meta_key in self.metadata:
				self.label = str(self.metadata[meta_key])
			else:
				self.label = ''
	
	def getNewick(self):
		nwk = []
		if len(self.children):
			nwk.append('(')
			for i,child in enumerate(self.children):
				if i > 0:
					nwk.append(',')
				nwk.append(child.getNewick())
			nwk.append(')')
		if self.label:
			nwk.append(self.label)
		if "branch_length" in self.metadata:
			nwk.append(':')
			nwk.append(str(self.metadata["branch_length"]))
		return ''.join(nwk)

	def getNewickWithCommentedMetadata(self):
		nwk = []
		# recursively get children
		if len(self.children):
			nwk.append('(')
			for i,child in enumerate(self.children):
				if i > 0:
					nwk.append(',')
				nwk.append(child.getNewickWithCommentedMetadata())
			nwk.append(')')
		# get this node's label (may be empty string, which is fine)
		if self.label:
			nwk.append(self.label)
		# get the branch length, if one is stored
		if "branch_length" in self.metadata:
			nwk.append(':')
			nwk.append(str(self.metadata["branch_length"]))
		# get the other metadata (everything except branch length), if present, in a comment
		if len(self.metadata):
			meta_keys = sorted(list(self.metadata.keys()))
			try:
				meta_keys.remove("branch_length")
			except ValueError:
				pass
			if len(self.metadata):
				nwk.append("[")
				for i,meta_key in enumerate(meta_keys):
					if i > 0:
						nwk.append(',')
					meta_value = self.metadata[meta_key]
					value_quote = ''
					try:
						float(meta_value)
					except ValueError:
						value_quote = '"'
					nwk.append(f"\"{meta_key}\"={value_quote}{meta_value}{value_quote}")
				nwk.append("]")
		return ''.join(nwk)

	def getJson(self):
		j = [f'{{"label":"{self.label}","metadata":{{']
		for i,k in enumerate(sorted(list(self.metadata.keys()))):
			if i > 0:
				j.append(',')
			j.append(f'"{k}":')
			if type(self.metadata[k]) is str:
				j.append(f'"{self.metadata[k]}"')
			else:
				j.append(f'{self.metadata[k]}')
		j.append('},"children":[')
		for i,child in enumerate(self.children):
			if i > 0:
				j.append(',')
			j.append(child.getJson())
		j.append(']}')
		return ''.join(j)

	def getPrettyJson(self, indent=0):
		tabs = ''.join(['\t'] * indent)
		
		# label
		j = [f'{tabs}{{\n{tabs}\t"label": "{self.label}",\n{tabs}\t"metadata":']

		# metadata
		if len(self.metadata):
			j.append(f'\n{tabs}\t\t{{\n')
			for i,k in enumerate(sorted(list(self.metadata.keys()))):
				if i > 0:
					j.append(',\n')
				j.append(f'{tabs}\t\t\t"{k}":')
				if type(self.metadata[k]) is str:
					j.append(f'"{self.metadata[k]}"')
				else:
					j.append(f'{self.metadata[k]}')
			j.append(f'\n{tabs}\t\t')
		else:
			j.append(" {")

		j.append(f'}},\n{tabs}\t"children":')

		# children
		if len(self.children):
			j.append(f'\n{tabs}\t\t[\n')
			for i,child in enumerate(self.children):
				if i > 0:
					j.append(',\n')
				j.append(child.getPrettyJson(indent=indent+3))
			j.append(f'\n{tabs}\t\t')
		else:
			j.append(" [")
		j.append(f']\n{tabs}}}')

		# return
		return ''.join(j)
	
	def getAscii(self, prefix="", children_prefix=""):
		output = prefix + self.label + '\n'
		if self.hasChildren():
			for i in range(0, len(self.children) - 1, 1):
				#next_prefix = children_prefix + "|-- "
				#next_children_prefix = children_prefix + "|   "
				#output += self.children[i].getAscii(prefix=next_prefix, children_prefix=next_children_prefix)
				output += self.children[i].getAscii(prefix=f"{children_prefix}|-- ", children_prefix=f"{children_prefix}|   ")
			#output += self.children[-1].getAscii(prefix="", children_prefix="")
			output += self.children[-1].getAscii(prefix=f"{children_prefix}'-- ", children_prefix=f"{children_prefix}    ")
		return output

	# make str(some_node) meaningful
	def __str__(self):
		return f'{{ label: "{self.label}", metadata: {str(self.metadata)}, children: {len(self.children)} }}'
	
	# make print(some_node) meaningful
	def __repr__(self):
		return "Node: " + self.__str__()
	

# ------------- MAIN ----------------------------- ||
if __name__ == "__main__":
	sys.stderr.write("ERROR: This is a module, it is meant to be imported -- not run directly!\n")
	sys.exit(1)

