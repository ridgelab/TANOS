#! /bin/env python3

__author__ = "Brandon Pickett"

# ----------- IMPORTS ---------------------------- ||
import sys

# ----------- CLASSES ---------------------------- ||
class MalformedNewickTree(Exception):
	"""
	Exception raised when a Newick-formatted tree string is malformed.
	"""
	pass

class Node:
	"""
		Represent a node in a phylogenetic tree.

		Each node may be either a leaf (taxon) or an internal clade.
		Nodes store references to their immediate children, a label, 
		and a dictionary of metadata.

		Attributes:
			children (list[Node]): Nodes immediately following current node.
			label (str): Taxon or clade name.
			metadata (dict[str, Any]): Dictionary of associated branch lengths or other annotations.

			Notes:
				Metadata values may include both numeric and string data. 
				A standard branch length can be accessed via self.metadata["branch_length"].
	"""

	# ----------- CONSTRUCTORS ------------------- ||
	def __init__(self):
		"""
		Initialize an empty Node in a phylogenetic tree.
		
		The node is created with no children, an empty label, and an empty metadata dictionary.
		"""
		self.children = []
		self.label = ""
		self.metadata = {}
	
	def initializeNode(self, newick, index=0):
		"""
		Recursively parse a Newick string and populate Node attributes.

		Args:
			newick (str): Newick-formatted tree string.
			index (int, optional): Current parsing position in string.

		Returns:
			int: Updated index after parsing node.

		"""
		# 1 - Conceptual str.lstrip() beginning at position index.
		index = self.__consumeNewickWhitespace__(newick, index=index)

		# 2 - Process children (recursively if neeeded).
		while index < len(newick) and ( newick[index] == '(' or newick[index] == ',' ):
			index += 1 # get past the recursive signal (left paren or comma)
			self.children.append(Node())
			index = self.children[-1].initializeNode(newick, index=index)
			index = self.__consumeNewickWhitespace__(newick, index=index)

		# 3 - Process label (quote or unquoted).
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

		# 4 - Process branch length.
		if index < len(newick):
			index = self.__getFromNewickAndPossiblySetBranchLength__(newick, index)

		else:
			raise MalformedNewickTree("Reached end of tree (while about to process a potential branch length) without encountering a semi-colon")

		# 5 - Process end of node (possibly of entire tree).
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
		
	# ----------- PRIVATE HELPER FUNCTIONS  --------------- ||
	def __consumeNewickWhitespace__(self, newick, index=0):
		"""
		Skip index past all whitespace in Newick string.

		Args:
			newick (str): Newick string.
			index (int): Current index.

		Returns:
			int: Index of first character after whitespace.
		"""
		while index < len(newick) and newick[index].isspace():
			index += 1
		return index

	def __getQuotedNewickLabel__(self, newick, index=0):
		"""
		Get a quoted label from a Newick string.

		Args:
			newick (str): Newick string.
			index (int): Starting index that points to a quote.

		Returns:
			tuple[str, int]: Label and updated index.
		"""
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
		"""
		Get an unquoted label from a Newick string.

		Args:
			newick (str): Newick string.
			index (int): Starting index.

		Returns:
			tuple[str, int]: Extracted label and updated index.
		"""
		starting_index = index
		index += 1 # skip first label character (already validated)
		while index < len(newick) and not ( newick[index].isspace() or newick[index] in "),:;" ):
			index += 1
		return newick[starting_index:index], index
	
	def __getFromNewickAndPossiblySetBranchLength__(self, newick, index=0):
		"""
		Parse and store branch length metadata if present.

		Args:
			newick (str): Newick string.
			index (int): Current index.

		Returns:
			int: Updated index.
		"""
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
				# or it is a right parenthesis, comma, semi-colon, or space
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

			# index is either a right parenthesis, comma, semi-colon, or space
			# now let's skip past spaces to ensure it is one of the other three
			index = self.__consumeNewickWhitespace__(newick, index=index)

		return index
	
	# ----------- PUBLIC MEMBER FUNCTIONS ----------- ||
	def isEqualBasedOnSetOfLeafLabels(self, other):
		""""
		Compare self to another node to determine if they have the same leaf labels.

		This method assumes that leaf labels within each tree have no duplicates.

		Args:
			other (Node): Another tree node.

		Returns:
			bool: True if both nodes have the same leaf labels, False otherwise.
		"""
		return self.isEqualBasedOnPreFetchedSetOfLeafLabels(sorted(other.getLeafLabels()))

	def isEqualBasedOnPreFetchedSetOfLeafLabels(self, leaf_labels):
		"""
		Compare self's leaf labels to pre-fetached, then sorted, leaf labels.

		Args:
			leaf_labels (list[str]): A sorted list of leaf labels.

		Returns:
			bool: True if this node's leaf labels match, False otherwise.
		"""
		return sorted(self.getLeafLabels()) == leaf_labels
	
	def isLeaf(self):
		"""
		Determine if node is a leaf.

		Returns:
			bool: True if leaf (node has no children), False otherwise. 
		"""
		return not self.hasChildren()
	
	def hasChildren(self):
		"""
		Check if node has children.

		Returns:
			bool: True if this node has children, False otherwise.
		"""
		return bool(len(self.children))
	
	def hasGrandChildren(self):
		"""
		Check if node has grandchildren:

		Returns:
			bool: True if any child node has children, False otherwise.
		"""
		for child in self.children:
			if child.hasChildren():
				return True
		return False
	
	def getLeafLabels(self):
		"""
		Get all leaf labels in subtree.

		Returns:
			list[str]: List of leaf labels.
		"""
		leaves = []
		if self.isLeaf():
			leaves.append(self.label)
		else:
			for child in self.children:
				leaves.extend(child.getLeafLabels())
		# assuming the tree does _not_ have leaves w/ identical labels
		return leaves
	
	def getEachSubTreeLeafLabelSets(self):
		"""
		Get leaf labels for all subtrees.

		Returns: 
			list[list[str]]: Nested list of leaf labels where each inner list is the labels of one subtree.
		"""
		leaves = []
		if self.isLeaf():
			leaves.append([self.label])
		else:
			for child in self.children:
				leaves.extend(child.getEachSubTreeLeafLabelSets())
			leaves.append(self.getLeafLabels())
		return leaves

	def getEachSubTreeLeafLabelSetStrs(self):
		"""
		Get leaf labels for all subtrees as concatonated strings.

		Returns:
			list[str]: List of leaf labels where each element in the list is a concatonated 
				string of all subtree leaf labels.
		"""
		leaves = []
		if self.isLeaf():
			leaves.append(self.label)
		else:
			for child in self.children:
				leaves.extend(child.getEachSubTreeLeafLabelSetStrs())
			leaves.append(''.join(self.getLeafLabels()))
		return leaves
	
	def containsSubtreeBasedOnSetOfLeafLabels(self, node):
		"""
		Check whether this node contains subtree with matching leaf labels as `node.`

		Args:
			node (Node): Node with leaf labels of interest.

		Returns:
			bool: True if a matching subtree exists, False otherwise. 
		"""
		subtree_of_interest = sorted(node.getLeafLabels())
		return self.containsSubtreeBasedOnPreFetchedSetOfLeafLabels(subtree_of_interest)
		
	def containsSubtreeBasedOnPreFetchedSetOfLeafLabels(self, leaf_labels):
		"""
		Check whether this node contains a subtree with a pre-featch set of leaf labels.

		Args:
			leaf_labels (list[str]): Sorted list of leaf labels of subtree.
		
		Returns:
			bool: True if a matching subtree exists, False otherwise. 
		"""
		for child in self.children:
			if child.containsSubtreeBasedOnPreFetchedSetOfLeafLabels(leaf_labels):
				return True
		return self.isEqualBasedOnPreFetchedSetOfLeafLabels(leaf_labels)
	
	def generateNodesViaDepthFirstTraversal(self):
		"""
		Generator to yield all nodes in subtree using *depth-first* traversal.

		Yields:
			Node: Each node in the subtree.
		"""
		for child in self.children:
			yield from child.generateNodesViaDepthFirstTraversal()
		yield self
	
	def scoreResiliency(self, taxa_x_trees, meaningful=True):
		"""
		Compute and store a taxa-resiliency score in the node's metadata.

		Where a taxa-resiliency score is a measure of how robust the subtree is to exclusion of individual taxa.

		Args:
			taxa_x_trees (dict[str, list[Node]]): Mapping of taxon names to lists of trees.
			meaningful (bool, optional): Whether the root of a given subtree has a meaningful resiliency score.
				Where nodes with no granchildren do not have a meaningful score, and are assigned a deafult score of 1.

		Sets self.metadata["taxa-resilency"] with a float score between 0 and 1.
		"""
		score = 0
		if meaningful:
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
		"""
		Replace branch length in metadata with another value.

		Args:
			meta_key (str): Metadata key to update the value of self.metadata["branch_length"].
		"""
		if meta_key in self.metadata:
			self.metadata["branch_length"] = self.metadata[meta_key]
	
	def replaceInternalLabelWithOtherValue(self, meta_key):
		"""
		Replace the label of an internal node with a metadata value.

		Args:
			meta_key (str): Metadata key to update the value of self.label.
		"""
		if self.hasChildren():
			if meta_key in self.metadata:
				self.label = str(self.metadata[meta_key])
			else:
				self.label = ''
	
	def getNewick(self):
		"""
		Generate a Newick string representation of the tree or subtree.

		Returns:
			str: Newick-formatted string for this subtree.
		"""
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
		"""
		Generate a Newick string with additional, commented metadata.

		Returns:
			str: Newick-formatted string with branch lengths and other metadata.
		"""
		nwk = []
		# recursively get children
		if len(self.children):
			nwk.append('(')
			for i,child in enumerate(self.children):
				if i > 0:
					nwk.append(',')
				nwk.append(child.getNewickWithCommentedMetadata())
			nwk.append(')')
		# get node label
		if self.label:
			nwk.append(self.label)
		# get the branch length, if stored
		if "branch_length" in self.metadata:
			nwk.append(':')
			nwk.append(str(self.metadata["branch_length"]))
		# get other metadata, if present, in a comment
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
		"""
		Generate compact JSON representation of tree or subtree.

		Returns:
			str: JSON string representing tree or subtree.
		"""
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
		"""
		Generate a human-readable, indented JSON representation of tree or subtree.

		Args:
			indent (int, optional): Number of tab indents.

		Returns:
			str: "Pretty" JSON string.
		"""
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

		return ''.join(j)
	
	def getAscii(self, prefix="", children_prefix=""):
		"""
		Generate an ASCII art representation of the tree or subtree.

		Args:
			prefix(str, optional): Prefix for this node's line.

		Returns:
			str: ASCII representation of the tree or subtree.
		"""
		output = prefix + self.label + '\n'
		if self.hasChildren():
			for i in range(0, len(self.children) - 1, 1):
				output += self.children[i].getAscii(prefix=f"{children_prefix}|-- ", children_prefix=f"{children_prefix}|   ")
			output += self.children[-1].getAscii(prefix=f"{children_prefix}'-- ", children_prefix=f"{children_prefix}    ")
		return output
	
# ------------- STRING REPRESENTATIONS ----------- ||
	def __str__(self):
		"""
		String representation of the node for print() and str().

		Returns:
			str: Label, metadata, and number of children of this node.
		"""
		return f'{{ label: "{self.label}", metadata: {str(self.metadata)}, children: {len(self.children)} }}'
	
	def __repr__(self):
		"""
		Official string representation of the node.

		Returns:
			str: Node description, also including label, metadata, and number of children of this node.
		"""
		return "Node: " + self.__str__()
	

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

