
import sys
import re

class MalformedPhyFile(Exception):
	pass

if __name__ == "__main__":
	
	ifn = sys.argv[1]
	odn = sys.argv[2]

	headers = []
	seqs = []
	num_tax = 0
	num_pos = 0

	# extract header and sequence data from original file
	with open(ifn, 'r') as ifd:
		fields = re.sub(r"\s+", r"\t", ifd.readline().rstrip('\n')).split('\t')
		num_tax = int(fields[0])
		num_pos = int(fields[1])

		i = 0
		line = ifd.readline()
		while line != '' and i < num_tax:
			fields = line.strip().split()
			header = fields[0]
			seq = ''.join(fields[1:])
			print(f"Processing input {header}...", file=sys.stderr)
			headers.append(header)
			seqs.append([seq])
			line = ifd.readline()
			i += 1

		round = 2
		while line != '':
			if not round % 100:
				print(f"Processing round {round}...", file=sys.stderr)
			if line != '\n':
				raise MalformedPhyFile("Expecting a blank line right now...why don't I have one?")
			line = ifd.readline() # skip blank line and get next line
			i = 0
			while line != '' and i < num_tax:
				seqs[i].append(re.sub(r"\s", '', line.rstrip('\n')))
				line = ifd.readline()
				i += 1
			round += 1
	
	# validate lengths
	if len(headers) != num_tax or len(seqs) != num_tax:
		raise MalformedPhyfile("Expected headers and seqs lists to be the same length.")

	# collapse seqs[i] lists into strings
	for i in range(num_tax):
		seqs[i] = ''.join(seqs[i])
	
	# validate sequence lengths
	for seq in seqs:
		if len(seq) != num_pos:
			raise MalformedPhyFile(f"Expected all sequences to be {num_pos} long.")
	
	# validate headers
	if len(headers) != len(frozenset(headers)):
		raise MalformedPhyFile("Headers were not all unique")
	
	# write out the new files. One file for each header. That file will have all sequences
	# except for the one belonging to the that header.
	for index,header in enumerate(headers):
		print(f"Writing output without {header}...", file=sys.stderr)
		with open(f"{odn}/{header}.fa", 'w') as ofd:
			i = 0
			while i < num_tax:
				if i != index:
					ofd.write(f">{headers[i]}\n{seqs[i]}\n")
				i += 1


