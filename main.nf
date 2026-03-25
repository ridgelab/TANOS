//need to add options to run locally or on slurm or on aws
process iqTreeModelTest {
	//container …
	
	input:
		path input_alignment
		val output_prefix
	
	output:
		MODEL="GTR+F+I+G4"

	script:
  """
    bash scripts/01-iqtreeModelTest.submit ${input_alignment} ${output_prefix}
  """
}


process iqTreeTree {
	//container …
	
	input:
    path input_alignment
    val model
	
	output:


	script:
  """
    bash scripts/02-iqtreeTree.submit ${input_alignment} ${model}
  """
}


process jackknifeAlignment {
	//container …
	
	input:
    path input_alignment
	
	output:

	script:
  """
    bash scripts/03-jackknifeAlignment.sh ${input_alignment}
  """
}


/*
* Pipeline parameters
*/
params {
   input: Path = 'data/orig/supermatrix_dna.phy'
   output_prefix: String = ‘data/modeltest’
   batch: String = 'batch'
}


workflow {


   main:
   // create a channel for inputs from a CSV file
   greeting_ch = channel.fromPath(params.input)


   // emit a greeting
   iqTreeModelTest(tree_ch)
   iqTreeTree(tree_ch, iqTreeModelTest.out)
   jackknifeAlignment(tree_ch)
   iqTreeJackknife(jackknifeAlignment.out)
   calcScore(iqTreeJackknife.out)


   publish:
   tree_file = iqTreeTree.out
   jackknife_tree_file = jackknifeAlignment.out
}


output {
   tree_file {
       path 'data/mainTree/tree'
       mode 'copy'
   }
   jackknife_tree_file {
       path 'data/jackknife/aln'
       mode 'copy'
   }
}
