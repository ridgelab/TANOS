# Will I have to change calcScore to make sure it doesn’t get messed up if there’s a label indicating a parameter but no actual parameter?
# How would I include a parameter for data directory

process iqTreeModelTest {
	container …
	
	input:
		path input_alignment
		val output_prefix
	
	output:
		val model

	script:
}


process iqTreeTree {
	container …
	
	input:
		path input_alignment
		val model
	
	output:


	script:
}


process jackknifeAlignment {
	container …
	
	input:
		path input_alignment
	
	output:


	script:
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
