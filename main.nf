//need to add options to run locally or on slurm or on aws
process model_test {
	//container …
	
	input:
    //val job_name
		path input_alignment
	
	output:
    path "data/modelTest", emit: files
    stdout, emit: model

	script:
  """
    bash scripts/01-iqtreeModelTest.submit -i ${input_alignment}
    awk -F': ' '/Best-fit model according to BIC:/ {print $2}' data/modelTest/modelTest.iqtree
  """
}


process main_tree {
	//container …
	
	input:
    //val job_name
    path input_alignment
    val model
	
	output:
    path "data/mainTree"

	script:
  """
    bash scripts/02-iqtreeTree.submit -i ${input_alignment} -m ${model}
  """
}


process jackknife_alignment {
	//container …
	
	input:
    path input_alignment
	
	output:
    path "data/jackknife/aln"

	script:
  """
    bash scripts/03-jackknifeAlignment.sh -i ${input_alignment}
  """
}

process jackknife_tree {
  input:
    path jackknife_aln
    val model

  output:
    path "data/jackknife/tree"
  
  script:
  """
    bash scripts/04-iqtreeJackknife.submit -m ${model}
  """

}


/*
* Pipeline parameters
*/
params {
   input_alignment: Path = 'data/orig/supermatrix_dna.phy'
   //job name

   //need to have options for local, slurm, or aws
}


workflow {


   main:
   greeting_ch = channel.fromPath(params.input_alignment)

   // emit a greeting
   model_test(tree_ch)
   main_tree(tree_ch, iqTreeModelTest.out.model)
   jackknife_alignment(tree_ch)
   jackknife_tree(jackknife_alignment.out)
   calcScore(iqTreeJackknife.out)


   publish:
   model_files = model_test.out.files
   tree_files = main_tree.out
   jackknife__alignment_files = jackknife_alignment.out
   jackknife_tree_files = jackknife_tree.out
}


output {
  model_files {
    path 'data/modelTest'
    mode 'copy'
  }
  tree_files {
    path 'data/mainTree'
    mode 'copy'
  }
  jackknife__alignment_files {
    path 'data/jackknife/aln'
    mode 'copy'
  }
  jackknife_tree_file {
    path 'data/jackknife/tree'
    mode 'copy'
  }
}
