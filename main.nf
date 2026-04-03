//need to add options to run locally or on slurm or on aws
process model_test {
	//container …
	
	input:
  val job_name
	path input_alignment
	
	output:
  path "modelTest", emit: files
  path "model.txt", emit: model

	script:
  """
  mkdir -p modelTest

  bash scripts/01-iqtreeModelTest.submit \
    -j ${job_name}_model_test \
    -i ${input_alignment}

  awk -F': ' '/Best-fit model according to BIC:/ {print \$2}' \
  modelTest/modelTest.iqtree > model.txt
  """
}

/*
* Pipeline parameters
*/
params {
   input_alignment = 'data/orig/supermatrix_dna.phy'
   job_name = 'tanos'

   //need to have options for local, slurm, or aws
}


workflow {
   main:
   tree_ch = Channel.fromPath(params.input_alignment, checkIfExists: true)

   model_test(params.job_name, tree_ch)
   model_test.out.model.view()

   publish:
   model_files = model_test.out.files
}


output {
  model_files {
    path 'results/modelTest'
    mode 'copy'
  }
}
