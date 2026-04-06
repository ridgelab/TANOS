include { MODEL_TEST } from './processes/01-iqtreeModelTest.nf'

/*
* Pipeline parameters
*/
params.job_name = 'tanos'
params.input_alignment = 'data/orig/supermatrix_dna.phy'

workflow {
   tree_file = file(params.input_alignment, checkIfExists: true)
   MODEL_TEST(params.job_name, tree_file)
   MODEL_TEST.out.model.view()
}
