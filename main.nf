include { MODEL_TEST } from './processes/01-iqtreeModelTest.nf'
include { MAIN_TREE } from './processes/02-iqtreeTree.nf'
include { JACKKNIFE_ALIGNMENT } from './processes/03-jackknifeAlignment.nf'

/*
* Pipeline parameters
*/
params.input_alignment = 'data/orig/supermatrix_dna.phy'
params.outdir = 'results'

workflow {
   tree_file = file(params.input_alignment, checkIfExists: true)
   MODEL_TEST(tree_file)
   MODEL_TEST.out.model.view()
   MAIN_TREE(tree_file, MODEL_TEST.out.model)
   JACKKNIFE_ALIGNMENT(tree_file)
}