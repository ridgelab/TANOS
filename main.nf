include { MODEL_TEST } from './processes/01-iqtreeModelTest.nf'
include { MAIN_TREE } from './processes/02-iqtreeTree.nf'
include { JACKKNIFE_ALIGNMENT } from './processes/03-jackknifeAlignment.nf'
include { IQTREE_JACKKNIFE } from './processes/04-iqtreeJackknife.nf'
include { JACKKNIFE_SUPPORT } from './processes/05-jackknifeSupport.nf'

/*
* Pipeline parameters
*/

workflow {

    tree_file = file(params.input_alignment, checkIfExists: true)

    MODEL_TEST(tree_file)

    MODEL_TEST.out.model.view()

    model_ch = MODEL_TEST.out.model.map { model_file ->
        model_file.text.trim()
    }

    MAIN_TREE(tree_file, model_ch)

    JACKKNIFE_ALIGNMENT(tree_file)
    jackknife_ch = JACKKNIFE_ALIGNMENT.out.files
        .flatten()
        .flatMap { fa ->
            (1..params.replicates).collect { rep ->
                tuple(fa, rep)
            }
        }

    IQTREE_JACKKNIFE(jackknife_ch, model_ch)
    jackknife_treefiles = IQTREE_JACKKNIFE.out.collect()

    JACKKNIFE_SUPPORT(tree_file, jackknife_treefiles, MAIN_TREE.out.treefile, model_ch)
}
