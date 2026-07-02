include { MODEL_TEST } from './processes/01-iqtreeModelTest.nf'
include { MAIN_TREE } from './processes/02-iqtreeTree.nf'
include { JACKKNIFE_ALIGNMENT } from './processes/03-jackknifeAlignment.nf'
include { IQTREE_JACKKNIFE } from './processes/04-iqtreeJackknife.nf'

/*
* Pipeline parameters
*/
params.input_alignment = 'data/orig/supermatrix_dna.phy'
params.outdir = 'results'

Channel
    .fromPath("${params.input_alignment}")
    .map { file ->
        def taxon = file.simpleName.replaceAll(/\.phy$/, "")
        tuple(taxon, file)
    }

workflow {

    tree_file = file(params.input_alignment, checkIfExists: true)

    MODEL_TEST(tree_file)

    MODEL_TEST.out.model.view()

    model_ch = MODEL_TEST.out.model.map {   model_file ->
        model_file.text.trim()
    }

    MAIN_TREE(tree_file, model_ch)

    JACKKNIFE_ALIGNMENT(tree_file)
        .map { fa ->
            def rep = fa.name.replaceAll(/\\.fa$/, "").toInteger()
            tuple(fa, rep)
        }
        | IQTREE_JACKKNIFE
}