include { MODEL_TEST } from './processes/01-iqtreeModelTest.nf'
include { MAIN_TREE } from './processes/02-iqtreeTree.nf'
include { JACKKNIFE_ALIGNMENT } from './processes/03-jackknifeAlignment.nf'
include { IQTREE_JACKKNIFE } from './processes/04-iqtreeJackknife.nf'
include { JACKKNIFE_SUPPORT } from './processes/05-jackknifeSupport.nf'
include { JACKKNIFE_CONVERGENCE } from './processes/06-jackknifeConvergence.nf'

/*
* Pipeline parameters
*/

workflow {

    if (params.replicates < 1) {
        error "params.replicates must be at least 1"
    }

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
    jackknife_records = IQTREE_JACKKNIFE.out.trees
    main_tree_value = MAIN_TREE.out.treefile.first()

    jackknife_runs = jackknife_records.collect()
        .flatMap { records ->
            (1..params.replicates).collect { replicate_count ->
                def replicate_treefiles = records
                    .findAll { record -> record[0] <= replicate_count }
                    .sort { a, b -> a[0] <=> b[0] }
                    .collect { record -> record[1] }
                tuple(replicate_count, replicate_treefiles)
            }
        }

    JACKKNIFE_SUPPORT(jackknife_runs, main_tree_value)

    support_trees = JACKKNIFE_SUPPORT.out.support_tree
        .map { record -> record[1] }
        .collect()

    JACKKNIFE_CONVERGENCE(
        support_trees,
        file("${projectDir}/analyze_jackknife_convergence.py")
    )
    JACKKNIFE_CONVERGENCE.out.summary.view()
}
