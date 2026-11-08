include { MODEL_TEST } from './processes/01-iqtreeModelTest.nf'
include { MAIN_TREE } from './processes/02-iqtreeTree.nf'
include { JACKKNIFE_ALIGNMENT } from './processes/03-jackknifeAlignment.nf'
include { IQTREE_JACKKNIFE } from './processes/04-iqtreeJackknife.nf'
include { JACKKNIFE_SUPPORT } from './processes/05-jackknifeSupport.nf'
include { JACKKNIFE_SUPPORT_COUNT } from './processes/06-jackknifeSupportCount.nf'
include { JACKKNIFE_CONVERGENCE } from './processes/07-jackknifeConvergence.nf'

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
    jackknife_treefiles = jackknife_records.map { record -> record[1] }.collect()

    JACKKNIFE_SUPPORT(tree_file, jackknife_treefiles, MAIN_TREE.out.treefile, model_ch)

    support_script = file("${projectDir}/processes/05-jackknifeSupport.py")
    jackknife_runs = jackknife_records.collect()
        .combine(MAIN_TREE.out.treefile.collect())
        .flatMap { records, main_trees ->
            (1..params.replicates).collect { replicate_count ->
                def replicate_treefiles = records
                    .findAll { record -> record[0] <= replicate_count }
                    .collect { record -> record[1] }
                tuple(replicate_count, replicate_treefiles, main_trees[0], support_script)
            }
        }

    JACKKNIFE_SUPPORT_COUNT(jackknife_runs)

    support_trees = JACKKNIFE_SUPPORT_COUNT.out.support_tree
        .map { record -> record[1] }
        .collect()

    JACKKNIFE_CONVERGENCE(
        support_trees,
        file("${projectDir}/analyze_jackknife_convergence.py")
    )
    JACKKNIFE_CONVERGENCE.out.summary.view()
}
