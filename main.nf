// Set PyProphet subsample ratio to 1/nr_samples
// subsample_ratio = (1 / file("${params.dia_folder}/*.mzXML").size()).round(3)

process irt_linear_filter_pqp {
    scratch true
    stageInMode 'copy'
    
    publishDir 'Results/openswath', mode:'link'
    
    input:
    file target_library from file(params.target_library)

    output:
    file "linear_irt.pqp" into irtLinearFilterOut
    
    """
    easypqp reduce --in=$target_library \
    --out=linear_irt.pqp \
    --bins=$params.irt_linear_filter_bins \
    --peptides=$params.irt_linear_filter_peptides
    """
}


process irt_nonlinear_filter_pqp {
    scratch true
    stageInMode 'copy'
    
    publishDir 'Results/openswath', mode:'link'
    
    input:
    file target_library from file(params.target_library)

    output:
    file "nonlinear_irt.pqp" into irtNonlinearFilterOut
    
    """
    easypqp reduce --in=$target_library \
    --out=nonlinear_irt.pqp \
    --bins=$params.irt_nonlinear_filter_bins \
    --peptides=$params.irt_nonlinear_filter_peptides
    """
}


process openswath {
    scratch 'ram-disk'
    stageInMode 'copy'
    cpus params.openswath_threads
    
    publishDir 'Results/openswath', mode:'link'

    input:
    file library from file(params.library)
    file mzxml from file("${params.dia_folder}/*.mzXML")
    file swath_windows from file("${params.openswath_swath_windows}")
    file irt_linear from irtLinearFilterOut
    file irt_nonlinear from irtNonlinearFilterOut

    output:
    file "*.osw" into openswathOut
    
    """
    OpenSwathWorkflow -in $mzxml \
    -tr $library \
    -tr_irt $irt_linear \
    -tr_irt_nonlinear $irt_nonlinear \
    -out_osw `basename $mzxml .mzXML`.osw \
    -threads $params.openswath_threads \
    -swath_windows_file $swath_windows \
    -min_upper_edge_dist 0 \
    -mz_extraction_window 30 \
    -mz_extraction_window_unit ppm \
    -mz_extraction_window_ms1 20 \
    -mz_extraction_window_ms1_unit ppm \
    -mz_correction_function regression_delta_ppm \
    -use_ms1_traces \
    -irt_mz_extraction_window 50 \
    -irt_mz_extraction_window_unit ppm \
    -rt_extraction_window 600 \
    -RTNormalization:estimateBestPeptides \
    -RTNormalization:alignmentMethod lowess \
    -RTNormalization:lowess:span 0.05 \
    -RTNormalization:outlierMethod none \
    -RTNormalization:estimateBestPeptides \
    -RTNormalization:MinBinsFilled 5 \
    -Scoring:stop_report_after_feature 5 \
    -Scoring:TransitionGroupPicker:compute_peak_quality false \
    -Scoring:Scores:use_ms1_mi \
    -Scoring:Scores:use_mi_score \
    -readOptions cache \
    -batchSize 1000 \
    -ms1_isotopes 3
    """
}


openswathOut.into{ openswathOut1; openswathOut2 }


process pyprophet_subsample {
    scratch true
    stageInMode 'copy'
    
    publishDir 'Results/PyProphet/RunSpecific', mode:'link'
    
    input:
    file osw from openswathOut1

    output:
    file "*.osws" into pyprophetSubsampleOut

    """
    pyprophet subsample --subsample_ratio=$params.pyprophet_subsample_subsample_ratio \
    --in=$osw \
    --out=${osw}s
    """
}


process pyprophet_merge {
    scratch true
    stageInMode 'copy'
    
    publishDir 'Results/PyProphet/RunSpecific', mode:'link'
    
    input:
    file library from file(params.library)
    file osws from pyprophetSubsampleOut.collect()

    output:
    file "subsampled.osw" into pyprophetMergeOut
    
    """
    pyprophet merge --template=$library \
    --out=subsampled.osw \
    $osws
    """
}


process pyprophet_learn {
    scratch true
    stageInMode 'copy'
    
    cpus params.pyprophet_learn_threads
    
    publishDir 'Results/PyProphet/RunSpecific', mode:'link'
    
    input:
    file merged_osws from pyprophetMergeOut

    output:
    file "model.osw" into pyprophetLearnOut
    file "model_ms1ms2_report.pdf"

    """
    pyprophet score --classifier=$params.pyprophet_classifier \
    --in $merged_osws \
    --out model.osw \
    --level=ms1ms2 \
    --xeval_num_iter=3 \
    --ss_initial_fdr=0.05 \
    --ss_iteration_fdr=0.01 \
    --pi0_method=$params.pyprophet_learn_pi0_method \
    --threads=$params.pyprophet_learn_threads
    """
}


process pyprophet_apply {
    scratch true
    stageInMode 'copy'
    
    publishDir 'Results/PyProphet/RunSpecific', mode:'link'
    
    input:
    file osw from openswathOut2
    file model from pyprophetLearnOut

    output:
    file "*.oswa" into pyprophetApplyOut
    file "*.oswr" into pyprophetApplyReducedOut
    file "*.pdf"
   
    """
    pyprophet score --in $osw\
    --out ${osw}a \
    --group_id=feature_id \
    --classifier=$params.pyprophet_classifier \
    --apply_weights=$model \
    --level=ms1ms2


    pyprophet reduce --in ${osw}a \
    --out ${osw}r
    """
}


process pyprophet_global {
    scratch true
    stageInMode 'copy'
    
    publishDir 'Results/PyProphet/Global', mode:'link'

    input:
    file library from file(params.library)
    file oswr from pyprophetApplyReducedOut.collect()

    output:
    file "model.oswr" into pyprophetGlobalOut
    file "*.pdf"

    """
    pyprophet merge --template $library \
    --out model.oswr $oswr

    pyprophet peptide --context=global \
    --in model.oswr 

    pyprophet protein --context=global \
    --in model.oswr
    """
}


process pyprophet_backpropagate {
    scratch true
    stageInMode 'copy'
    
    publishDir 'Results/PyProphet/Integrated', mode:'link'

    input:
    file oswa from pyprophetApplyOut
    file model from pyprophetGlobalOut

    output:
    file "*.osw" into pyprophetBackpropagateOut

    """
    pyprophet backpropagate --apply_scores $model \
    --in $oswa --out ${oswa.baseName}.osw
    """
}


process tric_prepare {
    scratch true
    stageInMode 'copy'

    publishDir 'Results/Tric', mode:'link'

    input:
    file osw from pyprophetBackpropagateOut

    output:
    file "*.tsv" into tricPrepareOut

    """
    pyprophet export --in $osw \
    --format=legacy_merged \
    --max_rs_peakgroup_qvalue=0.2 \
    --ipf_max_peptidoform_pep=1.0 \
    --max_global_peptide_qvalue=0.01 \
    --max_global_protein_qvalue=0.01
    """
}


process tric_feature_alignment {
    scratch true
    stageInMode 'copy'

    publishDir 'Results/Tric', mode:'link'

    input:
    file tsvs from tricPrepareOut.collect()

    output:
    file "feature_alignment.tsv"
    file "feature_alignment_matrix.tsv"

    """
    feature_alignment.py --in $tsvs \
    --out feature_alignment.tsv \
    --out_matrix feature_alignment_matrix.tsv \
    --method LocalMST \
    --realign_method lowess_cython \
    --max_rt_diff 60 \
    --mst:useRTCorrection True \
    --mst:Stdev_multiplier 3.0 \
    --fdr_cutoff 0.01 \
    --max_fdr_quality 0.05 \
    --alignment_score $params.alignment_score
    """
}
