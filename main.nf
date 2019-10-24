// Set PyProphet subsample ratio to 1/nr_samples
// subsample_ratio = (1 / file("${params.dia_folder}/*.mzXML").size()).round(3)

process irt_filter_pqp {
    scratch true
    stageInMode "copy"
    
    publishDir "Results/SpectraST"
    
    input:
    file spec_lib from file(params.spec_lib)
    val bins from params.irt_filter_bins

    output:
    file "irt*.pqp" into irtOut
    
    """
    Rscript /home/phrt/workflows/NF-OpenSwath/bin/hrirt.R $spec_lib irt_${bins}.pqp $bins $params.irt_filter_peptides
    """
}


process irt_filter_nonlinear_pqp {
    scratch true
    stageInMode "copy"
    
    publishDir "Results/SpectraST"
    
    input:
    file spec_lib from file(params.spec_lib)
    val bins from params.irt_filter_nonlinear_bins

    output:
    file "irt*.pqp" into irtNonlinearOut
    
    """
    Rscript /home/phrt/workflows/NF-OpenSwath/bin/hrirt.R $spec_lib irt_${bins}.pqp $bins $params.irt_filter_peptides
    """
}


process openswath {
    scratch "ram-disk"
    stageInMode "copy"
    cpus params.os_threads
    
    publishDir "Results/OSW"

    input:
    file spec_lib from file(params.spec_lib)
    file mzxml from file("${params.dia_folder}/*.mzXML")
    file swath_windows from file("${params.swath_windows}")
    file irt_lib from irtOut
    file irt_nonlinear_lib from irtNonlinearOut

    output:
    file "*.osw" into openswathOut
    
    """
    OpenSwathWorkflow -in $mzxml \
    -tr $spec_lib \
    -tr_irt_nonlinear $irt_nonlinear_lib \
    -out_osw `basename $mzxml .mzXML`.osw \
    -threads $params.os_threads \
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
    -RTNormalization:outlierMethod none \
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
    stageInMode "copy"
    
    publishDir "Results/PyProphet/RunSpecific"
    
    input:
    file osw from openswathOut1

    output:
    file "*.osws" into pypSubsampleOut

    """
    pyprophet subsample --subsample_ratio=$params.subsample_ratio \
    --in=$osw \
    --out=${osw}s
    """
}


process pyprophet_merge {
    scratch true
    stageInMode "copy"
    
    publishDir "Results/PyProphet/RunSpecific"
    
    input:
    file spec_lib from file(params.spec_lib)
    file osws from pypSubsampleOut.collect()

    output:
    file "subsampled.osw" into pypMergeOut
    
    """
    pyprophet merge --template=$spec_lib \
    --out=subsampled.osw $osws
    """
}

//TODO: Remove split if not used
pypMergeOut.into{ pypMergeOut1; pypMergeOut2 }

process pypropht_learn {
    scratch true
    stageInMode "copy"
    
    cpus params.pyp_threads
    
    publishDir "Results/PyProphet/RunSpecific"
    
    input:
    file merged from pypMergeOut1

    output:
    file "model.osw" into pypLearnOut
    file "*.pdf"

    """
    pyprophet score --classifier=$params.pyp_classifier \
    --in $merged \
    --out model.osw \
    --level=ms1ms2 \
    --ss_initial_fdr=0.15 \
    --ss_iteration_fdr=0.05 \
    --threads=$params.pyp_threads \
    --ss_num_iter=$params.ss_num_iter \
    --pi0_lambda $params.pi0_lambda 
    """
}

    // pyprophet score --classifier=XGBoost \
    // --in model.osw \
    // --out model.osw \
    // --level=ms1 \
    // --ss_initial_fdr=0.1 \
    // --ss_iteration_fdr=0.05 \
    // --pi0_lambda 0.05 0.2 0.05 \
    // --threads=$params.pyp_threads

process pyprophet_apply {
    scratch true
    stageInMode "copy"
    
    publishDir "Results/PyProphet/RunSpecific"
    
    input:
    file osw from openswathOut2
    file model from pypLearnOut

    output:
    file "*.oswa" into pyApplyOut
    file "*.oswr" into pyApplyRedOut
    file "*.pdf"
   
    """
    pyprophet score --in $osw\
    --out ${osw}a \
    --group_id=feature_id \
    --classifier=$params.pyp_classifier \
    --apply_weights=$model \
    --level=ms1ms2


    pyprophet reduce --in ${osw}a \
    --out ${osw}r
    """
}


    // pyprophet score --in $osw \
    // --out $osw \
    // --group_id=feature_id \
    // --classifier=XGBoost \
    // --apply_weights=$model \
    // --level=ms1


process pyprophet_global {
    scratch true
    stageInMode "copy"
    
    publishDir "Results/PyProphet/Global"

    input:
    file spec_lib from file(params.spec_lib)
    file oswr from pyApplyRedOut.collect()

    output:
    file "model.oswr" into pypGlobalOut
    file "*.pdf"

    """
    pyprophet merge --template $spec_lib \
    --out model.oswr $oswr

    pyprophet peptide --context=global \
    --in model.oswr 

    pyprophet protein --context=global \
    --in model.oswr
    """
}


process pyprophet_backpropagate {
    scratch true
    stageInMode "copy"
    
    publishDir "Results/PyProphet/Integrated"

    input:
    file osw from pyApplyOut
    file model from pypGlobalOut

    output:
    file "*.oswab" into pypBackpropagateOut

    """
    pyprophet backpropagate --apply_scores $model \
    --in $osw --out ${osw}b
    """
}


process tric_prepare {
    scratch true
    stageInMode "copy"

    publishDir "Results/Tric"

    input:
    file osw from pypBackpropagateOut

    output:
    file "*.tsv" into tricPrepareOut

    """
    pyprophet export --in $osw \
    --format=legacy_merged \
    --max_rs_peakgroup_qvalue=0.05 \
    --ipf_max_peptidoform_pep=1.0 \
    --max_global_peptide_qvalue=1.0 \
    --max_global_protein_qvalue=1.0
    """
}


process tric_feature_alignment {
    scratch true
    stageInMode "copy"

    publishDir "Results/Tric"

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
