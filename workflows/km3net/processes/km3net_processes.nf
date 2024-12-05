process ConvertFilesKM3NeT{
    errorStrategy 'ignore'

    input:
    // path input_file
    // val detectorname
    tuple path(input_file), val(detectorname), val(shower_reco_name), val(energy_threshold), val(energy_low), val(energy_high), path(ar_shower_file), path(ar_track_file), val(trackscore_threshold)

    output:
    // path "*.h5", emit: converted_file
    tuple path("*.h5"), val(detectorname), val(shower_reco_name), val(energy_threshold), val(energy_low), val(energy_high), path(ar_shower_file), path(ar_track_file), val(trackscore_threshold), emit: converted_file
    // tuple val(ratio), path("*.h5"), val(iteration), emit: converted_file

    publishDir "${params.output_dir}/input", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/ConvertKM3NeTFiles.py -i "${input_file}" -o "./" 
    """
}

process ClassifyEventsKM3NeT{
    input:
    tuple path(input_file), val(detectorname), val(shower_reco_name), val(energy_threshold), val(energy_low), val(energy_high), path(ar_shower_file), path(ar_track_file), val(trackscore_threshold)
    val parampid_folder
    
    output:
    tuple path("*classified.h5"), val(detectorname), val(shower_reco_name), val(energy_threshold), val(energy_low), val(energy_high), path(ar_shower_file), path(ar_track_file), val(trackscore_threshold), emit: classified_file
    
    publishDir "${params.output_dir}/input", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/ClassifyEventsKM3NeT.py -i "${input_file}" -p "${parampid_folder}" -o "./" --detector ${detectorname}
    """

}

process BlindDataKM3NET{
    input:
    tuple path(input_file), val(detectorname), val(shower_reco_name), val(energy_threshold), val(energy_low), val(energy_high), path(ar_shower_file), path(ar_track_file), val(trackscore_threshold)
    output:
    tuple path("*blinded.h5"), val(detectorname), val(shower_reco_name), val(energy_threshold), val(energy_low), val(energy_high), path(ar_shower_file), path(ar_track_file), val(trackscore_threshold), emit: blinded_file

    publishDir "${params.output_dir}/blinded_data", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/BlindDataKM3NeT.py -i "${input_file}" -o "./"
    """

}

process CreateEventListKM3NeT{
    input:
    tuple path(input_file), val(detectorname)
    // tuple val(ratio), path(input_file), val(iteration)
    path source_specs_file
    val dist
    val energy_threshold

    output:
    path "*eventlist.hdf5", emit: eventlist
    tuple val(ratio), path('*eventlist.hdf5'), val(iteration), emit: eventlist

    publishDir "${params.output_dir}/eventlists", mode: 'link', overwrite: true;

    """
    python3  /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/CreateEventListKM3NeT.py -i "${input_file}" -o "./" -s "${source_specs_file}" --energy_th ${energy_threshold} --dist ${dist} --detector ${detectorname}
    """
}

process CorrectEventListKM3NeT{
    input:
    path input_file
    // path input_file
    // tuple val(ratio), path(input_file), val(iteration)
    path source_specs_file
    
    output:
    path "*corrected.hdf5", emit: corrected_eventlist
    // tuple val(ratio), path("*corrected.hdf5"), val(iteration), emit: corrected_eventlist

    publishDir "${params.output_dir}/eventlists", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/CorrectEventListKM3NeT.py -i "${input_file}" -o "./" -s "${source_specs_file}"
    """
}

process InjectSignalKM3NeT{
    input:
    tuple val(ratio), path(input_file), val(iteration)
    val frequency
    val pulseshape
    val df
    val baseline
    val a
    val phi
    val kappa

    
    output:
    tuple val(ratio), path("*signal*hdf5"), val(iteration), emit: injected_signal

    publishDir "${params.output_dir}/signal", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/InjectSignalKM3NeT.py -i ${input_file} -o "./" --ratio ${ratio} --pulseshape ${pulseshape} --df ${df} --frequency ${frequency} --baseline ${baseline} --a ${a} --phi ${phi} --kappa ${kappa}
    """
}

process CombineEventListsKM3NeT{
    input:
    tuple val(ratio), path(input_files), val(iteration)
    //path input_files
    output:
    tuple val(ratio), path("*combined_eventlist.hdf5"), val(iteration), emit: combined_file
    //path "*combined_eventlist.hdf5", emit: combined_file

    publishDir "${params.output_dir}/combined_eventlists", mode: 'link', overwrite: true;

    script:
    def inputFilesString = input_files.collect { "-i ${it}" }.join(' ')
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/CombineEventListsKM3NeT.py ${inputFilesString} -o "./"
    """
}

// This worked: python3 /home/hpc/capn/capn107h/software/psr/src/scripts/CombineEventListsKM3NeT.py -i '*signal_mvm.hdf5' -o "./"

process EpochFoldingKM3NeT{
    input:
    tuple val(ratio), path(input_file), val(iteration)
    val frequency
    val number_of_testf
    val df
    val nbin
    val segment_size

    output:
    tuple val(ratio), path("*epochfolding_results.hdf5"), val(iteration), emit: hdf5;
    path "*.png", emit: plot;

    publishDir "${params.output_dir}/epoch_folding", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/EpochFoldingKM3NeT.py -i "${input_file}" -o "./" --frequency ${frequency} --number_of_testf ${number_of_testf} --df ${df} --nbin ${nbin} --ratio ${ratio} --iteration ${iteration} --segment_size ${segment_size}
    """
}

process Chi2HistogramKM3NeT{
    input:
    tuple val(ratio), path(input_files), val(iteration)
    val nhbins

    output:
    path "*maxchi2.hdf5", emit: hdf5
    path "*.png", emit: plot

    publishDir "${params.output_dir}/maxchi2", mode: 'link', overwrite: true

    script:
    def inputFilesString = input_files.collect { "-i ${it}" }.join(' ')
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/Chi2HistogramKM3NeT.py ${inputFilesString} -o "./" --ratio ${ratio} --nhbins ${nhbins}
    """
}

process SignalNoiseStatisticsKM3NeT{
    input:
    path input_files
    val nbin

    output:
    path "*StatisticOverSNR.hdf5", emit: hdf5
    path "*StatisticOverSNR_plotlin.png", emit: plot_lin
    path "*StatisticOverSNR_plotlog.png", emit: plot_log

    publishDir "${params.output_dir}/eff_statistic", mode: 'link', overwrite: true
    
    script:
    def inputFilesString = input_files.collect { "-i ${it}" }.join(' ')
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/SignalNoiseStatisticsKM3NeT.py ${inputFilesString} -o "./" --nbin ${nbin}
    """
}

process AngularResolutionKM3NeT{
    input:
    path input_file
    val detector
    val runtype
    val recotype
    val pltscale

    output: 
    path "*AngularResolutionOverEnergy*.hdf5", emit: hdf5
    path "*TestPlotAngularRes*.png", emit: plot
    path "*HistogramSeparations*.png", emit: histogram

    publishDir "${params.output_dir}/angular_resolutions", mode: 'link', overwrite: true

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/AngularResolutionKM3NeT.py -i "${input_file}" -o "./" --detector ${detector} --runtype ${runtype} --recotype ${recotype} --pltscale ${pltscale}
    """
}

process MultifileAngularResolutionKM3NeT{
    input:
    path input_files
    val detectorname
    val runtype
    val recotype
    val pltscale

    output: 
    path "*AngularResolutionOverEnergy*.hdf5", emit: hdf5
    path "*TestPlotAngularRes*.png", emit: plot 
    path "*HistogramSeparations*.png", emit: histogram

    publishDir "${params.output_dir}/angular_resolutions", mode: 'link', overwrite: true

    script:
    def inputFilesString = input_files.collect { "-i ${it}" }.join(' ')
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/MultifileAngularResolutionKM3NeT.py ${inputFilesString} -o "./" --detector ${detectorname} --runtype ${runtype} --recotype ${recotype} --pltscale ${pltscale}
    """
}

process CombineAngularResolutionKM3NeT{
    input:
    path input_files
    val detector
    val runtype
    val recotype
    val pltscale
    

    output: 
    path "*AngularResolutionOverEnergy*.hdf5", emit: hdf5
    path "*TestPlotAngularRes*.png", emit: plot

    publishDir "${params.output_dir}/angular_resolutions", mode: 'link', overwrite: true

    script:
    def inputFilesString = input_files.collect { "-i ${it}" }.join(' ')
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/CombineAngularResolutionKM3NeT.py ${inputFilesString} -o "./" --detector ${detector} --runtype ${runtype} --recotype ${recotype} --pltscale ${pltscale}
    """
}

process MultiRatioWorkflowKM3NeT{
    input:
    set files, ratio, run_id, workflow_output from workflowInputs

    output:
    path("${workflow_output}")

    script:
    """
    mkdir -p ${workflow_output}

    nextflow run \
        -process.queueSize 4 \
        -with-singularity ${params.singularity_image} \
        -entrypoint workflow_executor \
        --files ${files} \
        --ratio ${ratio} \
        --run_id ${run_id} \
        --workflow_output ${workflow_output}
    """
}

process CreateEventListKM3NeT_new{
    // errorStrategy 'ignore'
    input:
    tuple path(input_file), val(detectorname), val(shower_reco_name), val(energy_threshold), val(energy_low), val(energy_high), path(ar_shower_file), path(ar_track_file), val(trackscore_threshold)
    // tuple val(ratio), path(input_file), val(iteration)
    path source_specs_file
    val delta_search_min


    output:
    path "*eventlist_new.hdf5", emit: eventlist_new
    // tuple val(ratio), path('*eventlist.hdf5'), val(iteration), emit: eventlist

    publishDir "${params.output_dir}/eventlists", mode: 'link', overwrite: true;

    """
    python3  /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/CreateEventListKM3NeT_new.py -i "${input_file}" -o "./" -s "${source_specs_file}" -e "${ar_shower_file}" -t "${ar_track_file}" --energy_th ${energy_threshold} --trackscore_th ${trackscore_threshold} --detector ${detectorname} --energy_low ${energy_low} --energy_high ${energy_high} --shower_reco_name ${shower_reco_name} --delta_search_min ${delta_search_min}
    """
}


process ExtractHitFeatures{
    errorStrategy 'ignore'

    input:
    path(input_file)
    val(detectorname)

    output:
    // path "*.h5", emit: converted_file
    path "*.h5", emit: converted_file
    // tuple val(ratio), path("*.h5"), val(iteration), emit: converted_file

    publishDir "${params.output_dir}/extracted_hitfeatures", mode: 'link', overwrite: true;

    script:
    def scriptOptions = detectorname == 'ORCA' ? '-jsh -jg' : '-as -jg'
    """
    python3 /home/hpc/capn/capn107h/software/newhitfeatures/scripts/extractor_ARCA.py -f "${input_file}" -dir "./" ${scriptOptions}
    """
}

process ConcatFiles{

    input:
    path input_files

    output:
    // path "*.h5", emit: converted_file
    path "*concatenated.h5", emit: concatenated_file

    publishDir "${params.output_dir}/concatenated_hitfeatures", mode: 'link', overwrite: true;
    
    """
    python3 /home/hpc/capn/capn107h/software/newhitfeatures/scripts/concat.py -f ${input_files.join(' ')} -o "./extracthitfeatures_concatenated.h5"
    """
}

process TrainPID{

    input:
    path input_file
    path columntable_file
    val detectorname
    
    output:
    // path "*.h5", emit: converted_file
    path "*.rdf", emit: rd_file
    path "Accuracy*.png", emit: accuracy_plot;
    path "Separability*.png", emit: separability_plot;
    path "pid_track_score*.png", emit: trackscore_plot;

    publishDir "${params.output_dir}/trained_rdfiles", mode: 'link', overwrite: true;

    """
    python3 /home/hpc/capn/capn107h/software/parampid/APC_PID.py -i "${input_file}" -w "./" -d ${detectorname} -c ${columntable_file} -t track -s track_score
    """

}

process ApplyClassifierPID{

    input:
    path input_file
    path rd_file
    path columntable_file

    output:
    path "*.h5", emit: output_file

    publishDir "${params.output_dir}/classified_files", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/parampid/scripts/applyRDF.py -i ${input_file} -r ${rd_file} -c ${columntable_file} -o ${input_file.toString().replace('.h5', '_scored.h5')} -s track_score
    """
}