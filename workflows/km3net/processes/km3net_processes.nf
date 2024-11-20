process ConvertFilesKM3NeT{
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
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/ConvertKM3NeTFiles.py -i "${input_file}" -o "./" 
    """
}

process ClassifyEventsKM3NeT{
    input:
    tuple path(input_file), val(detectorname), val(shower_reco_name), val(energy_threshold), val(energy_low), val(energy_high), path(ar_shower_file), path(ar_track_file), val(trackscore_threshold)

    output:
    tuple path("*classified.h5"), val(detectorname), val(shower_reco_name), val(energy_threshold), val(energy_low), val(energy_high), path(ar_shower_file), path(ar_track_file), val(trackscore_threshold), emit: classified_file

    publishDir "${params.output_dir}/input", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/ClassifyEventsKM3NeT.py -i "${input_file}" -o "./" 
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
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/BlindDataKM3NeT.py -i "${input_file}" -o "./"
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
    // tuple val(ratio), path('*eventlist.hdf5'), val(iteration), emit: eventlist

    publishDir "${params.output_dir}/eventlists", mode: 'link', overwrite: true;

    """
    python3  /home/hpc/capn/capn107h/software/psr/src/scripts/CreateEventListKM3NeT.py -i "${input_file}" -o "./" -s "${source_specs_file}" --energy_th ${energy_threshold} --dist ${dist} --detector ${detectorname}
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
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/CorrectEventListKM3NeT.py -i "${input_file}" -o "./" -s "${source_specs_file}"
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
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/InjectSignalKM3NeT.py -i ${input_file} -o "./" --ratio ${ratio} --pulseshape ${pulseshape} --df ${df} --frequency ${frequency} --baseline ${baseline} --a ${a} --phi ${phi} --kappa ${kappa}
    """
}

process CombineEventListsKM3NeT{
    input:
    // tuple val(ratio), path(input_files), val(iteration)
    path input_files
    output:
    // tuple val(ratio), path("*combined_eventlist.hdf5"), val(iteration), emit: combined_file
    path "*combined_eventlist.hdf5", emit: combined_file

    publishDir "${params.output_dir}/combined_eventlists", mode: 'link', overwrite: true;

    script:
    def inputFilesString = input_files.collect { "-i ${it}" }.join(' ')
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/CombineEventListsKM3NeT.py ${inputFilesString} -o "./"
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

    output:
    tuple val(ratio), path("*epochfolding_results.hdf5"), val(iteration), emit: hdf5;
    path "*.png", emit: plot;

    publishDir "${params.output_dir}/epoch_folding", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/EpochFoldingKM3NeT.py -i "${input_file}" -o "./" --frequency ${frequency} --number_of_testf ${number_of_testf} --df ${df} --nbin ${nbin} --ratio ${ratio} --iteration ${iteration}
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
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/Chi2HistogramKM3NeT.py ${inputFilesString} -o "./" --ratio ${ratio} --nhbins ${nhbins}
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
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/SignalNoiseStatisticsKM3NeT.py ${inputFilesString} -o "./" --nbin ${nbin}
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
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/AngularResolutionKM3NeT.py -i "${input_file}" -o "./" --detector ${detector} --runtype ${runtype} --recotype ${recotype} --pltscale ${pltscale}
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
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/MultifileAngularResolutionKM3NeT.py ${inputFilesString} -o "./" --detector ${detectorname} --runtype ${runtype} --recotype ${recotype} --pltscale ${pltscale}
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
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/CombineAngularResolutionKM3NeT.py ${inputFilesString} -o "./" --detector ${detector} --runtype ${runtype} --recotype ${recotype} --pltscale ${pltscale}
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
    input:
    tuple path(input_file), val(detectorname), val(shower_reco_name), val(energy_threshold), val(energy_low), val(energy_high), path(ar_shower_file), path(ar_track_file), val(trackscore_threshold)
    // tuple val(ratio), path(input_file), val(iteration)
    path source_specs_file
    // path ar_shower_file
    // path ar_track_file
    // val trackscore_threshold
    // val energy_threshold
    // val energy_low
    // val energy_high
    // val shower_reco_name

    output:
    path "*eventlist_new.hdf5", emit: eventlist_new
    // tuple val(ratio), path('*eventlist.hdf5'), val(iteration), emit: eventlist

    publishDir "${params.output_dir}/eventlists", mode: 'link', overwrite: true;

    """
    python3  /home/hpc/capn/capn107h/software/psr/src/scripts/CreateEventListKM3NeT_new.py -i "${input_file}" -o "./" -s "${source_specs_file}" -e "${ar_shower_file}" -t "${ar_track_file}" --energy_th ${energy_threshold} --trackscore_th ${trackscore_threshold} --detector ${detectorname} --energy_low ${energy_low} --energy_high ${energy_high} --shower_reco_name ${shower_reco_name}
    """
}
