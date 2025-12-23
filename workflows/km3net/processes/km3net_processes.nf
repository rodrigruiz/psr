process ConvertFilesKM3NeT{
    errorStrategy 'ignore'

    input:
    // path input_file
    // val detectorname
    tuple path(input_file), val(detectorname)

    output:
    // path "*.h5", emit: converted_file
    tuple path("*.h5"), val(detectorname), emit: converted_file
    // tuple val(ratio), path("*.h5"), val(iteration), emit: converted_file

    publishDir "${params.output_dir}/input", mode: 'link', overwrite: true //, enabled: false
    //publishDir = false

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/ConvertKM3NeTFiles.py -i "${input_file}" -o "./" --detector ${detectorname}
    """
}

process AddTrackScoreKM3NeT{
    errorStrategy 'ignore'
    input:
    tuple path(input_file), val(detectorname)
    val parampid_folder
    
    output:
    tuple path("*classified.h5"), val(detectorname), emit: classified_file
    
    publishDir "${params.output_dir}/input", mode: 'link', overwrite: true //, enabled: false

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/AddTrackScoreKM3NeT.py -i "${input_file}" -p "${parampid_folder}" -o "./" --detector ${detectorname} --random_track_score
    # python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/AddTrackScoreKM3NeT.py -i "${input_file}" -p "${parampid_folder}" -o "./" --detector ${detectorname}

    """

}

process BlindDataKM3NET{
    errorStrategy 'ignore'
    input:
    tuple path(input_file), val(detectorname)
    output:
    tuple path("*blinded.h5"), val(detectorname), emit: blinded_file

    publishDir "${params.output_dir}/blinded_data", mode: 'link', overwrite: true, enabled: false

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

    publishDir "${params.output_dir}/eventlists", mode: 'link', overwrite: true, enabled: false

    """
    python3  /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/CreateEventListKM3NeT.py -i "${input_file}" -o "./" -s "${source_specs_file}" --energy_th ${energy_threshold} --dist ${dist} --detector ${detectorname}
    """
}

process CorrectEventListKM3NeT{
    errorStrategy 'ignore'
    input:
    tuple path(input_file), val(detectorname)
    // path input_file
    // tuple val(ratio), path(input_file), val(iteration)
    path source_specs_file
    
    output:
    tuple path("*corrected.hdf5"), val(detectorname),  emit: corrected_eventlist
    // tuple val(ratio), path("*corrected.hdf5"), val(iteration), emit: corrected_eventlist

    publishDir "${params.output_dir}/eventlists_preprocessed", mode: 'link', overwrite: true //, enabled: false

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/CorrectEventListKM3NeT.py -i "${input_file}" -o "./" -s "${source_specs_file}"
    """
}

process InjectSignalKM3NeT_old{
    errorStrategy 'ignore'
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

process InjectSignalKM3NeT_backup{
    //errorStrategy 'ignore'
    input:
    tuple val(ratio), path(input_file), val(detectorname), path(total_events_file), val(iteration), val(energy_min), val(gamma)
    val frequency
    val pulseshape
    val df
    val baseline
    val a
    val phi
    val kappa
    val method


    
    output:
    tuple val(ratio), path("*signal*hdf5"), val(detectorname), val(iteration), val(energy_min), val(gamma), emit: injected_signal
    path("*.png"), emit: plot , optional : true
    path("*total_events2.txt"), emit: total_events , optional : true //gets created for each injection, so not only once...

    publishDir "${params.output_dir}/signal", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/InjectSignalKM3NeT_new.py -i ${input_file} -o "./" --total_events_file ${total_events_file} --ratio ${ratio} --pulseshape ${pulseshape} --df ${df} --frequency ${frequency} --baseline ${baseline} --a ${a} --phi ${phi} --kappa ${kappa} --method ${method} --E_min ${energy_min} --gamma ${gamma} #--plot True
    """
}

process InjectSignalKM3NeT{
    //errorStrategy 'ignore'
    input:
    tuple val(rate), path(input_file), val(detectorname), val(iteration), val(energy_min), val(gamma)
    val frequency
    val pulseshape
    val df
    val baseline
    val a
    val phi
    val kappa
    val method
    path energy_response_file


    
    output:
    tuple val(rate), path("*signal*hdf5"), val(detectorname), val(iteration), val(energy_min), val(gamma), emit: injected_signal
    path("*.png"), emit: plot , optional : true
    path("*total_events2.txt"), emit: total_events , optional : true //gets created for each injection, so not only once...

    publishDir "${params.output_dir}/signal", mode: 'link', overwrite: true //, enabled: false

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/InjectSignalKM3NeT.py -i ${input_file} -o "./" -e ${energy_response_file} --rate ${rate} --pulseshape ${pulseshape} --df ${df} --frequency ${frequency} --baseline ${baseline} --a ${a} --phi ${phi} --kappa ${kappa} --method ${method} --E_min ${energy_min} --gamma ${gamma} #--plot True
    """
}

process CombineEventListsRuns{
    input:
    //tuple val(ratio), path(input_files), val(iteration)
    //tuple val(input_files), val(source_specs_file), val(delta_search_min), val(filestype), val(detector)
    tuple path(input_files), val(detector)
    path source_specs_file
    val delta_search_min
    val filestype
    //val detector

    output:
    //tuple val(ratio), path("*combined_eventlist*.hdf5"), val(iteration), emit: combined_file
    tuple path("*_combined_*.hdf5"), val(detector),  emit: combined_file
    path("*lightcurve.png"), emit: lightcurve_plot, optional: true
    path("*total_events.txt"), emit: total_events, optional: true 
    path("*zenith_over_time.png"), emit: zenith_over_time_plot, optional: true
    path("*energy_check.txt"), emit: energy_check, optional: true
    path("*gtis.pkl"), emit: gti_file, optional: true

    publishDir "${params.output_dir}/combined_eventlists", mode: 'link', overwrite: true;

    script:
    def inputFilesString = input_files.collect { "-i ${it}" }.join(' ')
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/CombineEventListsKM3NeT.py ${inputFilesString} -o "./" -s ${source_specs_file} --delta_search_min ${delta_search_min} --filestype ${filestype} --detector ${detector}
    """
}

process CombineEventListsDetectors{
    input:
    //tuple val(ratio), path(input_files), val(iteration)
    //tuple val(input_files), val(source_specs_file), val(delta_search_min), val(filestype), val(detector)
    tuple val(ratio), path(input_files), val(detectorname), val(iteration), val(energy_min), val(gamma)
    path source_specs_file
    val delta_search_min
    val filestype
    //val detector

    output:
    //tuple val(ratio), path("*combined_eventlist*.hdf5"), val(iteration), emit: combined_file
    tuple val(ratio), path("*combinedet*.hdf5"), val(iteration), val(energy_min), val(gamma), emit: combined_file
    path("*lightcurve.png"), emit: lightcurve_plot, optional: true
    path("*total_events.txt"), emit: total_events, optional: true 
    path("*zenith_over_time.png"), emit: zenith_over_time_plot, optional: true
    path("*energy_check.txt"), emit: energy_check, optional: true

    publishDir "${params.output_dir}/combined_eventlists", mode: 'link', overwrite: true;

    script:
    def inputFilesString = input_files.collect { "-i ${it}" }.join(' ')
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/CombineEventListsKM3NeT.py ${inputFilesString} -o "./" -s ${source_specs_file} --delta_search_min ${delta_search_min} --filestype ${filestype} --detector 'arcaorca' --combinedet 'True'
    """
}

// This worked: python3 /home/hpc/capn/capn107h/software/psr/src/scripts/CombineEventListsKM3NeT.py -i '*signal_mvm.hdf5' -o "./"

process EpochFoldingKM3NeT{
    input:
    tuple val(ratio), path(input_file), path(gti_file), val(iteration), val(energy_min), val(gamma), val(length)
    val frequency
    val number_of_testf
    val df
    val nbin
    val segment_size
    // path gti_file

    output:
    tuple val(ratio), path("*ef.hdf5"), val(iteration), val(energy_min), val(gamma), val(length), emit: hdf5;
    path "*.png", emit: plot;

    publishDir "${params.output_dir}/epoch_folding", mode: 'link', overwrite: true, enabled: false

    script:
    """
    echo "Running Epoch Folding for ratio=${ratio}, iteration=${iteration}"
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/EpochFoldingKM3NeT.py -i "${input_file}" -o "./" --gti_files ${gti_file} --frequency ${frequency} --number_of_testf ${number_of_testf} --df ${df} --nbin ${nbin} --ratio ${ratio} --iteration ${iteration} --segment_size ${segment_size}
    """
}

process Chi2HistogramKM3NeT{
    input:
    tuple val(ratio), path(input_files), val(iteration), val(energy_min), val(gamma), val(length)
    val nhbins

    output:
    tuple path("*maxchi2.hdf5"), val(energy_min), val(gamma), val(length), emit: hdf5
    path "*.png", emit: plot

    publishDir "${params.output_dir}/maxchi2", mode: 'link', overwrite: true, enabled: false

    script:
    def inputFilesString = input_files.collect { "-i ${it}" }.join(' ')
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/Chi2HistogramKM3NeT.py ${inputFilesString} -o "./" --ratio ${ratio} --nhbins ${nhbins}
    """
}

process SignalNoiseStatisticsKM3NeT{
    input:
    tuple path(input_files), val(energy_min), val(gamma), val(length)
    val nbin
    val frequency
    val angle
    path total_events_file

    output:
    tuple path("*StatisticOverRate.hdf5"), val(energy_min), val(gamma), val(length), emit: hdf5
    path "*StatisticOverRate_plotlin.png", emit: plot_lin
    path "*StatisticOverRate_plotlog.png", emit: plot_log

    publishDir "${params.output_dir}/eff_statistic", mode: 'link', overwrite: true
    
    script:
    def inputFilesString = input_files.collect { "-i ${it}" }.join(' ')
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/SignalNoiseStatisticsKM3NeT.py ${inputFilesString} -o "./" --total_events_file ${total_events_file} --nbin ${nbin} --frequency ${frequency} --angle ${angle} --E_min ${energy_min} --gamma ${gamma} --length ${length}
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
    val energy_low
    val energy_high
    val nbins

    output: 
    path "*AngularResolutionOverEnergy*.hdf5", emit: hdf5
    path "*TestPlotAngularRes*.png", emit: plot 
    path "*HistogramSeparations*.png", emit: histogram

    publishDir "${params.output_dir}/angular_resolutions", mode: 'link', overwrite: true

    script:
    def inputFilesString = input_files.collect { "-i ${it}" }.join(' ')
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/MultifileAngularResolutionKM3NeT.py ${inputFilesString} -o "./" --detector ${detectorname} --runtype ${runtype} --recotype ${recotype} --pltscale ${pltscale} --energy_low ${energy_low} --energy_high ${energy_high} --nbins ${nbins}
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

process CreateEventListKM3NeT_new_old{
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

process CreateEventListKM3NeT_new {
    errorStrategy 'ignore'

    input:
    tuple path(input_file), val(detectorname)
    path source_specs_file
    val delta_search_min

    path ar_shower_file_arca
    path ar_track_file_arca
    path ar_shower_file_orca
    path ar_track_file_orca    

    val energy_low_arca
    val energy_high_arca
    val energy_low_orca
    val energy_high_orca

    val energy_threshold_arca
    val trackscore_threshold_arca
    val muonscore_threshold_arca
    val energy_threshold_orca
    val trackscore_threshold_orca
    val muonscore_threshold_orca

    val cone_all


    output:
    tuple path("*eventlist_new.hdf5"), val(detectorname), emit: eventlist_new

    publishDir "${params.output_dir}/eventlists", mode: 'link', overwrite: true, enabled: false

    script:

    """
    # Select the correct files based on the detectorname
    

    if [[ "${detectorname}" == "arca" ]]; then
        python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/CreateEventListKM3NeT_new.py \
            -i "${input_file}" \
            -o "./" \
            -s "${source_specs_file}" \
            -e "${ar_shower_file_arca}" \
            -t "${ar_track_file_arca}" \
            --energy_th ${energy_threshold_arca} \
            --trackscore_th ${trackscore_threshold_arca} \
            --muonscore_th ${muonscore_threshold_arca} \
            --detector ${detectorname} \
            --energy_low ${energy_low_arca} \
            --energy_high ${energy_high_arca} \
            --delta_search_min ${delta_search_min} \
            --cone_all ${cone_all}
    elif [[ "${detectorname}" == "orca" ]]; then
        python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/CreateEventListKM3NeT_new.py \
            -i "${input_file}" \
            -o "./" \
            -s "${source_specs_file}" \
            -e "${ar_shower_file_orca}" \
            -t "${ar_track_file_orca}" \
            --energy_th ${energy_threshold_orca} \
            --trackscore_th ${trackscore_threshold_orca} \
            --muonscore_th ${muonscore_threshold_orca} \
            --detector ${detectorname} \
            --energy_low ${energy_low_orca} \
            --energy_high ${energy_high_orca} \
            --delta_search_min ${delta_search_min} \
            --cone_all ${cone_all}
    else
        echo "Error: Unknown detector name: ${detectorname}"
        exit 1
    fi

    # Run the Python script with the selected files

    """
}



process ExtractHitFeatures {
    errorStrategy 'ignore'
    //label 'short'

    input:
    path input_file
    val detectorname
    val filetype

    output:
    path "*.h5", emit: converted_file

    publishDir "${params.output_dir}/extracted_hitfeatures", mode: 'link', overwrite: true

    script:
    // Determine which script to use based on the detector name
    def scriptPath = detectorname == 'ORCA' ? 'extractor.py' : 'extractor_ARCA.py'
    
    // Set the appropriate script options
    def scriptOptions = detectorname == 'ORCA' ? '-jsh -jg' : '-as -jg'
    //def scriptOptions = detectorname == 'ORCA' ? '-as -jg' : '-as -jg'
    // Allow option for data files
    def filetypeOption = filetype == 'data' ? '-dat' : ''
    def detxOption = detectorname == 'ORCA' ? '-d "/home/hpc/capn/capn107h/software/psr/workflows/km3net/km3net_parampid/calibration_00000148_00016545.detx"' : ''

    """
    python3 /newhitfeatures/scripts/${scriptPath} -f "${input_file}" -dir "./" ${detxOption} ${scriptOptions} ${filetypeOption}
    #python3 /home/hpc/capn/capn107h/software/newhitfeatures/scripts/extractor_ARCA.py -f "${input_file}" -dir "./" ${scriptOptions} ${filetypeOption}
    """
}

process ConcatFiles{
    label 'short'

    input:
    path input_files

    output:
    // path "*.h5", emit: converted_file
    path "*concatenated.h5", emit: concatenated_file

    publishDir "${params.output_dir}/concatenated_hitfeatures", mode: 'link', overwrite: true;
    
    """
    python3 /newhitfeatures/scripts/concat.py -f ${input_files.join(' ')} -o "./extracthitfeatures_concatenated.h5"
    """
}

process TrainPID{

    label 'long'

    input:
    path input_file
    path columntable_file
    val detectorname
    val classtype
    val colname
    
    output:
    path "*.h5", emit: pid_output
    path "*.rdf", emit: rd_file
    path "Accuracy*.png", emit: accuracy_plot, optional: true
    path "Separability*.png", emit: separability_plot, optional: true
    path "pid_track_score*.png", emit: trackscore_plot, optional: true
    path "muon_score*.png", emit: muonscore_plot, optional: true

    publishDir "${params.output_dir}/trained_rdfiles", mode: 'link', overwrite: true;

    """
    python3 /parampid/APC_PID.py -i "${input_file}" -w "./" -d ${detectorname} -c ${columntable_file} -t ${classtype} -sn ${colname}
    """

}

process ApplyClassifierPID{
    label 'long'

    input:
    path input_file
    path rd_file
    path columntable_file
    val colname
    //val muon_score_threshold = null

    output:
    path "*.h5", emit: output_file

    publishDir "${params.output_dir}/classified_files", mode: 'link', overwrite: true;

    script:
    """
    python3 /parampid/scripts/applyRDF.py -i ${input_file} -r ${rd_file} -c ${columntable_file} -o ${input_file.toString().replace('.h5', '_scored.h5')} -s ${colname}
    """
}

process TestProcess{

    input:
    path input_file

    output:
    path "*.txt"

    script:
    def filename = input_file.getBaseName()
    def txt_file = "${filename}.txt"
    """
    echo "$filename" > $txt_file
    """

}

process TestProcess2{

    input:
    val number

    output:
    path "*.txt"

    script:
    def txt_file = "testfile${number}.txt"
    """
    echo "$number" > $txt_file
    """
}

process FilterMuonScore{
    label 'short'

    input:
    path input_file
    val muon_score_threshold

    output:
    path "*.h5", emit: filtered_file

    publishDir "${params.output_dir}/filtered_files", mode: 'link', overwrite: true;

    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/FilterMuonScore_ParamPID.py --input "${input_file}" --threshold ${muon_score_threshold} -o "./"
    """
}


process PlotSkyMapsKM3NeT{
    errorStrategy 'ignore'
    input:
    tuple path(input_file), val(detectorname)
    //path input_file
    //val detectorname
    path source_specs_file
    val output_folder
    val radius
    val num_bins


    
    output:
    path("*_SkyMap.png"), emit: events_plot , optional : true
    path("*_HistoSkyMap.png"), emit: histo_plot , optional : true 
    path("*_RA_Dec_2DHist.png"), emit: radec_histo_plot , optional: true
    path("*_Theta_Phi_2DHist.png"), emit: thetaphi_histo_plot , optional: true
    path("*_EnergyHistogram.png"), emit: energy_hist , optional: true
    path("*_RA_Dec_Plot.png"), emit: radec_plot , optional: true

    publishDir "${params.output_dir}/${output_folder}", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/PlotSkyMapsKM3NeT.py -i ${input_file} -o "./" -s "${source_specs_file}" --radius ${radius} --num_bins ${num_bins} --detectorname ${detectorname}
    """
}

process PlotSkyMapsKM3NeTCombined{
    errorStrategy 'ignore'
    input:
    tuple val(ratio), path(input_file), val(iteration)
    val detectorname
    path source_specs_file
    val output_folder
    val radius
    val num_bins


    
    output:
    path("*_SkyMap.png"), emit: events_plot , optional : true
    path("*_HistoSkyMap.png"), emit: histo_plot , optional : true 
    path("*_RA_Dec_2DHist.png"), emit: radec_histo_plot , optional: true
    path("*_Theta_Phi_2DHist.png"), emit: thetaphi_histo_plot , optional: true
    path("*_EnergyHistogram.png"), emit: energy_hist , optional: true
    path("*_RA_Dec_Plot.png"), emit: radec_plot , optional: true
    path("*_PhiProjection.png"), emit: phi_projection_plot, optional: true
    path("*_ThetaProjection.png"), emit: theta_projection_plot, optional: true
    path("*_RAProjection.png"), emit: ra_projection_plot, optional: true
    path("*_DecProjection.png"), emit: dec_projection_plot, optional: true

    publishDir "${params.output_dir}/${output_folder}", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/PlotSkyMapsKM3NeT.py -i ${input_file} -o "./" -s "${source_specs_file}" --radius ${radius} --num_bins ${num_bins} --detectorname ${detectorname}
    """
}

process SummarizeEventFilesKM3NeT{
    input:
    val trigger

    script:
    """
    echo "Running final summary..."
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/SummarizeEventFilesKM3NeT.py
    """
}

process FindGTIsKM3NeTSingle{
    input:
    tuple path(input_file), val(detectorname), path(total_events_file)
    val dt_gtis

    output:
    path("*_gtis.pkl") ,  emit: gti_file

    publishDir "${params.output_dir}/gtis", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/FindGTIsKM3NeT.py -i ${input_file} -o "./" --df ${dt_gtis}
    """
}



process FindGTIsKM3NeTCombined{
    input:
    // tuple val(ratio), path(input_files), val(iteration), val(energy_min), val(gamma)
    path(input_files)
    val dt_gtis

    output:
    path("*combined_gtis.pkl") ,  emit: gti_file
    //tuple val(ratio), path(input_files), val(iteration), val(energy_min), val(gamma), path("*combined_gtis.pkl"), emit: gti_file_tuple
    path("*GTIs.png"), emit: gti_lightcurve_plot, optional: true

    publishDir "${params.output_dir}/gtis", mode: 'link', overwrite: true;

    script:
    def inputFilesString = input_files.collect { "-i ${it}" }.join(' ')
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/FindGTIsKM3NeT.py ${inputFilesString} -o "./" --df ${dt_gtis} --combine
    """
}

process ApplyCutsKM3NeT{
    input:
    tuple val(ratio), path(input_file), val(detectorname), val(iteration), val(energy_min_inj), val(gamma)
    val energy_min
    val energy_max
    val nhits_tr_min
    val nhits_sh_min
    val lik_tr_min
    val lik_sh_min
    val beta0_tr_max
    val zenith_min
    val trackscore_low_max
    val trackscore_high_min

    output:
    tuple val(ratio), path("*cutsapplied.hdf5"), val(detectorname), val(iteration), val(energy_min_inj), val(gamma), emit: cut_file

    publishDir "${params.output_dir}/cutfiles", mode: 'link', overwrite: true, enabled: false

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/ApplyCutsKM3NeT.py -i ${input_file} -o "./" --detector ${detectorname} \
        --energy_min ${energy_min} --energy_max ${energy_max} --nhits_tr_min ${nhits_tr_min} --nhits_sh_min ${nhits_sh_min} \
        --lik_tr_min ${lik_tr_min} --lik_sh_min ${lik_sh_min} --beta0_tr_max ${beta0_tr_max} --zenith_min ${zenith_min} \
        --trackscore_low_max ${trackscore_low_max} --trackscore_high_min ${trackscore_high_min}
                                                                                       
    """
}

process CalculateWeightsKM3NeT{
    input:
    tuple path(input_files), val(detectorname)

    output:
    path("*weights.hdf5"), emit: weights

    publishDir "${params.output_dir}/weights", mode: 'link', overwrite: true;

    script:
    def inputFilesString = input_files.collect { "-i ${it}" }.join(' ')
    """
    #python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/AtmosphericNeutrinoWeightsKM3NeT.py -i "${input_files}" -o "./"
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/AtmosphericNeutrinoWeightsKM3NeT.py ${inputFilesString} -o "./" --detector ${detectorname}
    """
}

process AddWeightsKM3NeT_old {
    errorStrategy 'ignore'

    input:
    tuple path(input_file), val(detectorname)
    val weights_folder

    output:
    tuple path("*with_weights.h5"), val(detectorname), emit: weighted_file

    publishDir "${params.output_dir}/input", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/AddWeightsKM3NeT.py -i "${input_file}" -p "${weights_folder}" -o "./" --detector ${detectorname}
    """
}

process AddWeightsKM3NeT {
    // errorStrategy 'ignore'

    input:
    tuple path(input_file), val(detectorname)
    path orca_weights_file
    path arca_weights_file

    output:
    tuple path("*with_weights.h5"), val(detectorname), emit: weighted_file

    publishDir "${params.output_dir}/input", mode: 'link', overwrite: true, enabled: false

    script:
    def weight_file = detectorname == 'orca' ? orca_weights_file : arca_weights_file
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/AddWeightsKM3NeT_new.py \
        -i "${input_file}" \
        -w "${weight_file}" \
        -o "./" \
        --detector ${detectorname}
    """
}

process SelectSubsetKM3NeT {

    input:
    tuple val(ratio), path(input_file), val(iteration), val(energy_min_inj), val(gamma), val(length) //, path(gti_file)
    path(gti_file)
    // val fraction
    
    output:
    tuple val(ratio), path("*days.hdf5"), path("*subset_gtis.pkl"), val(iteration), val(energy_min_inj), val(gamma), val(length), emit: subset_file

    publishDir "${params.output_dir}/subsets", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/SelectDataSetLengthKM3NeT.py -i ${input_file} -o "./" --gti_files ${gti_file} --length ${length}
    """

}


process InjectCutSelectFoldKM3NeT {

    input:
    tuple val(rate), path(input_file), val(detectorname), val(iteration), val(energy_min_inj), val(gamma), val(length)
    path(gti_file)

    val frequency
    val pulseshape
    val df
    val baseline
    val a
    val phi
    val kappa
    val method
    path energy_response_file

    val energy_min
    val energy_max
    val nhits_tr_min
    val nhits_sh_min
    val lik_tr_min
    val lik_sh_min
    val beta0_tr_max
    val zenith_min
    val trackscore_low_max
    val trackscore_high_min

    val number_of_testf
    val testdf
    val nbin
    val segment_size

    val run_name


    output:
    tuple val(rate), path("*ef.hdf5"), val(iteration), val(energy_min_inj), val(gamma), val(length), emit: efhdf5;
    path "*total_events.txt", emit: total_events_file
    path "*.png", emit: plot, optional: true

    publishDir "${params.output_dir}/epoch_folding", mode: 'link', overwrite: true, enabled: false

    script:

    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/InjectCutSelectFoldKM3NeT.py \
        -i "${input_file}" -o "./" -e ${energy_response_file} --rate ${rate} \
        --pulseshape ${pulseshape} --df ${df} --frequency ${frequency} --baseline ${baseline} --a ${a} --phi ${phi} --kappa ${kappa} \
        --method ${method} --E_min ${energy_min_inj} --gamma ${gamma} --detector ${detectorname} \
        --energy_min ${energy_min} --energy_max ${energy_max} --nhits_tr_min ${nhits_tr_min} --nhits_sh_min ${nhits_sh_min} \
        --lik_tr_min ${lik_tr_min} --lik_sh_min ${lik_sh_min} --beta0_tr_max ${beta0_tr_max} --zenith_min ${zenith_min} \
        --trackscore_low_max ${trackscore_low_max} --trackscore_high_min ${trackscore_high_min} \
        --gti_files ${gti_file} --length ${length} \
        --number_of_testf ${number_of_testf} --testdf ${testdf} --nbin ${nbin} --iteration ${iteration} --segment_size ${segment_size} \
        --run_name ${run_name}

    """
}

process RateOverGammaSensitivityPlotKM3NeT{

    input:
    tuple path(input_files), val(energy_min), val(gamma), val(length)

    output:
    path("*.png"), emit: sensitivity_plot
    path("*Gamma.hdf5"), emit: sensitivity_data


    publishDir "${params.output_dir}/sensitivity", mode: 'link', overwrite: true

    script:
    def inputFilesString = input_files.collect { "-i ${it}" }.join(' ')
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/RateOverGammaSensitivityPlotKM3NeT.py ${inputFilesString} -o "./" 
    """
}

process RateOverEminSensitivityPlotKM3NeT{

    input:
    tuple path(input_files), val(energy_min), val(gamma), val(length)

    output:
    path("*.png"), emit: sensitivity_plot
    path("*E_min.hdf5"), emit: sensitivity_data

    publishDir "${params.output_dir}/sensitivity", mode: 'link', overwrite: true

    script:
    def inputFilesString = input_files.collect { "-i ${it}" }.join(' ')
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/km3net/RateOverEminSensitivityPlotKM3NeT.py ${inputFilesString} -o "./" 
    """
}