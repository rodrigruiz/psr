process ConvertFilesKM3NeT{
    input:
    tuple val(ratio), path(input_file)

    output:
    tuple val(ratio), path("*.h5"), emit: converted_file

    publishDir "${params.output_dir}/input", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/ConvertKM3NeTFiles.py -i "${input_file}" -o "./"
    """

}

process CreateEventListKM3NeT{
    input:
    tuple val(ratio), path(input_file)
    path source_specs_file
    val dist
    val energy_threshold

    output:
    tuple val(ratio), path('*eventlist.hdf5'), emit: eventlist

    publishDir "${params.output_dir}/eventlists", mode: 'link', overwrite: true;

    """
    python3  /home/hpc/capn/capn107h/software/psr/src/scripts/CreateEventListKM3NeT.py -i "${input_file}" -o "./" -s "${source_specs_file}" --energy_th ${energy_threshold} --dist ${dist}
    """
}

process CorrectEventListKM3NeT{
    input:
    tuple val(ratio), path(input_file)
    path source_specs_file
    
    output:
    tuple val(ratio), path("*corrected.hdf5"), emit: corrected_eventlist

    publishDir "${params.output_dir}/eventlists", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/CorrectEventListKM3NeT.py -i "${input_file}" -o "./" -s "${source_specs_file}"
    """
}

process InjectSignalKM3NeT{
    input:
    tuple val(ratio), path(input_file)
    val frequency
    val pulseshape
    val df
    val baseline
    val a
    val phi
    val kappa

    
    output:
    tuple val(ratio), path("*signal*hdf5"), emit: injected_signal

    publishDir "${params.output_dir}/signal", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/InjectSignalKM3NeT.py -i ${input_file} -o "./" --ratio ${ratio} --pulseshape ${pulseshape} --df ${df} --frequency ${frequency} --baseline ${baseline} --a ${a} --phi ${phi} --kappa ${kappa}
    """
}

process CombineEventListsKM3NeT{
    input:
    tuple val(ratio), path(input_files)
    
    output:
    tuple val(ratio), path("*combined_eventlist.hdf5"), emit: combined_file

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
    tuple val(ratio), path(input_file)
    val frequency
    val number_of_testf
    val df
    val nbin

    output:
    path "*epochfolding_results.hdf5", emit: hdf5;
    path "*.png", emit: plot;

    publishDir "${params.output_dir}/epoch_folding", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/EpochFoldingKM3NeT.py -i "${input_file}" -o "./" --frequency ${frequency} --number_of_testf ${number_of_testf} --df ${df} --nbin ${nbin} --ratio ${ratio}
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

