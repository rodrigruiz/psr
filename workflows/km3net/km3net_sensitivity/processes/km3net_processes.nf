process ConvertFilesKM3NeT{
 input:
 path input_files;

 output:
 path "*.h5";

 publishDir "${params.output_dir}/input", mode: 'link', overwrite: true;

 script:
 """
 python3 /home/hpc/capn/capn107h/software/psr/src/scripts/ConvertKM3NeTFiles.py -i "${input_files}" -o "./"
 """

}

process CreateEventListKM3NeT{
    input:
    path input_files;
    path source_specs_file;
    val dist;
    val energy_threshold;

    output:
    path '*eventlist.hdf5';

    publishDir "${params.output_dir}/eventlists", mode: 'link', overwrite: true;

    """
    python3  /home/hpc/capn/capn107h/software/psr/src/scripts/CreateEventListKM3NeT.py -i "${input_files}" -o "./" -s "${source_specs_file}" --energy_th ${energy_threshold} --dist ${dist}
    """
}

process CorrectEventListKM3NeT{
    input:
    path input_files;
    path source_specs_file;
    
    output:
    path "*corrected.hdf5";

    publishDir "${params.output_dir}/eventlists", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/CorrectEventListKM3NeT.py -i "${input_files}" -o "./" -s "${source_specs_file}"
    """
}

process InjectSignalKM3NeT{
    input:
    path input_files;
    val frequency;
    val ratio;
    val pulseshape;
    val df;
    val baseline;
    val a;
    val phi;
    val kappa;

    
    output:
    path "*signal*hdf5";

    publishDir "${params.output_dir}/signal", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/InjectSignalKM3NeT.py -i ${input_files} -o "./" --ratio ${ratio} --pulseshape ${pulseshape} --df ${df} --frequency ${frequency} --baseline ${baseline} --a ${a} --phi ${phi} --kappa ${kappa}
    """
}

process CombineEventListsKM3NeT{
    input:
    path input_files;
    
    output:
    path "*combined_eventlist.hdf5";

    publishDir "${params.output_dir}/combined_eventlists", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/CombineEventListsKM3NeT.py -i '*signal_mvm.hdf5' -o "./"
    """
}

process EpochFoldingKM3NeT{
    input:
    path input_file;
    val frequency;
    val number_of_testf;
    val df;
    val nbin;

    output:
    path "*epochfolding_results.hdf5", emit: hdf5;
    path "*.png", emit: plot;

    publishDir "${params.output_dir}/epoch_folding", mode: 'link', overwrite: true;

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/EpochFoldingKM3NeT.py -i "${input_file}" -o "./" --frequency ${frequency} --number_of_testf ${number_of_testf} --df ${df} --nbin ${nbin}
    """
}

