process ConvertFilesKM3NeT{
    input:
    path input_files;

    output:
    path output_dir;

    script:
    """
    python3 ConvertFiles.py -i "${input_files}" -o "./"
    """

}

process CreateEventListKM3NeT{
    input:
    path input_files;
    path source_specs_file;
    val dist;
    val energy_threshold;

    output:
    path output_dir;

    """
    python3  CreateEventListKM3NeT.py -i "${input_files}" -o "${output_dir}" -s "${source_specs_file}" --energy_th ${energy_threshold} --dist ${dist}
    """
}

process CorrectEventListKM3NeT{
    input:
    path input_files;
    path source_specs_file;
    
    output:
    path output_dir;

    script:
    """
    python3 CorrectEventListKM3NeT.py -i "${input_files}" -o "${output_dir}" -s "${source_specs_file}"
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
    path output_dir;

    script:
    """
    python3 InjectSignalKM3NeT.py -i "${input_files}" -o "${output_dir}" --ratio ${ratio} --pulseshape ${pulseshape} --df ${df} --frequency ${frequency} --baseline ${baseline} --a ${a} --phi ${phi} --kappa ${kappa}
    """
}

process CombineEventListsKM3NeT{
    input:
    path input_files;
    
    output:
    path output_dir;

    script:
    """
    python3 CombineEventListsKM3NeT.py -i "${input_files}" -o "${output_dir}"
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
    path output_dir;

    script:
    """
    python3 EpochFoldingKM3NeT.py -i "${input_files}" -o "${output_dir}" --frequency ${frequency} --number_of_testf ${number_of_testf} --df ${df} --nbin ${nbin}
    """
}
