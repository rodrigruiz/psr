nextflow.enable.dsl = 2

workflow {
    run_python_script()
}

process run_python_script {
   

    script:
    """
    python3 /home/hpc/capn/capn107h/software/psr/src/scripts/SensitvityStudyKM3NeT.py
    """
}
