input = [source_file : "${WORKFLOW_DIR}/inputs/sources/Vela_X-1.h5",
	     km3net_arca_files : "${WORKFLOW_DIR}/../mc_files_arcav9.txt",
         km3net_orca_files : "${WORKFLOW_DIR}/../mc_files_orcav9.txt",
         km3net_files: "${WORKFLOW_DIR}/inputs/orca_runs_anue_test.txt",
         detector : 'arca',
         runtype : 'anue',
         recotype : 'jshower',
         pltscale : 'linear',
         energy_low : 1e-1,
         energy_high : 1e3,
         n_energy_bins : 50,
         ]

         // /home/hpc/capn/capn107h/software/psr/workflows/km3net/km3net_sensitivity/inputs/runs.txt