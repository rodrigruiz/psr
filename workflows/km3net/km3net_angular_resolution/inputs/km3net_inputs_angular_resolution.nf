input = [source_file : "${WORKFLOW_DIR}/inputs/sources/Vela_X-1.h5",
	     km3net_arca_files : "${WORKFLOW_DIR}/inputs/runs_muons_test.txt",
         km3net_orca_files : "${WORKFLOW_DIR}/inputs/orca_runs_anue_test.txt",
         detector : 'orca',
         runtype : 'anue',
         recotype : 'jshower',
         pltscale : 'linear',
         energy_low : 1,
         energy_high : 1e2,
         n_energy_bins : 37,
         ]

         // /home/hpc/capn/capn107h/software/psr/workflows/km3net/km3net_sensitivity/inputs/runs.txt