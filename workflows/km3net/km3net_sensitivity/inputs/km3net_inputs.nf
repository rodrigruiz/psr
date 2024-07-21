input = [dist : 70.0,
         energy_threshold : 0,
         bins_per_file  : 200,
         pulseshape : "mvm",
         frequency : 0.003,
         df : 1e-5,
         testf_df: 1e-5,
         ratio_min : 0.01, // Lower limit of the ratio
         ratio_max : 0.2, // Upper limit of the ratio
         ratio_step : 0.01, // Step for the ratio
         repetitions : 100, // Number of repetitions for each ratio
         number_of_testf : 200,
         nbin : 32,
         kappa : 5.0,
         a : 1.0,
         baseline : 0.0,
         phi : 0.0,
         nhbins : 13,
         source_file : "${WORKFLOW_DIR}/inputs/sources/Vela_X-1.h5",
	     km3net_files : "${WORKFLOW_DIR}/inputs/runs_test.txt",
         ]

         // /home/hpc/capn/capn107h/software/psr/workflows/km3net/km3net_sensitivity/inputs/runs.txt