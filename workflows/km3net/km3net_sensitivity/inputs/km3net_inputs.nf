input = [dist : 3.0,
         energy_threshold : 1e2,
         trackscore_threshold : 0.3,
         bins_per_file  : 200,
         pulseshape : "mvm",
         frequency : 0.003,
         df : 1e-5,
         testf_df: 1e-5,
         ratio_min : 0.01, // Lower limit of the ratio
         ratio_max : 0.76, // Upper limit of the ratio
         ratio_step : 0.05, // Step for the ratio
         repetitions : 30, // Number of repetitions for each ratio
         number_of_testf : 200,
         nbin : 32,
         kappa : 5.0,
         a : 1.0,
         baseline : 0.0,
         phi : 0.0,
         nhbins : 5,
         source_file : "${WORKFLOW_DIR}/inputs/sources/Vela_X-1.h5",
	     km3net_arca_numu_files : "${WORKFLOW_DIR}/inputs/runs_muons_test.txt",
         km3net_arca_anue_files : "${WORKFLOW_DIR}/inputs/runs_test.txt",
         km3net_orca_files : "${WORKFLOW_DIR}/inputs/orca_runs_muons_test.txt",
         ar_shower_file : "/home/hpc/capn/capn107h/software/nextflow_output/angular_resolutions/AngularResolutionOverEnergy_aashower_mcv8.1.gsg_anue-CCHEDIS_1e2-1e8GeV.sirene.jterbr00013.hdf5",
         ar_track_file : "/home/hpc/capn/capn107h/software/nextflow_output/angular_resolutions/AngularResolutionOverEnergy_jmuon_mcv8.1.gsg_numu-CCHEDIS_1e2-1e8GeV.sirene.jterbr000132.hdf5",
         ]

         // /home/hpc/capn/capn107h/software/psr/workflows/km3net/km3net_sensitivity/inputs/runs.txt