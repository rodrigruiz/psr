input = [km3net_arca_energy_threshold : 1e2,
         km3net_arca_energy_low : 2,
         km3net_arca_energy_high : 8,
         km3net_arca_trackscore_threshold : 0.3,
         km3net_orca_energy_threshold : 0,
         km3net_orca_energy_low : 0,
         km3net_orca_energy_high : 2,
         km3net_orca_trackscore_threshold : 0.3,
         bins_per_file  : 200,
         pulseshape : "mvm",
         frequency : 0.003,
         df : 1e-5,
         testf_df: 1e-5,
         ratio_min : 0.01, // Lower limit of the ratio
         ratio_max : 0.46, // Upper limit of the ratio
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
	     km3net_arca_numu_files : "${WORKFLOW_DIR}/inputs/arca_runs_numu_test.txt",
         km3net_arca_anue_files : "${WORKFLOW_DIR}/inputs/arca_runs_anue_test.txt",
         km3net_orca_numu_files : "${WORKFLOW_DIR}/inputs/orca_runs_numu_test.txt",
         km3net_orca_anue_files : "${WORKFLOW_DIR}/inputs/orca_runs_anue_test.txt",
         ar_shower_file_arca : "${WORKFLOW_DIR}/inputs/angular_res_files/AngularResolutionOverEnergy_aashower_mcv8.1.gsg_anue-CCHEDIS_1e2-1e8GeV.sirene.jterbr00013.hdf5",
         ar_track_file_arca : "${WORKFLOW_DIR}/inputs/angular_res_files/AngularResolutionOverEnergy_jmuon_mcv8.1.gsg_numu-CCHEDIS_1e2-1e8GeV.sirene.jterbr000132.hdf5",
         ar_shower_file_orca : "${WORKFLOW_DIR}/inputs/angular_res_files/AngularResolutionOverEnergy_jshower_gsg_elec-CC_1.0-100.0GeV_orca.hdf5",
         ar_track_file_orca : "${WORKFLOW_DIR}/inputs/angular_res_files/AngularResolutionOverEnergy_jmuon_gsg_muon-CC_1.0-100.0GeV_orca.hdf5",
         ]

         // /home/hpc/capn/capn107h/software/psr/workflows/km3net/km3net_sensitivity/inputs/runs.txt