nextflow.enable.dsl = 2

include{
    ConvertFilesKM3NeT;
    BlindDataKM3NET;
    AddTrackScoreKM3NeT;
    CreateEventListKM3NeT_new as CreateSelectedEventList;
    CreateEventListKM3NeT_new as CreateFullSkyEventList;
    CorrectEventListKM3NeT;
    InjectSignalKM3NeT;
    CombineEventListsRuns as CombineEventListsRunsSelected;
    CombineEventListsRuns as CombineEventListsRunsFullSky;
    CombineEventListsDetectors;
    EpochFoldingKM3NeT;
    Chi2HistogramKM3NeT;
    SignalNoiseStatisticsKM3NeT;
    PlotSkyMapsKM3NeT as PlotSkyMapsKM3NeTSelected;
    PlotSkyMapsKM3NeT as PlotSkyMapsKM3NeTFullSky;
    PlotSkyMapsKM3NeTCombined;
    SummarizeEventFilesKM3NeT as SummarizeEventFiles;
    FindGTIsKM3NeTSingle;
    FindGTIsKM3NeTCombined;
    ApplyCutsKM3NeT;
    AddWeightsKM3NeT;
    SelectSubsetKM3NeT;
} from '../processes/km3net_processes.nf'

evaluate(new File(params.input_file))
input.delta_search_min = params.delta_search_min ?: input.delta_search_min
input.filestype        = params.filestype ?: input.filestype
input.source_file      = params.source_file ?: input.source_file
input.cone_all         = params.cone_all ?: input.cone_all

def generateList(start, end, step) {
    def list = []
    for (def i = start; i <= end; i += step) {
        list.add(i)
    }
    return list
}


workflow{
    
    snr_list = generateList(input.ratio_min, input.ratio_max, input.ratio_step)
    rate_list = [ 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]
    // rate_list = [5.0, 10.0, 25.0, 50.0]
    // rate_list = [5.0]
    iteration_list = (1..input.repetitions).toList()
    //e_min_list = [5e4, 8e4, 1e5, 2e5, 5e5, 6e5, 1e6, 2e6]
    //gamma_list = [1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6]
    
    e_min_list = [1e3]
    //gamma_list = [1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0]
    gamma_list = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0] 
    //length_list = [5.0, 10.0, 15.0, 20.0, 25.0]
    //length_list = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0]
    length_list = [60.0]

    /*      
    Channel
    .fromPath(input.km3net_arca_numu_files)
    .splitText(by: 1)
    .combine(Channel.of('arca'))
    //.combine(Channel.of('aashower'))
    //.combine(Channel.of(input.km3net_arca_energy_threshold))
    //.combine(Channel.of(input.km3net_arca_energy_low))
    //.combine(Channel.of(input.km3net_arca_energy_high))
    //.combine(Channel.of(input.ar_shower_file_arca))
    //.combine(Channel.of(input.ar_track_file_arca))
    //.combine(Channel.of(input.km3net_arca_trackscore_threshold))
    .set {ARCA_numu_Files_Channel}

    Channel
    .fromPath(input.km3net_arca_anue_files)
    .splitText(by: 1)
    .combine(Channel.of('arca'))
    //.combine(Channel.of('aashower'))
    //.combine(Channel.of(input.km3net_arca_energy_threshold))
    //.combine(Channel.of(input.km3net_arca_energy_low))
    //.combine(Channel.of(input.km3net_arca_energy_high))
    //.combine(Channel.of(input.ar_shower_file_arca))
    //.combine(Channel.of(input.ar_track_file_arca))
    //.combine(Channel.of(input.km3net_arca_trackscore_threshold))
    .mix(ARCA_numu_Files_Channel)
    .set {ARCA_Files_Channel}
    
    Channel
    .fromPath(input.km3net_orca_numu_files)
    .splitText(by: 1)
    .combine(Channel.of('orca'))
    .combine(Channel.of('jshower'))
    .combine(Channel.of(input.km3net_orca_energy_threshold))
    .combine(Channel.of(input.km3net_orca_energy_low))
    .combine(Channel.of(input.km3net_orca_energy_high))
    //.combine(Channel.of(input.ar_shower_file_orca))
    //.combine(Channel.of(input.ar_track_file_orca))
    .combine(Channel.of(input.km3net_orca_trackscore_threshold))   
    .set {ORCA_numu_Files_Channel}

    Channel
    .fromPath(input.km3net_orca_anue_files)
    .splitText(by: 1)
    .combine(Channel.of('orca'))
    .combine(Channel.of('jshower'))
    .combine(Channel.of(input.km3net_orca_energy_threshold))
    .combine(Channel.of(input.km3net_orca_energy_low))
    .combine(Channel.of(input.km3net_orca_energy_high))
    //.combine(Channel.of(input.ar_shower_file_orca))
    //.combine(Channel.of(input.ar_track_file_orca))
    .combine(Channel.of(input.km3net_orca_trackscore_threshold)) 
    .mix(ORCA_numu_Files_Channel)
    .mix(ARCA_Files_Channel)
    .set {Files_Channel}
    */

    /* WORKS 2025-04-09
    Channel
    .fromPath(input.km3net_arca_data_files)
    .splitText(by: 1)
    .combine(Channel.of('arca'))
    .set {ARCA_Files_Channel}
    */
    /*
    Channel
    .fromPath(input.km3net_orca_data_files)
    .splitText(by: 1)
    .combine(Channel.of('arca'))
    .mix(ARCA_Files_Channel)
    .set {Files_Channel}
    */



    if (input.filestype == 'data') {
        Channel
            .fromPath(input.km3net_arca_data_files)
            .splitText(by: 1)
            .combine(Channel.of('arca'))
            .set { ARCA_Files_Channel }
    } else if (input.filestype == 'mc') {
        Channel
            .fromPath(input.km3net_arca_mc_files)
            .splitText(by: 1)
            .combine(Channel.of('arca'))
            .set { ARCA_Files_Channel }
    } else {
        error "Unsupported input.filestype: ${input.filestype}"
    }

    if (input.filestype == 'data') {
        Channel
            .fromPath(input.km3net_orca_data_files)
            .splitText(by: 1)
            .combine(Channel.of('orca'))
            .mix(ARCA_Files_Channel)
            .set { Files_Channel }
    } else if (input.filestype == 'mc') {
        Channel
            .fromPath(input.km3net_orca_mc_files)
            .splitText(by: 1)
            .combine(Channel.of('orca'))
            .mix(ARCA_Files_Channel)
            .set { Files_Channel }
    } else {
        error "Unsupported input.filestype: ${input.filestype}"
    }

    // SNR_Channel = Channel.fromList(snr_list)
    SNR_Channel = Channel.fromList(rate_list)
    Iteration_Channel = Channel.fromList(iteration_list)
    Emin_Channel = Channel.fromList(e_min_list)
    Gamma_Channel = Channel.fromList(gamma_list)
    // Convert files and add track score
    
    // ConvertFilesKM3NeT(ARCA_Files_Channel) //.view()
    // AddTrackScoreKM3NeT(ConvertFilesKM3NeT.out, input.parampid_folder)

    AddTrackScoreKM3NeT(ARCA_Files_Channel, input.parampid_folder)

    // Add weights for MC only
    if (input.filestype == 'mc') {
        AddWeightsKM3NeT(AddTrackScoreKM3NeT.out, input.orca_weights_file, input.arca_weights_file)
        eventlist_input = AddWeightsKM3NeT.out
    } else {
        eventlist_input = AddTrackScoreKM3NeT.out
    }

    CreateFullSkyEventList(
    eventlist_input,
    input.source_file,
    180.0,
    input.ar_shower_file_arca,
    input.ar_track_file_arca,
    input.ar_shower_file_orca,
    input.ar_track_file_orca,
    input.energy_low_arca,
    input.energy_high_arca,
    input.energy_low_orca,
    input.energy_high_orca,
    input.energy_threshold_arca,
    input.trackscore_threshold_arca,
    input.muonscore_threshold_arca,
    input.energy_threshold_orca,
    input.trackscore_threshold_orca,
    input.muonscore_threshold_orca,
    'True',
    )

    // Blind the data
    BlindDataKM3NET(eventlist_input)
    
    // Create Event Lists
    //CreateEventListKM3NeT_new(BlindDataKM3NET.out, input.source_file, input.delta_search_min, input.ar_shower_file_arca, input.ar_track_file_arca, input.ar_shower_file_orca, input.ar_track_file_orca)
    CreateSelectedEventList(
    BlindDataKM3NET.out,
    //AddTrackScoreKM3NeT.out,
    input.source_file,
    input.delta_search_min,
    input.ar_shower_file_arca,
    input.ar_track_file_arca,
    input.ar_shower_file_orca,
    input.ar_track_file_orca,
    input.energy_low_arca,
    input.energy_high_arca,
    input.energy_low_orca,
    input.energy_high_orca,
    input.energy_threshold_arca,
    input.trackscore_threshold_arca,
    input.muonscore_threshold_arca,
    input.energy_threshold_orca,
    input.trackscore_threshold_orca,
    input.muonscore_threshold_orca,
    input.cone_all,
    )

    CombineEventListsRunsFullSky(CreateFullSkyEventList.out.groupTuple(by: 1),input.source_file, input.delta_search_min, input.filestype)
    if (input.plot_all == 'True') {
        PlotSkyMapsKM3NeTFullSky(CombineEventListsRunsFullSky.out.combined_file, input.source_file, "combined_eventlists", input.delta_search_min, 120)
    }
    FindGTIsKM3NeTSingle(CombineEventListsRunsFullSky.out.combined_file, 600)
    //CombineEventListsRunsFullSky.out.combined_file.view()
    //CombineEventListsRunsFullSky.out.combined_file.collect().view()
    CombineEventListsRunsFullSky.out.combined_file
        .map { it instanceof List || it instanceof Tuple ? it[0] : it }
        .collect()
        .set{ CombinedEventlistsOutput}
    // CombinedEventlistsOutput.view()
    FindGTIsKM3NeTCombined(CombinedEventlistsOutput,600)
    
    // Apply correction to the combined event lists
    //CorrectEventListKM3NeT(CreateEventListKM3NeT_new.out, input.source_file)  

        // Group eventlists by detectorname
    //CreateEventListKM3NeT_new.out.groupTuple(by: 1).set { EventList_Grouped_By_Detector }

    //EventList_Grouped_By_Detector.map { eventlist_files, detectorname  -> 
    //tuple(eventlist_files, input.source_file, input.delta_search_min, input.filestype, detectorname)
    //}.set { CombinedInputChannel }
    //CombineEventListsKM3NeT(CombinedInputChannel)
    //CombineEventListsKM3NeT(CreateEventListKM3NeT_new.out.collect(),input.source_file, input.delta_sesarch_min, input.filestype, input.detector)
    //CreateEventListKM3NeT_new.out.groupTuple(by: 1) //.view()
    CombineEventListsRunsSelected(CreateSelectedEventList.out.groupTuple(by: 1),input.source_file, input.delta_search_min, input.filestype)
    PlotSkyMapsKM3NeTSelected(CombineEventListsRunsSelected.out.combined_file, input.source_file, "combined_eventlists", input.delta_search_min, 120)
    CorrectEventListKM3NeT(CombineEventListsRunsSelected.out.combined_file, input.source_file) 

    Combined_Channel = SNR_Channel.combine(CorrectEventListKM3NeT.out).combine(Iteration_Channel).combine(Emin_Channel).combine(Gamma_Channel)
    // Combined_Channel.view()
    /*
    Combined_Channel = SNR_Channel
        .combine(CorrectEventListKM3NeT.out)              // (ratio, file, det)
        .map { ratio, file, det -> tuple(ratio, file, det) }
        .combine(Iteration_Channel)                       // ((ratio, file, det), iter)
        .map { it -> tuple(it[0][0], it[0][1], it[0][2], it[1]) } // ✅ Fix: wrap in `tuple(...)`
        .join(CombineEventListsRuns.out.total_events) { a, b -> a[2] == b[0] }  // match det
        .map { ratio, file, det, iter, _, total_file -> 
            tuple(ratio, file, det, iter, total_file) 
        }
    */
    //Combined_Channel.view()
    
    //FindGTIsSingle(CorrectEventListKM3NeT.out, 600)
    // Maybe new channel has to be created first which is then further processed?

    InjectSignalKM3NeT(Combined_Channel, input.frequency, input.pulseshape, input.df, input.baseline, input.a, input.phi, input.kappa, input.inject_method)
    //PlotSkyMapsKM3NeTInjected(InjectSignalKM3NeT.out.injected_signal, input.source_file, "signal", input.delta_search_min, 120)
    // Skymap plot fails because of output tuple from injectsignals, expects only single filepath, check for epochfolding process inputs
    // Grouped_Channel = InjectSignalKM3NeT.out.injected_signal.groupTuple(ratio, injected_file, detectorname, iteration -> tuple(ratio, iteration) ).view()
    


    /*
    InjectSignalKM3NeT.out
        .map { tuple(ratio, injected_file, detectorname, iteration) -> tuple(tuple(ratio, iteration), injected_file) }
        .groupTuple()
        .map { key, files ->
            def ratio = key[0]
            def iteration = key[1]
            def combined_file = files.flatten()  // all files for that (ratio, iteration)

            // Since you have multiple files, you probably want to pass a *list of files* to EpochFolding
            tuple(ratio, combined_file, iteration)
        }
        .set { Injected_Signals_Channel }

    CombineEventListsDetectors(Injected_Signals_Channel, input.source_file, input.delta_search_min, input.filestype)
    EpochFoldingKM3NeT(CombineEventListsDetectors, input.frequency, input.number_of_testf, input.testf_df, input.nbin, input.folding_segment_size)
    */
    // InjectSignalKM3NeT.out.injected_signal.view()
    //InjectSignalKM3NeT.out.injected_signal.groupTuple(by:[0,3]).view()
    ApplyCutsKM3NeT(InjectSignalKM3NeT.out.injected_signal,0,1e10,1,1,0,0,1,0,0.5,0.5)
    //ApplyCutsKM3NeT.out.cut_file.view()
    //ApplyCutsKM3NeT.out.cut_file.groupTuple(by:[0,3,4,5])
    CombineEventListsDetectors(ApplyCutsKM3NeT.out.cut_file.groupTuple(by:[0,3,4,5]),input.source_file, input.delta_search_min, input.filestype)
    // CombineEventListsDetectors(InjectSignalKM3NeT.out.injected_signal.groupTuple(by:[0,3]),input.source_file, input.delta_search_min, input.filestype)
    //PlotSkyMapsKM3NeTCombined(CombineEventListsDetectors.out.combined_file,'arcaorca' , input.source_file, "combined_eventlists", input.delta_search_min, 120)
    //PlotSkyMapsKM3NeTCombined(CombineEventListsDetectors.out.combined_file, 'combined' ,input.source_file, "combined_eventlists", input.delta_search_min, 120)
    // FindGTIsKM3NeTCombined(CombinedEventListsDetectors.out.combined_file,600)

    // 
    // INSERT SPLIT INTO DIFFERENT DATASET LENGTHS HERE:
    // TAKE OUTPUT OF COMBINED EVENTLISTS AND MAKE NEW 'LENGTH' CHANNEL
    // NEW SCRIPT TAKES THE CONBINED EVENTLIST AND LENGTH AS INPUT AND CUTS SUBSET OF CORRESPONDING LENGTH OUT OF SET
    // AND FEEDS IT INTO THE REMAINING WORKFLOW, NAME HAS TO STORE THE USED LENGTH
    // ALLOWS TO INVESTIGATE HOW SENSITIVITY CHANGES WITH DATASET LENGTH FOR GIVEN RATE -> 3D PLOT
    //
    Length_Channel = Channel.fromList(length_list)
    SubsetChannel =  CombineEventListsDetectors.out.combined_file.combine(Length_Channel) //.combine(FindGTIsKM3NeTCombined.out.gti_file)
    SelectSubsetKM3NeT(SubsetChannel, FindGTIsKM3NeTCombined.out.gti_file)
    EpochFoldingKM3NeT(SelectSubsetKM3NeT.out.subset_file, input.frequency, input.number_of_testf, input.testf_df, input.nbin, input.folding_segment_size) //, FindGTIsKM3NeTCombined.out.gti_file)
    // EpochFoldingKM3NeT(CombineEventListsDetectors.out.combined_file, input.frequency, input.number_of_testf, input.testf_df, input.nbin, input.folding_segment_size) //, FindGTIsKM3NeTCombined.out.gti_file)
    //Combined_Channel = SNR_Channel.combine(CorrectEventListKM3NeT.out).combine(Iteration_Channel)
    //InjectSignalKM3NeT(Combined_Channel, input.frequency, input.pulseshape, input.df, input.baseline, input.a, input.phi, input.kappa)
    //CombineEventListsKM3NeT(InjectSignalKM3NeT.out.groupTuple(by: [0,2]))
    //EpochFoldingKM3NeT(CombineEventListsKM3NeT.out, input.frequency, input.number_of_testf, input.testf_df, input.nbin, input.folding_segment_size)
    //EpochFoldingKM3NeT.out.hdf5.groupTuple(by: [0, 3]).view()
    Chi2HistogramKM3NeT(EpochFoldingKM3NeT.out.hdf5.groupTuple(by: [0, 3, 4, 5]), input.nhbins)

    SignalNoiseStatisticsKM3NeT(Chi2HistogramKM3NeT.out.hdf5.groupTuple(by: [1, 2, 3]), input.nbin, input.frequency, input.delta_search_min, InjectSignalKM3NeT.out.total_events.last() )
    // CHANGE HERE HOW THE TOTAL EVENTS FILE FROM INJECTSIGNAL IS PARSED
    // NOT DOING WHAT IT SHOULD
    // MAYBE PASS THROUGH EPOCH FOLDING AND SO ON
    def summary_trigger = SignalNoiseStatisticsKM3NeT.out[0].collect().map { "done" }

    SummarizeEventFiles(summary_trigger)
    //SummarizeEventFiles()
}
workflow.onComplete = {
  println "Pipeline complete"
  println "Command line: $workflow.commandLine"
  def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

    sendMail(to: 'hwarnhofer@km3net.de', subject: 'My pipeline execution', body: msg)
}

workflow.onError = {
  println "Oops... something went wrong"
}