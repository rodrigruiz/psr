nextflow.enable.dsl = 2

include{
    ConvertFilesKM3NeT;
    BlindDataKM3NET;
    AddTrackScoreKM3NeT;
    CreateEventListKM3NeT_new;
    CorrectEventListKM3NeT;
    InjectSignalKM3NeT;
    CombineEventListsKM3NeT;
    EpochFoldingKM3NeT;
    Chi2HistogramKM3NeT;
    SignalNoiseStatisticsKM3NeT;
} from '../processes/km3net_processes.nf'

evaluate(new File(params.input_file))

def generateList(start, end, step) {
    def list = []
    for (def i = start; i <= end; i += step) {
        list.add(i)
    }
    return list
}


workflow{
    
    snr_list = generateList(input.ratio_min, input.ratio_max, input.ratio_step)
    iteration_list = (1..input.repetitions).toList()


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


    Channel
    .fromPath(input.km3net_arca_data_files)
    .splitText(by: 1)
    .combine(Channel.of('arca'))
    .set {ARCA_Files_Channel}

    /*
    Channel
    .fromPath(input.km3net_orca_data_files)
    .splitText(by: 1)
    .combine(Channel.of('arca'))
    .mix(ARCA_Files_Channel)
    .set {Files_Channel}
    */

    SNR_Channel = Channel.fromList(snr_list)
    Iteration_Channel = Channel.fromList(iteration_list)

    // Convert files and add track score
    ConvertFilesKM3NeT(ARCA_Files_Channel)
    AddTrackScoreKM3NeT(ConvertFilesKM3NeT.out, input.parampid_folder)

    // Blind the data
    // BlindDataKM3NET(AddTrackScoreKM3NeT.out)
    
    // Create Event Lists
    //CreateEventListKM3NeT_new(BlindDataKM3NET.out, input.source_file, input.delta_search_min, input.ar_shower_file_arca, input.ar_track_file_arca, input.ar_shower_file_orca, input.ar_track_file_orca)
    CreateEventListKM3NeT_new(
    AddTrackScoreKM3NeT.out,
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
    input.energy_threshold_orca,
    input.energy_threshold_orca
    )

    // Apply correction to the combined event lists
    CorrectEventListKM3NeT(CreateEventListKM3NeT_new.out, input.source_file)  
    CombineEventListsKM3NeT(CorrectEventListKM3NeT.out.collect(),input.source_file)

    // Create a combined channel (SNR, Corrected Event List, Iteration info)
    Combined_Channel = SNR_Channel.combine(CombineEventListsKM3NeT.out.combined_file).combine(Iteration_Channel)
    
    // Inject signal into the combined channel
    // InjectSignalKM3NeT(Combined_Channel, input.frequency, input.pulseshape, input.df, input.baseline, input.a, input.phi, input.kappa)
    
    // Combine the event lists after signal injection
    // CombineEventListsKM3NeT(InjectSignalKM3NeT.out.groupTuple(by: [0,2]))

    // Epoch folding and Chi2 histogram
    EpochFoldingKM3NeT(Combined_Channel, input.frequency, input.number_of_testf, input.testf_df, input.nbin, input.folding_segment_size)
    //Chi2HistogramKM3NeT(EpochFoldingKM3NeT.out.hdf5.groupTuple(by: 0), input.nhbins)
    //SignalNoiseStatisticsKM3NeT(Chi2HistogramKM3NeT.out.hdf5.collect(), input.nbin)
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