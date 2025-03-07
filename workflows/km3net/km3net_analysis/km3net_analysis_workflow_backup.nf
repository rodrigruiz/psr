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

    Channel
    .fromPath(input.km3net_arca_numu_files)
    .splitText(by: 1)
    .combine(Channel.of('arca'))
    .combine(Channel.of('aashower'))
    .combine(Channel.of(input.km3net_arca_energy_threshold))
    .combine(Channel.of(input.km3net_arca_energy_low))
    .combine(Channel.of(input.km3net_arca_energy_high))
    //.combine(Channel.of(input.ar_shower_file_arca))
    //.combine(Channel.of(input.ar_track_file_arca))
    .combine(Channel.of(input.km3net_arca_trackscore_threshold))
    .set {ARCA_numu_Files_Channel}

    Channel
    .fromPath(input.km3net_arca_anue_files)
    .splitText(by: 1)
    .combine(Channel.of('arca'))
    .combine(Channel.of('aashower'))
    .combine(Channel.of(input.km3net_arca_energy_threshold))
    .combine(Channel.of(input.km3net_arca_energy_low))
    .combine(Channel.of(input.km3net_arca_energy_high))
    //.combine(Channel.of(input.ar_shower_file_arca))
    //.combine(Channel.of(input.ar_track_file_arca))
    .combine(Channel.of(input.km3net_arca_trackscore_threshold))
    .mix(ARCA_numu_Files_Channel)
    .set {ARCA_Files_Channel}
    
    /*
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


    SNR_Channel = Channel.fromList(snr_list)
    Iteration_Channel = Channel.fromList(iteration_list)

    ConvertFilesKM3NeT(ARCA_Files_Channel)

    AddTrackScoreKM3NeT(ConvertFilesKM3NeT.out, input.parampid_folder)

    BlindDataKM3NET(AddTrackScoreKM3NeT.out)
    CreateEventListKM3NeT_new(BlindDataKM3NET.out, input.source_file, input.delta_search_min, input.ar_shower_file_arca, input.ar_track_file_arca, input.ar_shower_file_orca, input.ar_track_file_orca)
    CorrectEventListKM3NeT(CreateEventListKM3NeT_new.out, input.source_file)  

    Combined_Channel = SNR_Channel.combine(CorrectEventListKM3NeT.out).combine(Iteration_Channel)
    InjectSignalKM3NeT(Combined_Channel, input.frequency, input.pulseshape, input.df, input.baseline, input.a, input.phi, input.kappa)
    CombineEventListsKM3NeT(InjectSignalKM3NeT.out.groupTuple(by: [0,2]))

    EpochFoldingKM3NeT(CombineEventListsKM3NeT.out, input.frequency, input.number_of_testf, input.testf_df, input.nbin, input.folding_segment_size)
    Chi2HistogramKM3NeT(EpochFoldingKM3NeT.out.hdf5.groupTuple(by: 0), input.nhbins)
    SignalNoiseStatisticsKM3NeT(Chi2HistogramKM3NeT.out.hdf5.collect(), input.nbin)

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