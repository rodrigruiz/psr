nextflow.enable.dsl = 2

include{
    ConvertFilesKM3NeT;
    CreateEventListKM3NeT;
    CorrectEventListKM3NeT;
    InjectSignalKM3NeT;
    CombineEventListsKM3NeT;
    EpochFoldingKM3NeT;
    Chi2HistogramKM3NeT;
    SignalNoiseStatisticsKM3NeT;
} from './processes/km3net_processes.nf'

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
    // def snr_list = [0.05,0.2,0.4] //(input.ratio_min..input.ratio_max).step(input.ratio_step).toList()

    Channel
    .fromPath(input.km3net_files)
    .splitText(by: 1)
    .set {Files_Channel}

    SNR_Channel = Channel.fromList(snr_list)
    Iteration_Channel = Channel.fromList(iteration_list)
    //.view()

    // Combined_Channel = SNR_Channel.combine(Files_Channel).combine(Iteration_Channel)

    ConvertFilesKM3NeT(Files_Channel)
    CreateEventListKM3NeT(ConvertFilesKM3NeT.out, input.source_file, input.dist, input.energy_threshold)
    CorrectEventListKM3NeT(CreateEventListKM3NeT.out, input.source_file)

    Combined_Channel = SNR_Channel.combine(CorrectEventListKM3NeT.out).combine(Iteration_Channel)
    InjectSignalKM3NeT(Combined_Channel, input.frequency, input.pulseshape, input.df, input.baseline, input.a, input.phi, input.kappa)

    // InjectSignalKM3NeT(CorrectEventListKM3NeT.out, input.frequency, input.pulseshape, input.df, input.baseline, input.a, input.phi, input.kappa) //.out.injected_signal.groupTuple(by: 0).view().set { File_Collection }
    CombineEventListsKM3NeT(InjectSignalKM3NeT.out.groupTuple(by: [0,2]))
    EpochFoldingKM3NeT(CombineEventListsKM3NeT.out, input.frequency, input.number_of_testf, input.testf_df, input.nbin)
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