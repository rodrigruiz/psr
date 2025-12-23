nextflow.enable.dsl = 2

include {
    InjectSignalKM3NeT;
    CombineEventListsDetectors;
    EpochFoldingKM3NeT;
    Chi2HistogramKM3NeT;
    SignalNoiseStatisticsKM3NeT;
    ApplyCutsKM3NeT;
    SummarizeEventFilesKM3NeT as SummarizeEventFiles;
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

    // def param_file = file(params.parameter_file)

    Channel
        .fromPath("${params.output_dir}/*corrected.hdf5")
        .map { file -> tuple(file, file.name.contains("orca") ? "orca" : "arca", file("${params.output_dir}/${file.baseName.replace('.hdf5','')}_total.txt")) }
        .set { CorrectedFiles_Channel }

    /*
    Channel
        .from(param_file.readLines())
        .map { line ->
            def (ratio, iteration) = line.tokenize()
            tuple(ratio.toDouble(), iteration.toInteger())
        }
        .set { ParamChannel }
    */

    SNR_Channel = Channel.fromList(snr_list)
    Iteration_Channel = Channel.fromList(iteration_list)
    Combined_Channel = SNR_Channel.combine(CorrectedFiles_Channel).combine(Iteration_Channel)

    InjectSignalKM3NeT(Combined_Channel, input.frequency, input.pulseshape, input.df, input.baseline, input.a, input.phi, input.kappa, input.inject_method)
    ApplyCutsKM3NeT(InjectSignalKM3NeT.out.injected_signal, 0, 1e10, 1, 1, 0, 0, 1, 0, 0.5, 0.5)

    CombineEventListsDetectors(ApplyCutsKM3NeT.out.cut_file.groupTuple(by: [0, 3]), input.source_file, input.delta_search_min, input.filestype)

    EpochFoldingKM3NeT(
        CombineEventListsDetectors.out.combined_file,
        input.frequency,
        input.number_of_testf,
        input.testf_df,
        input.nbin,
        input.folding_segment_size,
        file(input.gti_file)
    )

    Chi2HistogramKM3NeT(EpochFoldingKM3NeT.out.hdf5.groupTuple(by: 0), input.nhbins)

    SignalNoiseStatisticsKM3NeT(
        Chi2HistogramKM3NeT.out.hdf5.collect(),
        input.nbin,
        input.frequency,
        input.delta_search_min,
        CombineEventListsDetectors.out.total_events.last()
    )

    def summary_trigger = SignalNoiseStatisticsKM3NeT.out[0].collect().map { "done" }
    SummarizeEventFiles(summary_trigger)
}
