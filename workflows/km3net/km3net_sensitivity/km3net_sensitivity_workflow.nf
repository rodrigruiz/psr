nextflow.enable.dsl = 2

include {
    ConvertFilesKM3NeT
    CreateEventListKM3NeT
    CorrectEventListKM3NeT
    InjectSignalKM3NeT
    CombineEventListsKM3NeT
    EpochFoldingKM3NeT
} from './processes/km3net_processes.nf'

evaluate(new File(params.input_file))

workflow runForRatio {
    take:
    file, ratio

    main:
    convert = ConvertFilesKM3NeT(file)
    event_list = CreateEventListKM3NeT(convert.out, input.source_file, input.dist, input.energy_threshold)
    corrected_list = CorrectEventListKM3NeT(event_list.out, input.source_file)
    injected_signal = InjectSignalKM3NeT(corrected_list.out, input.frequency, ratio, input.pulseshape, input.df, input.baseline, input.a, input.phi, input.kappa)
    combined_list = CombineEventListsKM3NeT(injected_signal.out.collect())
    epoch_folding = EpochFoldingKM3NeT(combined_list.out, input.frequency, input.number_of_testf, input.testf_df, input.nbin)

    emit:
    epoch_folding.out
}

workflow {
    def ratios = [0.1, 0.2, 0.3, 0.4]  // example ratio values

    Channel
    .fromPath(input.km3net_files)
    .splitText(by: 1)
    .set { files }
    
    files
    .combine(Channel.fromList(ratios))
    .map { file, ratio -> tuple(file, ratio) }
    .flatMap { params -> 
        runForRatio(file: params[0], ratio: params[1])
    }
}