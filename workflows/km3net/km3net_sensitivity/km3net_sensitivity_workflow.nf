nextflow.enable.dsl = 2

include{
    ConvertFilesKM3NeT;
    CreateEventListKM3NeT;
    CorrectEventListKM3NeT;
    InjectSignalKM3NeT;
    CombineEventListsKM3NeT;
    EpochFoldingKM3NeT;
} from './processes/km3net_processes.nf'

evaluate(new File(params.input_file))

workflow{

    def snr_list = [0.05,0.1,0.2,0.3,0.4] //(input.ratio_min..input.ratio_max).step(input.ratio_step).toList()

    Channel
    .fromPath(input.km3net_files)
    .splitText(by: 1)
    .set {Files_Channel}

    SNR_Channel = Channel.fromList(snr_list)

    Combined_Channel = SNR_Channel.combine(Files_Channel)
    // .view()
    
    ConvertFilesKM3NeT(Combined_Channel);

    CreateEventListKM3NeT(ConvertFilesKM3NeT.out, input.source_file, input.dist, input.energy_threshold)
    CorrectEventListKM3NeT(CreateEventListKM3NeT.out, input.source_file)
    InjectSignalKM3NeT(CorrectEventListKM3NeT.out, input.frequency, input.pulseshape, input.df, input.baseline, input.a, input.phi, input.kappa) //.out.injected_signal.groupTuple(by: 0).view().set { File_Collection }
    //CombineEventListsKM3NeT(File_Collection)
    //InjectSignalKM3NeT.out.groupTuple(by: 0).view()
    CombineEventListsKM3NeT(InjectSignalKM3NeT.out.groupTuple(by: 0))
    EpochFoldingKM3NeT(CombineEventListsKM3NeT.out, input.frequency, input.number_of_testf, input.testf_df, input.nbin)

} 