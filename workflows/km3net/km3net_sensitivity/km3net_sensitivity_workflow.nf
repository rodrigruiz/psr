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

    Channel
    .fromPath(input.km3net_files)
    .splitText(by: 1)
    .set {files}
    
    ConvertFilesKM3NeT(files);
    CreateEventListKM3NeT(ConvertFilesKM3NeT.out, input.source_file, input.dist, input.energy_threshold);
    CorrectEventListKM3NeT(CreateEventListKM3NeT.out, input.source_file);
    InjectSignalKM3NeT(CorrectEventListKM3NeT.out, input.frequency, input.ratio, input.pulseshape, input.df, input.baseline, input.a, input.phi, input.kappa);
    CombineEventListsKM3NeT(InjectSignalKM3NeT.out.collect());
    // EpochFoldingKM3NeT(CombineEventListsKM3NeT.out, input.frequency, input.number_of_testf, input.testf_df, input.nbin);

}
