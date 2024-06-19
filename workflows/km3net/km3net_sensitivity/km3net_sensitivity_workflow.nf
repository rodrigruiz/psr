nextflow.enable.dsl = 2

params.input_files = '/home/hpc/capn/mppi104h/wecapstor3/out/ARCA/KM3NeT_00000133/v8.1/reco/*11.root'
params.source_file = '/home/hpc/capn/capn107h/software/hdf5SourceFiles/Vela_X-1.h5'
params.output_dir = '.'

include{
    ConvertFilesKM3NeT;
    CreateEventListKM3NeT;
    CorrectEventListKM3NeT;
    InjectSignalKM3NeT;
    CombineEventListsKM3NeT;
    EpochFoldingKM3NeT;
} from '../km3net_processes.nf'

evaluate(new File("../inputs/km3net_inputs.nf"))

workflow{
    ConvertFilesKM3NeT(params.input_files);
    CreateEventListKM3NeT(ConvertFilesKM3NeT.out, params.source_file, input.dist, input.energy_threshold);
    CorrectEventListKM3NeT(CreateEventListKM3NeT.out, params.source_file);
    InjectSignalKM3NeT(CorrectEventListKM3NeT.out, input.frequency, input.ratio, input.pulseshape, input.df, input.baseline, input.a, input.phi, input.kappa);
    CombineEventListsKM3NeT(InjectSignalKM3NeT.out);
    EpochFoldingKM3NeT(CombineEventListsKM3NeT.out, input.frequency, input.number_of_testf, input.df, input.nbin);

}