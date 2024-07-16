nextflow.enable.dsl = 2

params.source_file = '/home/hpc/capn/capn107h/software/hdf5SourceFiles/Vela_X-1.h5'

include{
    ConvertFilesKM3NeT;
    CreateEventListKM3NeT;
    CorrectEventListKM3NeT;
    InjectSignalKM3NeT;
    CombineEventListsKM3NeT;
    EpochFoldingKM3NeT;
} from '/home/hpc/capn/capn107h/software/psr/workflows/km3net/km3net_processes.nf'

evaluate(new File("/home/hpc/capn/capn107h/software/psr/workflows/km3net/inputs/km3net_inputs.nf"))

workflow{

    files = channel.of('/home/hpc/capn/mppi104h/wecapstor3/out/ARCA/KM3NeT_00000133/v8.1/reco/mcv8.1.gsg_anue-CCHEDIS_1e2-1e8GeV.sirene.jterbr00013288.jchain.aashower.1.root',
    '/home/hpc/capn/mppi104h/wecapstor3/out/ARCA/KM3NeT_00000133/v8.1/reco/mcv8.1.gsg_anue-CCHEDIS_1e2-1e8GeV.sirene.jterbr00013289.jchain.aashower.2.root',
    '/home/hpc/capn/mppi104h/wecapstor3/out/ARCA/KM3NeT_00000133/v8.1/reco/mcv8.1.gsg_anue-CCHEDIS_1e2-1e8GeV.sirene.jterbr00013290.jchain.aashower.3.root',
    '/home/hpc/capn/mppi104h/wecapstor3/out/ARCA/KM3NeT_00000133/v8.1/reco/mcv8.1.gsg_anue-CCHEDIS_1e2-1e8GeV.sirene.jterbr00013291.jchain.aashower.4.root',
    '/home/hpc/capn/mppi104h/wecapstor3/out/ARCA/KM3NeT_00000133/v8.1/reco/mcv8.1.gsg_anue-CCHEDIS_1e2-1e8GeV.sirene.jterbr00013292.jchain.aashower.5.root',
    '/home/hpc/capn/mppi104h/wecapstor3/out/ARCA/KM3NeT_00000133/v8.1/reco/mcv8.1.gsg_anue-CCHEDIS_1e2-1e8GeV.sirene.jterbr00013293.jchain.aashower.6.root',
    '/home/hpc/capn/mppi104h/wecapstor3/out/ARCA/KM3NeT_00000133/v8.1/reco/mcv8.1.gsg_anue-CCHEDIS_1e2-1e8GeV.sirene.jterbr00013294.jchain.aashower.7.root',
    '/home/hpc/capn/mppi104h/wecapstor3/out/ARCA/KM3NeT_00000133/v8.1/reco/mcv8.1.gsg_anue-CCHEDIS_1e2-1e8GeV.sirene.jterbr00013295.jchain.aashower.8.root',
    '/home/hpc/capn/mppi104h/wecapstor3/out/ARCA/KM3NeT_00000133/v8.1/reco/mcv8.1.gsg_anue-CCHEDIS_1e2-1e8GeV.sirene.jterbr00013296.jchain.aashower.9.root',
    '/home/hpc/capn/mppi104h/wecapstor3/out/ARCA/KM3NeT_00000133/v8.1/reco/mcv8.1.gsg_anue-CCHEDIS_1e2-1e8GeV.sirene.jterbr00013297.jchain.aashower.10.root',
    '/home/hpc/capn/mppi104h/wecapstor3/out/ARCA/KM3NeT_00000133/v8.1/reco/mcv8.1.gsg_anue-CCHEDIS_1e2-1e8GeV.sirene.jterbr00013298.jchain.aashower.11.root');
    
    ConvertFilesKM3NeT(files);
    CreateEventListKM3NeT(ConvertFilesKM3NeT.out, input.source_file, input.dist, input.energy_threshold);
    CorrectEventListKM3NeT(CreateEventListKM3NeT.out, input.source_file);
    InjectSignalKM3NeT(CorrectEventListKM3NeT.out, input.frequency, input.ratio, input.pulseshape, input.df, input.baseline, input.a, input.phi, input.kappa);
    CombineEventListsKM3NeT(InjectSignalKM3NeT.out.collect());
    EpochFoldingKM3NeT(CombineEventListsKM3NeT.out, input.frequency, input.number_of_testf, input.testf_df, input.nbin);

}