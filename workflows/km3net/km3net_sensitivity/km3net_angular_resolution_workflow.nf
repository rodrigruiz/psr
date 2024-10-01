nextflow.enable.dsl = 2

include{
    ConvertFilesKM3NeT;
    AngularResolutionKM3NeT;
    CombineAngularResolutionKM3NeT;
} from './processes/km3net_processes.nf'

evaluate(new File(params.input_file))

workflow{

    println "Input file path: ${input.km3net_arca_files}"

    Channel.fromPath('/home/hpc/capn/capn107h/software/psr/workflows/km3net/km3net_sensitivity/inputs/runs_test.txt').view()
    //.fromPath(input.km3net_arca_files)
    .splitText(by: 1)
    .set {ARCAFiles_Channel}
    .view()
    .combine(Channel.of('arca'))

    ConvertFilesKM3NeT(ARCAFiles_Channel)
    AngularResolutionKM3NeT(ConvertFilesKM3NeT.out,input.detector,input.runtype,input.recotype,input.pltscale)
    CombineAngularResolutionKM3NeT(AngularResolutionKM3NeT.out.hdf5.collect(),input.detector,input.runtype,input.recotype,input.pltscale)
}