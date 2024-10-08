nextflow.enable.dsl = 2

include{
    ConvertFilesKM3NeT;
    AngularResolutionKM3NeT;
    MultifileAngularResolutionKM3NeT;
    CombineAngularResolutionKM3NeT;
} from './processes/km3net_processes.nf'

evaluate(new File(params.input_file))

workflow{

    println "Input file path: ${input.km3net_arca_files}"

    Channel
    .fromPath(input.km3net_arca_files)
    .splitText(by: 1)
    .combine(Channel.of('arca'))
    .set {Files_Channel}

    ConvertFilesKM3NeT(Files_Channel)
    // AngularResolutionKM3NeT(ConvertFilesKM3NeT.out,input.detector,input.runtype,input.recotype,input.pltscale)
    // CombineAngularResolutionKM3NeT(AngularResolutionKM3NeT.out.hdf5.collect(),input.detector,input.runtype,input.recotype,input.pltscale)
    def converted_file_paths = ConvertFilesKM3NeT.out.map { file, detectorname -> file }
    
    MultifileAngularResolutionKM3NeT(converted_file_paths.collect(),input.detector,input.runtype,input.recotype,input.pltscale)


}