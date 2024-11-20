nextflow.enable.dsl = 2

include{
    ConvertFilesKM3NeT;
    AngularResolutionKM3NeT;
    MultifileAngularResolutionKM3NeT;
    CombineAngularResolutionKM3NeT;
} from './processes/km3net_processes.nf'

evaluate(new File(params.input_file))

workflow{

    println "Input file path: ${input.km3net_orca_files}"

    Channel
    .fromPath(input.km3net_orca_files)
    .splitText(by: 1)
    .combine(Channel.of('orca'))
    .set {Files_Channel}

    ConvertFilesKM3NeT(Files_Channel)

    def converted_file_paths = ConvertFilesKM3NeT.out.map { file, detectorname -> file }
    
    MultifileAngularResolutionKM3NeT(converted_file_paths.collect(),input.detector,input.runtype,input.recotype,input.pltscale,input.energy_low,input.energy_high,input.n_energy_bins)


}