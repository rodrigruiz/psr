nextflow.enable.dsl = 2

include{
    CalculateWeightsKM3NeT
} from '../processes/km3net_processes.nf'

evaluate(new File(params.input_file))

workflow{

    println "Input file path: ${input.km3net_arca_files}"

    Channel
    .fromPath(input.km3net_arca_files)
    .splitText(by: 1)
    .combine(Channel.of('arca'))
    .set {Files_Channel}

    CalculateWeightsKM3NeT(Files_Channel.groupTuple(by: 1))
}