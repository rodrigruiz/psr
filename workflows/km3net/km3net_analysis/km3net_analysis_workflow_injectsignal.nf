nextflow.enable.dsl = 2

include{
    InjectSignalKM3NeT;
} from '../processes/km3net_processes.nf'

evaluate(new File(params.input_file))
input.delta_search_min = params.delta_search_min ?: input.delta_search_min
input.filestype        = params.filestype ?: input.filestype
input.source_file      = params.source_file ?: input.source_file
input.cone_all         = params.cone_all ?: input.cone_all

def generateList(start, end, step) {
    def list = []
    for (def i = start; i <= end; i += step) {
        list.add(i)
    }
    return list
}


workflow{
    
    snr_list = generateList(input.ratio_min, input.ratio_max, input.ratio_step)
    rate_list = [ 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]
    // rate_list = [5.0, 10.0, 25.0, 50.0]
    // rate_list = [5.0]
    iteration_list = (1..input.repetitions).toList()
    //e_min_list = [5e4, 8e4, 1e5, 2e5, 5e5, 6e5, 1e6, 2e6]
    //gamma_list = [1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6]
    
    e_min_list = [1e2]
    // gamma_list = [1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0]
    gamma_list = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0] 
    //length_list = [5.0, 10.0, 15.0, 20.0, 25.0]
    //length_list = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0]
    length_list = [90.0]


    SNR_Channel = Channel.fromList(rate_list)
    Iteration_Channel = Channel.fromList(iteration_list)
    Emin_Channel = Channel.fromList(e_min_list)
    Gamma_Channel = Channel.fromList(gamma_list)
   
    CorrectedEventList_Channel = Channel.fromPath(input.corrected_eventlist)
    Combined_Channel = SNR_Channel.combine(CorrectedEventList_Channel).combine(Channel.of('arca')).combine(Iteration_Channel).combine(Emin_Channel).combine(Gamma_Channel)
    

    InjectSignalKM3NeT(Combined_Channel, input.frequency, input.pulseshape, input.df, input.baseline, input.a, input.phi, input.kappa, input.inject_method, input.arca_energy_response_file)
}
workflow.onComplete = {
  println "Pipeline complete"
  println "Command line: $workflow.commandLine"
  def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

    sendMail(to: 'hwarnhofer@km3net.de', subject: 'My pipeline execution', body: msg)
}

workflow.onError = {
  println "Oops... something went wrong"
}