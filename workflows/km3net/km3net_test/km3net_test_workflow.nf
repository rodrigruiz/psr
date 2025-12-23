nextflow.enable.dsl = 2

include{
    TestProcess2;
} from '../processes/km3net_processes.nf'

evaluate(new File(params.input_file))


workflow{

    Channel
    .fromPath(input.km3net_arca_anue_files)
    .splitText(by: 1)
    .collate(10)
    //.buffer(size: 10, remainder: true)
    .set {ARCA_Files_Channel}
    //.buffer(size: 10, remainder: true)
    //.collate(10)


    Channel.of(1,2,3,4,5,6,7,8,9)
    .buffer(size: 3)
    .view()
    .set {NumberChannel}


    TestProcess2(NumberChannel)

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