nextflow.enable.dsl = 2

include{
    ExtractHitFeatures as ExtractHitFeaturesTraining;
    ExtractHitFeatures as ExtractHitFeaturesClassification;
    ConcatFiles;
    TrainPID;
    ApplyClassifierPID;
} from '../processes/km3net_processes.nf'

evaluate(new File(params.input_file))



workflow{
    

    Channel
    .fromPath(input.mc_files_arca)
    .splitText(by: 1)
    .set {TrainFiles_Channel}
    
    
    ExtractHitFeaturesTraining(TrainFiles_Channel, input.detectorname, 'mc')
    ConcatFiles(ExtractHitFeaturesTraining.out.converted_file.collect())
    TrainPID(ConcatFiles.out.concatenated_file,input.column_table_file, input.detectorname)

    Channel
    .fromPath(input.data_files_arca)
    .splitText(by: 1)
    .set {Files_Channel}

    ExtractHitFeaturesClassification(Files_Channel, input.detectorname, 'data')
    ApplyClassifierPID(ExtractHitFeaturesClassification.out.converted_file, TrainPID.out.rd_file, input.column_table_file)



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