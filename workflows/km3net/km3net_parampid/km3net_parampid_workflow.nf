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
    .fromPath(input.train_files_arca_numu)
    .splitText(by: 1)
    .set {Numu_TrainFiles_Channel}

    Channel
    .fromPath(input.train_files_arca_anue)
    .splitText(by: 1)
    .mix(Numu_TrainFiles_Channel)
    .set {TrainFiles_Channel}
    
    
    ExtractHitFeaturesTraining(TrainFiles_Channel, input.detectorname)
    ConcatFiles(ExtractHitFeaturesTraining.out.converted_file.collect())
    TrainPID(ConcatFiles.out.concatenated_file,input.column_table_file, input.detectorname)

    Channel
    .fromPath(input.files_arca_numu)
    .splitText(by: 1)
    .set {Numu_Files_Channel}

    Channel
    .fromPath(input.files_arca_anue)
    .splitText(by: 1)
    .mix(Numu_Files_Channel)
    .set {Files_Channel}

    ExtractHitFeaturesClassification(Files_Channel, input.detectorname)
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