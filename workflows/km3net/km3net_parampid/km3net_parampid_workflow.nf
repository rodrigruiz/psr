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
    .set {ARCA_numu_TrainFiles_Channel}

    Channel
    .fromPath(input.train_files_arca_anue)
    .splitText(by: 1)
    .mix(ARCA_numu_TrainFiles_Channel)
    .set {ARCA_TrainFiles_Channel}
    
    
    ExtractHitFeaturesTraining(ARCA_TrainFiles_Channel,'ARCA')
    ConcatFiles(ExtractHitFeaturesTraining.out.converted_file.collect())
    TrainPID(ConcatFiles.out.concatenated_file,input.column_table_file,'ARCA')

    Channel
    .fromPath(input.files_arca_numu)
    .splitText(by: 1)
    .set {ARCA_numu_Files_Channel}

    Channel
    .fromPath(input.files_arca_anue)
    .splitText(by: 1)
    .mix(ARCA_numu_Files_Channel)
    .set {ARCA_Files_Channel}

    ExtractHitFeaturesClassification(ARCA_Files_Channel, 'ARCA')
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