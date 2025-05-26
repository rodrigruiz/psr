nextflow.enable.dsl = 2

include{
    ExtractHitFeatures as ExtractHitFeaturesTraining;
    ExtractHitFeatures as ExtractHitFeaturesClassification;
    ConcatFiles;
    TrainPID as TrainPIDTrack;
    TrainPID as TrainPIDMuon1;
    TrainPID as TrainPIDMuon2;
    ApplyClassifierPID as ApplyClassifierPIDTrack;
    ApplyClassifierPID as ApplyClassifierPIDMuon;
    FilterMuonScore as FilterMuonScore;
} from '../processes/km3net_processes.nf'

evaluate(new File(params.input_file))



workflow{
    

    Channel
    .fromPath(input.mc_files)
    .splitText(by: 1)
    .set {TrainFiles_Channel}
    
    
    ExtractHitFeaturesTraining(TrainFiles_Channel, input.detectorname, 'mc')
    ConcatFiles(ExtractHitFeaturesTraining.out.converted_file.collect())

    TrainPIDTrack(ConcatFiles.out.concatenated_file,input.column_table_file, input.detectorname,'track','track_score')

    //TrainPIDMuon1(ConcatFiles.out.concatenated_file,input.column_table_file, input.detectorname,'muon','muon_score1')
    //FilterMuonScore(TrainPIDMuon1.out.pid_output , 0.8)
    //TrainPIDMuon2(FilterMuonScore.out.filtered_file,input.column_table_file, input.detectorname,'muon','muon_score2')
    
    Channel
    .fromPath(input.data_files)
    .splitText(by: 1)
    .set {Files_Channel}

    ExtractHitFeaturesClassification(Files_Channel, input.detectorname, 'data')

    ApplyClassifierPIDTrack(ExtractHitFeaturesClassification.out.converted_file, TrainPIDTrack.out.rd_file, input.column_table_file, 'track_score')
    //ApplyClassifierPIDMuon(ApplyClassifierPIDTrack.out.output_file, TrainPIDMuon2.out.rd_file,  input.column_table_file, 'muon_score')


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