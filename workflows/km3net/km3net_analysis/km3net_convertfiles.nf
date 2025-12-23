nextflow.enable.dsl = 2

include{
    ConvertFilesKM3NeT;
} from '../processes/km3net_processes.nf'

evaluate(new File(params.input_file))
input.filestype        = params.filestype ?: input.filestype
workflow{

    if (input.filestype == 'data') {
        Channel
            .fromPath(input.km3net_arca_data_files)
            .splitText(by: 1)
            .combine(Channel.of('arca'))
            .set { ARCA_Files_Channel }
    } else if (input.filestype == 'mc') {
        Channel
            .fromPath(input.km3net_arca_mc_files)
            .splitText(by: 1)
            .combine(Channel.of('arca'))
            .set { ARCA_Files_Channel }
            
    } else if (input.filestype == 'mcv') {
        Channel
            .fromPath("/home/hpc/capn/capn107h/software/psr/workflows/km3net/mc_files_arca.txt")
            .splitText(by: 1)
            .combine(Channel.of('arca'))
            .set { ARCA_Files_Channel }

        ConvertFilesKM3NeT(ARCA_Files_Channel)
    } else if (input.filestype == 'mc_orca') {
        Channel
            .fromPath("/home/hpc/capn/capn107h/software/psr/workflows/km3net/mc_files_orcav9.txt")
            .splitText(by: 1)
            .combine(Channel.of('orca'))
            .set { ORCA_Files_Channel }

        ConvertFilesKM3NeT(ORCA_Files_Channel)
    } else if (input.filestype == 'data_orca') {
        Channel
            .fromPath("/home/hpc/capn/capn107h/software/psr/workflows/km3net/data_files_orca_v9.2_backup2.txt")
            .splitText(by: 1)
            .combine(Channel.of('orca'))
            .set { ORCA_Files_Channel }

        ConvertFilesKM3NeT(ORCA_Files_Channel)
    } else if (input.filestype == 'data_arca') {
        Channel
            .fromPath("/home/hpc/capn/capn107h/software/psr/workflows/km3net/data_files_arca_v9_backup2.txt")
            .splitText(by: 1)
            .combine(Channel.of('arca'))
            .set { ORCA_Files_Channel }

        ConvertFilesKM3NeT(ORCA_Files_Channel)
    } else {
        error "Unsupported input.filestype: ${input.filestype}"
    }


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