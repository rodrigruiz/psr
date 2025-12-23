nextflow.enable.dsl = 2

include{
    ConvertFilesKM3NeT;
    AddTrackScoreKM3NeT;
    CreateEventListKM3NeT_new as CreateFullSkyEventList;
    CombineEventListsRuns as CombineEventListsRunsSelected;
    CombineEventListsRuns as CombineEventListsRunsFullSky;
    CombineEventListsDetectors;
    PlotSkyMapsKM3NeT as PlotSkyMapsKM3NeTFullSky;
    FindGTIsKM3NeTCombined;
    AddWeightsKM3NeT;
} from '../processes/km3net_processes.nf'

evaluate(new File(params.input_file))
input.delta_search_min = params.delta_search_min ?: input.delta_search_min
input.filestype        = params.filestype ?: input.filestype
input.source_file      = params.source_file ?: input.source_file
input.cone_all         = params.cone_all ?: input.cone_all
input.convert          = params.convert ?: input.convert

def generateList(start, end, step) {
    def list = []
    for (def i = start; i <= end; i += step) {
        list.add(i)
    }
    return list
}


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
        AddTrackScoreKM3NeT(ConvertFilesKM3NeT.out, input.parampid_folder)
        eventlist_input_mcv = AddTrackScoreKM3NeT.out
    } else if (input.filestype == 'mc_orca') {
        Channel
            .fromPath("/home/hpc/capn/capn107h/software/psr/workflows/km3net/mc_files_orcav9.txt")
            .splitText(by: 1)
            .combine(Channel.of('orca'))
            .set { ORCA_Files_Channel }

        ConvertFilesKM3NeT(ORCA_Files_Channel)
        AddTrackScoreKM3NeT(ConvertFilesKM3NeT.out, input.parampid_folder)
        eventlist_input_mc_orca = AddTrackScoreKM3NeT.out
    } else if (input.filestype == 'data_orca') {
        Channel
            .fromPath("/home/hpc/capn/capn107h/software/psr/workflows/km3net/data_files_orca_v9.2_backup.txt")
            .splitText(by: 1)
            .combine(Channel.of('orca'))
            .set { ORCA_Files_Channel }

        ConvertFilesKM3NeT(ORCA_Files_Channel)
        AddTrackScoreKM3NeT(ConvertFilesKM3NeT.out, input.parampid_folder)
        eventlist_input_data_orca = AddTrackScoreKM3NeT.out
    } else if (input.filestype == 'data_arca') {
        Channel
            .fromPath("/home/hpc/capn/capn107h/software/psr/workflows/km3net/data_files_arca_v9_backup.txt")
            .splitText(by: 1)
            .combine(Channel.of('arca'))
            .set { ARCA_Files_Channel }

        //ConvertFilesKM3NeT(ORCA_Files_Channel)
        AddTrackScoreKM3NeT(ARCA_Files_Channel, input.parampid_folder)
        eventlist_input_data_orca = AddTrackScoreKM3NeT.out
    } else {
        error "Unsupported input.filestype: ${input.filestype}"
    }

    if (input.filestype == 'data') {
        Channel
            .fromPath(input.km3net_orca_data_files)
            .splitText(by: 1)
            .combine(Channel.of('orca'))
            .mix(ARCA_Files_Channel)
            .set { Files_Channel }
    } else if (input.filestype == 'mc') {
        Channel
            .fromPath(input.km3net_orca_mc_files)
            .splitText(by: 1)
            .combine(Channel.of('orca'))
            .mix(ARCA_Files_Channel)
            .set { Files_Channel }
    } else {
        //error "Unsupported input.filestype: ${input.filestype}"
    }


    /*
    if (input.convert == 'True'){
        ConvertFilesKM3NeT(ARCA_Files_Channel) //.view()
        AddTrackScoreKM3NeT(ConvertFilesKM3NeT.out, input.parampid_folder)
    }
    else{
        AddTrackScoreKM3NeT(ARCA_Files_Channel, input.parampid_folder)
    }

    // Add weights for MC only
    if (input.filestype == 'mc') {
        AddWeightsKM3NeT(AddTrackScoreKM3NeT.out, input.orca_weights_file, input.arca_weights_file)
        eventlist_input = AddWeightsKM3NeT.out
    } else {
        eventlist_input = AddTrackScoreKM3NeT.out
    }
    */

    if (input.add_trackscore == 'True'){
        if (input.convert == 'True'){
            ConvertFilesKM3NeT(ARCA_Files_Channel) //.view()
            AddTrackScoreKM3NeT(ConvertFilesKM3NeT.out, input.parampid_folder)
        }
        else{
            AddTrackScoreKM3NeT(ARCA_Files_Channel, input.parampid_folder)
        }
        if (input.filestype == 'mc') {
            AddWeightsKM3NeT(AddTrackScoreKM3NeT.out, input.orca_weights_file, input.arca_weights_file)
            eventlist_input = AddWeightsKM3NeT.out
        } else if (input.filestype == 'data') {
            eventlist_input = AddTrackScoreKM3NeT.out
        }
    }
    else if (input.filestype == 'data' | input.filestype == 'mc') {
        eventlist_input = ARCA_Files_Channel

    }

    if (input.filestype == 'mcv') {
        eventlist_input = eventlist_input_mcv
    } else if (input.filestype == 'mc_orca') {
        eventlist_input = eventlist_input_mc_orca
    } else if (input.filestype == 'data_orca') {
        eventlist_input = eventlist_input_data_orca
    }

    CreateFullSkyEventList(
    eventlist_input,
    input.source_file,
    180.0,
    input.ar_shower_file_arca,
    input.ar_track_file_arca,
    input.ar_shower_file_orca,
    input.ar_track_file_orca,
    input.energy_low_arca,
    input.energy_high_arca,
    input.energy_low_orca,
    input.energy_high_orca,
    input.energy_threshold_arca,
    input.trackscore_threshold_arca,
    input.muonscore_threshold_arca,
    input.energy_threshold_orca,
    input.trackscore_threshold_orca,
    input.muonscore_threshold_orca,
    'True',
    )
        
    CombineEventListsRunsFullSky(CreateFullSkyEventList.out.groupTuple(by: 1),input.source_file, input.delta_search_min, input.filestype)
    if (input.plot_all == 'True') {
        PlotSkyMapsKM3NeTFullSky(CombineEventListsRunsFullSky.out.combined_file, input.source_file, "combined_eventlists", input.delta_search_min, 120)
    }

    CombineEventListsRunsFullSky.out.combined_file
        .map { it instanceof List || it instanceof Tuple ? it[0] : it }
        .collect()
        .set{ CombinedEventlistsOutput}
    // CombinedEventlistsOutput.view()
    FindGTIsKM3NeTCombined(CombinedEventlistsOutput,600)

    CreateSelectedEventList(
    eventlist_input,
    //AddTrackScoreKM3NeT.out,
    input.source_file,
    input.delta_search_min,
    input.ar_shower_file_arca,
    input.ar_track_file_arca,
    input.ar_shower_file_orca,
    input.ar_track_file_orca,
    input.energy_low_arca,
    input.energy_high_arca,
    input.energy_low_orca,
    input.energy_high_orca,
    input.energy_threshold_arca,
    input.trackscore_threshold_arca,
    input.muonscore_threshold_arca,
    input.energy_threshold_orca,
    input.trackscore_threshold_orca,
    input.muonscore_threshold_orca,
    input.cone_all,
    )
    


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