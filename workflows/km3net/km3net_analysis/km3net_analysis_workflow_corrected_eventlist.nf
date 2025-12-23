nextflow.enable.dsl = 2

include{
    ConvertFilesKM3NeT;
    BlindDataKM3NET;
    AddTrackScoreKM3NeT;
    CreateEventListKM3NeT_new as CreateSelectedEventList;
    CorrectEventListKM3NeT;
    CombineEventListsRuns as CombineEventListsRunsSelected;
    PlotSkyMapsKM3NeT as PlotSkyMapsKM3NeTSelected;
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
        error "Unsupported input.filestype: ${input.filestype}"
    }

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
        } else {
            eventlist_input = AddTrackScoreKM3NeT.out
        }
        BlindDataKM3NET(eventlist_input)
    }
    else {
        BlindDataKM3NET(ARCA_Files_Channel)

    }
    //AddTrackScoreKM3NeT(ARCA_Files_Channel, input.parampid_folder)

    // Add weights for MC only



    // Blind the data
    // BlindDataKM3NET(eventlist_input)
    
    // Create Event Lists
    //CreateEventListKM3NeT_new(BlindDataKM3NET.out, input.source_file, input.delta_search_min, input.ar_shower_file_arca, input.ar_track_file_arca, input.ar_shower_file_orca, input.ar_track_file_orca)
    CreateSelectedEventList(
    BlindDataKM3NET.out,
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

    CombineEventListsRunsSelected(CreateSelectedEventList.out.groupTuple(by: 1),input.source_file, input.delta_search_min, input.filestype)
    PlotSkyMapsKM3NeTSelected(CombineEventListsRunsSelected.out.combined_file, input.source_file, "combined_eventlists", input.delta_search_min, 120)
    CorrectEventListKM3NeT(CombineEventListsRunsSelected.out.combined_file, input.source_file) 


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