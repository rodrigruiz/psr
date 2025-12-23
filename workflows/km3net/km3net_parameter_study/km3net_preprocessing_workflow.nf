nextflow.enable.dsl = 2

include {
    ConvertFilesKM3NeT;
    AddTrackScoreKM3NeT;
    AddWeightsKM3NeT;
    CreateEventListKM3NeT_new as CreateSelectedEventList;
    CreateEventListKM3NeT_new as CreateFullSkyEventList;
    BlindDataKM3NET;
    CombineEventListsRuns as CombineEventListsRunsFullSky;
    CombineEventListsRuns as CombineEventListsRunsSelected;
    PlotSkyMapsKM3NeT as PlotSkyMapsKM3NeTFullSky;
    PlotSkyMapsKM3NeT as PlotSkyMapsKM3NeTSelected;
    CorrectEventListKM3NeT;
    FindGTIsKM3NeTSingle;
    FindGTIsKM3NeTCombined;
} from '../processes/km3net_processes.nf'


evaluate(new File(params.input_file))
input.delta_search_min  = params.delta_search_min ?: input.delta_search_min
input.filestype         = params.filestype ?: input.filestype
input.source_file       = params.source_file ?: input.source_file
input.filestype         = params.filestype ?: input.filestype
input.plot_all          = params.plot_all ?: input.plot_all

workflow {

    Channel
        .fromPath(input.km3net_arca_data_files)
        .splitText(by: 1)
        .combine(Channel.of('arca'))
        .set { ARCA_Files_Channel }

    Channel
        .fromPath(input.km3net_orca_data_files)
        .splitText(by: 1)
        .combine(Channel.of('orca'))
        .mix(ARCA_Files_Channel)
        .set { Files_Channel }

    if (input.filestype == 'mc') {
        Channel
            .fromPath(input.km3net_arca_mc_files)
            .splitText(by: 1)
            .combine(Channel.of('arca'))
            .set { ARCA_Files_Channel }

        Channel
            .fromPath(input.km3net_orca_mc_files)
            .splitText(by: 1)
            .combine(Channel.of('orca'))
            .mix(ARCA_Files_Channel)
            .set { Files_Channel }
    }

    ConvertFilesKM3NeT(Files_Channel)
    AddTrackScoreKM3NeT(ConvertFilesKM3NeT.out, input.parampid_folder)

    def eventlist_input = (input.filestype == 'mc')
        ? AddWeightsKM3NeT(AddTrackScoreKM3NeT.out, input.orca_weights_file, input.arca_weights_file).out
        : AddTrackScoreKM3NeT.out

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

    BlindDataKM3NET(eventlist_input)

    CreateSelectedEventList(
        BlindDataKM3NET.out,
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

    CombineEventListsRunsFullSky(CreateFullSkyEventList.out.groupTuple(by: 1), input.source_file, input.delta_search_min, input.filestype)
    
    if (input.plot_all == 'True') {
        PlotSkyMapsKM3NeTFullSky(CombineEventListsRunsFullSky.out.combined_file, input.source_file, "combined_eventlists", input.delta_search_min, 120)
    }
    
    CombineEventListsRunsSelected(CreateSelectedEventList.out.groupTuple(by: 1), input.source_file, input.delta_search_min, input.filestype)
    PlotSkyMapsKM3NeTSelected(CombineEventListsRunsSelected.out.combined_file, input.source_file, "combined_eventlists", input.delta_search_min, 120)

    FindGTIsKM3NeTSingle(CombineEventListsRunsFullSky.out.combined_file, 600)

    CombineEventListsRunsFullSky.out.combined_file
        .map { it instanceof List || it instanceof Tuple ? it[0] : it }
        .collect()
        .set{ CombinedEventlistsOutput}
    CombinedEventlistsOutput.view()
    FindGTIsKM3NeTCombined(CombinedEventlistsOutput,600)
    

    CorrectEventListKM3NeT(CombineEventListsRunsSelected.out.combined_file, input.source_file)
}
