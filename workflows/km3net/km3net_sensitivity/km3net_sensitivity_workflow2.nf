nextflow.enable.dsl = 2

include {
    ConvertFilesKM3NeT;
    CreateEventListKM3NeT;
    CorrectEventListKM3NeT;
    InjectSignalKM3NeT;
    CombineEventListsKM3NeT;
    EpochFoldingKM3NeT;
} from './processes/km3net_processes.nf'

evaluate(new File(params.input_file))

workflow {
    def ratiolist = [0.1, 0.2, 0.3, 0.4]  
    def num_iterations = 3  

    
    Channel
        .fromPath(input.km3net_files)
        .splitText(by: 1)
        .set { files }

    
    ratios = Channel.fromList(ratiolist)

    
    combinedChannel = ratios.combine(files).flatMap { ratio, file -> 
        (1..num_iterations).flatMap { [ratio, file] }
    }

    // Split the combined channel into separate channels for each element
    ratioChannel = combinedChannel.map { it[0] }
    fileChannel = combinedChannel.map { it[1] }
    iterationChannel = combinedChannel.map { it[2] }

    // View the combined channel for debugging purposes
    combinedChannel.view()

    // Pass only the file paths to ConvertFilesKM3NeT
    ConvertFilesKM3NeT(fileChannel)

    // Pass the output of ConvertFilesKM3NeT and the ratios to CreateEventListKM3NeT
    CreateEventListKM3NeT(ConvertFilesKM3NeT.out, input.source_file, input.dist, input.energy_threshold)

    // Pass the output of CreateEventListKM3NeT to CorrectEventListKM3NeT
    CorrectEventListKM3NeT(CreateEventListKM3NeT.out, input.source_file)

    // Pass the output of CorrectEventListKM3NeT, the ratios, and the iteration number to InjectSignalKM3NeT
    InjectSignalKM3NeT(CorrectEventListKM3NeT.out, input.frequency, ratioChannel, input.pulseshape, input.df, input.baseline, input.a, input.phi, input.kappa)

    // Collect the output of InjectSignalKM3NeT and pass it to CombineEventListsKM3NeT
    CombineEventListsKM3NeT(InjectSignalKM3NeT.out.collect())

    // Pass the output of CombineEventListsKM3NeT to EpochFoldingKM3NeT
    EpochFoldingKM3NeT(CombineEventListsKM3NeT.out, input.frequency, input.number_of_testf, input.testf_df, input.nbin, iterationChannel)
}
