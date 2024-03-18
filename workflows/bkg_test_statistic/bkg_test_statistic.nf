
/*
 * author: rgruiz
 */

/*
 * run this workflow with the following command:
 * nextflow run km3buu.nf -c profiles.config -profile woody -with-dag flowchart.png -with-timeline -with-trace -with-report
 */

/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Import processes
 */

include {
  GetFiles;
  LoadTimeSeries;
  FindGTIs;
  SplitTimeSeries;
  CreateEventList;
  CorrectEventList;
  CalculateChi2;
} from './processes/analysis.nf'


/*
 * A way to include configuration parameters is to evaluate a separate script that includes the definition of the parameters as groovy variables.
 * In this file, the parameters are defined as [key,value] pairs in associative arrays. Their elements can be accessed here as <array_name>.key.
 * The array containing the input parameters is called "inputs".
 */
evaluate(new File("./inputs/inputs.nf"))

/*
 * Workflow definition.
 */
workflow {
    files = channel.of( "Antares_053140_total_rates.hdf5",
			"Antares_053150_total_rates.hdf5" );
    
    GetFiles(files);

    LoadTimeSeries(GetFiles.out);  

    SplitTimeSeries(LoadTimeSeries.out);

    CreateEventList(SplitTimeSeries.out.flatten());
    
    CorrectEventList(CreateEventList.out);

    FindGTIs(LoadTimeSeries.out);

    CalculateChi2(CorrectEventList.out,FindGTIs.out);
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

    sendMail(to: 'rgracia@km3net.de', subject: 'My pipeline execution', body: msg)
}

workflow.onError = {
  println "Oops... something went wrong"
}
