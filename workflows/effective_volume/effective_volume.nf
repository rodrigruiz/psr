
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
  Km3buuSingleEnergy_cylinder;
} from './processes/km3buu.nf'

include {
  Km3sim;
} from './processes/light.nf'

include {
  JaanetPreprocessor;
  EffectiveVolume;
} from './processes/analysis.nf'

include {
  CrossSection;
} from './processes/plotting.nf'

/*
 * A way to include configuration parameters is to evaluate a separate script that includes the definition of the parameters as groovy variables.
 * In this file, the parameters are defined as [key,value] pairs in associative arrays. Their elements can be accessed here as <array_name>.key.
 * The array containing the input parameters is called "inputs".
 */
evaluate(new File("./parameters/inputs.nf"))

/*
 * Workflow definition.
 */
workflow {
    energies = channel.of(5.0, 10.0, 15.0);

    Km3buuSingleEnergy_cylinder(input.n_events,
		                        input.n_runs,
			                    input.interaction,
			                    input.flavor ,
			                    input.target_z,
			                    input.target_a,
			                    energies,
			                    input.radius);

  // Km3sim(Km3buuSingleEnergy.out.km3buu_root,
  // 	 input.detector,
  // 	 input.pmt_parameters,
  // 	 input.schema);

  // JaanetPreprocessor(Km3sim.out);

  // EffectiveVolume(energies.collect(),
  // 		  input.radius,
  // 		  Km3sim.out.collect(),
  // 		  JaanetPreprocessor.out.collect(),
  // 		  input.target_z,
  // 		  input.target_a,
  // 		  input.flavor,
  // 		  input.interaction);

  // CrossSection(Km3buuSingleEnergy.out.km3buu_xs.collect());
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
