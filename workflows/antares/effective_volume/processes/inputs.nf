process jcreate_gdml_file {

 input:
  path  detector_file;
  env   GDML_SCHEMA_DIR;
  env   CAN_MARGIN_M;
  val   debug;

 output:
  path "*.gdml", emit: gdml_file;

 publishDir "${params.output_dir}/inputs", mode: 'link', overwrite: true;

 script:
  """
  basename="\$(basename ${detector_file})";

  output_file="\${basename%.*}.can_margin_\${CAN_MARGIN_M}m.gdml";

  version="\$(JPrintDetector -a ${detector_file} -O "version" -d 2 --)";

  JConvertDetectorFormat \\
    -a ${detector_file}  \\
    -o \${output_file}   \\
    -V \${version^^}     \\
    -d ${debug} --
  """
}

process jcreate_pmt_parameters_file {

 input:
  path detector_file;
  val  debug;

 output:
  path "*.txt", emit: pmt_parameters_file;

 publishDir "${params.output_dir}/inputs", mode: 'link', overwrite: true;

 script:
  """
  output_file="pmt_parameters.txt";

  JEditPMTParameters \\
    -a ${detector_file}         \\
    -M "-1 -1 set TTS_ns -38.0" \\
    -o \${output_file}          \\
    -d ${debug} --
  """
};
