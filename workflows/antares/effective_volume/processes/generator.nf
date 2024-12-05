
process Km3buuSingleEnergy_cylinder{

 input:
    val n_events;
    val n_runs;
    val interaction;
    val flavor;
    val target_z;
    val target_a;
    val energy;
    val radius;

 output:
    path "km3buu*.root", emit: km3buu_root;
    path "xs*.dat", emit: km3buu_xs;

 publishDir "${params.output_dir}/generator", mode: 'link', overwrite: true;

 script:
 """
 #!/usr/bin/python3
 import glob
 import re
 from shutil import copyfile
 from tempfile import TemporaryDirectory
 from km3buu.jobcard import generate_neutrino_jobcard
 from km3buu.jobcard import XSECTIONMODE_LOOKUP, PROCESS_LOOKUP, FLAVOR_LOOKUP
 from km3buu.ctrl import run_jobcard
 from km3buu.geometry import CANVolume, SphericalVolume
 from km3buu.output import GiBUUOutput, write_detector_file

 #jc = generate_neutrino_jobcard(${n_events}, ${n_runs}, "${interaction}", "${flavor}", ${energy}, (${target_z},${target_a}),fluxfile=None)

 jc = generate_neutrino_jobcard(${n_events}, ${n_runs}, "${interaction}", "${flavor}", ${energy}, (${target_a},${target_z}),fluxfile=None)

 outdir = "."

 output_filename = "km3buu_${flavor}_${interaction}_${energy}_Z${target_z}A${target_a}.root"

 run_jobcard(jc, outdir)

 fobj = GiBUUOutput(outdir)

 volume = CANVolume(radius=118.15, zmin=0.0, zmax=498.0)

 write_detector_file(fobj, geometry=volume, ofile=output_filename)

 output_xs = "xs_${flavor}_${interaction}_${energy}_Z${target_z}A${target_a}_"

 regexp = re.compile('neutrino_absorption_cross_section_(.+).dat')
 for file in glob.glob("*.dat"):
     if (regexp.match(file)!=None):
         copyfile(file,output_xs+regexp.split(file)[1]+'.dat')
 """
}
