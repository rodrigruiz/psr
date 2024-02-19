
process JaanetPreprocessor{

 input:
    path input_file;

 output:
    path "jaanet_preprocessor*.root";

 script:
 """
 KM3SIM_FILE=${input_file}
 OUTPUT_FILE=\${KM3SIM_FILE/km3sim/jaanet_preprocessor}
 JAAnetPreprocessor -f ${input_file} -o \${OUTPUT_FILE} -H
 """
}

process EffectiveVolume{

 input:
  val  energy;
  val  radius;
  path unfiltered_file;
  path filtered_file;
  val target_a;
  val target_z;
  val flavor;
  val interaction;

 output:
  path "*.png", emit: plt;
  path "*.txt", emit: txt;

 script:
 """
#!/usr/bin/python3
import km3io
import sys
import re
import numpy as np
import matplotlib.pyplot as plt


def atof(text):
    try:
        retval=float(text)
    except ValueError:
        retval=text
    return retval

def natural_keys(text):
    return [atof(c) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', text)]

output_ascii = "effective_volume_${flavor}_${interaction}_Z${target_z}A${target_a}.txt"
output_plt = "effective_volume_${flavor}_${interaction}_Z${target_z}A${target_a}.png"

energies = sorted(${energy})
jaanet_files = sorted("${filtered_file}".split(), key=natural_keys)
km3sim_files = sorted("${unfiltered_file}".split(), key=natural_keys)

generated_events=[]
detected_events =[]
n = len(energies)

if (len(km3sim_files) != n or len(jaanet_files) != n):
    sys.exit("Not all the lists have the same length")

regexp = re.compile('(\\w+)_(.+)_(\\w+)(.root)')

for i in range(n):
    if (regexp.match(km3sim_files[i])!=None and regexp.match(jaanet_files[i])!=None):
        e_gen = regexp.split(km3sim_files[i])[2]
        e_det = regexp.split(jaanet_files[i])[2]

        if (float(e_gen) == energies[i] and float(e_det) == energies[i]):
            generated_events.append(len(km3io.OfflineReader(km3sim_files[i]).events))
            detected_events .append(len(km3io.OfflineReader(jaanet_files[i]).events))
        else:
            print('Energies arenot the same'+str(float(e_gen))+'!='+str(float(e_gen))+'!='+str(energies[i]))
    else:
        sys.exit("File names don't match regexp")

v = 4*np.pi*np.power(${radius},3)/3
p_hat = [i/j for i,j in zip(detected_events,generated_events)]
err = [v*np.sqrt(p*(1-p)/n_gen) for p,n_gen in zip(p_hat,generated_events)]
eff_vol = [v*i/j for i,j in zip(detected_events,generated_events)]
print(len(energies))
print(len(eff_vol))
print(len(err))
np.savetxt(output_ascii, np.column_stack([energies,eff_vol,err]),delimiter=';')

fig, ax = plt.subplots( nrows=1, ncols=1 )
ax.errorbar(energies, eff_vol, yerr=err, color='blue')
ax.set_ylabel('Effective volume [m^3]')
ax.set_xlabel('Energy [GeV]')
ax.set_yscale('log')
ax.grid(linestyle='--', which='both')

fig.savefig(output_plt)
plt.close(fig)

 """
}
