process GetFiles{

 input:
    val input_file;
 output:
    path "*.hdf5";  
 script:
 """
 cp /home/hpc/capn/mppi104h/saturn/data/ANTARES/rates/hdf5/${input_file} ${input_file} 
 """

}

process LoadTimeSeries{

 input:
    path input_file;

 output:
    path "*combined.hdf5";

 script:
 """
 load-timeseries -i ${input_file} -o ./ --blind
 """
}

process FindGTIs{

  input:
    path input_file;

  output:
    path "*gtis*.pkl"

  script:
  """
  find-gtis ${input_file} -o ./ --filepattern 'Antares_(\\d*)_total_rates_combined.hdf5' --combine
  """
}

process SplitTimeSeries{
  
  input:
    path input_file;
    val bins_per_file;

  output:
    path "*total_rates_[!combined]*.hdf5";

  script:
  """
  split-timeseries -i ${input_file} -o ./ --bins_per_file ${bins_per_file}
  """
}

process CreateEventList{
  
  input:
    path input_file;
    val scaling;

  output:
    path "*eventlist*.hdf5";

  script:
  """
  create-eventlist -i ${input_file} -o ./ --filepattern 'Antares_(\\d*)_total_rates_(\\d*).hdf5' --scaling ${scaling}
  """
}

process CorrectEventList{

  input:
    path input_file;
    val correction;
    val rajd;
    val decjd;
    val Porb;
    val axsini;
    val e;
    val omega;
    val Tpi2;

  output:
    path "*corrected.hdf5";

  script:
  """
  correct-eventlist -i ${input_file} -o ./ --filepattern 'Antares_(\\d*)_eventlist_(\\d*).hdf5' --correction ${correction} --rajd ${rajd} --decjd ${decjd} --Porb ${Porb} --axsini ${axsini} --e ${e} --omega ${omega} --Tpi2 ${Tpi2}
  """
}

process CalculateChi2{

  input:
    path input_file;
   
  output:
    path "*chi2*.hdf5";

  script:
  """
  calculate-chi2s ${input_file} -o ./ --filepattern 'Antares_(\\d*)_eventlist_(\\d*)_corrected.hdf5'
  """

}

process InjectSignal{
  input:
    path input_file;
    val a;
    val pulseshape;
    val frequency;
  output:
    path "*signal*.hdf5";

  script:
  """
  inject-signal ${input_file} -o ./ --filepattern 'Antares_(\\d*)_eventlist_(\\d*)_corrected.hdf5' --a ${a} --pulseshape '${pulseshape}' --frequency ${frequency}
  """ 
}

FoldEventlist{

  input:
    path input_file;
    val frequency;
    val signal_strength;
  output:
    path "*folded*mvm.hdf5";

  script:
  """
  fold-eventlist ${input_file} -o ./ --filpattern 'Antares_(\\d*)_eventlist_(\\d*)_signal_${signal_strength}_mvm.hdf5' --frequency ${frequency} --signal_strength ${signal_strength}
  """
}

process EFAddedProfile{
  input:
    val frequency;
    val signal_strength;
  output:
    path "*added.hdf5";
 
  script:
  """
  ef-add-profile -i ./  -o ./ --filepattern 'Antares_*_folded_*_${signal_strength}' --signal_strength ${signal_strength} --frequency ${frequency}
  """
}
