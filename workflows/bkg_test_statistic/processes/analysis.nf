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
 load-timeseries -i ${input_file} -o ./
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

  output:
    path "*total_rates_[!combined]*.hdf5";

  script:
  """
  split-timeseries -i ${input_file} -o ./ --bins_per_file 200
  """
}

process CreateEventList{
  
  input:
    path input_file;

  output:
    path "*eventlist*.hdf5";

  script:
  """
  create-eventlist -i ${input_file} -o ./ --filepattern 'Antares_(\\d*)_total_rates_(\\d*).hdf5' --scaling 1e10
  """
}

process CorrectEventList{

  input:
    path input_file;

  output:
    path "*corrected.hdf5";

  script:
  """
  correct-eventlist -i ${input_file} -o ./ --filepattern 'Antares_(\\d*)_eventlist_(\\d*).hdf5' --correction bary+bin --rajd 02h43m40.4252869512s --decjd +61d26m03.757456824s --Porb 27.6943 --axsini 115.531 --e 0.1029 --omega -74.05 --Tpi2 2458116.09700
  """
}

process CalculateChi2{

  input:
    path input_file;
    path gti_file;
  output:
    path "*chi2*.hdf5";

  script:
  """
  calculate-chi2s ${input_file} -o ./ --filepattern 'Antares_(\\d*)_eventlist_(\\d*)_corrected.hdf5' --expocorr --gti_files '${gti_file}'
  """

}
