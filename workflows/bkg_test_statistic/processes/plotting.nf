process PlotHist{
  
  input:
    path input_file;
  output:
    path output_file;

  script:
  """
  plot-hist ${input_file} -o ./
  """
   
}
