process PlotHist{
  
  input:
    path input_file;
  output:
    path "*.png";

  script:
  """
  plot-hist ${input_file} -o ./
  """
   
}
