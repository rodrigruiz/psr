process Km3sim{

 input:
    path input_file;
    val detector_file;
    val parameters_file;
    val schema;

 output:
    path "km3sim*.root";

 script:
 """
 KM3BUU_FILE=${input_file}
 OUTPUT_FILE=\${KM3BUU_FILE/km3buu/km3sim}
 KM3Sim -s 1234 -o \${OUTPUT_FILE} -d ${detector_file} -p ${parameters_file} -i ${input_file}
 """
}
