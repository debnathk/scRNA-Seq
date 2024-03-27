<?php

// Read mapping file
$hashmap_file = "miRNA_int_hashmap_alcohol_data.txt";
$file_handle = fopen($hashmap_file,'r');
$file_contents = explode("\n",file_get_contents($hashmap_file)); //read a file into a string broken into each line
$mirna_list = array();

// Loop to store all the miRNAs in an internal list
for ($i = 0; $i<count($file_contents); $i++)
{
  $eachline = explode("\t",$file_contents[$i]); //Break each line into 2 words
  $mirna_list[$i] = $eachline[1]; //Store the miRNA name into an array
}


// Open output file
$output_file= "output.txt";
$fp_output_file = fopen($output_file, 'w+') or die("Could not open output.txt file");

// Open miRNA-miRNA network file
$mirna_network = "consensus.txt";
$fp_mirna_network = fopen($mirna_network,'r');
$limit = 1; // cutting off at 0.9000 score

if($fp_mirna_network)
 {
  while(($edge = fgets($fp_mirna_network))!==false and $limit < 165000)
   {
     //Break the edge into strings 
      $edge_data = str_getcsv($edge,"\t");

      // find their respective integers 
      for($i=0; $i<count($mirna_list); $i++)
      {
       if($edge_data[0] == $mirna_list[$i]) // Loop though the miRNA list
         { $int_m1 = $i+1; }
      }
  

      for($i=0; $i<count($mirna_list); $i++)
      {
       if($edge_data[1] == $mirna_list[$i])
         { $int_m2 = $i+1; } 
      }
  
      $converted_edge = $int_m1." ".$int_m2." ".$edge_data[2]."\r\n"; // Make the converted edge
      fwrite($fp_output_file,$converted_edge); // Print it to output file
      $limit++;
   }
  fclose($fp_mirna_network);
 }
else echo "Could not open the miRNA network file";   

fclose($file_handle);
fclose($fp_output_file);

?>

