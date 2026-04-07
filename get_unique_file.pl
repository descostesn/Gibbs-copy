use Getopt::Std;

getopt('piol');

if(length($opt_i) == 0)
{
     printf "Warning: input file not provided!\n";
}
else
{

$path = $opt_p;
if(length($opt_p) == 0)
{
   printf "The directory (-p) of the file is not specified\n";
   printf "The default directory is set as the current folder.\n";
   $path = "./";
}
$bowtie_file = $opt_i;
$screen_result_file = (length($opt_o)>0)?($opt_o):("unique_screen_result.bed");
$region_len = (length($opt_l)>0)?($opt_l):(147);

$command = "mkdir ".$path."/unique_profile";
system($command);

##### get unique tags
$in_file = $path."/".$bowtie_file;
$out_file = $path."/unique_profile/unique_tag.bed";

open(IN,"$in_file");
@line = <IN>;
close(IN);

open(OUT,">$out_file");
$prev = select(OUT);

for($i=0;$i<=$#line;$i++)
{
  chop($line[$i]);
  @temp = split(/\t/,$line[$i]);
  @temp2 = split(/,/,$temp[1]);
  
  if($#temp2 == 0)
  {
     printf "%s\t%s\n",$temp2[0],$temp[0];
  }
  
}
close(OUT);
select($prev);


############## split unique tag
@chr = ();
$chrom_file = $path."/chromosome_index.txt";
open(IN,"$chrom_file");
@chromosome = <IN>;
close(IN);
for($i=0;$i<=$#chromosome;$i++)
{
   chop($chromosome[$i]);
   $chr[$i+1] = $chromosome[$i];
}
$chr_num = $#chromosome + 1;

$in_file = $path."/unique_profile/unique_tag.bed";
open(IN,"$in_file");
@line = <IN>;
close(IN);

for($i=1;$i<=$chr_num;$i++)
{
   $out_file = $path."/unique_profile/split_unique_".$chr[$i].".bed";
   open(OUT,">$out_file");
   $prev = select(OUT);
   
   for($j=0;$j<=$#line;$j++)
   {
      @temp = split(/\s+/,$line[$j]);
      @temp2 = split(/\>/,$temp[0]);
      
      if($temp2[0] eq $chr[$i])
      {
         printf "%d\n",$temp2[1];
      }
      
   }
   
   
   close(OUT);
   select($prev);
}

###############sort unique tags
for($i=1;$i<=$chr_num;$i++)
{
   $in = $path."/unique_profile/split_unique_".$chr[$i].".bed";
   $out = $path."/unique_profile/sort_".$chr[$i].".bed";
   
   open(IN,"$in");
   @line = <IN>;
   close(IN);
   
   open(OUT,">$out");
   $prev = select(OUT);
   
   for($j=0;$j<=$#line;$j++)
   {
     chop($line[$j]);
   }
            
   @order = sort {$a <=> $b} @line;
   
   for($j=0;$j<=$#order;$j++)
   {
     printf "%d\n",$order[$j];
   }
   
   
   close(OUT);
   select($prev);
}

#################slide window
for($i=1;$i<=$chr_num;$i++)
{
           $in = $path."/unique_profile/sort_".$chr[$i].".bed";
           $out = $path."/unique_profile/window_tag_".$chr[$i].".bed";
           
           open(IN,"$in");
           @line = <IN>;
           close(IN);
           
           open(OUT,">$out");
           $prev = select(OUT);
           
           $low_bound = 0;
           $high_bound = 1;
           
           $window_start = $line[$low_bound]-$region_len+1;
           
             while($high_bound <= $#line)
             {
                while(($line[$high_bound] < ($line[$high_bound-1] + $region_len)) && ($high_bound <=$#line))
                {
                   $high_bound ++;
                }
                if($high_bound < $#line-1)
                {
                  $high_bound --;
                }
              
                for($window_start = $line[$low_bound]-$region_len+1;$window_start <= $line[$high_bound];$window_start=$window_start + 50)
                {
                   $tag = 0;
                   for($j=$low_bound;$j<=$high_bound;$j++)
                   {
                      if(($line[$j] >= $window_start) && ($line[$j]<=($window_start+146)))
                      {
                        $tag ++;
                      }
                   }
                   printf "%d\t%d\t%d\n",$window_start,$window_start+$region_len-1,$tag;
                }
                
                $low_bound = $high_bound + 1;
                $high_bound = $low_bound + 1;
              
             }
           
           
           
           
           close(OUT);
           select($prev);
           
}

############ delete_duplicate
for($i=1;$i<=$chr_num;$i++)
{
           $in = $path."/unique_profile/window_tag_".$chr[$i].".bed";
           $out = $path."/unique_profile/shrink_windowtag_".$chr[$i].".bed";
           
           open(IN,"$in");
           @line = <IN>;
           close(IN);
           
           open(OUT,">$out");
           $prev = select(OUT);
           
             for($j=0;$j<=$#line;$j++)
             {
               chop($line[$j]);
             }
           
             @temp = split(/\s+/,$line[0]);
             $max_tag = $temp[2];
             $max_window_index = 0;
           
             $index1 = 0;
             $index2 = 1;
           
             while($index2 <= $#line)
             {
               @temp1 = split(/\s+/,$line[$index1]);
               @temp2 = split(/\s+/,$line[$index2]);
             
               if($temp2[0] <= $temp1[1])
               {
                 if($temp2[2] >= $max_tag)
                 {
                   $max_window_index = $index2;
                   $max_tag = $temp2[2];
                 }
               }
             
               else
               {
                 printf "%s\n",$line[$max_window_index];
                 $max_tag = $temp2[2];
                 $max_window_index = $index2;
               }
               
               $index1 ++;
               $index2 ++;
             
             }
           
           
           close(OUT);
           select($prev);
           
}
      
#################### output

      $out = $path."/".$screen_result_file;
      open(OUT,">$out");
      $prev = select(OUT);
      
      printf "track name=\"unique_tag_site\" description=\"unique_tag_site\"\n";
      
      $thresh = 3;
      
      
      for($i=1;$i<=$chr_num;$i++)
      {
           $in = $path."/unique_profile/shrink_windowtag_".$chr[$i].".bed";
           
           open(IN,"$in");
           @line = <IN>;
           close(IN);
           
             for($j=0;$j<=$#line;$j++)
             {
               chop($line[$j]);
               @temp = split(/\s+/,$line[$j]);
               
               if($temp[2] > $thresh)
               {
                 printf "%s\t%s\n",$chr[$i],$line[$j];
               }
               
             }
           
      }
      
      
      close(OUT);
      select($prev);
      
      

}




$command = "rm -r ".$path."/unique_profile";
system($command);