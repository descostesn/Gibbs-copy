use Getopt::Std;

$first = $ARGV[0];
getopt('pfuolramth');



### for help ###

if(($first eq "-h") || ($first eq "-help") || ($first eq "help"))                     
{
  printf "Help:\n";
  printf "-p: directory of the formatted initial mapping file and the file of sites found by unique tags alone\n";
  printf "-f: name of the initial mapping file\n";
  printf "-u: name of the file of sites found by unique tags alone\n";
  printf "-o: name of output file, default value: final_mapping.bed\n";
  printf "-l: length of the adjacent genomic region used to count co-located tags, default value: 147bp\n";
  printf "-r: the maximum tag count for calculating likelihood ratio, default value: 50\n";
  printf "-a: the relative confidence of ambiguous tags, default value: 0.2\n";
  printf "-m: the maximum number of iterations, default value: 5\n";
}



### required parameter warnings ###

elsif((length($opt_p) == 0) || (length($opt_f) == 0) || (length($opt_u) == 0))         
{
   if(length($opt_p) == 0)
   {
     printf "Warning: directory of the files need to be provided!\n";
     printf "Readme: -p\n";
   }
   if(length($opt_f) == 0)
   {
     printf "Warning: the initial mapping file need to be provided!\n";
     printf "Readme: -f\n";
   }
   if(length($opt_u) == 0)
   {
     printf "Warning: the file of sites found by unique tags need to be provided!\n";
     printf "Readme: -u\n";
   }
}



### core of the script ###

else                                                                                  
{

    ##################### read parameters ######################
    
    $path = $opt_p;
    $mapping_file = $opt_f;
    $unique_file = $opt_u;
    $output_name = (length($opt_o) > 0)? ($opt_o):("final_mapping.bed");
    $region_length = (length($opt_l) > 0)? ($opt_l):147;
    $max_tag_thresh = (length($opt_r) > 0)? ($opt_r):50;
    $ambiguous_conf = (length($opt_a) > 0)? ($opt_a):0.2;
    $max_iter_num = (length($opt_m) > 0)? ($opt_m):5;
    $threshold = (length($opt_t) > 0)? ($opt_t):4;
  
    printf "path: %s\nmapping_file: %s\nunique_file: %s\noutput_file_name: %s\nregion_length: %d\nmax_tag_thresh: %d\nambiguous_conf: %f\nmax_iteration: %d\n",$path,$mapping_file, $unique_file,$output_name, $region_length,$max_tag_thresh,$ambiguous_conf,$max_iter_num;
    ### show all the parameters used ###
    

    ##################### start Gibbs sampling #########################
    
    srand();
    
    $command_temp = "mkdir ".$path."/temp";
    system($command_temp);
    printf "directory \"temp\" built\n";

    gibbs_signal_score($mapping_file,$unique_file,$region_length,$max_tag_thresh,$path,$ambiguous_conf);
    printf "parameterized\n";  ### parameterize the likelihood ratios ###

    gibbs_initial_pos($mapping_file,$path);
    printf "location initialized\n"; ### initialize the mapping positions of tags ###


    for($iteration=1;$iteration<=$max_iter_num;$iteration++)    ######  iteration loop  ######
    {
        gibbs_split_pos($path);
        printf "%d\tsplited\n",$iteration;
  
        gibbs_sort_pos($path);
        printf "%d\tsorted\n",$iteration;
   
        gibbs_count_adjacent($region_length,$path);
        printf "%d\tcounted\n",$iteration;
   
        gibbs_count_score($ambiguous_conf,$max_tag_thresh,$path);
        printf "%d\tscored\n",$iteration;
  
        if($iteration == $max_iter_num)
        {
            gibbs_core($mapping_file,$path,1);
        }
        else
        {
            gibbs_core($mapping_file,$path,0);
        }
        
        printf "iteration %d finished\n",$iteration;
  
    }

    ###################### make output file #########################
    
    format_output($path,$mapping_file,$output_name);
    printf "finished\n";
    
    $command_final = "rm -r ".$path."/temp";
    system($command_final);

}#end of else




#####################################################################################################################################################################################




####################**************************************   functions  ************************************##########################
####################****************************************************************************************##########################
####################****************************************************************************************##########################
####################****************************************************************************************##########################
####################****************************************************************************************##########################
####################**************************************   functions  ************************************##########################




###############################################  gibbs_signal_score   ##############################################

sub gibbs_signal_score
{
   #my ($mapping_file,$unique_file,$region_length,$max_tag_thresh,$path,$ambiguous_conf) = @_;
   
   my $mapping_file_pre = $_[0];
   my $unique_file_pre = $_[1];
   my $region_length = $_[2];
   
   my $max_tag = $_[3];
   
   my $path = $_[4];
   
   my $ambiguous_conf = $_[5];
   
   
   my $mapping_file = $path."/".$mapping_file_pre;
   my $unique_file = $path."/".$unique_file_pre;
   
   my $out = $path."/temp/LR_score.bed";
   open(OUT,">$out");
   $prev = select(OUT);
   
   open(IN,"$unique_file");
   my $line = <IN>;
   
   my $total = 0;
   my $region = 0;
   my @temp = ();
   
   while(length($line) > 1)
   {
     chop($line);
     @temp = split(/\s+/,$line);
     $total = $total + $temp[$#temp];
     $region ++;
     $line = <IN>;
   }
   my $ave = $total/($region);
   
   close(IN);
   
   my $sd = 0;
   open(IN,"$unique_file");
   $line = <IN>;
   
   while(length($line) > 1)
   {
     @temp = split(/\s+/,$line);
     $sd = $sd + ($temp[$#temp] - $ave)**2;
     $line = <IN>;
   }
   $sd = sqrt($sd/($region-1));
   
   close(IN);
   
   
   open(IN2,"$mapping_file");
   $line = <IN2>;
   
   
   my $total_tag = 0;
   my $genome_length = 3107677273;
   
   while(length($line) > 1)
   {
      $total_tag ++;
      $line = <IN2>;
   }
   
   my $ave_2 = ($region_length*$total_tag)/$genome_length;
   
   $ave = $ave + $ave*$ambiguous_conf*(($total_tag-$total)/$total);
   
   close(IN2);
   
   my $LLR=0;
   
   for($i=1;$i<=$max_tag;$i++)
   {
      #$LLR = $i*(log($ave) - log($ave_2)) - ($ave - $ave_2);
      $LLR = -0.5*log(2*3.1416*$sd*$sd) - 0.5*($i-$ave)**2/($sd*$sd) - $i*log($ave_2) + $ave_2 +$i*log($i) - $i;
      printf "%f\n",exp($LLR);
   }
   
   
   close(OUT);
   select($prev);
   
}






###############################################  gibbs_initial_pos   ##############################################

sub gibbs_initial_pos
{
   my $mapping_file_pre = $_[0];
   my $path = $_[1];
   
   my $mapping_file = $path."/".$mapping_file_pre;
   
   srand();
   
   my $in_file = $mapping_file;
   my $out_file = $path."/temp/position.bed";
   
   open(IN,"$in_file");
   my $line = <IN>;
   
   
   open(OUT,">$out_file");
   $prev = select(OUT);
   
   my @temp = ();
   my @temp2 = ();
   my @chr_pos = ();
   
   while(length($line) > 1)
   {
     chop($line);
     @temp = split(/\t/,$line);
     @temp2 = split(/,/,$temp[1]);
     
     if($#temp2 == 0)
     {
        @chr_pos = split(/\>/,$temp2[0]);
        printf "%s\t%d\tu\n",$chr_pos[0],$chr_pos[1];
     }
     else
     {
        for($j=0;$j<=$#temp2;$j++)
        {
          @chr_pos = split(/\>/,$temp2[$j]);
          printf "%s\t%d\ta\n",$chr_pos[0],$chr_pos[1];
        }
     }
     
     $line = <IN>;
   }
   
   close(IN);
   
   close(OUT);
   select($prev);
   
}





###############################################  gibbs_split_pos   ##############################################

sub gibbs_split_pos
{
   
   my $path = $_[0];
   
   my $chrom_file = $path."/chromosome_index.txt";
   
   my @chr = ();
   my @chromosome = ();
   open(IN,"$chrom_file");
   @chromosome = <IN>;
   close(IN);
   for($i=0;$i<=$#chromosome;$i++)
   {
      chop($chromosome[$i]);
      $chr[$i+1] = $chromosome[$i];
   }
   my $chr_num = $#chromosome + 1;
   
   my %chr_hash = ();
   for($i=1;$i<=$chr_num;$i++)
   {
      $chr_hash{$chr[$i]} = $i;
   }

   
   my $in_file = $path."/temp/position.bed";
   open(IN,"$in_file");
   my $line = <IN>;
   
   
   my @pos_array = ();
   my @index_array = ();
   
   for($i=1;$i<=$chr_num;$i++)
   {
     $index_array[$i] = 0;
   }
   
   my @temp = ();
   my $chrom = 0;
   my $insert = "";
   
   while(length($line) > 1)
   {
      @temp = split(/\s+/,$line);
      
      $chrom = $chr_hash{$temp[0]};
      
      $insert = $temp[1]."\t".$temp[2];
      
      $pos_array[$chrom][$index_array[$chrom]] = $insert;
      $index_array[$chrom] ++;
      
      $line = <IN>;
   }
   
   
   for($i=1;$i<=$chr_num;$i++)
   {
      my $out_file = $path."/temp/split_".$chr[$i].".bed";
      open(OUT,">$out_file");
      $prev = select(OUT);
      
      for($k=0;$k<$index_array[$i];$k++)
      {
         printf "%s\n",$pos_array[$i][$k];
      }
      
      
      close(OUT);
      select($prev);
      
   }
   
   close(IN);
   
   @pos_array = ();
  
}





###############################################  gibbs_sort_pos   ##############################################

sub gibbs_sort_pos
{
   my $path = $_[0];
   
   my $chrom_file = $path."/chromosome_index.txt";
   
   my @chr = ();
   my @chromosome = ();
   
   open(IN,"$chrom_file");
   @chromosome = <IN>;
   close(IN);
   for($i=0;$i<=$#chromosome;$i++)
   {
      chop($chromosome[$i]);
      $chr[$i+1] = $chromosome[$i];
   }
   my $chr_num = $#chromosome + 1;
   
   
   for($i=1;$i<=$chr_num;$i++)
   {
      my $in = $path."/temp/split_".$chr[$i].".bed";
      my $out = $path."/temp/sort_".$chr[$i].".bed";
      
      open(IN,"$in");
      my @line = <IN>;
      close(IN);
      
      open(OUT,">$out");
      $prev = select(OUT);
      
      my %hash = ();
      my @ref = ();
      my @order = ();
      
      my @temp = ();
      
      for($j=0;$j<=$#line;$j++)
      {
        chop($line[$j]);
        @temp = split(/\s+/,$line[$j]);
        $ref[$j] = $temp[0];
        $hash{$temp[0]} = $temp[1];
      }
               
      @order = sort {$a <=> $b} @ref;
      
      for($j=0;$j<=$#order;$j++)
      {
        printf "%d\t%s\n",$order[$j],$hash{$order[$j]};
      }
      
      
      close(OUT);
      select($prev);
      
   }
   
}




###############################################  gibbs_count_adjacent   ##############################################

sub gibbs_count_adjacent
{
   
   my $region_length = $_[0];
   my $path = $_[1];
  
   my $chrom_file = $path."/chromosome_index.txt";
   
   my @chr = ();
   my @chromosome = ();
   
   open(IN,"$chrom_file");
   @chromosome = <IN>;
   close(IN);
   for($i=0;$i<=$#chromosome;$i++)
   {
      chop($chromosome[$i]);
      $chr[$i+1] = $chromosome[$i];
   }
   my $chr_num = $#chromosome + 1;
   
   
   for($i=1;$i<=$chr_num;$i++)
   {
      my $in = $path."/temp/sort_".$chr[$i].".bed";
      my $out = $path."/temp/count_".$chr[$i].".bed";
      
      open(IN,"$in");
      my @line = <IN>;
      close(IN);
      
      open(OUT,">$out");
      $prev = select(OUT);
      
         #for($j=0;$j<=$#line;$j++)
         #{
          # chop($line[$j]);
         #}
      
         for($j=0;$j<=$#line;$j++)
         {
            my $count_u = 0;
            my $count_a = 0;
            
            my @temp = ();
            
            @temp = split(/\s+/,$line[$j]);
            my $position = $temp[0];
            my $ua_self = $temp[1];
            
            my $min = 0;
            my $max = 0;
            
            if($j-20 > 0)
            {
              $min = $j-20;
            }
            else
            {
              $min = 0;
            }
            if($#line - $j > 20)
            {
              $max = $j + 20;
            }
            else
            {
              $max = $#line;
            }
         
            for($k=$min;$k<=$max;$k++)
            {
               @temp = split(/\s+/,$line[$k]);
               
               if(abs($temp[0] - $position) <= $region_length)
               {
                 if($temp[1] eq 'u')
                 {
                    $count_u ++;
                 }
                 else
                 {
                    $count_a ++;
                 }
               }
            }
            
            printf "%d\t%d\t%d\n",$position,$count_u,$count_a;
         
         }
      
      
      close(OUT);
      select($prev);
      
   }
     
}



###############################################  gibbs_count_score   ##############################################

sub gibbs_count_score
{
   my $ambiguous_conf = $_[0];
   my $max_tag = $_[1];
   my $path = $_[2];
   
   my $parameter_file = $path."/temp/LR_score.bed";
   open(IN,"$parameter_file");
   my @score = <IN>;
   close(IN);
   
   for($i=0;$i<=$#score;$i++)
   {
    chop($score[$i]);
   }
   
 
   my $chrom_file = $path."/chromosome_index.txt";
   
   my @chr = ();
   my @chromosome = ();
   
   open(IN,"$chrom_file");
   @chromosome = <IN>;
   close(IN);
   for($i=0;$i<=$#chromosome;$i++)
   {
      chop($chromosome[$i]);
      $chr[$i+1] = $chromosome[$i];
   }
   my $chr_num = $#chromosome + 1;
   
   
   
   
   for($i=1;$i<=$chr_num;$i++)
   {
       my $in = $path."/temp/count_".$chr[$i].".bed";
       my $out = $path."/temp/score_".$chr[$i].".bed";
       
       open(IN,"$in");
       my @line = <IN>;
       close(IN);
       
       open(OUT,">$out");
       $prev = select(OUT);
       
       my @temp = ();
       my $number = 0;
       
       for($j=0;$j<=$#line;$j++)
       {
          chop($line[$j]);
          @temp = split(/\s+/,$line[$j]);
          
          $number = int($temp[1] + $ambiguous_conf*$temp[2]);
          
          if($number <= $max_tag)
          {
             if($number >= 1)
             {
               $u_index = $number-1;
               printf "%d\t%f\n",$temp[0],$score[$u_index];
             }
             else
             {
               printf "%d\t%f\n",$temp[0],exp($temp[2]- 1-(1/$ambiguous_conf))*$score[0];
             }
            
             
          }
          else
          {
             printf "%d\t%f\n",$temp[0],$score[$#score];
          }
          
       }
       
       
       close(OUT);
       select($prev);
      
   }   
   
   
}




###############################################  gibbs_core   ##############################################

sub gibbs_core
{
    
      my $mapping_file_pre = $_[0];
      my $path = $_[1];
      my $strand = $_[2];
      
      my $mapping_file = $path."/".$mapping_file_pre;
      
      my $chrom_file = $path."/chromosome_index.txt";
      
      my @chr = ();
      my @chromosome = ();
      
      open(IN,"$chrom_file");
      @chromosome = <IN>;
      close(IN);
      for($i=0;$i<=$#chromosome;$i++)
      {
         chop($chromosome[$i]);
         $chr[$i+1] = $chromosome[$i];
      }
      my $chr_num = $#chromosome + 1;
      
      my %chr_hash = ();
      for($i=1;$i<=$#chr_num;$i++)
      {
         $chr_hash{$chr[$i]} = $i;
      }
      
      
      srand();
      
      
      
      
      my %pre_hash;
      my %post_hash;
      my @ref;
      my $k = 0;
      
      my @temp = ();
      my $index = "";
      my $count = 0;
      
      for($i=1;$i<=$chr_num;$i++)
      {
         my $hash_in = $path."/temp/score_".$chr[$i].".bed";
         open(IN,"$hash_in");
         my @line = <IN>;
         close(IN);
         
         for($j=0;$j<=$#line;$j++)
         {
            chop($line[$j]);
            @temp = split(/\s+/,$line[$j]);
            $index = $chr[$i].">".$temp[0];
            $count = $temp[1];
            $pre_hash{$index} = $count;
            #$post_hash{$index} = 0;
            
            $ref[$k] = $index;
            $k ++;
         }
         
      }
      
      my $total_position = $k;
      
      printf "hash constructed\n";
      #printf "test_hash:%f\n",$pre_hash{'chr1>555898'};
      #printf "test_hash:%f\n",$pre_hash{'chr1>556025'};
      
      
      
      my $map_file = $mapping_file;
      open(IN,"$map_file");
      my $record = <IN>;
      
      my $out = $path."/temp/position.bed";
      open(OUT,">$out");
      $prev = select(OUT);
      
      
      while(length($record) > 2)
      {
         chop($record);
         @temp = split(/\t/,$record);
         my @pos = split(/,/,$temp[1]);
         
         my $total = 0;
         my @tag;
         my @pure_pos;
         my @relative_tag;
         
         for($x=0;$x<=$#pos;$x++)
         {
           $tag[$x] = 0.000002;
         }
         
         my @tempo = ();
         my $pure_pos = "";
         my $max_tag_num = 0;
         my $adjacent = 0;
         my $adja_temp_pos = "";
         
         for($j=0;$j<=$#pos;$j++)
         {
           @tempo = split(/\>/,$pos[$j]);
           if($#tempo == 2)
           {
              $pure_pos[$j] = $tempo[0].">".$tempo[1];
           
              $max_tag_num = 0;
              
              no warnings 'uninitialized';
              
              for($adjacent = $tempo[1] - 10;$adjacent <= $tempo[1] + 10;$adjacent ++)
              {
                 $adja_temp_pos = $tempo[0].">".$adjacent;
                 if($pre_hash{$adja_temp_pos} > $max_tag_num)
                 {
                   $max_tag_num = $pre_hash{$adja_temp_pos};
                 }
              }
              
              $tag[$j] += $max_tag_num;
              #$tag[$j] += $pre_hash{$pure_pos[$j]};
              $total += $tag[$j];
           }
           else
           {
              $tag[$j] = 0;
           }
         }
         
         my $seed = 0;
         my $look = 0;
         my $y = 0;
         my @select = ();
         
         if($total > 0)
         {
            for($j=0;$j<=$#pos;$j++)
            {
              $relative_tag[$j] = $tag[$j]/$total;
              #$post_hash{$pure_pos[$j]} += $relative_tag;
            }
         
            $seed = rand(1);
            $look = $relative_tag[0];
            $y = 0;
         
            while(($seed > $look) && ($y <=$#pos))
            {
              $y ++;
              $look += $relative_tag[$y];
            }
         
            @select = split(/\>/,$pos[$y]);
            
            #printf "%f\t%f\t%f\t%d\t%d\t",$seed,$look,$relative_tag[$y],$y,$#pos;  
            
            if($strand == 0)
            {
              if($#pos == 0)
              {
                printf "%s\t%d\tu\n",$select[0],$select[1];
              }
              else
              {
                printf "%s\t%d\ta\n",$select[0],$select[1];
              }
            }
            else
            {
              if($#pos == 0)
              {
                printf "%s\t%d\t%s\tu\n",$select[0],$select[1],$select[2];
              }
              else
              {
                printf "%s\t%d\t%s\ta\n",$select[0],$select[1],$select[2];
              }
            }
            
        
         }
         
         $record = <IN>;
      }
      
      
      close(IN);
      
      close(OUT);
      select($prev);
      
 
 }






###############################################  format_output   ##############################################

sub format_output
{
    
    my $path = $_[0];
    my $mapping_file = $_[1];
    my $output_name = $_[2];
    
    my $in = $path."/".$mapping_file;
    open(IN,"$in");
    my $tag = <IN>;
    
    my $in2 = $path."/temp/position.bed";
    open(IN2,"$in2");
    my $line = <IN2>;
    
    my $out = $path."/".$output_name;
    open(OUT,">$out");
    $prev = select(OUT);
    
    my @temp = ();
    
    while(length($tag) > 1)
    {
       @temp = split(/\t/,$tag);
       printf "%s\t",$temp[0];
       
       chop($line);
       @temp = split(/\s+/,$line);
       printf "%s\t%d\t%s\t%s\n",$temp[0],$temp[1],$temp[2],$temp[3];
       
       $tag = <IN>;
       $line = <IN2>;
    }
    
    close(OUT);
    select($prev);
    close(IN);
    close(IN2);
    
    
}













































