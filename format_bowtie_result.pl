use Getopt::Std;
getopt('pio');

if(length($opt_i) == 0)
{
  
     printf "Warning: input file not provided!\n";
 
}

else
{

  $path = $opt_p;
  if(length($opt_p) == 0)
  {
       printf "The directory (-p) of the file is not specified.\n";
       printf "The default directory is set as the current folder.\n";
       $path = "./";
  }
  $bowtie_file = $opt_i;
  $change_file = (length($opt_o)>0)?($opt_o):("mapping_result.bed");

  $in_file = $path."/".$bowtie_file;
  $out_file = $path."/".$change_file;
  
  $out_chrom_file = $path."/chromosome_index.txt";
  
  @chrom = ();
  $pointer = 0;

  open(IN,"$in_file");
  $line = <IN>;

  open(OUT,">$out_file");
  $prev = select(OUT);

    @temp = split(/\t/,$line);
    printf "%s\t%s>%d>%s,",$temp[0],$temp[2],$temp[3]+1,$temp[1];
    $tag_name = $temp[0];
    $first_pos = $temp[3]+1;
    $chr = $temp[2];

    $chrom[0] = $chr;
    $pointer ++;

    $line = <IN>;
    while(length($line) > 2)
    {
        @temp = split(/\t/,$line);
        
           $flag = 0;
           for($k=0;$k<$pointer;$k++)
           {
              if($chrom[$k] eq $temp[2])
              {
                 $flag = 1;
              }
           }
           if($flag == 0)
           {
              $chrom[$pointer] = $temp[2];
              $pointer ++;
           }
        

        if($temp[0] eq $tag_name)
        {
           if((!(($temp[3]+1) == $first_pos)) || ($temp[2] ne $chr))
           {
            printf "%s>%d>%s,",$temp[2],$temp[3]+1,$temp[1];
           }
        }
        else
        {
           printf "\n";
           printf "%s\t%s>%d>%s,",$temp[0],$temp[2],$temp[3]+1,$temp[1];
           $tag_name = $temp[0];
           $first_pos = $temp[3] + 1;
           $chr = $temp[2];
        }
   
   
        $line = <IN>;
    }

  close(OUT);
  select($prev);
  close(IN);
  
  
  open(OUT,">$out_chrom_file");
  $prev = select(OUT);
  
  for($x=0;$x<$pointer;$x++)
  {
     printf "%s\n",$chrom[$x];
  }
  
  close(OUT);
  select($prev);
  
  
  

}