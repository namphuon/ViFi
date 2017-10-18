use threads;
use threads::shared;
use phylo;
use File::Spec::Functions qw(rel2abs);
use File::stat;
use File::Basename;
use strict;
use File::Path;
use Data::Dumper;
use LWP::Simple;    #use Math::Random::Zipf;
use phylo;
use Getopt::Long qw(:config no_ignore_case);

my $dir = "";

my $Usage = "Usage: $0 -p prefix -b bamfile -t threads -d directory -H hmm_list\n";
my %hmms = ();
my %opt = (bamfile => "",
  thread => 24,
  thread_per_job => 1,  
  dir => "",
  hmm_list => "",
  prefix => "test" );
  GetOptions(\%opt, "bamfile|b=s", "dir|d:s", "thread|t:i", 'prefix|p=s', 'thread_per_job|j:i', 'hmm_list|H=s','debug:s');

run_pipeline();

sub run_pipeline {
  read_hmm_file();
  if (not -e "$opt{dir}/logs") {
    `mkdir -p $opt{dir}/logs`;
    `mkdir -p $opt{dir}/temp`;    
  }
  
  if (not -e "$opt{dir}/temp/unmapped.fas" or -s "$opt{dir}/temp/unmapped.fas" == 0) {
    prepare_unmapped_sequences();
  }  
  
  #Search for any viral reads on a subset of the data
  job_pool(build_hmm_jobs(),new Worker({'function',\&hmm_worker,'name','HmmSearch','threads_per_job',$opt{thread_per_job}}));
  if (not -e "$opt{dir}/temp/reduced.csv") {
    print "Extracting Fas\n";  
    process_results();  
  }  
  if (not -e "$opt{dir}/temp/$opt{prefix}_R2.fastq.gz") {
    print "Extracting FastQ\n";
    extract_reads_unmapped();
  }
  if (not -e "$opt{dir}/temp/viral.sam") {
    map_all_viral_reads();
  }
  
  return;
  
  if (not -e "$opt{dir}/temp/$opt{prefix}_R2.fastq.gz") {
    print "Extracting FastQ\n";
    extract_reads_fastq_gz();
  }
  find_strain();    
    
  if ($opt{debug} eq "true") {
  }
  return;  
  my $scores = read_results();
  
  
  if (not -e "$opt{dir}/temp/reduced.fasta") {
    process_results($scores);  
  }

  if (not -e "$opt{dir}/temp/$opt{prefix}_R2.fastq.gz") {
    extract_reads_fastq_gz();
  }
  return;
}             

sub map_all_viral_reads {
  if (not -e "$opt{dir}/temp/$opt{prefix}_R2.fastq.gz") {
    return;
  }
  my $start_time = time();
  print "Mapping viral to the human genome: $start_time\n";
  `bwa mem -t $opt{thread} $opt{human_reference} $opt{dir}/temp/$opt{prefix}_R1.fastq.gz $opt{dir}/temp/$opt{prefix}_R2.fastq.gz > $opt{dir}/temp/human.viral.sam`;
  `bwa mem -t $opt{thread} $opt{viral_reference} $opt{dir}/temp/$opt{prefix}_R1.fastq.gz $opt{dir}/temp/$opt{prefix}_R2.fastq.gz > $opt{dir}/temp/viral.sam`;  
  my $end_time = time()-$start_time;    
  print "Finished mapping viral reads to the references genome: $end_time\n";  
  `echo "$end_time" > $opt{dir}/logs/bwa.human.time`;
}

#Computes the number of randomly selected reads necessary in order to reach the 
#probability of not having any integrated viral reads
sub compute_reads {
  return log($opt{min_percent})/(log(1-($opt{virus_base}/($opt{virus_base}+$opt{human_base}))));
}

sub read_hmm_file {
  open(INPUT, "$opt{hmm_list}");
  my $idx = 0;
  while (my $line = Phylo::trim(<INPUT>."")) {
    $hmms{$idx} = $line;
    $idx++;    
  }
  close(INPUT);
}

sub fix_reads {
  my $file = $_[0];
  if (-e "$file.bak") {
    return;
  }
  open(OUTPUT, ">$file.fixed");
  open(INPUT, "$file");
  while (my $seq = <INPUT>) {
    my $line = <INPUT>;
    #print "$seq\n";
    if ($line =~ m/([^\s+ATCGN])/) {
      #print "delete $seq\n";
    } else {
      print OUTPUT "$seq$line";
    }
  }
  close(INPUT);
  close(OUTPUT);
  `mv $file $file.bak`;
  `mv $file.fixed $file`;
}

sub run_hmms {
  read_hmm_file();
  foreach my $key (sort keys %hmms) {
    if (not -e "$opt{dir}/temp/$key.hmm_search") {
      my $temp_file = Phylo::get_temp_file();    
      `/usr/bin/time -v -o $opt{dir}/logs/hmmsearch.$key.time nhmmer -o $temp_file --cpu $opt{thread} $hmms{$key} $opt{dir}/$opt{prefix}.fasta >& $opt{dir}/logs/hmmsearch.$key.log`;    
      `mv $temp_file $opt{dir}/temp/$key.hmm_search`;
    }
  }
}

sub convert_search_into_scores {
  my $input_file = $_[0];
  my $output_file = $_[1];
  my $prefix = $_[2];
  
  open(INPUT, "$input_file");
  my %results = ();
  while (my $line = <INPUT>) {      
    while ($line !~ m/Query:/ and eof INPUT != 1) {
      $line = Phylo::trim("".<INPUT>."");
      if (not defined $line) {
        last;
      }
    }
    $line =~ m/:\s+([^\s]+)/;
    my $hmm = $1;
    my $start_match = 0;
    my $quit = 0;
    while (defined $line and $quit == 0) {
      #print $line;
      $line = Phylo::trim($line);
      if ($line =~ m/inclusion threshold/) {
        $line = <INPUT>;
        next;
      }
      if ($line =~ m/Scores for complete hits/ and not $start_match) {
        $line = <INPUT>;
        $line = <INPUT>;
        $start_match = 1;
      } elsif ($start_match and $line ne "") {
        my @res = split(/\s+/, $line);
        my $name = $res[3];
        if (not defined $name) {
          print "Reading $input_file\n";
          exit;
        }
        my $direction = "forward";
        if ($res[4]-$res[5] > 0) {
          $direction = "reverse";
        }        
        if (not defined $results{$hmm}->{$name} or ($results{$hmm}->{$name}->[1] < $res[1])) {
          $results{$hmm}->{$name} = [$res[0], $res[1], $direction, $res[4], $res[5]];
        }
      } elsif ($start_match) {
        $quit = 1;
      }
      $line = <INPUT>;
    }
  }  
  close(INPUT);
  
  my $temp_file = Phylo::get_temp_file();  
  open(OUTPUT, ">$temp_file");  
  foreach my $key (keys %results) {
    my $hmm_name = $key;
    if (defined $prefix) {
      $hmm_name = $prefix;
    }
    foreach my $query (keys %{$results{$key}}) {
      my $new_name = $query;
      print OUTPUT "$new_name,$results{$key}->{$query}->[0],$results{$key}->{$query}->[1],$results{$key}->{$query}->[2],$hmm_name,$results{$key}->{$query}->[3],$results{$key}->{$query}->[4]\n";
    }
  }
  close(OUTPUT);  
  `mv $temp_file $output_file`;  
  return \%results;  
}

sub map_to_reference {
  my $start_time = time();
  print "Mapping to the human genome: $start_time\n";
  if ($opt{experiment} eq "wgs") {
    `bwa mem -t $opt{thread} $opt{human_reference} $opt{forward} $opt{reverse} | samtools view -b -o $opt{dir}/temp/temp.bam > $opt{dir}/logs/bwa.human.log 2>&1`;
    `mv $opt{dir}/temp/temp.bam $opt{dir}/temp/human.bwa.bam`;
  } elsif ($opt{experiment} eq "rnaseq") {
    #`hisat2 -x $opt{human_reference} -1 $opt{forward} -2 $opt{reverse} -p  $opt{thread} | samtools view -b -o $opt{dir}/temp/human.bwa.bam > $opt{dir}/logs/bwa.human.log 2>&1`;
  }
  my $end_time = time()-$start_time;    
  print "Finished mapping to the human genome: $end_time\n";  
  `echo "$end_time" > $opt{dir}/logs/bwa.human.time`;
}

sub prepare_unmapped_sequences {
  my $start_time = time();
  
  my $counter = 0;
  open CMD,'-|',"samtools fastq -f 4 $opt{bamfile}" or die $@;
  my $line;
  open(OUTPUT, ">$opt{dir}/temp/temp.fas");
  open(MAP, ">$opt{dir}/temp/temp.map");  
  while (defined($line=<CMD>)) {
    my $seqname = $line;
    chomp $seqname;
    my $name = $opt{prefix} . "_" . $counter;
    $counter++;
    
    my $seq = <CMD>;
    chomp $seq;
    my $plus = <CMD>;
    my $qual = <CMD>;
    print MAP "$seqname\t$name\n";
    print OUTPUT ">$name\n$seq\n";
  }
  close(MAP);  
  close(OUTPUT);
  close(CMD);
  `mv $opt{dir}/temp/temp.fas $opt{dir}/temp/unmapped.fas`;
  `mv $opt{dir}/temp/temp.map $opt{dir}/temp/unmapped.map`;  
  my $end_time = time() - $start_time;
}

sub join_fastq_file {  
  my $start_time = time();
  my $sample_number = int(compute_reads()+1);

  `mkdir -p $opt{dir}/splits/`;
  my %map = ("a","t","t","a","c","g","g","c", "A","T","T","A","C","G","G","C","N","N");  
  open MAP,">:gzip", "$opt{dir}/$opt{prefix}.map.gz" || die $!;            
  my $idx = 0;  
  open FAS,">", "$opt{dir}/splits/$idx.fas" || die $!;            
  open LOG,">", "$opt{dir}/logs/join.log" || die $!;            
  
  print LOG "Building $opt{dir}/splits/$idx.fas:\t0\n";
  
  open FORWARD,"<:gzip",$opt{forward} || open IN,$opt{forward} || die $!;      
  open REVERSE,"<:gzip",$opt{reverse} || open IN,$opt{reverse} || die $!;        
  
  my $counter = 0;
  my $ns = "N" x 10;
  while(my $forward_name = <FORWARD>){
    chomp $forward_name;
    my $forward_seq = <FORWARD>;
    chomp $forward_seq;
    
    my $forward_plus = <FORWARD>;
    my $forward_quality = <FORWARD>;

    my $reverse_name = <REVERSE>;
    chomp $reverse_name;
    my $reverse_seq = <REVERSE>;
    chomp $reverse_seq;
    
    $reverse_seq =~ tr/ACGTacgt/TGCAtgca/;
    
    my $reverse_plus = <REVERSE>;
    my $reverse_quality = <REVERSE>;

    my $seq = $forward_seq.$ns.reverse $reverse_seq;
    my $rev_comp = reverse $seq;
    $rev_comp =~ tr/ACGTacgt/TGCAtgca/;
    
    print FAS ">$opt{prefix}_$counter\_UF\n$seq\n";        
    print FAS ">$opt{prefix}_$counter\_UR\n$rev_comp\n";
    print MAP "$opt{prefix}_$counter\t$forward_name$reverse_name\n";
    
    $counter++;
    if ($counter % $sample_number ==  ($sample_number-1)) {
      print "$counter merging\n";
      $idx++;
      close(FAS);
      open FAS,">", "$opt{dir}/splits/$idx.fas" || die $!;                  
      my $time = time()-$start_time;
      print LOG "Building $opt{dir}/splits/$idx.fas:\t$time\n";      
    }            
  }
  close(FAS);  
  close(FORWARD);
  close(MAP);
  close(REVERSE);
  my $time = time()-$start_time;
  print LOG "Finished:\t$time\n";  
}

sub subsample_and_split {
  my $sample_number = int(compute_reads()+1)*2;
  my $idx = 0;
  my $counter = 0;
  `mkdir -p $opt{dir}/splits/`;
  open(INPUT, "$opt{dir}/$opt{prefix}.fasta");
  open(OUTPUT, "$opt{dir}/splits/$idx.fas");  
  while (my $name = <INPUT>) {
    my $line = <INPUT>;
    print OUTPUT "$name$line";
    if ($counter % $sample_number == $sample_number-1) {
      close(OUTPUT);
      $idx++;
      open(OUTPUT, "$opt{dir}/splits/$idx.fas");        
    }
    $counter++;
  }
  close(OUTPUT);
}

sub select_hmms {
  if (not -e "$opt{dir}/temp/sample_result") {
    `cat $opt{dir}/temp/*.0.*.hmm_search > $opt{dir}/temp/sample_result`;
  }
  if (not -e "$opt{dir}/temp/sample_result.csv") {
    convert_search_into_scores("$opt{dir}/temp/sample_result","$opt{dir}/temp/sample_result.csv");
  }  
  my %scores = read_scores("$opt{dir}/temp/hmm_search.csv", "single");  

}

sub process_results {
  my @files = <$opt{dir}/temp/*.hmm_search>;
  my %jobs = ();
  foreach my $file (@files) {
    $file =~ m/(\d+)\.hmm_search/;
    my $hmm = $1;
    $jobs{"$hmm"} = ["hmm_search", $file, "hmm_csv", "$file.csv", "subset", 0, "hmm", $hmm];
  }
  if (not -e "$opt{dir}/temp/reduced.csv") {
    print "Building csv\n";
    my %reduced :shared;  
    job_pool(\%jobs,new Worker({'function',\&read_hmms_worker,'name','read_hmms','threads_per_job', 1, "process", \&reduce_results, "results", \%reduced}));
    open(OUTPUT, ">$opt{dir}/temp/reduced.csv");
    local $" = ",";
    foreach my $key (keys %reduced) {
      print OUTPUT "@{$reduced{$key}}\n";
    }
    close(OUTPUT);
  }
  
  if (not -e "$opt{dir}/temp/reduced.fas") {
    print "Reading csv\n";  
    my %scores :shared;
    my %temp_scores = %{read_scores("$opt{dir}/temp/reduced.csv", "combined")};
    foreach my $key (keys %temp_scores) {
      $scores{$temp_scores{$key}->[0]} = "";
      delete $temp_scores{$key};
    }
    #my %reads :shared;
    my @files = <$opt{dir}/splits/*.fas>;
    my %join_jobs = ();
    foreach my $file (@files) {
      $file =~ m/(\d+)\.fas/;
      my $subset = $1;
      if (not -e "$opt{dir}/temp/subset.$subset") {
        $join_jobs{$subset} = ["subset_file", $file, "scores", \%scores, "subset", $subset];
      }
    }
    print "Running jobpool\n";      
    #job_pool(\%join_jobs,new Worker({'function',\&find_reads_worker,'name','find_reads','threads_per_job', 1, "process", \&reduce_reads, "results", \%reads}));    
    job_pool(\%join_jobs,new Worker({'function',\&find_reads_worker, 'name','find_reads','threads_per_job', 1}));
    `cat $opt{dir}/temp/subset.* > $opt{dir}/temp/reduced.fas`;
    `rm $opt{dir}/temp/subset.*`;
  }
}

sub reduce_results {
  my $hmm_results = $_[0];
  my $shared_results :shared;
  $shared_results = $_[1];
  {
    lock($shared_results);
    foreach my $name (keys %{$hmm_results}) {
      if (not defined $shared_results->{$name} or $shared_results->{$name}->[2] < $hmm_results->{$name}->[2]) {
        $shared_results->{$name} = $hmm_results->{$name};
      }    
    }
  }
}

sub combined_map {
  my $keep = $_[0];
  my $shared_results :shared;
  $shared_results = $_[1];
  {
    lock($shared_results);
    if (not defined $shared_results->{1}) {
      $shared_results->{1} = {};
      $shared_results->{2} = {};
    }
    foreach my $name (keys %{$keep}) {
      $keep->{$name} =~ m/(@[^@]+)(@[^@]+)/;
      $shared_results->{1}->{$1}="";
      $shared_results->{2}->{$2}="";
    }
  }
}

sub combined_fastq {
  my $keep = $_[0];  
  my $shared_results :shared;
  $shared_results = $_[1];
  {
    lock($shared_results);
    foreach my $read (keys %{$keep}) {
      foreach my $name (keys %{$keep->{$read}}) {
        $shared_results->{$read}->{$name}=$keep->{$read}->{$name};      
      }
    }
  }
}


sub read_scores {
  my $input = $_[0];
  my $type = $_[1];
  my $evalue = 10;
  
  my %scores = ();
  open(INPUT, "$input");
  #my $line = <INPUT>;
  while (my $line = Phylo::trim(<INPUT>."")) {
    my @results = split(/,/, $line);
    my $name = $results[0];
    if ($results[1] > $evalue) {
      next;
    }
    if (not defined $scores{$name} or $scores{$name}->[2] < $results[2]) {
      $scores{$name} = \@results;
    }
  }  
  close(INPUT);  
  return \%scores;
}


use Thread::Queue;

sub extract_reads_from_fastq {
  if (not -e "$opt{dir}/temp/reduced.fas" or -s "$opt{dir}/temp/reduced.fas" == 0) {
    return;
  }
  my %keep :shared;
  my %other = %{Phylo::read_fasta_file("$opt{dir}/temp/reduced.fas",0)};
  foreach my $key (keys %other) {
    my $read = $key;
    $read =~ s/_[^\_]+$//;
    delete $other{$key};
    $keep{$read} = "";
  }
  
  my %map :shared;
  my %read1 :shared;
  my %read2 :shared;
  
  $map{1} = \%read1;
  $map{2} = \%read2;
}

sub extract_reads_unmapped {
  if (not -e "$opt{dir}/temp/reduced.fas" or -s "$opt{dir}/temp/reduced.fas" == 0) {
    return;
  }
  my %keep = ();
  my %other = %{Phylo::read_fasta_file("$opt{dir}/temp/reduced.fas",0)};
  foreach my $key (keys %other) {
    my $read = $key;
    $read =~ s/_[^\_]+$//;
    delete $other{$key};
    $keep{$read} = "";
  }
  
  my %map = ();
  
  open(INPUT, "$opt{dir}/temp/unmapped.map");
  while (my $line = <INPUT>) {
    chomp $line;
    my @results = split(/\s+/, $line);
    if (defined $keep{$results[1]}) {
      $map{$results[0]} = "";
    } else {
      delete $keep{$results[1]};
    }
  }
  close(INPUT);
  
  open(OUT1, ">:gzip", "$opt{dir}/temp/$opt{prefix}_R1.fastq.gz.temp");    
  open(OUT2, ">:gzip", "$opt{dir}/temp/$opt{prefix}_R2.fastq.gz.temp");    
  
  my $fh = undef;
  open($fh, "<","$opt{dir}/temp/unmapped.fastq");
  while (my $line = <$fh>) {
    chomp($line);
    $line =~ m/([^\s]+)\/(\d+)/;
    my ($name1,$read1) = ($1,$2);
    my $seq1 = <$fh>;
    my $comment1 = <$fh>;
    my $quality1 = <$fh>;

    $line = <$fh>;
    chomp($line);
    $line =~ m/([^\s]+)\/(\d+)/;
    my ($name2,$read2) = ($1,$2);
    my $seq2 = <$fh>;
    my $comment2 = <$fh>;
    my $quality2 = <$fh>;
    
    if (($name1 ne $name2) or ($read1 ne 1) or ($read2 ne 2)) {
      die "Failed to parse unmapped correctly\n";
    }    
    if (defined $map{$name1}) {
      print OUT1 "$name1$\/$read1\n$seq1$comment1$quality1";
      print OUT2 "$name2$\/$read2\n$seq2$comment2$quality2";      
    }
  }
  close(OUT1);
  close(OUT2);
  `mv $opt{dir}/temp/$opt{prefix}_R1.fastq.gz.temp $opt{dir}/temp/$opt{prefix}_R1.fastq.gz`;
  `mv $opt{dir}/temp/$opt{prefix}_R2.fastq.gz.temp $opt{dir}/temp/$opt{prefix}_R2.fastq.gz`;
}


sub extract_reads_fastq_gz {
  if (not -e "$opt{dir}/temp/reduced.fas" or -s "$opt{dir}/temp/reduced.fas" == 0) {
    return;
  }
  my %keep :shared;
  my %other = %{Phylo::read_fasta_file("$opt{dir}/temp/reduced.fas",0)};
  foreach my $key (keys %other) {
    my $read = $key;
    $read =~ s/_[^\_]+$//;
    delete $other{$key};
    $keep{$read} = "";
  }
  
  my %map :shared;
  my %read1 :shared;
  my %read2 :shared;
  
  $map{1} = \%read1;
  $map{2} = \%read2;
  
  open(INPUT, "$opt{dir}/temp/unmapped.map");
  while (my $line = <INPUT>) {
    chomp $line;
    my @results = split(/\s+/, $line);
    if (defined $keep{$results[1]}) {
      $map{1}->{$results[0]} = "";
      $map{2}->{$results[0]} = "";      
    } else {
      delete $keep{$results[1]};
    }
  }
  close(INPUT);
    
  job_pool(build_read_fastq_jobs(\%map),new Worker({'function',\&find_fastq_worker,'name','ParallelFastQRead', 'threads_per_job', $opt{thread}/2}));      
}

#Code example taken from http://www.perlmonks.org/index.pl?node_id=735923
sub worker {
  my $tid = threads->tid;
  my( $Qwork, $Qresults ) = @_;
  while( my $work = $Qwork->dequeue ) {
    my $result;
    my %command = @{$work};    
    my $temp_file = Phylo::get_temp_file();    
    print "Running /usr/bin/time -v -o $opt{dir}/logs/hmmsearch.$command{key}.time nhmmer -o $temp_file --cpu $command{thread} $command{hmm} $opt{dir}/$opt{prefix}.fasta > $opt{dir}/logs/hmmsearch.$command{key}.log 2>&1\n";
    `/usr/bin/time -v -o $opt{dir}/logs/hmmsearch.$command{key}.time nhmmer -o $temp_file --cpu $command{thread} $command{hmm} $opt{dir}/$opt{prefix}.fasta > $opt{dir}/logs/hmmsearch.$command{key}.log 2>&1`;        
    `mv $temp_file $opt{dir}/temp/$command{key}.hmm_search`;

    ## Process $work to produce $result ##    
    $result = "$tid : result for workitem $work";
    $Qresults->enqueue( $result );
  }
  $Qresults->enqueue( undef ); ## Signal this thread is finished
}

sub build_split_jobs {
  `mkdir -p $opt{dir}/splits/`;

  my $sample_size = int(compute_reads()+1)*2;  
  my $file_size = stat("$opt{dir}/$opt{prefix}.fasta")->size;
  my $chunk_size = $file_size/$opt{thread};
  my @array = map($_*$chunk_size,(0..($opt{thread})));
  my %jobs = ();
  foreach my $key (0..$opt{thread}-1) {
    $jobs{$key} = ["idx",$key,"location",\@array, "sample_size", $sample_size];
  }
  split_function($jobs{0});
  return \%jobs;    
}

sub split_function {
  my %commands = @{$_[0]};    
  my $idx = $commands{"idx"};
  my $start_pos =  $commands{"location"}->[$idx];
  my $end_pos = $commands{"location"}->[$idx+1];
  my $sample_size = $commands{"sample_size"};
  my $fh = undef;
  open($fh, "$opt{dir}/$opt{prefix}.fasta");
  seek ($fh, $commands{"location"}->[$commands{"idx"}], 0);

  #Now check to see if you're a start of a line, if not, trash it till you get to the next
  #valid sequence

  my $name = <$fh>;
  if ($name !~ m/^>/) {
    $name = <$fh>;
  }

  if ($name !~ m/^>/) {
    $name = <$fh>;
  }

  if ($name !~ m/^>/) {    
    print "WTF?\n";
    exit;
  }
  my $line = <$fh>;
  my $counter = 0;
  my $c_idx = 0;  
  my $temp_file = Phylo::get_temp_file(); 
  
  print "Building $opt{dir}/splits/$idx.$c_idx.fas\n";
  #print sprintf("Start Idx: %d, Current: %d, End: %d\n", $idx,tell $fh,$end_pos);  
  my $out = undef;
  open($out, ">$temp_file");
  print $out "$name$line";

  while ($name = <$fh> and ((tell $fh) < $end_pos or not defined $end_pos)) {
    $line = <$fh>;
    print $out "$name$line";  
    if ($counter % $sample_size == $sample_size-1) {
      close($out);
      `mv $temp_file $opt{dir}/splits/$idx.$c_idx.fas`;
      $c_idx++;
      print "Building $opt{dir}/splits/$idx.$c_idx.fas\n";       
      #print sprintf("End Idx: %d, Current: %d, End: %d\n", $idx,tell $fh,$end_pos);      
      open($out, ">$temp_file");  
    }
    $counter++;    
  }
  #print sprintf("End Idx: %d, Current: %d, End: %d\n", $idx,tell $fh,$end_pos);

  #Now check to see if the current line is in the middle of a sequence name, if so, print it out
  if ($name =~ m/^>/) {
    $line = <$fh>;
    if (Phylo::trim($line) ne "") {
      print $out "$name$line";    
    }
  }
  close($out);
  `mv $temp_file $opt{dir}/splits/$idx.$c_idx.fas`;  
}

sub build_hmm_jobs {
  my $reduced = $_[0];
  
  my %jobs = ();
  foreach my $hmm (sort keys %hmms) {
    my $key = "$hmm";
    if (not -e "$opt{dir}/temp/$key.hmm_search") {
      $jobs{$key} = ['hmm', "$hmms{$hmm}", 'key', $key, "result", "$opt{dir}/temp/$key.hmm_search", "file", "$opt{dir}/temp/unmapped.fas",'thread', $opt{thread_per_job}];
    }
  }
  return \%jobs;    
}


#Code example taken from http://www.perlmonks.org/index.pl?node_id=735923
sub job_pool {
  my %jobs = %{$_[0]};
  my $worker = $_[1];  #worker function
  our $THREADS = int($opt{thread}/$worker->{'threads_per_job'});
  my $Qwork = new Thread::Queue;
  my $Qresults = new Thread::Queue;

  ## Create the pool of workers
  my @pool = map{
      threads->create( $worker->{'function'}, $Qwork, $Qresults )
  } 1 .. $THREADS;

  ## Get the work items (from somewhere)
  ## and queue them up for the workers
  my @order = keys %jobs;
  #   if (scalar @order != 0) {
  #     Phylo::fisher_yates_shuffle(\@order);
  #   }
  foreach my $key (@order) {
    $Qwork->enqueue($jobs{$key});
  }
  
  ## Tell the workers there are no more work items
  $Qwork->enqueue( (undef) x $THREADS );

  ## Process the results as they become available
  ## until all the workers say they are finished.
  for ( 1 .. $THREADS ) {
      while( my $result = $Qresults->dequeue ) {
          ## Do something with the result ##
          $worker->{'process'}($result,$worker->{"results"});
      }
  }
  ## Clean up the threads
  $_->join for @pool;  
}

sub hmm_worker {
  my $tid = threads->tid;
  my( $Qwork, $Qresults ) = @_;
  while( my $work = $Qwork->dequeue ) {
    my $result;
    my %command = @{$work};    
    my $temp_file = Phylo::get_temp_file();
    if (not -e "$command{result}" or -s "$command{result}" == 0) {
      print "Running /usr/bin/time -v -o $opt{dir}/logs/hmmsearch.$command{key}.time nhmmer -o $temp_file --cpu $command{thread} $command{hmm} $command{file} > $opt{dir}/logs/hmmsearch.$command{key}.log 2>&1\n";
      `/usr/bin/time -v -o $opt{dir}/logs/hmmsearch.$command{key}.time nhmmer -o $temp_file --cpu $command{thread} $command{hmm} $command{file} > $opt{dir}/logs/hmmsearch.$command{key}.log 2>&1`;        
      `mv $temp_file $command{result}`;
    }

    ## Process $work to produce $result ##    
    $result = "$tid : result for workitem $work";
    $Qresults->enqueue( $result );
  }
  $Qresults->enqueue( undef ); ## Signal this thread is finished
}

sub split_worker {  
  my $tid = threads->tid;
  my( $Qwork, $Qresults ) = @_;
  while( my $work = $Qwork->dequeue ) {
    my $result;
    split_function($work);
    ## Process $work to produce $result ##    
    $result = "$tid : result for workitem $work";
    $Qresults->enqueue( $result );
  }
  $Qresults->enqueue( undef ); ## Signal this thread is finished
}

sub read_hmms_worker {  
  my $tid = threads->tid;
  my( $Qwork, $Qresults ) = @_;
  while( my $work = $Qwork->dequeue ) {
    my $result;
    my %command = @{$work};      
    if (not -e $command{hmm_csv} and -e $command{hmm_search}) {
      convert_search_into_scores($command{hmm_search},$command{hmm_csv},$command{hmm});
    }
    my $scores = read_scores($command{hmm_csv}, "single");    
    
    ## Process $work to produce $result ##    
    $result = "$tid : result for workitem $work";
    $Qresults->enqueue( $scores );
  }
  $Qresults->enqueue( undef ); ## Signal this thread is finished
}

sub build_read_map_jobs {
  my $keep = $_[0];
  my $file_size = stat("$opt{dir}/$opt{prefix}.map.gz")->size;
  my $thread = 1;
  my $chunk_size = $file_size/$thread;
  my @array = map(int($_*$chunk_size),(0..($thread)));
  my %jobs = ();
  foreach my $key (0..$thread-1) {
    $jobs{$key} = ["idx", $key,"location",\@array, "keep", $keep, "prefix", $opt{prefix}];
  }
  return \%jobs;    
}

sub build_read_fastq_jobs {
  my $keep = $_[0];
  my %jobs = ();
  foreach my $read (1,2) {
    my $file_size = "";
    my $file = "";
    if ($read == 1) {
      $file = $opt{forward};
      $file_size = stat("$opt{forward}")->size; 
    }
    else {
      $file = $opt{reverse};
      $file_size = stat("$opt{reverse}")->size;
    }
    my $thread = 1;
    my $chunk_size = $file_size/$thread;
    my @array = map(int($_*$chunk_size),(0..($thread)));
    foreach my $key (0..$thread-1) {
      $jobs{"$key.$read"} = ["outfile", "$opt{dir}/temp/$opt{prefix}_R$read.fastq.gz", "idx", $key,"location",\@array, "file", $file, "keep", $keep, "read", $read, "prefix", $opt{prefix}];
    }    
  }  
  return \%jobs;    
}



sub parallel_read_map {
  my %commands = @{$_[0]};    
  my $idx = $commands{"idx"};
  my $start_pos =  $commands{"location"}->[$idx];
  my $end_pos = $commands{"location"}->[$idx+1];
  $end_pos = undef;
  my $keep = $commands{"keep"};
  my $prefix = $commands{"prefix"};  
  my $fh = undef;
  my %map = ();
  open($fh, "<:gzip","$opt{dir}/$opt{prefix}.map.gz");
  seek ($fh, $commands{"location"}->[$commands{"idx"}], 0);
  my $line;
  my $count = 0;
  my $total = 0;
  my $to_find = scalar keys %{$keep};
  #while ($line = <$fh> and ((tell $fh) < $end_pos or not defined $end_pos)) {  
  while ($line = <$fh>) {
    #if ($line !~ m/^$prefix/) {
    #  next;
    #}
    chomp($line);
    my @results = split(/\t/, $line);
    if (defined $keep->{$results[0]}) {
      $map{$results[0]} = $results[1];
      $total++;
    }
    if ($count % 1000000 == 999999) {
      print sprintf("%d; Found %d; to go %d\n", $idx, scalar keys %map, $count);
    }
    if ($total ==  $to_find) {
      last;
    }
    $count++;
  }
  if ($line ne "") {
    chomp($line);
    my @results = split(/\t/, $line);
    if (defined $keep->{$results[0]}) {
      $map{$results[0]} = $results[1];
    }  
  }
  return \%map;
}

sub write_read_fastq {
  my %commands = @{$_[0]};    
  my $idx = $commands{"idx"};
  my $keep = $commands{"keep"};
  open(OUTPUT, ">:gzip", $commands{"outfile"});    
  my $fh = undef;
  open($fh, "<:gzip","$commands{file}") or open($fh, "<","$commands{file}") or die "Cannot open $commands{file}\n";
  seek ($fh, $commands{"location"}->[$commands{"idx"}], 0);
  my $found = 0;
  my $total = scalar keys %{$keep->{$commands{read}}};
  my $count = 0;
  while (my $line = <$fh>) {
    chomp($line);
    $line =~ m/([^\s]+)/;
    my $name = $1;
    if (defined $keep->{$commands{read}}->{$name}) {
      my $seq = <$fh>;
      my $comment = <$fh>;
      my $quality = <$fh>;
      print OUTPUT "$line\n$seq$comment$quality";
      $found++;
    }
    $count++;
    if ($count % 1000000 == 999999) {
      print sprintf("Found %d; %d lines\n",$found,$count);
    }
    if ($found == $total) {
      print sprintf("Found %d; %d lines\n",$found,$count);    
      last;
    }
  }
  close(OUTPUT);
}


sub find_map_worker {  
  my $tid = threads->tid;
  my( $Qwork, $Qresults ) = @_;
  while( my $work = $Qwork->dequeue ) {
    my $result;
    print "Starting $tid find_map_worker\n";    
    my $map = parallel_read_map($work);
    
    ## Process $work to produce $result ##    
    $result = "$tid : result for workitem $work";
    print "Ending $tid find_map_worker\n";    
    $Qresults->enqueue( $map );    
  }
  $Qresults->enqueue( undef ); ## Signal this thread is finished
}

sub find_fastq_worker {  
  my $tid = threads->tid;
  my( $Qwork, $Qresults ) = @_;
  while( my $work = $Qwork->dequeue ) {
    my $result;
    print "Starting $tid find_fastq_worker\n";    
    #my $map = parallel_read_fastq($work);
    my $map = write_read_fastq($work);
    
    ## Process $work to produce $result ##    
    $result = "$tid : result for workitem $work";
    print "Ending $tid find_fastq_worker\n";    
    #$Qresults->enqueue( $map );
    $Qresults->enqueue( $result );    
  }
  $Qresults->enqueue( undef ); ## Signal this thread is finished
}


sub find_reads_worker {  
  my $tid = threads->tid;
  my( $Qwork, $Qresults ) = @_;
  while( my $work = $Qwork->dequeue ) {
    my $result;
    my %keep = ();
    my %command = @{$work};      
    
    if (not -e "$opt{dir}/temp/subset.$command{subset}") {
      print "Starting $tid find_reads_worker for $command{subset_file}\n";    
      my $scores = $command{scores};    
      my $fid;
      open($fid, "$command{subset_file}");
      while (my $name = <$fid>) {
        chomp($name);
        my $sequence = <$fid>;
        chomp($sequence);
        $name =~ m/>(.*)/;
        my $name = $1;
        if (defined $scores->{$name}) {
          $keep{$name} = $sequence;        
        }
      }
      close($fid);
      Phylo::write_alignment(\%keep, "$opt{dir}/temp/subset.$command{subset}");
    }
    ## Process $work to produce $result ##    
    $result = "$tid : result for workitem $work";
    print "Ending $tid find_reads_worker\n";    
    $Qresults->enqueue($result);
    
  }
  $Qresults->enqueue( undef ); ## Signal this thread is finished
}

sub reduce_reads {
  my $reads = $_[0];
  my $shared_results :shared;
  $shared_results = $_[1];
  {
    lock($shared_results);
    foreach my $name (keys %{$reads}) {
      $shared_results->{$name} = $reads->{$name};
    }
  }
}


{
package Worker;
  sub new  { 
    my $class = shift; 
    my $args = shift;
    my $self = {'threads_per_job', 1, 'function', undef, 'process', \&process, 'results', undef};
    foreach my $key (keys %{$args}) {
      $self->{$key} = $args->{$key};
    }
    bless ($self, $class);
    return $self;
  }
  
  sub process {
    print "Done! $_[0]\n";
    return;
  }
  1;
}