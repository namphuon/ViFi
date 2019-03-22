#!/usr/bin/perl

package Phylo;

#
#This perl script performs lots of modular tasks.  It is used to generate statistical scores, eventually mask alignments, estimate phylogenetic trees, and score them.
########

#cat green_genes/SUBPROBLEMconservativetree | sed -e 's/[0-9][0-9]*/\!/g' -e 's/[^\!]//g'|tr -d '\n'|wc -c
use strict;
use warnings;
use File::Spec::Functions qw(rel2abs);
use File::Basename;
use File::stat;
use Time::localtime;
use PerlIO::gzip;
use POSIX qw(ceil floor);
use Scalar::Util qw(looks_like_number);
use Data::Dumper;
no strict 'refs';

my $temp_dir = $ENV{'TMPDIR'}; 

#Returns temp file
sub get_temp_file {
  my $curr_temp_dir = $_[0];
  if (not defined $curr_temp_dir) {
    $curr_temp_dir = $temp_dir;
    if (not defined $curr_temp_dir) {
      $curr_temp_dir = '/tmp';
    }
  }

  my $temp_file = $curr_temp_dir . "/". Phylo::get_random_name($curr_temp_dir);
  while (-e $temp_file) {
    $temp_file = $curr_temp_dir . "/".Phylo::get_random_name($curr_temp_dir);
  }
  return $temp_file;  
}

sub get_random_name {
  my $target_dir = $_[0];
    
  if (not defined $target_dir) {
    $target_dir = "";
  }
  local $"="";
  my @name_str = map { ("a".."z", "A".."Z", 0..9)[rand 62] } 1..20;
  my $name = "temp@name_str";
  while (-e "$target_dir$name") {
    @name_str = map { ("a".."z", 0..9)[rand 36] } 1..20;
    $name = "temp@name_str";
  }
  return $name;
}

sub read_fasta_file {
  my $input_file = $_[0];
  my $verify = $_[1];
  my %name_map = ();
  
  if (not defined $verify) {
    $verify = 1;
  }
  
  if (not -e $input_file) {
    print "$input_file does not exist!\n";
    exit;
  }
  open(FILE, "<$input_file");  
  my $line = "";
  $line = <FILE>;
  my $length = -1;  
  while (defined $line and $line !~ />([^\n\r\f]*)/) {
    $line = <FILE>;
  }
  while (defined $line and $line =~ />([^\n\r\f]*)/) {
    my $sequence = "";
    my $name = $1;
    $line = <FILE>;    
    while (defined $line and not $line =~ />([^\n\r\f]*)/) {
      $line =~ s/\s//g;
      $sequence = $sequence . trim($line);
      $line = <FILE>;
    }
    if ($length == -1 and $verify) {
      $length = length($sequence);
    } else {
      if ($length != length($sequence) and $verify) {
        print "Error in sequence length $name \n";
        return "Error in sequence length $name\n";
      }
    }
    if (length(trim($sequence)) == 0 and $verify) {
      print "Error empty alignment\n";
      return "Error empty alignment\n";
    
    }         
    if (not $sequence =~ /[^-?]/ and $verify) {      
      print "Error all indel characters $name\n";
      return "Error all indel characters $name\n";      
    }
    $name_map{$name} = $sequence;
  }
  close(FILE);
  return \%name_map;
}


#Trims whitespace from strings
sub trim {
  my $string = $_[0];
  if (not defined($string)) {
      return undef;
  }
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}

#Given an alignment map and output file, write to file
sub write_alignment {
  my %aln = %{$_[0]};
  my $output_file = $_[1];
    
  open(OUTPUT, ">$output_file");
  foreach my $key (sort keys %aln) {
    print OUTPUT ">$key\n$aln{$key}\n";
  }
  close(OUTPUT);
}
1;
