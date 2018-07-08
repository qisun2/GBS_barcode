#!/usr/local/bin/perl -w
use strict;
use Getopt::Std;

my $sizelimit = 40;


my %opts;

unless (getopts("c:b", \%opts))
{
    print "Command: GBS_barcode.pl <options> yourBarCodeFile yourEnzymeFile\n";
	print "Options: -b keep barcode in the output files\n";
	print "Options: -c check for 'N' in this number of base pair at the 5' region.\n";
	print "Note: Put your barcode file, enzyme file and all Illumina qseq or fastq files in one directory. \"cd\" to the directory. Then execute the command.\n";
	print "You can use the template file /shared_data/qisun/GBS_barcode_template.txt and and enzyme_template.txt\n";
	exit;
}

# constants
my ($barcodefile) = $ARGV[0];  #default $qseqORfastq
my ($enzymefile) = $ARGV[1];

my $keepbarcode = 0;
if (exists $opts{"b"}) 
{
	$keepbarcode = 1;
}

my $baseToCheckN =0;
if (exists $opts{"c"}) 
{
	$baseToCheckN =$opts{"c"} ;
}


our %compnt = ("A"=>"T", 
			"T"=>"A",
			"G"=>"C",
			"C"=>"G",
			"N"=>"N",
			);
unless ($barcodefile) 
{
	print "Command: GBS_barcode.pl <options> yourBarCodeFile yourEnzymeFile\n";
	print "Options: -b keep barcode in the output files\n";
	print "Note: Put your barcode file, enzyme file and all Illumina qseq or fastq files in one directory. \"cd\" to the directory. Then execute the command.\n";
	print "You can use the template file /shared_data/qisun/GBS_barcode_template.txt and enzyme_template.txt\n";
	exit;
}

unless ($enzymefile) 
{
	print "Command: GBS_barcode.pl <options> yourBarCodeFile yourEnzymeFile\n";
	print "Options: -b keep barcode in the output files\n";
	print "Note: Put your barcode file, enzyme file and all Illumina qseq or fastq files in one directory. \"cd\" to the directory. Then execute the command.\n";
	print "You can use the template file /shared_data/qisun/GBS_barcode_template.txt and and enzyme_template.txt\n";
	exit;
}


unless (-e $barcodefile)
{
	print "File $barcodefile does not exits!\n";
	exit;
}

unless (-e $enzymefile)
{
	print "File $enzymefile does not exits!\n";
	exit;
}

my $workingdir = "./";
if ($barcodefile=~/\//) 
{
	 $workingdir= $barcodefile;
	 $workingdir=~s/\/[^\/]+$//;
}
open (STAT, ">$workingdir/stat.txt") || die "Can not create statistic file. permission denied!";

#open OUT, ">tmp";
my ($enzyme, $finalsize, $ends, $restrictionsize, $EnzymeEndSize) = parse_enzyme_file($enzymefile);

		$enzyme=~s/^\s+//; $enzyme=~s/\s+$//;
		$enzyme=uc $enzyme;


$ends = uc $ends;
my @likelyends = $ends=~/(\w+)/g;
my ($barcoderef, $minBarLen, $maxBarLen) = get_barcode($barcodefile, $enzyme, $finalsize, $restrictionsize);
my %file2barcodes = %{$barcoderef};

my $line1;
my $line2;
my $line3;
my $line4;
my $line5;
my $line6;
my $line7;
my $line8;

my %ends2count = ();
my $barcode;
my $filehandle;
my ($pass_count_lane, $read_count_lane);
my ($sstr, $qstr, $sstr2, $qstr2, $linename, $samplename, $codesize, $paired, $codestr);
my %linecount;
#my %line;
my $id = 0;

my ($read_count,$pass_count) = 0;

##check files
foreach my $file (keys %file2barcodes)
{
	unless (-e "$workingdir/$file") 
	{
		print "File $workingdir/$file does not exists!!\n";
		exit;
	}
}
my %allsamples = ();
foreach my $file (keys %file2barcodes) 
{ 
	print "start read file $file\n";
  %linecount  = ();
  ($pass_count_lane, $read_count_lane) = (0, 0);

  my %tag2sample = %{$file2barcodes{$file}};
  my %filehandles = ();
  my $paired_flag = "s";
  FHhandle:foreach(values %tag2sample)
  {
		
		($samplename, $paired) = split "\t", $_;
		
		next FHhandle if (exists $filehandles{$samplename});
		local *FILE;
		local *FILE2;
		my $overwrite = ">";
		if (exists $allsamples{$samplename}) 
		{
			$overwrite = ">>";
		}
		if ($paired =~/p/i) 
		{
			
			open(FILE, "${overwrite}$workingdir/${samplename}_1.fastq") || die "Error: cannot write to file $!\n";
			#print "xxx${_}xxx\n";
			$filehandles{$samplename} = *FILE;
			print "output file ${samplename}_1.fastq ${samplename}_2.fastq\n";
			open(FILE2, "${overwrite}$workingdir/${samplename}_2.fastq") || die "Error: cannot write to file $!\n";
			#print "xxx${_}xxx\n";
			$filehandles{"${samplename}_2"} = *FILE2;
			$paired_flag = "p";
		}
		else
		{
			open(FILE, "${overwrite}$workingdir/${samplename}.fastq") || die "Error: cannot write to file $!\n";
			print "output file ${samplename}\n";
			#print "xxx${_}xxx\n";
			$filehandles{"${samplename}"} = *FILE;
		}
		$allsamples{$samplename} = "";
   }

	#exit;
	if ($file=~/\.gz$/) 
	{
		open IN, "gunzip -c $workingdir/$file |" or die "ERROR: Cannot open input file: $!\n\n\n";
	}
	else
	{
		open IN, "<$workingdir/$file" or die "ERROR: Cannot open input file: $!\n\n\n";
	}
	print  "\n\nReading data from the input file:   $file...\n";

	if ( $paired_flag eq "p") 
	{
		my $file2 = $file;
		$file2=~s/(.+)_1/${1}_2/;
		if ($file2 eq $file) 
		{
			print "Error: $file is not the name of the paired end file\n";
			exit;
		}
		if ($file=~/\.gz$/) 
		{
			open IN2, "gunzip -c $workingdir/$file2 |" or die "ERROR: Cannot open input file: $!\n\n\n";
		}
		else
		{
			open IN2, "<$workingdir/$file2" or die "ERROR: Cannot open input file: $!\n\n\n";
		}

		#open IN2, "<$workingdir/$file2" or die "ERROR: Cannot open input file: $!\n\n\n";
		print "  --- paired-end files ...\n";
	}

	##test qseq or fastq
	my $testline1 = <IN>; <IN>; my $testline3= <IN>; close IN;
	my $fileformat = "qseq";
	if (($testline1 =~ /^\@/) && ($testline3 =~ /^\+/))
	{
		$fileformat = "fastq";
	}
	elsif (($testline1 =~ /\t/) && ($testline3 =~ /\t/))
	{
		$fileformat = "qseq";
	}
	else
	{
		print "Error: File not qseq or fastq file\n";
		exit;
	}
	close IN;

	#reopen the file
	if ($file=~/\.gz$/) 
	{
		open IN, "gunzip -c $workingdir/$file |" or die "ERROR: Cannot open input file: $!\n\n\n";
	}
	else
	{
		open IN, "<$workingdir/$file" or die "ERROR: Cannot open input file: $!\n\n\n";
	}

	my $totalrecords = 0;
	FILELOOP:while ($line1 = <IN>) {
	#my $line1_whole = $_;
	$totalrecords ++;
	chomp $line1;
	#last FILELOOP if ($totalrecords > 10000);
	if ($fileformat=~/qseq/) 
	{	
		my @data = split "\t", $line1;
		($line2, $line4) = @data[8, 9];
	}
	else
	{
		$line2 =<IN>;
		chomp $line2;
		$line3 =<IN>;
		chomp $line3;
		$line4 =<IN>;
		chomp $line4;
	}

	($sstr, $qstr, $sstr2, $qstr2, $codestr, $linename, $barcode, $samplename,  $codesize) = ();
	#$readsize=0;
				$sstr = $line2;
				$qstr = $line4;
				$codestr = "";
            TAGLOOP:for (my $j=$minBarLen;$j<=$maxBarLen;$j++){

                my $code = substr($sstr,0,($j+$restrictionsize));
                if (exists $tag2sample{$code}){
                    ($samplename, $paired, $codesize) = split "\t", $tag2sample{$code};
					$sstr = substr ($sstr, $codesize);
					$qstr = substr ($qstr, $codesize);
					$codestr = $code;
					last TAGLOOP;
                }
            }


#	TAGLOOP:foreach my $code(keys %tag2sample) 
#	{
#		if ($line2 =~ /^$code/)
#		{
#			$sstr = $1;
#			$linename = $tag2sample{$code};
#			($samplename, $paired, $readsize, $codesize) = split "\t", $linename;
#			$barcode = $code;
#			last TAGLOOP;
#		}
#	}
	#my $line2_whole = "";

	if ($paired_flag eq "p")
	{
		if ($fileformat=~/qseq/) 
		{
			$line5 = <IN2>;
			chomp $line5;
			my @data = split "\t", $line5;
			($line6, $line8) = @data[8, 9];
		}
		else
		{
			$line5 = <IN2>;
			chomp $line5;
			$line6 =<IN2>;
			chomp $line6;
			$line7 = <IN2>;
			chomp $line7;
			$line8 =<IN2>;
			chomp $line8;
		}

		$sstr2 = $line6;

		$qstr2 = $line8;

	}

	
	next FILELOOP unless ($samplename);
	my ($seqsize1, $seqsize2) = (0, 0);
	if ($keepbarcode ==0) 
	{
		my $checkNstr = $sstr;
		if ($baseToCheckN>0) 
		{
			$checkNstr = substr($sstr, 0, $baseToCheckN);
		}

		next FILELOOP if ($checkNstr=~/\.|N/) ;

		
			ELOOP1:for (my $k=0;$k<@likelyends;$k++){
				if ($sstr =~ /$likelyends[$k]/){
						
				  #my $sqlen = length($sstr);
				  my $pos = index($sstr,$likelyends[$k]);
				  if ($pos<$sizelimit) 
				  {
					  next FILELOOP;
				  }
				#substr ($sstr, ($pos + $EnzymeEndSize)) = "N" x 100;
				#substr ($qstr, ($pos + $EnzymeEndSize)) = "!" x 100;
				$seqsize1 = $pos + $EnzymeEndSize;
				$ends2count{$likelyends[$k]} ++;
				last ELOOP1
				}
			 }
			 if ( $paired_flag eq "p")  
			 {
				my $endstr = revcom_adaptor($codestr);
				if ($sstr2 =~ /$endstr/){
						
				  #my $sqlen = length($sstr2);
				   my $pos = index($sstr2,$endstr);
				   #substr ($sstr2, ($pos + $EnzymeEndSize)) = "N" x 100;
				   #substr ($qstr2, ($pos + $EnzymeEndSize)) = "!" x 100;
				   $seqsize2 = $pos + $EnzymeEndSize;
				   #$ends2count{$likelyends[$k]} ++;
				}


			 }


		
		  $sstr=substr($sstr, 0, $finalsize);

		  $qstr=substr($qstr, 0, $finalsize);


	}



	  #print OUT $line1_whole, $line2_whole, $barcode, "\n";
	  if ( $paired_flag eq "p")  
	  {
		    if ($seqsize1 == $seqsize2) 
			{		  	$id ++;
				$linecount{$samplename} ++;
		  		$filehandle = $filehandles{$samplename};
				if ($keepbarcode ==0) 
				{
					print $filehandle "\@${id}/1\n";
					print $filehandle "$sstr\n";
					print $filehandle "+${id}/1\n";
					print $filehandle "$qstr\n";
				}
				else
				{
					print $filehandle "$line1\n";
					print $filehandle "$line2\n";
					print $filehandle "$line3\n";
					print $filehandle "$line4\n";
				}

		  		$sstr2=substr($sstr2, 0, $finalsize);
		  		$qstr2=substr($qstr2, 0, $finalsize);
				$filehandle = $filehandles{"${samplename}_2"};
				#print $filehandle "$barcode\n";
				#print $filehandle "$line2_whole";
				if ($keepbarcode ==0) 
				{
					print $filehandle "\@${id}/2\n";
					print $filehandle "$sstr2\n";
					print $filehandle "+${id}/2\n";
					print $filehandle "$qstr2\n";
				}
				else
				{
					print $filehandle "$line5\n";
					print $filehandle "$line6\n";
					print $filehandle "$line7\n";
					print $filehandle "$line8\n";
				}
		    }


	  }
	  else
	 {

		 $id ++;
		$linecount{$samplename} ++;
		$filehandle = $filehandles{$samplename};
		if ($keepbarcode ==0) 
		{
			print $filehandle "\@$id\n";
			print $filehandle "$sstr\n";
			print $filehandle "+$id\n";
			print $filehandle "$qstr\n";
		}
		else
		{
			print $filehandle "$line1\n";
			print $filehandle "$line2\n";
			print $filehandle "$line3\n";
			print $filehandle "$line4\n";
		}


	 }

	}


  close IN;


   foreach(values %filehandles)
  {
    close $_;
   }
    print STAT "File: $file\tTotal reads: $totalrecords\n";
	my $passed  =0;
	my %allsamples = ();
	foreach my $code(keys %tag2sample) 
	{
		my ($samplename, $paired, $codesize) = split "\t", $tag2sample{$code};
		$allsamples{$samplename} = 0
	}
	
	foreach my $samplename (keys %allsamples) {
		my $count =0;
		if (exists $linecount{$samplename}) 
		{
			$count = $linecount{$samplename};
		}
		print STAT $samplename,  "\t", $count,"\t",  "\n";
		$passed += $count;
	}
	print STAT "-------------\n";
	print STAT "Combined: $passed reads passed filter\n";
	print STAT "-------------\n\n\n";

	print STAT "Incomplete cuts or short fragments\n";
	foreach  (keys %ends2count) 
	{
		print STAT $_, "\t", $ends2count{$_}, "\n";
	}
}


#print STAT "Reads per Line:";
#foreach  (keys %line) 
#{
#	print STAT $_, "\t", $line{$_}, "\n";
#}
#
#
#print STAT "\n\nFinished!\n\n" 
#           . "  In total, $pass_count reads passed\n"
#           . "  out of    $read_count examined\n\nSee the following file for output:\n\t$outfile\n";
#
#print STAT "\n\n\nFare thee well, Panzean!!\n\n\n";

sub get_barcode
{
	my ($barcodefile, $enzyme, $finalsize, $restrictionsize) = @_;
	open (IN, $barcodefile) || die "Error: cannot open bar code file!!!\n";
	my $line = <IN>;
	chomp $line;
	$line=~s/"//g;
	my @data = split "\t", $line;
	unless ($data[0]=~/SampleName/i) 
	{
		print "Error: bar code file not in right format! first column title should be  'SampleName'\n";
		exit;
	}
	unless ($data[3]=~/file/i) 
	{
		print "Error: bar code file not in right format! 4th column title should be  'file'\n";
		exit;
	}

	my %file2code; 
	my %allsamplenames = ();
	my ( $minBarLen, $maxBarLen) = (100, 0);
	LOOP:while (<IN>) 
	{
		chomp;
		s/"//g;
		my ($sample, $code, $paired, $file) = split "\t";
		next LOOP unless ($file);
		next LOOP unless ($sample=~/\w/);


		$file=~s/\s//g;
		$file=~s/.+\///;


		
		$sample=~s/\s//g;

		$code=~s/^\s+//; $code=~s/\s+$//;
		$code = uc $code;

		my $len = length($code);
        if ($len > $maxBarLen){
            $maxBarLen = $len;
        }
        if ($len < $minBarLen){
            $minBarLen = $len;
        }
		
		#if (exists $allsamplenames{$sample}) 
		#{
		#	print "Error: sample $sample duplicated!!!\n";
		#	exit;
		#}
		$allsamplenames{$sample} = "";

		my $enzymeforsize = $enzyme;
		$enzymeforsize=~s/\[(\w+)\]/N/g;
		my $variations = $1;
		if ($variations) 
		{
			my @variations =  $variations=~/([ACGT])/g;
			#my $restsize = $finalsize - $restrictionsize;
		
			foreach  (@variations) 
			{
				my $enzymevariation = $enzymeforsize;
				$enzymevariation=~s/N/$_/;
				${$file2code{$file}}{$code."$enzymevariation"}  = $sample. "\t". $paired. "\t".  (length $code) ;
			}
		}
		else
		{
			${$file2code{$file}}{$code.$enzyme}  = $sample. "\t". $paired. "\t".  (length $code) ;
		}
	
	}

	if (($minBarLen==100) ||  ($maxBarLen==0)) 
	{
		print "no valid info in barcode file!"; exit;
	}
	return (\%file2code, $minBarLen, $maxBarLen) ; 
}

sub parse_enzyme_file
{
	my $enzymefile = shift;
	open (IN, $enzymefile) || die "Error: cannot open bar code file!!!\n";
	my ($enzyme, $size, $ends ) = ("", 0, "");
	my $count =0;
	while (my $line=<IN>) 
	{
		if ($line=~/^Enzyme:/) 
		{
			$line=~/^Enzyme:(.*)/;
			$enzyme = $1;
			$enzyme =~s/^\s+//;
			$enzyme =~s/\s+$//;
			$count ++;
		}
		elsif ($line=~/^FinalSize:/)
		{
			$line=~/^FinalSize:(.*)/;
			$size = $1;
			$size =~s/^\s+//;
			$size =~s/\s+$//;
			$count ++;
		}
		elsif ($line=~/^Ends:/)
		{
			$line=~/^Ends:(.*)/;
			$ends = $1;
			$ends=~s/^\s+//;
			$ends=~s/\s+$//;
			$count ++;
		}
		elsif ($line=~/^EnzymeEndSize:/)
		{
			$line=~/^EnzymeEndSize:(.*)/;
			$EnzymeEndSize = $1;
			$EnzymeEndSize=~s/^\s+//;
			$EnzymeEndSize=~s/\s+$//;
			$count ++;
		}
	}
	unless ($count == 4) 
	{
		print "Enzyme file either does not exist or there is a problem!\n";
		exit;
	}

	my $enzymeforsize = $enzyme;
	$enzymeforsize=~s/\[(\w+)\]/N/g;
	my $restrictionsize = length $enzymeforsize;


	return ($enzyme, $size, $ends, $restrictionsize, $EnzymeEndSize );
}


sub revcom_adaptor
{
	my $input = shift;
	my $output = "";
	my @nt = reverse ($input=~/(\w)/g);
		foreach  (@nt) 
		{
			$output .=$compnt{$_};
		}
		$output .="AGATCGGA";
		$output = substr($output, 0, 8);
	return $output;
}
