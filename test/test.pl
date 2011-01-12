#!/usr/bin/perl -w
#
#########################################################################
#
# Shell script for running DALTON perl-based tests 
#
# The tests are defined in the file "checklist".
#
# First version by Luca Frediani, 2006
#
#########################################################################
#
#
# Revised Sep-2007 Hans Joergen Aa. Jensen
#   - ignore case when searching for PATTERN (as if "grep -i PATTERN"  instead of "grep PATTERN")
#
# todo
# override mechanism is primitive: list-type overrides need to be 
# "overridden" with lists. It should not be necessary
#
# Setting default command line options
#
my $options = "";
my $all = 0;
my $checkref = 0;
my $help = '';
my $keep = '';
my $explist ='';
my $check = '';
my $keyword = '';
my $oldout = '';
my $quiet = '';
my $biglog = 'perl_test.log';
my $bigerr = 'perl_test.err';
chomp(my $tstdir =`pwd`);
my $dalton   = $tstdir."/../bin/dalton";
$verbose = 0;
%verb = ('test_inp'  => 0,
	 'fetch_out' => 1,
	 'check_res' => 2,
	 'fetch_inp' => 3,
	 'check_det' => 4,
	 'all'       => 5);
use Getopt::Long;
use Pod::Usage;
use Tie::File;
#use Getopt::Long 2.34;
parse_input();
##############################################################################

############################################################################
# define useful tools, main code at the bottom.
############################################################################

############################################################################
sub parse_input {    
    my $result = GetOptions ("dalton=s" => \$dalton,  
			     "options=s"=> \$options, 
			     "keyword=s"=> \$keyword, 
			     "check=s"  => \$check, 
			     "err=s"    => \$bigerr, 
			     "log=s"    => \$biglog,
			     "list=s"   => \$explist, 
			     "verbose=i"=> \$verbose,
			     "all"      => \$all,
			     "checkref" => \$checkref,
			     "keep"     => \$keep,
			     "quiet"    => \$quiet,
			     "oldout"   => \$oldout,
			     "tstdir=s" => \$tstdir,
			     "help"     => \$help);
    pod2usage(1) if $help;
    pod2usage(3) if $result != 1;
}
############################################################################
sub load_def_check_param ($$$) {
#
# Load defaults values for the check parameters
# Defaults are loaded from a corresponding default check in checklist
#
    my %check = %{$_[0]};
    my $tstdir = $_[1];
    my %files = %{$_[$#_]};
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});
    my $result = 0;

    if ((not defined $check{type}) or $check{type} < 0) {
	print $log "The present check does not have a valid type\n";
	print $log "Check type set to zero\n";
	$check{type} = 0;
    }

    my $type = $check{type};
    print $log "Setting default values for present check. Type $type\n" if $verbose > $verb{check_det};
    $def_name="default_".$type;

    my %def_check = load_check($def_name,$tstdir,\%files);

    if ($verbose > $verb{all}) {
	print $log "Check list before setting default values:\n";
	print_hash(\%check,$log);
	print $log "Check list containing default values:\n";
	print_hash(\%def_check,$log);
    } 

    foreach $def_key (keys %def_check) {
	if ((not defined $check{$def_key})) {
	    print $log "Setting default value for parameter $def_key\n" if $verbose > $verb{check_det};
	    $check{$def_key} = $def_check{$def_key};
	}
    }
    if ($verbose > $verb{all}) {
	print $log "Check list after setting default values:\n";
	print_hash(\%check,$log);
    }
    %{$_[0]} = %check;
}
############################################################################
sub control_check_list ($$) {
#
# Check the integrity if each loaded check. 
# If the check is incomplete or not well defined it is skipped
#
    my %check = %{$_[0]};
    my %files = %{$_[$#_]};
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});
    my $result = 0;
    if ((not defined $check{type}) or $check{type} < 0) {
	print $log "The present check does not have a valid type and will be skipped\n";
	print $err "The present check does not have a valid type and will be skipped\n";
	return $result;
    }
    if (not defined $check{name}) {
	print $log "The present check does not have a name and will be skipped\n";
	print $err "The present check does not have a name and will be skipped\n";
	return $result;
    }
    my $type = $check{type};
    my $name = $check{name};
    print {$log} "Integrity verification of check $name. Type $type\n" if $verbose > $verb{check_det};

    my %key_list =();
  SWITCH: {
      if ($type == 0) {
	  %key_list=(string => 'string',
		     name   => 'string',
		     type   => 'integer',
		     pos    => 'integer',
		     abs    => 'integer',
		     thr    => 'float');
	  %max_val=(abs => 1);
	  %min_val=(abs => 0,
		    pos => 1,
		    thr => 1.0e-20);
	  last SWITCH;
      } 
      if ($type == 1) {
	  %key_list=(string => 'string',
		     name   => 'string',
		     type   => 'integer',
		     pos    => 'integer',
		     maxout => 'integer',
		     maxget => 'integer',
		     abs    => 'integer',
		     thr    => 'float');
	  %max_val=(abs => 1);
	  %min_val=(abs => 0,
		    pos => 1,
		    maxout => -1,
		    maxget => -1,
		    thr => 1.0e-20);
	  last SWITCH;
      } 
      if ($type == 2) {
	  %key_list=(start_string   => 'string',
		     end_string     => 'string',
		     name           => 'string',
		     start_offset   => 'integer',
		     end_offset     => 'integer',
		     type           => 'integer',
		     pos            => 'integer',
		     line           => 'integer',
		     abs            => 'integer',
		     thr            => 'float');
	  %max_val=(abs => 1);
	  %min_val=(abs => 0,
		    start_offset => 0,
		    pos => 1,
		    line => 1,
		    thr => 1.0e-20);
	  last SWITCH;
      } 
      if ($type == 3) {
	  %key_list=(string   => 'string',
		     name     => 'string',
		     type     => 'integer',
		     maxout   => 'integer',
		     maxget   => 'integer',
		     pos      => 'integer_list',
		     abs      => 'integer_list',
		     thr      => 'float_list');
	  %max_val=(abs => 1);
	  %min_val=(abs => 0,
		    pos => 1,
		    maxout => -1,
		    maxget => -1,
		    thr => 1.0e-20);
	  last SWITCH;
      } 
      if ($type == 4) {
	  %key_list=(name   => 'string',
		     string => 'string',
		     type   => 'integer',
		     maxout => 'integer',
		     post   => 'integer',
		     rembeg => 'integer',
		     remend => 'integer');
	  %max_val=(post   => 5);
	  %min_val=(rembeg => 0,
		    remend => 0,
		    post   => 0,
		    maxout => -1,
		    thr => 1.0e-20);
	  last SWITCH;
      } 
      if ($type == 5) {
	  %key_list=(  name         => 'string',
		       type         => 'integer',
		       start_string => 'string',
		       end_string   => 'string',
		       start_offset => 'integer',
		       end_offset   => 'integer',
		       post         => 'integer');
	  %max_val=(post         => 0);
	  %min_val=(start_offset => 0,
		    post         => 0);
	  last SWITCH;
      } 
      if ($type == 6) {
	  %key_list=(start_string => 'string',
		     end_string   => 'string',
		     name         => 'string',
		     type         => 'integer',
		     start_offset => 'integer',
		     end_offset   => 'integer',
		     pos          => 'integer_list',
		     rel          => 'integer_list',
		     abs          => 'integer_list',
		     thr          => 'float_list');
	  %max_val=(abs => 1, 
	            rel => 1);
	  %min_val=(abs => 0,
		    rel => 0,
		    pos => 1,
		    start_offset => 0,
		    thr => 1.0e-13);
	  last SWITCH;
      } 
  }
    $result = verify_hash_element_types(\%check,
					\%key_list,
					\%min_val,
					\%max_val,
					\%files);
    return $result;
}
############################################################################
sub check_hash_element($$$$$$) {
    my $key   = $_[0];
    my $val   = $_[1];
    my $type  = $_[2];
    my $min   = $_[3];
    my $max   = $_[4];
    my %files = %{$_[$#_]};
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});
    my $result = 0;
    my $error = 0;
    if ($type != check_num($val)) {
	print $log "wrong type for parameter $key\n" if $verbose > $verb{check_det};
	print $err "wrong type for parameter $key\n";
	$error++;
    } 
    if (defined $min && $val < $min) {
	print $log "value of parameter $key is not allowed. Value=$val Min=$min\n" if $verbose > $verb{check_det};
	print $err "value of parameter $key is not allowed. Value=$val Min=$min\n";
	$error++;
    } 
    if (defined $max && $val > $max) {
	print $log "value of parameter $key is not allowed. Value=$val Max=$max\n" if $verbose > $verb{check_det};
	print $err "value of parameter $key is not allowed. Value=$val Max=$max\n";
	$error++;
    }
    $result = 1 if ($error == 0);
    return $result;
}
############################################################################
sub verify_hash_element_types($$$$$) {
#
# This routine verifies that the values of the checklist hash
# are existing and contain what they should.
#
    my %right_type=('string'       => 0,
		    'integer'      => 1,
		    'float'        => 2,
		    'integer_list' => 1,
		    'float_list'   => 2);
    
    my %check    = %{$_[0]};
    my %key_list = %{$_[1]};
    my %min = %{$_[2]};
    my %max = %{$_[3]};
    my %files    = %{$_[$#_]};
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});
    my $n_right = 0;
    my $n_wrong = 0;
    my $result = 0;
    my $type = "";
    my $value = "";
    my $ref_type = "";
    
    foreach $key (keys(%key_list)) {
	if (defined $check{$key}) {
	    $type  = $right_type{$key_list{$key}};
	    $value = $check{$key};
	    $ref_type = ref($value);
	    if ($ref_type eq 'ARRAY') {
		my $n_array_right = 0;
		my $n_array_wrong = 0;
		foreach $elm (@{$value}) {
		    if(check_hash_element($key,$elm,$type,$min{$key},$max{$key},\%files)) {
			$n_array_right++;
		    } else { 
			$n_array_wrong++;
		    }
		}
		if (($n_array_right == @{$value}) and ($n_array_wrong == 0)) {
		    $n_right++;
		    if ($verbose > $verb{all}) {
			print $log "Array \"$key\" has been checked\n";
			print $log "Content \"@{$value}\" is correct\n";
		    }
		} else {
		    $n_wrong++;
		    print $err "Array \"$key\" has been checked\n";
		    print $err "Content \"@{$value}\" is NOT correct\n";
		}
	    } else {
		if(check_hash_element($key,$value,$type,$min{$key},$max{$key},\%files)) {
		    $n_right++;
		} else { 
		    $n_wrong++;
		}
	    }
	} else {
	    $n_wrong++;
	    print $err "Element \"$key\" has been checked\n";
	    print $err "Content is NOT defined\n";
	}
    }
    $result = 1 if (($n_right == (keys(%key_list))) and ($n_wrong == 0));
    return $result;	    
}
############################################################################
sub load_check($$$) {
#
# Load the check defined by the first argument from the library of available checks
#
    my $check  = $_[0];
    my $tstdir = $_[1];
    my %files = %{$_[$#_]};
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});
    my $file = "$tstdir/checklist";
    my $params = "";
    open OUT, $file or die;
    my @list = get_portion("$check start","$check end",1,1,\%files);
    close OUT;
    if (@list) {
	my $string = '';
	foreach (@list) {$string = $string.$_}
	if ($verbose > $verb{check_det}) {
	    print $log "----$check start----------------\n";
	    print $log "$string\n";
	    print $log "----$check end------------------\n";
	}
	eval $string;
    } else {
	$params={name => "$check"};
	print {$log} "Warning: Check $check not found on checklist file\n" if $verbose > $verb{test_inp};
	print {$err} "Warning: Check $check not found on checklist file\n";
    } 
    print {$log} "Check $check has been loaded\n" if $verbose > $verb{test_inp};
    return %$params;
}
###########################################################################
sub chomp_and_remove_spaces {
    foreach $_ (@_) {
	chomp;
#remove all spaces
	s/(\s*)//;
    }
}
###########################################################################
sub chomp_and_shorten_spaces {
    foreach $_ (@_) {
	chomp;
#remove all leading spaces
	s/^(\s*)//;
#reduce non-leading spaces to one space
	s/(\s+)/ /;
    }
}
###########################################################################
sub parse_test_file($$$) {
    my $name   = $_[0];
    my $tstdir = $_[1];
    my %files = %{$_[$#_]};
    my ($log, $err, $dal, $mol, $ref, $pot) = ($_[$#_]{log}, 
					 $_[$#_]{err}, 
					 $_[$#_]{dal}, 
					 $_[$#_]{mol}, 
					 $_[$#_]{ref},
					 $_[$#_]{pot});
    my @list=();
    my @check_list=();
    my $filename=$tstdir."/".$name.".tst";
    open OUT, $filename or die;
    @list = get_portion("START DALINP","END DALINP",1,1,\%files);
    close OUT;
    foreach $item (@list) {
	print $dal $item;
    }

    open OUT, $filename or die;
    @list = get_portion("START MOLINP","END MOLINP",1,1,\%files);
    close OUT;
    foreach $item (@list) {
	print $mol $item;
    }

    open OUT, $filename or die;
    @list = get_portion("START REFOUT","END REFOUT",1,1,\%files);
    close OUT;
    foreach $item (@list) {
	print $ref $item;
    }

    open OUT, $filename or die;
    @list = get_portion("START POTINP","END POTINP",1,1,\%files);
    close OUT;
    foreach $item (@list) {
	print $pot $item;
    }

    open OUT, $filename or die;
    @list = get_portion("START CHECKLIST","END CHECKLIST",1,1,\%files);
    close OUT;
    chomp_and_shorten_spaces(@list);
    my $idx = 0;
    my @allchecks = ();
    my $check = "";
    foreach $check (@list) {
		$idx++;
		%check = load_check($check, $tstdir, \%files);
		while (($idx <= $#list) and ($list[$idx]=~/OVERRIDE/)) {
			my $override = splice(@list,$idx,1);
			my @override = split(/\s+/,$override);
			if ($#override == 2) {
				$check{$override[1]} = $override[2];
			}
			elsif ($#override > 2) {
				$check{$override[1]} = [splice @override, 2];
			}
			else {
				print $log "override --@override-- has a wrong number of arguments and will not be included\n";
			}
		}
		push @allchecks,{%check}
    }
    close OUT;
    return @allchecks;
}

##################################################################################
sub check_num($)
# check whether the argument is numeric and returns an integer
#  0 non numeric
#  1 integer
#  2 real
{	      
    $_ = shift;
    if (/^\d+$/)                                               { return 1 }
    elsif (/^-?\d+$/)                                          { return 1 }
    elsif (/^[+-]?\d+$/)                                       { return 1 }
    elsif (/^-?\d+\.?\d*$/)                                    { return 2 }
    elsif (/^-?(?:\d+(?:\.\d*)?|\.\d+)$/)                      { return 2 }
    elsif (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/) { return 2 }
    else                                                       { return 0 }
}


############################################################################
sub compare($$$$$$) {
# 
# Compare two numbers up to the required precision
#
    my ($test_value, $ref_value, $thr, $abs, $rel) = ($_[0], $_[1], $_[2], $_[3], $_[4]);
    my %files = %{$_[$#_]};
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});

    my $check = check_num($test_value) * check_num($ref_value) * check_num($thr) * check_num($abs);    
    my $result = 0;
    if ($check == 0) {
	print $err "compare(): non numeric value(s) as input, test skipped\n";;
    } else {
	if ($abs == 1) {
	    $test_value = abs($test_value);
	    $ref_value  = abs($ref_value);
	}
	my $diff = 0;
	if ($rel) {
	    $diff = abs(($test_value - $ref_value)/$ref_value);
	} else {
	    $diff = abs($test_value - $ref_value);
	}
	if($diff <= $thr) {
	    print $log 
		"Numeric comparison successful\n" 
		if $verbose > $verb{check_det};
	    print $log 
		"Result: $test_value  Reference: $ref_value  Required accuracy: $thr  Difference: $diff\n" 
		if $verbose > $verb{all};
	    $result = 1;
	} else {
	    print $err "Numeric comparison failed!\n";
	    print $err "Result: $test_value  Reference: $ref_value  Required accuracy: $thr  Difference: $diff\n";
	}
    }
    return $result;
}
############################################################################
sub line_compare($$$$$$)
{
# 
# Compare two lines after some preprocessing (preprocessing needs to be implemented better)
#
    my ($test_value, $ref_value, $prepro, $n, $m) = ($_[0],$_[1],$_[2],$_[3],$_[4]);
    my %files = %{$_[$#_]};
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});
    my $result = 0;
    
  PREPRO:
    {
	if ($prepro == 0) 
	{
	    last PREPRO;
	}
	if ($prepro == 1) 
	{
	    $test_value=~s/^\s+//;
	    $ref_value=~s/^\s+//;
	    last PREPRO;
	}
	if ($prepro == 2) 
	{
	    while ($test_value=~s/\s//){};
	    while ($ref_value=~s/\s//){};
	    last PREPRO;
	}
	print $err "line_compare: Only preprocessing 1 is implemented.\n";
    }

    if($test_value eq $ref_value)
    {
	print $log "line_compare() successful\n" if $verbose > $verb{check_det};
	print $log "Value: $test_value\n  Reference: $ref_value\n" if $verbose > $verb{all};
	$result = 1;
    }
    else
    {
	print $err "line_compare(): failed comparison\n";
	print $err "input:\n test: $test_value\n  ref: $ref_value\n";
    }
    return $result;
}


###############################################################################
sub check_files
{
#
# Check if all necessary files are readable by the user
# input: list of file names, hash with needed filehandles
#
    my %files = %{pop(@_)};
    my ($log, $err, $biglog) = ($files{log}, $files{err}, $files{biglog});
    my $result = 0;
    foreach $file (@_) {
       next if (-f $file);
       $result++;
	print $err "File $file is not readable!\n";
    }
    if ($result) {
       print $err "One or more file are missing. Test will be skipped\n";
       print $log "One or more file are missing. Test will be skipped\n";
       print $log "Check error file for details\n";
       print $biglog "One or more file are missing. Test will be skipped\n";
       print $biglog "Check error file for details\n";
    }
    return $result;
}
###############################################################################
sub rm_empty_lines
{
    my @result=();
    foreach (@_)
    {
	push @result, $_ unless ($_ eq "");
    }
    return @result;
}
###############################################################################
sub post_coord {
#
# Input geometry check. Compare the coordinates from the test files with the ones obtained from the reference one
# assume it is called with the following input:
# number of elements in the test list  
# number of elements in the ref list  
# the test list 
# the ref list  
#
# assume the first line of each list contains the row "Total number of coordinates: XX" from the DALTON.OUT file
# assume both test and reference contain the coordinates in the same order 
# assume the coordinate value is the last in the line 
#
    my @test = @{$_[0]};
    my @ref  = @{$_[1]};
    my %files = %{$_[$#_]};
    my ($log, $err) = ($files{log},$files{err});
    my $test_n_coor = shift @test;
    my $ref_n_coor  = shift @ref;
    my $threshold = 1.0e-8;
    my $result = 0;
    my $n_right = 0;
    my $n_wrong = 0;
    print $log "post_coord(): input parameters:\n test:\n @test\n  ref:\n @ref\n" if $verbose > $verb{check_det};
#
# 
#
    @test = rm_empty_lines(@test);
    @ref  = rm_empty_lines(@ref);
    if($ncoor = check_n_coor($test_n_coor,$ref_n_coor))
    {
	for (my $idx = 0; $idx < $ncoor; $idx++)
	{
	    my @t_list = split(/\s+/,$test[$idx]);
	    my @r_list = split(/\s+/,$ref[$idx]);
	    if(($#t_list == $#r_list) && $#t_list >= 0)
	    {
		$n_right += compare($t_list[$#t_list],$r_list[$#r_list],$threshold,0,1.0e-5,\%files);
	    }
	    else {
		print $log "post_coord(): empty or not matching lists\n" if $verbose > $verb{check_det};
		print $log "post_coord():\n @t_list\n @r_list\n" if $verbose > $verb{all};
		print $err "post_coord(): empty or not matching lists\n";
		print $err "post_coord():\n @t_list\n @r_list\n";
		$n_wrong++
	    }
	}
    }
    else
    {
	print $log "post_coord(): N. of coordinates do not match\n" if $verbose > $verb{check_det};
	print $log "test:\n @test\n ref:\n @ref\n" if $verbose > $verb{all};
	print $err "post_coord(): N. of coordinates do not match\n";
	print $err "test:\n @test\n ref:\n @ref\n";
	$n_wrong++
    }
    $result = 1 if ($n_right == $ncoor) and ($n_wrong == 0);
    return $result;
}
###############################################################################
sub check_n_coor
{
#
# Compare the number of coordinates in test and ref
# return the number of coordinates if true
# return zero if false
#
    my $result = 0;
    my @t_list = split(/\s+/,$_[0]);
    my @r_list = split(/\s+/,$_[1]);
    my $nt = $t_list[$#t_list];
    my $nr = $r_list[$#r_list];
    $result = $nt if (check_num($nt) && check_num($nr) && ($nt == $nr));
    return $result;
}
#########################################################################################################
sub title_print
{
    my $FH = shift;
    my $string= shift;
    print {$FH} "--------------------------------------------\n".$string.
              "\n--------------------------------------------\n";
    
}
#################################################################################################
# test drivers
#################################################################################################
sub test_line {
#
# fetches the first line of text matching a string, extracts a number and compares it with the reference 
# 
    my ($string, $pos, $thr, $name, $abs, $rel, $num_type) = 
       ($_[0]{string}, $_[0]{pos}, $_[0]{thr}, $_[0]{name}, $_[0]{abs}, $_[0]{rel}, $_[0]{num_type});    
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});
    my %files = %{$_[$#_]};

    my @lll = %{$_[0]};
    print $log "Parameters: @lll\n\n" if $verbose > $verb{test_inp};
    my $result = 0;
     
    open OUT, $files{out_f} or die;
    $tst_line = get_line($string, \%files);
    close OUT;

    open OUT, $files{ref_f} or die;
    $ref_line = get_line($string, \%files);
    close OUT;

    @tst=();
    @ref=();
    @tst = extract_numbers($tst_line, $num_type, \%files) if $tst_line;
    @ref = extract_numbers($ref_line, $num_type, \%files) if $ref_line;

    my $lpos = $pos - 1; #we want to be user friendly!! 
    if (@tst > 0 and @ref > 0) {
	print $log "Checking values.....\n test value: $tst[$lpos]\n  ref value: $ref[$lpos]\n" if $verbose > $verb{check_res};
 	$result = compare($tst[$lpos], $ref[$lpos], $thr, $abs, $rel, \%files);
    } else {
	print $log "ERROR: empty test or reference list, check $name skipped!\n";
	print $err "ERROR: empty test or reference list, check $name skipped!\n";
	print $err "from output file: @tst\n from ref file: @ref\n";
    }
    return $result;
}
############################################################################
sub test_lines {
#
# fetches some/all lines of text matching a string, extracts a number and compares with the reference 
#
    my (     $name,       $string,       $maxout,       $maxget,       $pos,       $thr,       $abs,       $rel) = 
       ($_[0]{name}, $_[0]{string}, $_[0]{maxout}, $_[0]{maxget}, $_[0]{pos}, $_[0]{thr}, $_[0]{abs}, $_[0]{rel});
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});
    my %files = %{$_[$#_]};

    my @lll = %{$_[0]};
    print {$log} "Parameters: @lll\n\n" if $verbose > $verb{test_inp};
    my $result = 0;

    my $matching = 0;

    open OUT, $files{out_f} or die;
    my @test_lines = get_lines($string, $maxout, $maxget, \%files);
    print {$log} "Fetched list from output\n @test_lines\n\n" if $verbose > $verb{fetch_out};
    close OUT;
    open OUT, $files{ref_f} or die;
    my @ref_lines =  get_lines($string, $maxout, $maxget, \%files);
    print {$log} "Fetched list from reference\n @ref_lines\n\n" if $verbose > $verb{fetch_out};
    close OUT;
    
    if (not @test_lines) {
	print $err "Check $name:\n No matching lines found in output file.\n Check aborted\n";
    } elsif (not @ref_lines) {
	print $err "Check $name:\n No matching lines found in reference file.\n Check aborted\n";
    } elsif ($#ref_lines != $#test_lines) {
	print $err "Check $name:\n Line count in output does not match the reference.\n Check aborted\n";
    } else {
	print $log "Checking values.....\n" if $verbose > $verb{check_res};
	print $log "n.   test value    ref value\n" if $verbose > $verb{check_res};
	for ($i=0; $i<= $#test_lines; $i++) {
	    my $test = get_substring($test_lines[$i],$pos,\%files);
	    my $ref  = get_substring($ref_lines[$i], $pos,\%files);
	    print $log ($i+1)." $test $ref\n" if $verbose > $verb{check_res};
	    $matching = $matching + compare($test, $ref, $thr, $abs, $rel, \%files);
	}
	$result = 1 if $matching == $#test_lines + 1;
    }
    return $result;
}
#####################################################################
sub test_portion ($$) {
#
# fetches a portion of text, go to a specific line, extract a number and compares it to the reference
# fore more info about the input see test type 1
#
    my ($start_s,    $start_offset, $end_s, 
	$end_offset, $line_pos,     $name, 
	$abs,        $rel,          $pos, 
	$thr) = 
      ($_[0]{start_string}, $_[0]{start_offset}, $_[0]{end_string}, 
       $_[0]{end_offset},   $_[0]{line},         $_[0]{name}, 
       $_[0]{abs},          $_[0]{rel},          $_[0]{pos}, 
       $_[0]{thr});          
    my %files = %{$_[$#_]};
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});

    my @lll = %{$_[0]};
    print {$log} "Parameters: @lll\n\n" if $verbose > $verb{test_inp};
    my $result = 0;

    open OUT, $files{out_f} or die "cannot open dalton output file $files{out_f}";
    my @test_lines = get_portion($start_s, $end_s, $start_offset, $end_offset, \%files);
    close OUT;
    print {$log} "Fetched list from test output\n @test_lines\n\n" if $verbose > $verb{fetch_out};

    open OUT, $files{ref_f} or die "cannot open the reference output file $files{ref_f}.";
    my @ref_lines = get_portion($start_s, $end_s, $start_offset, $end_offset, \%files);
    close OUT;
    print {$log} "Fetched list from reference output\n @ref_lines\n\n" if $verbose > $verb{fetch_out};

    chomp_and_shorten_spaces(@test_lines,@ref_lines);

    if (not @test_lines) {
	print $err "Check $name:\n".
	    "No matching lines found in output file.\n".
	    "Check aborted\n";
    } elsif (not @ref_lines) {
	print $err "Check $name:\n".
	    "No matching lines found in reference file.\n".
	    "Check aborted\n";
    } elsif ($line_pos > $#test_lines + 1) {
	print $err "Check $name:\n".
	    "Not enough lines fetched from output file.\n".
	    "Check aborted\n";
    } elsif ($line_pos > $#ref_lines + 1) {
	print $err "Check $name:\n".
	    "Not enough lines fetched from reference file.\n".
	    "Check aborted\n";
    } else {
	print $log "Checking values.....\n" if $verbose > $verb{check_res};
	my $test = $test_lines[$line_pos-1];
	my $ref  = $ref_lines[$line_pos-1];	
	print $log "test line: $test\n ref line: $ref\n" if $verbose > $verb{check_res};
	my $t_val = get_substring($test,$pos,\%files);
	my $r_val = get_substring($ref, $pos,\%files);
	print $log "test value: $t_val\n ref value: $r_val\n" if $verbose > $verb{check_res};
	$result = compare($t_val,$r_val,$thr, $abs, $rel, \%files);
    }
    return $result;
}
############################################################################
sub test_lines_mult {
    my ($name, $string, $maxout, $maxget, $num_type) = 
       ($_[0]{name}, $_[0]{string}, $_[0]{maxout}, $_[0]{maxget}, $_[0]{num_type});
    my @pos = @{$_[0]{pos}};
    my @abs = @{$_[0]{abs}};
    my @thr = @{$_[0]{thr}};
    my @rel = @{$_[0]{rel}};
    
    my %files = %{$_[$#_]};
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});
    my @lll = %{$_[0]};
    print $log "Parameters: @lll\n\n" if $verbose > $verb{test_inp};
    my $result = 0;

    my $matching = 0;

    while ($#thr < $#pos) {
	push @thr, $thr[$#thr];
    }
    while ($#abs < $#pos) {
	push @abs, $abs[$#abs];
    }
    while ($#rel < $#pos) {
	push @rel, $rel[$#rel];
    }
    open OUT, $files{out_f};
    my @tst_lines = get_lines($string, $maxout, $maxget, \%files);
    print $log "Fetched list from output\n @tst_lines\n\n" if $verbose > $verb{fetch_out};
    close OUT;
    open OUT, $files{ref_f};
    my @ref_lines = get_lines($string, $maxout, $maxget, \%files); 
    print $log "Fetched list from reference\n @ref_lines\n\n" if $verbose > $verb{fetch_out};
    close OUT;

    if (not @tst_lines) {
	print $err "Check $name:\n".
	    "No matching lines found in output file.\n".
	    "Check aborted\n";
    } elsif (not @ref_lines) {
	print $err "Check $name:\n".
	    "No matching lines found in reference file.\n".
	    "Check aborted\n";
    }
    elsif ($#ref_lines != $#tst_lines) {
	print $err 
             "Number of lines in test file does not match the reference file 
              Test will be skipped: check results manually\n";
    }
    else {
	print $log "Checking values.....\n" if $verbose > $verb{check_res};
	for ($i = 0; $i <= $#tst_lines; $i++) {
	    my @tst = extract_numbers($tst_lines[$i], $num_type, \%files);
	    my @ref = extract_numbers($ref_lines[$i], $num_type, \%files);
	    for ($j=0; $j <= $#pos; $j++) {
		$tst_val = $tst[$pos[$j]-1];
		$ref_val = $ref[$pos[$j]-1];
		print $log "test value: $tst_val\n ref value: $ref_val\n" if $verbose > $verb{check_res};
		$matching = $matching + compare($tst_val, $ref_val, $thr[$j], $abs[$j], $rel[$j], \%files);
	    }
	}
    }
    $result = 1 if $matching == (@tst_lines) * (@pos);
    return $result;
}
############################################################################
sub test_line_fetch { 
#
# fetches a line, process it according to the defined options, then compares it with the reference
#
    my (     $name,       $string,       $maxout,       $type,       $m,            $n      ) = 
       ($_[0]{name}, $_[0]{string}, $_[0]{maxout}, $_[0]{post}, $_[0]{rembeg}, $_[0]{remend});
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});
    my %files = %{$_[$#_]};

    my @lll = %{$_[0]};
    print {$log} "Parameters: @lll\n\n" if $verbose > $verb{test_inp};
    my $result = 0;

    open OUT, $files{out_f} or die;
    my @test = get_line_old($string, $maxout, -1, \%files);
    print {$log} "Fetched list from output\n @test\n\n" if $verbose > $verb{fetch_out};
    close OUT;

    open OUT, $files{ref_f} or die;
    my @ref = get_line_old($string, $maxout, -1, \%files);
    print {$log} "Fetched list from reference\n @ref\n\n" if $verbose > $verb{fetch_out};
    close OUT;

    if (not @test) {
	 print $err "Warning: string $string not found ".
            "in testfile for testcase $name\n";
    } elsif (not @ref) {
	 print $err "Warning: string $string not found ".
            "in ref.file for testcase $name\n";
    } else {
	print $log "Checking values.....\n" if $verbose > $verb{check_res};
	$result = line_compare($test[0],$ref[0],$type,$m,$n,\%files);
	print $log "test value: $test[0]\n ref value: $ref[0]\n" if $verbose > $verb{check_res};
    }
    return $result;
}
############################################################################
sub test_portion_post {
#
# fetches a portion of text, then performs a specific test on it
#
    my ($name,            $start_s,                  $end_s,            $start_offset,      $end_offset,       $post) =     
       ($_[0]{name}, $_[0]{start_string}, $_[0]{end_string}, $_[0]{start_offset}, $_[0]{end_offset},      $_[0]{post});
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});
    my %files = %{$_[$#_]};

    my @lll = %{$_[0]};
    print {$log} "Parameters: @lll\n\n" if $verbose > $verb{test_inp};
    my $result = 0;

    open OUT, $files{out_f} or die;
    my @test_lines = get_portion($start_s, $end_s, $start_offset, $end_offset, \%files);
    print {$log} "Fetched list from output\n @test_lines\n\n" if $verbose > $verb{fetch_out};
    close OUT;

    open OUT, $files{ref_f} or die;
    my @ref_lines = get_portion($start_s, $end_s, $start_offset, $end_offset,\%files);
    print {$log} "Fetched list from reference\n @ref_lines\n\n" if $verbose > $verb{fetch_out};
    close OUT;
    chomp_and_shorten_spaces(@test_lines,@ref_lines);

  POSTPRO:
    {
	print $log "...sending lists to post processing n. $post...\n" if $verbose > $verb{check_res};
	if ($post == 0) 
	{
	    $result = post_coord(\@test_lines,\@ref_lines,\%files);
	    last POSTPRO;
	}
	print $err "test_portion_post: Only postprocessing 0 is implemented. Check $name skipped\n";
    }
    return $result;
}
############################################################################
# text parsing routines
############################################################################
sub get_portion ($$$$$) {
#
# search for $start_key and $end_key and returns the portion of the file 
# which is included there. The portion can be shrinked by $start_offset and $end_offset
# $end_offset can be negative (portion enlarged) whereas $start_offset can be zero or positive.
# a negative $start_offset will be considered as zero
#
    my ($start_key, $end_key, $start_offset, $end_offset) = ($_[0],$_[1],$_[2],$_[3]);
    my %files = %{$_[$#_]};
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});
    print {$log} "get_portion() called with parameters:\n" if $verbose > $verb{fetch_inp};
    print {$log} "start_key: $start_key, end_key: $end_key, start_offset: $start_offset, end_offset: $end_offset:\n" if $verbose > $verb{fetch_inp};
    
    my @lines =();
    my $get_more = 0;
    $start_offset = 0 if ($start_offset < 0);
    $get_more = - $end_offset if ($end_offset < 0);
    
    my $start_found = 0;
    my $end_found = 0;
    
    while(<OUT>) {
	if(/$start_key/) {
	    $lines[$#lines+1] = $_;
	    print $log "get_portion(): found start string\n" if $verbose > $verb{all};
	    last;
	}
    }
    while(<OUT>) {
	$lines[$#lines+1] = $_;
	if(($end_found == 0) && /$end_key/) {
	    $end_found = 1;
	    print $log "get_portion(): found end string\n" if $verbose > $verb{all};
	}
	$get_more--  if ($end_found == 1);
	last         if ($get_more < 0);
    }
    print $log "get_portion(): lines fetched from file:\n@lines\n" if $verbose > $verb{all};
    
    my @garbage = ();
    if(($start_offset > 0) and ($#lines >= $start_offset)) {
	@garbage = splice (@lines, 0, $start_offset);
	print $log "get_portion(): lines removed at the beginning:\n@garbage" if $verbose > $verb{all};
    } elsif ($start_offset > 0) {
	print $log "get_portion(): WARNING: no lines removed at the beginning.\n";
    }
    if (($end_offset > 0) and ($#lines >= $end_offset)) {
	@garbage = splice(@lines,-$end_offset) if($end_offset > 0);
	print {$log} "get_portion(): lines removed at the end:\n@garbage" if $verbose > $verb{all};
    } elsif ($end_offset > 0) {
	print $log "get_portion(): WARNING: no lines removed at the end.\n";
    }
    return @lines;
}
############################################################################
sub get_line_old($$$$) {
#
# search for $keyword in input and returns a list composed by the line
# itself followed by the splitted line (spaces removed). Returns an
# empty list if $keyword is not found. 
#
    my ($keyword, $position, $maxread) = ($_[0],$_[1],$_[2]);
    my %files = %{$_[$#_]};
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});
    print {$log} "get_line_old() called wit the following input:\n $keyword, $position, $maxread\n" if $verbose > $verb{fetch_inp};
    while(<OUT>) {
	$maxread--;
	if(/$keyword/) {
	    s/^\s+//; # remove leading spaces	    
	    @line=split(/\s+/,$_);
	    if ($position > 0)  {
		push @line, $line[$position-1];
	    }
	    print $log "get_line_old() successful.\n"  if $verbose > $verb{all};
	    print $log "returning the following list:\n $_ @line"  if $verbose > $verb{all};
	    return $_,@line;
	}
	last if ($maxread == 0);
    }
    print $err "get_line_old() could not find a matching line.\n";
    return ();
}
############################################################################
sub get_line($$) {
#
# search for $keyword in input and returns the first matching line in the output
# the file is parsed up to maxread lines (whole file if maxread is negative)
#
    my ($keyword, $maxread) = ($_[0], $_[1]);
    my %files = %{$_[$#_]};
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});
    print {$log} "get_line() called wit the following input:\n $keyword, $maxread\n" if $verbose > $verb{fetch_inp};
    $keyword = lc $keyword;
    while(<OUT>) {
	$_ = lc;
	$maxread--;
	if(/$keyword/) {
	    print $log "get_line_old() successful.\n"  if $verbose > $verb{all};
	    print $log "returning the following line:\n $_"  if $verbose > $verb{all};
	    return $_;
	}
	last if ($maxread == 0);
    }
    print $err "get_line() could not find a matching line.\n";
    return ();
}
##############################################################################
sub get_substring ($$$) {
    my @line=();
    my ($line,$position) = ($_[0], $_[1]);
    my %files = %{$_[$#_]};
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});
    print {$log} "get_substring(), parameters:\n line $line, position $position.\n" if $verbose > $verb{fetch_inp};
    my $result = '';
    $line=~s/^\s+//;
    @line=split(/\s+/,$line);
    if ($position > $#line + 1)
    {
	print $err "get_substring(), position $position not found in the following line\n"."$line\n";
	return $result;
    } else {
	$result = $line[$position-1];
	print $log "get_substring(), substring found: $result\n" if $ verbose > $verb{all};
	return $result;
    }
}
############################################################################
sub make_dirname() {
#
# create a unique directory name where files are placed
# based on date, time and process ID
#
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    $mon++; 
    $year = $year + 1900;
    my $name = "perl-pid.".$$."__".$year."_".$mon."_".$mday."__".$hour.".".$min;
    return $name;
}
############################################################################
sub get_lines ($$$$) {
#
# search for $keyword in input and returns a list of all the matching lines
# the input is parsed for a maximum $maxread lines and up to $maxget matches of $keyword
# negative values of $maxget and $maxread mean no corresponding limits
#
    

    my @lines   = ();
    my $keyword = $_[0];
    my $maxread = $_[1];
    my $maxget  = $_[2];
    my %files = %{$_[$#_]};
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});
    print {$log} "get_lines() called with the following input:\n $keyword, $maxread, $maxget\n" if $verbose > $verb{fetch_inp};
    while(<OUT>) {
	if(/$keyword/) {
	    s/^\s+//;
	    push @lines, $_;
	    $maxget--;
	    print {$log} "get_lines(): matching line found:\n$_\n" if $verbose > $verb{all};
	}
	$maxread--;
	last if ($maxread == 0 || $maxget == 0);
    }
    return @lines;
}
#########################################################################################################
sub make_testlist($$$$$$) {
#
# Prepares the list of tests to be performed
#
    my $name ="";
    my $file ="";
    my $all = $_[0];
    my @checks   = split(/\s+/,$_[1]);
    my @keywords = split(/\s+/,$_[2]);
    my @explist  = split(/\s+/,$_[3]);
    my @tstdir   = $_[4];
    my %files = %{$_[$#_]};
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});
    my %list = ();
    my @list = ();
    chomp (@allfiles = `ls $tstdir/*.tst`);
    print $log "available tests are:\n@allfiles\n" if $verbose > $verb{all};
#
#add tests asked for explicitly
#
  FILE:foreach $file (@allfiles) {
      $name = $file;
      $name =~ s/\.tst//;
      $name =~ s/$tstdir\///;
      if($all) {
	  $list{$name} = 1;
	  next FILE;
      }
      open OUT, $file or die;
      my @checklist = get_portion("START CHECKLIST","END CHECKLIST",1,1,\%files);
      close OUT;
      open OUT, $file or die;
      my @keylist = get_lines("KEYWORDS",-1,5,\%files);
      close OUT;
      foreach $_ (@checklist) {
	  foreach $check (@checks) {
	      if (/$check/) {
		  $list{$name} = 1;
		  next FILE;
	      }
	  }
      }
      foreach $_ (@keylist) {
	  foreach $key (@keywords) {
	      if (/$key/) {
		  $list{$name} = 1;
		  next FILE;
	      }
	  }
      }
  }
    foreach $name (@explist) {
	$list{$name} = 1;
    }
    @list = keys(%list);
    print $log "The following tests will be performed:\n";
    print $log "@list\n\n";
    return @list;
}
###########################################################################################################
sub preprocess_files($) {
#
# we need to have "clean" output and reference files before performing any checks
# things to be cleaned:
# 1) all real defined with D/d need to be changed in E 
#    since perl does not recognize xxxx.xxD+xxx and similar
#    as numbers
#
    my %files = %{$_[$#_]};
    my $file="";
    my $array=();
    foreach $file ($files{out_f},$files{ref_f})
    {    
	tie @array, 'Tie::File', $file or die "Cannot tie file $file";
	for (@array) {
	    s/([.\d])[DEd]([+-]\d)/$1e$2/g;
	}
	untie @array;
    }
}
###########################################################################################################
sub print_hash($$)
{
    my %hash_to_print = %{$_[0]};
    my $filehand = $_[1];
    print $filehand "keyword    value\n";
    foreach my $k (sort keys %hash_to_print) {
	print $filehand "$k => $hash_to_print{$k}\n";
    }
}
##############################################################################
sub extract_numbers($$$) {
    my %files = %{$_[$#_]};
    my ($line, $kind) = ($_[0], $_[1]);
    my ($log, $err) = ($_[$#_]{log}, $_[$#_]{err});
    print {$log} "extract_numbers() called with the following input:\n $line, $kind\n" if $verbose > $verb{fetch_inp};
    chomp $line;
    my $piece = "";
    my @strings  = ();
    my @integers = ();
    my @floats   = ();
    my @numbers  = ();
    foreach $piece (split(/\s+/,$line)) {
	my $check = check_num($piece);
	push @strings,  $piece if ($check == 0);
	push @integers, $piece if ($check == 1);
	push @floats,   $piece if ($check == 2);
	push @numbers,  $piece if ($check == 1 or $check == 2);
	if ($check < 0 or $check > 2) {
	    print $log "Argument $piece was found neither a number nor string.\n Please check manually\n ";
	    print $err "Argument $piece was found neither a number nor string.\n Please check manually\n ";
	}
    }
    if ($verbose >= $verb{fetch_out}) {
	print "line\n".$line."\n has been split in the following way:\n";
	print "Strings:\n".@strings."\n";
	print "Integers:\n".@integers."\n";
	print "Floats:\n".@floats."\n";
	print "Numbers:\n".@numbers."\n";
    }
    return @strings  if ($kind == 0);
    return @integers if ($kind == 1);
    return @floats   if ($kind == 2);
    return @numbers;
}
###########################################################################################################
#
# And finally the main code!!!!!
#
# Style notes
#
# 1) Always pass the reference to %files as the last argument to routines
#    the called routine should then have: my %files = %{$_[$#_]}
#                                         my ($log, $err) = ($_[$#_]{log}, $_[$_#]{err});
#
# 
#
# check test list definitions (to be implemented.....)
#
#check_test_cases();

my @passed_tests=();
my @failed_tests=();

# filehandles variables

my $BIGLOG  = "BIGLOG";
my $BIGERR  = "BIGERR";
my $LOGFILE = "LOGFILE";
my $ERRFILE = "ERRFILE";
my $DALINP  = "DALINP";
my $MOLINP  = "MOLINP";
my $POTINP  = "POTINP";
my $REFOUT  = "REFOUT";
my $DALOUT  = "DALOUT";

# list of global filehandles

my %bigfiles = (log    => $BIGLOG,
		err    => $BIGERR,
		biglog => $BIGLOG,
		bigerr => $BIGERR);

# list of file-specific filehandles

my %files = (log    => $LOGFILE,
	     err    => $ERRFILE,
	     dal    => $DALINP,
	     mol    => $MOLINP,
	     ref    => $REFOUT,
	     pot    => $POTINP,
	     out    => $DALOUT,
	     biglog => $BIGLOG,
	     bigerr => $BIGERR);

my $TOT_ERR = 0;

#
# create and enter test directory
#
my $dir = make_dirname();
mkdir $dir || die "Fatal: Unable to create perl-tests directory!";
open $BIGLOG, ">$biglog" or die;
open $BIGERR, ">$bigerr" or die;
chdir $dir || die "Fatal: Unable to enter perl-tests directory!"; 
my @testlist = make_testlist($all,$check,$keyword,$explist,$tstdir,\%bigfiles);
#
# Main loop
#
my $dalinp  = '';
my $molinp  = '';
my $dalout  = '';
my $refout  = '';
my $logdal  = '';
my $errdal  = '';
my $logfile = '';
my $errfile = '';
my $tgzfile = '';

print $BIGLOG "Perl testsuite started.....\n";
print $BIGLOG "Global logfiles will be stored in $biglog and $bigerr\n";
print $BIGLOG "Test-specific logfiles will be kept only for failed tests\n" unless $keep;
print $BIGLOG "Test-specific logfiles will be kept for all tests\n" if $keep;
unless ($quiet) {
    print "Perl testsuite started.....\n";
    print "Global logfiles will be stored in $biglog and $bigerr\n";
    print "Test-specific logfiles will be kept only for failed tests\n" unless $keep;
    print "Test-specific logfiles will be kept for all tests\n" if $keep;
}

TEST: foreach $test (@testlist) { #CAREFUL: this loop has also a "continue" block

    $dalinp  = $test.".dal";
    $molinp  = $test.".mol";
    $potinp  = $test.".pot";
    $dalout  = $test.".out";
    $refout  = $test.".ref";
    $logdal  = $test.".log";
    $errdal  = $test.".err";
    $logfile = $test.".lgf";
    $errfile = $test.".erf";
    $tgzfile = $test.".tar.gz";

    $files{log_f}    = $logfile;
    $files{err_f}    = $errfile;
    $files{dal_f}    = $dalinp;
    $files{mol_f}    = $molinp;
    $files{pot_f}    = $potinp;
    $files{ref_f}    = $refout;
    $files{out_f}    = $dalout;
    $files{biglog_f} = $biglog;
    $files{bigerr_f} = $bigerr;

    open $LOGFILE,">$logfile" or die;
    open $ERRFILE,">$errfile" or die;
    open $DALINP, ">$dalinp"  or die;
    open $MOLINP, ">$molinp"  or die;
    open $POTINP, ">$potinp"  or die;
    open $REFOUT, ">$refout " or die;

    print $BIGLOG "Now performing test $test ..............\n";
    print "Now performing test $test ..............\n" unless $quiet;
    title_print($LOGFILE,"Log information for test $test");
    title_print($ERRFILE,"Error information for test $test");
    $TOT_ERR = 0;
    if(-r "$tstdir/$test.tst") {
	@performed_checks = parse_test_file($test,$tstdir,\%files);
    } else {
	$TOT_ERR++;
	print $BIGLOG "File $test.tst does not exist or is not readable: test skipped!\n";
	print "File $test.tst does not exist or is not readable: test skipped!\n" unless $quiet;
	next TEST;
    }
#
# Test whether input files are corrupted
#
    if (check_files($dalinp, $molinp, $refout, \%files)) {
	title_print($BIGLOG, "Corrupted input files. Test $test will be skipped\n");
	title_print($LOGFILE,"Corrupted input files. Test $test will be skipped\n");
	title_print($ERRFILE,"Corrupted input files. Test $test will be skipped\n");
	$TOT_ERR++;
	next TEST;
    }
#
# run the test or copy the reference output
#
    if ($checkref) {
	system("cp $refout $dalout");
    } elsif ($oldout) {
	system("cp $tstdir/$dalout $dalout");
    } else { 
	system("$dalton $options $test 1> $logdal 2> $errdal");
    }
    if (check_files($dalout,\%files)) {
	title_print($BIGLOG ,"Corrupted output file. Test $test will be skipped\n");
	title_print($LOGFILE,"Corrupted output file. Test $test will be skipped\n");
	title_print($ERRFILE,"Corrupted output file. Test $test will be skipped\n");
	$TOT_ERR++;
	next TEST;
    }
    close $DALINP  or die;
    close $MOLINP  or die;
    close $REFOUT  or die;
#
# Preprocessing of output and reference files
# in order to be sure we have "clean" files
#
    preprocess_files(\%files);
#
#
# Main check loop
#
    my $ok = 0;
    my $idx = 0;
    my $n_test = $#performed_checks + 1;

  CHECK:while ($idx < $n_test) { #CAREFUL: this loop has also a "continue" block
      load_def_check_param($performed_checks[$idx],$tstdir,\%files);
      $name="UNDEF";
      $type=-1;
      $name = $performed_checks[$idx]{name} if defined $performed_checks[$idx]{name};
      $type = $performed_checks[$idx]{type} if defined $performed_checks[$idx]{type};
      unless (control_check_list($performed_checks[$idx],\%files)) {
	  print $LOGFILE "\n----------Check: $name-------------\n";
	  print $LOGFILE "Error in input check. Check number $idx will be skipped\n";
	  print $LOGFILE "Check name not yet defined, increase verbosity and check error files for more details\n";
	  print $BIGLOG  "Error in check input. Check number $idx will be skipped\n";
	  print $BIGLOG  "Check name not yet defined, increase verbosity and check error files for more details\n";
	  print $ERRFILE "Error in input check!\n";
	  print $ERRFILE "Increase verbosity for more details\n" if $verbose < $verb{check_res};
	  next CHECK;
	  
      }	  
      print $LOGFILE "\n----------Check: $name-------------\n";
      goto ("line","lines","portion","lines_mult","line_fetch","portion_post")[$type];
      
    line:         $ok =         test_line($performed_checks[$idx],\%files); next CHECK;
    lines:        $ok =        test_lines($performed_checks[$idx],\%files); next CHECK;
    portion:      $ok =      test_portion($performed_checks[$idx],\%files); next CHECK;
    lines_mult:   $ok =   test_lines_mult($performed_checks[$idx],\%files); next CHECK;
    line_fetch:   $ok =   test_line_fetch($performed_checks[$idx],\%files); next CHECK;
    portion_post: $ok = test_portion_post($performed_checks[$idx],\%files); next CHECK;
  }
#
# now checking for errors.....
#
    continue {
	if (not $ok) {
	    print "Error in check $name! See $test.lgf and $test.erf for details\n" unless $quiet;
	    print $BIGLOG  "Error in check $name! See $test.lgf and $test.erf for details\n";
	    print $LOGFILE "Error in check $name!\n----------------------------------------------\n\n";
	    print $ERRFILE "Error in check $name!\n----------------------------------------------\n\n";
	    $TOT_ERR++;
	} else {
	    print $LOGFILE "Check $name performed successfully.\n----------------------------------------------\n\n";
	}
	$idx++;
	$ok = 0;
    }
}
continue {
    close $LOGFILE or die;
    close $ERRFILE or die;
    if ($TOT_ERR > 0) {
	push @failed_tests, $test;
    } else {
	push @passed_tests, $test;
	foreach my $file ($dalinp, $molinp, $potinp, $dalout, $refout, 
			  $logdal, $errdal, $logfile, $errfile, $tgzfile) {
	    system("rm $file") if ((-f $file) and (not $keep));
	}
    }
}

if ((@passed_tests) or (@failed_tests)) {
    print "Tests finished!\n\n" unless $quiet;
    print $BIGLOG "Tests finished!\n\n";
    if (@passed_tests) {
	print "The following perl tests have been computed successfully:\n" unless $quiet;
	print "@passed_tests\n\n" unless $quiet;
	print $BIGLOG "The following perl tests have been computed successfully:\n";
	print $BIGLOG "@passed_tests\n\n";
    }
    if (@failed_tests) {
	print "The following perl tests have failed:\n" unless $quiet;
	print "@failed_tests\n" unless $quiet;
	print $BIGLOG "The following perl tests have failed:\n";
	print $BIGLOG "@failed_tests\n";
	print $BIGLOG "THERE IS A PROBLEM\n\n";
	print $BIGERR "The following perl tests have failed:\n";
	print $BIGERR "@failed_tests\n";
	print $BIGERR "THERE IS A PROBLEM\n\n";
        exit 1;
    } else {
	print "All tests completed successfully\n" unless $quiet;
	print $BIGLOG "All tests completed successfully\n";
    }
} else {
    print "WARNING! No tests have been performed!\n" unless $quiet;
    print "use test.pl --help for more info on how to choose the perl tests\n" unless $quiet;
    print $BIGLOG "WARNING! No tests have been performed!\n";
    print $BIGLOG "use test.pl --help for more info on how to choose the perl tests\n";
}

my $exit_status = 0;
if (($#passed_tests != $#testlist) or ( $#failed_tests >= 0 )) {
    my $exit_status = 1;
}
exit $exit_status;

##############################################################################
#
# TO DO
#
# 1) a better mechanism for parsing through the test and reference 
#    files: e.g. use Tie::File
# 2) check whether a hash element which is supposed to be a list is correctly
#    initialized as a list and viceversa. Currently we siliently check the
#    elemnts of a list if a list is given but there is no check to determine
#    if that is correct.
# 3) Better string search mechanisnm
#

__END__

=head1 Dalton 2.0 test suite

=head1 SYNOPSIS

test.pl [options] 

Options:
  
    keyword  value     explanation  
 
   --dalton  'string'  dalton script
   --options 'string'  option passed to the dalton script
   --keyword 'wrd1 wrd2 ..'  all tests with wrdX as a keyword will be performed
   --check   'wrd1 wrd2 ..'  all tests with wrdX as a check   will be performed
   --list    'tst1 tst2 ..'  all tests of the list will be performed
   --verbose [0-5]     verbosity level: 0 (default) 5 (print everything)
   --all               perform all tests (*.tst)
   --checkref          check reference files (consistency check)
   --quiet             output on STDOUT is suppressed
   --keep              keep all log files
   --oldout            checks a previously generated output
   --help              print detailed instructions

=head1 OPTIONS

=over 8

=item B<--dalton> (default `pwd`/../../bin/dalton)

Shell script to run Dalton (absolute path)

=item B<--options "dalton options">

Options passed to the dalton shell script

=item B<--log file>

Log file of the test script

=item B<--err file>

Error file of the test script

=item B<--keyword "wrd1 wrd2 ...">

The list of executed tests is made based on the given words.
The tests containing one of such words in the KEYWORDS line will be executed

Example:  --keyword "sym HF"

all tests having either "sym" or "HF" among the KEYWORDS will be executed.
NOTE: this filter is a simple REGEXP match therefore it is a bit sloppy now. Needs to be improved

Options passed to the dalton shell script (see dalton script for more info, default none)

=item B<--check "wrd1 wrd2 ..">

The list of executed tests is made based on the given words.
The tests containing one of the words in the list of checks will be executed

Example:  --check "dipole enehf"

all tests having either "dipole" or "enehf" among the list of checks will be executed 
NOTE: this filter is a simple REGEXP match therefore it is a bit sloppy now. Needs to be improved

=item B<--list "tst1 tst2 ..">

all tests given in the list are executed. Names are without the .tst extension

Example: --list "lf_dipole tpa_pcm_sym"

=item B<--verbose>

Set the verbosity level from 0 to 5
  0: standard output (default)
  0: some additional info from the main test routines
  1: print input of the main test routines
  2: print input of inner routines
  3: print details of performed checks
  4: print (almost) everything 

=item B<--all>

Perform all tests in the current directory based on `ls *.tst`

=item B<--checkref>

no real test is performed the script is run using the reference file twice. 
Useful for debugging and checking consistency of new tests

=item B<--quiet>

The output to STDOUT is suppressed. Used within the TEST script since there
the test.log file is attached to TESTLOG.

=item B<--keep>

All output files are kept. The default behavior is to keep the output files 
from failed tests only.

=item B<--oldout>

Checks previously generated outputs. Output files should be present in the 
same directory as the test script and have the same name as the .tst file
Mostly useful to debug the test script.

=item B<--help>

Print this help message to standard output and exits.

=head1 AUTHOR

Luca Frediani <luca@chem.uit.no>

=head1 COPYRIGHT

This programs is Copyright 2006, Luca Frediani

=cut

