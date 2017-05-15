#!/usr/bin/perl
#
# extract the names and the arguments of subroutines and funciones and
# what they call
#
# in this version (20021212) I generate a manual with frames in a directory called manual_html
#

# create de direcorty and open files
system('mkdir manual_html');

open(F,">manual_html/index.html");
open(T,">manual_html/top.html");
open(L,">manual_html/left.html");
open(M,">manual_html/main.html");

# print the needed HTML header on all 3 files
header(T);
header(L);
header(M);

# print the index.htm file
print F <<EOTEST ;
<html>
<head></head>
<frameset rows='10%,*'>
<frame src='top.html' name='top' >
  <frameset cols='20%,*'>
  <frame src='left.html' name='left' >
  <frame src='main.html' name='main' >
  </frameset>
</frameset>
</html>

EOTEST


print T <<EOTEST ;

<b>Andrea Ligabue, Copenaghen, December 2002</b>
<table border=1>
<tr> 
<td><a href='left.html#Files' target='left'>List of files and Soubroutines/Functions</a></td>
<td><a href='left.html#Subroutines' target='left'>Alfabetical list of the Soubroutines/Functions</a></td>
</tr>
</table>
<table border=1>
<tr>
<td><a href='left.html#lettera_A' target='left'>A</a></td>
<td><a href='left.html#lettera_B' target='left'>B</a></td>
<td><a href='left.html#lettera_C' target='left'>C</a></td>
<td><a href='left.html#lettera_D' target='left'>D</a></td>
<td><a href='left.html#lettera_E' target='left'>E</a></td>
<td><a href='left.html#lettera_F' target='left'>F</a></td>
<td><a href='left.html#lettera_G' target='left'>G</a></td>
<td><a href='left.html#lettera_H' target='left'>H</a></td>
<td><a href='left.html#lettera_I' target='left'>I</a></td>
<td><a href='left.html#lettera_J' target='left'>J</a></td>
<td><a href='left.html#lettera_K' target='left'>K</a></td>
<td><a href='left.html#lettera_L' target='left'>L</a></td>
<td><a href='left.html#lettera_M' target='left'>M</a></td>
<td><a href='left.html#lettera_N' target='left'>N</a></td>
<td><a href='left.html#lettera_O' target='left'>O</a></td>
<td><a href='left.html#lettera_P' target='left'>P</a></td>
<td><a href='left.html#lettera_Q' target='left'>Q</a></td>
<td><a href='left.html#lettera_R' target='left'>R</a></td>
<td><a href='left.html#lettera_S' target='left'>S</a></td>
<td><a href='left.html#lettera_T' target='left'>T</a></td>
<td><a href='left.html#lettera_U' target='left'>U</a></td>
<td><a href='left.html#lettera_V' target='left'>V</a></td>
<td><a href='left.html#lettera_W' target='left'>W</a></td>
<td><a href='left.html#lettera_X' target='left'>X</a></td>
<td><a href='left.html#lettera_Y' target='left'>Y</a></td>
<td><a href='left.html#lettera_Z' target='left'>Z</a></td>
</tr>
</table>

EOTEST

# get all files to read
@allfiles = @ARGV;

# open the log file
open(LOG,">struct.log") or die "Non riesco ad aprire il file struct.log";

undef %file;
undef %arg;
undef %line;
undef %call;
undef %called_by;

# loop over all files
foreach $file (@allfiles)
{ open(R,"<$file") or die "Non riesco ad aprire il file $in";


  $line = 0;

  while(<R>)
  { chomp();
    $line++;

    #  print "\n$line $_" ;
    # look for SUBROUTINES or FUNCTIONS ignoring comment lines
    if(/^\s{6,}SUBROUTINE\s+(.*)\((.*)/i or /^\s{6,}FUNCTION\s+(.*)\((.*)/i)
    { my $arg = $2; 
      $name = $1;

      my $trueline = $line;

      #  print "\narg: $arg name $name";
      # all in UPPERCASE
      $name =~ tr/a-z/A-Z/;
      $name =~ s/\s//g;
      $name =~ tr/A-Z_0-9//cd;

      # look if the argument list is more than one line long
      while($arg !~ /\)/)
      { $tmp = <R>;
      #  print "\n$line tmp: $tmp";
        chomp($tmp);
        # skip comment lines
        unless($tmp =~ /^\S/)
        { $arg .= $tmp; }
        $line++;

      }

      $arg =~ s/\)//g;
      $arg =~ s/\s//g;
      $arg =~ s/[^_,[:alnum:]]//g;
      # all in UPPERCASE
      $arg =~ tr/a-z/A-Z/;
 
      # check if the SUB or FUN was already defined elsewhere
      if(defined($file{$name}))
      { print LOG "\nSubroutine $name in file $file was already defined in file $file{$name}";
      }
      $file{$name} = $file;
      $arg{$name} = $arg;
      $line{$name} = $trueline;

      $tmp = <R>;
      $line++;
      while($tmp =~ /^[cC\*]/)
      { $commento{$name} .= $tmp;
        $tmp = <R>;
	$line++;
      } 
    }
    elsif(/^\s{6,}SUBROUTINE\s+(.*)/i or /^\s{6,}FUNCTION\s+(.*)/i)
    { 
      $name = $1;
      my $trueline = $line;

      #  print "\narg: $arg name $name";
      # all in UPPERCASE
      $name =~ tr/a-z/A-Z/;
      $name =~ s/\s//g;
      $name =~ tr/A-Z_0-9//cd;

      # check if the SUB or FUN was already defined elsewhere
      if(defined($file{$name}))
      { print LOG "\nSubroutine $name in file $file was already defined in file $file{$name}";
      }
      $file{$name} = $file;
      $arg{$name} = "";
      $line{$name} = $trueline;

      $tmp = <R>;
      $line++;
      while($tmp =~ /^[cC\*]/)
      { $commento{$name} .= $tmp;
        $tmp = <R>;
	$line++;
      } 
    }
    elsif(/^\s{6,}PROGRAM\s+(.*)/i)
    { 
      $name = $1;
      my $trueline = $line;

      #  print "\narg: $arg name $name";
      # all in UPPERCASE
      $name =~ tr/a-z/A-Z/;
      $name =~ s/\s//g;
      $name =~ tr/A-Z_0-9//cd;

      # check if the SUB or FUN was already defined elsewhere
      if(defined($file{$name}))
      { print LOG "\nSubroutine $name in file $file was already defined in file $file{$name}";
      }
      $file{$name} = $file;
      $arg{$name} = "";
      $line{$name} = $trueline;

      $tmp = <R>;
      $line++;
      while($tmp =~ /^[cC\*]/)
      { $commento{$name} .= $tmp;
        $tmp = <R>;
	$line++;
      } 
    }

    if((/^\s/) and /\s+CALL\s*(.*)/i)
    { my $sub = $1;
    
      $sub =~ s/(.*?)\(.*/$1/;
      $sub =~ s/\s//g;
      $sub =~ tr/a-z/A-Z/;
      $sub =~ tr/A-Z_0-9//cd;
      unless($sub =~ /^$/)
      { $call{$name}{$sub} .= "$line ";
        print "\n$name $sub"; 
        # called by section
	unless($called_by{$sub} =~ / $name /)
	{ $called_by{$sub} .= " $name ";  }
      }
    }
  }
  close(R);
}

foreach $file (@allfiles)
{
  # write the filename
  title(M,$file,1);

  foreach $sub (sort keys %file)
  { #print "\n$sub $file";
    if($file{$sub} =~ /^$file$/)
    { 
      # title
      title(M,"$sub",2);
      print M "line: $line{$sub} file: $file{$sub}";

      #argument list
      print M "<p>\n";
      my @arg = split(/,/,$arg{$sub});
      
      my $i=1;
      foreach $arg (@arg)
      { print M "($i) $arg " ; $i++; }

      
      # convert special characters
      $commento{$sub} =~ s/&/&amp;/g;
      $commento{$sub} =~ s/</&lt;/g;
      $commento{$sub} =~ s/>/&gt;/g;
      $commento{$sub} =~ s/"/&quot;/g;
      $commento{$sub} =~ s/\n/<br>\n/g;

      print M "<p>\n$commento{$sub}<p>\n";

      # called by
      print M "\n<p>\ncalled by:<p>";
      my @called_by = split(" ",$called_by{$sub});
      foreach $csub (sort {$a cmp $b} @called_by)
      { unless($myold eq $csub)
        { print M "<a href='main.html#$csub' target='main'>$csub</a>\n"; 
          $myold = $csub;
        }
      }

      # calling to
      print M "<p>\ncall to:";
      foreach $called_sub (sort  keys %{ $call{$sub} })
      { 
        print M "<br>\n<a href='main.html#$called_sub' target='main'>$called_sub</a> at line $call{$sub}{$called_sub}";
      }
    }
  }

}

# final infomration

title(L,"Files",1);

foreach $file (@allfiles)
{ print L "<a href='main.html#$file' target='main'><h2>$file</h2></a>";
  foreach $sub (sort keys %file)
  {
    if($file{$sub} =~ /^$file$/)
    { print L "<br>\n<a href='main.html#$sub' target='main'>$sub</a>";
     
    }
  }
}

# alfabaetical list of the subrotuines

title(L,"Subroutines",1);
$letold = 'Z';
foreach $sub (sort keys %file)
{ ($lett) = $sub =~ /^(.)/ ;
  unless($lett =~ /$letold/)
  { print L "\n<a name='lettera_$lett'><h3>$lett</h3></a>\n";
    $letold = $lett;
  }
  print L "\n<br><a href='main.html#$sub' target='main'>$sub</a>"; 
}


footer(T);
footer(L);
footer(M);

#########################################################
  
sub header()
{ my $file = $_[0];
  print $file "<html>\n<head></head><body>\n"; }

sub footer()
{ my $file = $_[0];
  print $file "</body></html>"; }

sub title()
{ my $name=$_[1];
  my $h = $_[2];
  my $file = $_[0];
  print $file "<hr><h$h><a name='$name'>$name</a></h$h>\n";
}
