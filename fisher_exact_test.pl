#!/usr/bin/perl -wT

  use strict;
  use POSIX;
  use IO::Handle;
  use DBI;
  use constant {
    DATABASE => 'sybase',
    DEBUG => 0,
    INTEGER_TOP => 2147483647,  # 2^31 - 1
  };
  use vars qw($database $dbh $db_data);
  use lib ("/home/fly8/werner/lib", ".", "..");
      
  {
    my ($in_file, $line, $accession, @value, $a, $b, $c, $d, $p, 
      $out_file, $seq_length);

#    print $0," @ARGV\n";

    $#ARGV==3 ||  die("usage: fisher_exact_test a b c d\n");
    ($a ,$b, $c, $d) = @ARGV;
    $p=&fisher_exact($a, $b, $c, $d);
    print "$p\n";
}
  sub fisher_exact {
    my (@value, $a, $b, $c, $d, $n, $p, @numerator, @denomenator, $myp);


    ($a, $b, $c, $d) = @_;
    $p = 0;
    $n = $a + $b + $c + $d;
#    print "$a, $b, $c, $d\n";
    while ($a >= 0 && $b >= 0 && $c >= 0 && $d >= 0) {

#      ($a, $b, $c, $d) = matrix(@value);
      @numerator = ($a + $b, $c + $d, $a + $c, $b + $d);
#      $a = 1 if ($a == 0);		# 0!
      @denomenator = ($n, $a, $b, $c, $d);
#      print "@numerator / @denomenator\n";
      $myp = solve_factorial_fraction(\@numerator, \@denomenator);
#      print "$myp\n";
      $p += $myp;
      ++$a;
      --$b;
      --$c;
      ++$d;
#      print "$a, $b, $c, $d\n";
  }
    $p = sprintf("%1.6g", $p); # %f would destroy accuracy for small p
    return $p;
  }


  sub solve_factorial_fraction {
    my( $ra_num, $ra_den ) = @_;                # arrarys for numerator and denominator
    
    my( $dum1, $dum2, $num_f, $den_f, $prime, $count, $num_m, $den_m );
    my( $i, $j, $ans );
    
    my @primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97);
    
    @$ra_num  = sort { $a <=> $b } @$ra_num;
    @$ra_den  = sort { $a <=> $b } @$ra_den;
    
    my @m_num = ();             # whats left over to multiply   (numerator)
    my @m_den = ();             #                               (denomenator)
    
    
    # Step one, cancel out 'factorial tails' and set up lists of
    # multiplicative integers on the numerator and denomenator
    
    while( $dum1 = @$ra_num  and  $dum2 = @$ra_den ) {
      $num_f = pop @$ra_num;
      $den_f = pop @$ra_den;
        
      while( $num_f > $den_f ) {  push @m_num, $num_f;  $num_f--;  }
      while( $den_f > $num_f ) {  push @m_den, $den_f;  $den_f--;  }
    }
    
    # Also include whats left as full factorials
    
    foreach $num_f (@$ra_num) { 
      while( $num_f > 1 ) {  push @m_num, $num_f;  $num_f--;  }
    }
    foreach $den_f (@$ra_den) {
      while( $den_f > 1 ) {  push @m_den, $den_f;  $den_f--;  }
    }

    # Step two is to use the first ten primes to simplify
    # as much as possible via factorisation and cancelation.
    
    foreach $prime (@primes) {
      $count = 0;
      foreach $num_m (@m_num) {
        while( $num_m%$prime == 0 ) { $num_m /= $prime; $count++; }
      }
      foreach $den_m (@m_den) {
        while( $den_m%$prime == 0 ) { $den_m /= $prime; $count--; }
      }
        
      if( $count > 0 ) { push @m_num, $prime**$count;    }
      if( $count < 0 ) { push @m_den, $prime**(-$count); }
    }
    
    
    # Step three is to look for any direct cancelations
    
    @m_num = sort {$b <=> $a} @m_num;
    @m_den = sort {$b <=> $a} @m_den;

    while( scalar(@m_num)>1 && $m_num[-1] == 1 ) { pop @m_num; }    # strip out the ones

    while( scalar(@m_den)>1 && $m_den[-1] == 1 ) { pop @m_den; }
    
    $i = 0;
    $j = 0;
    for( $i=0; $i<@m_num; $i++ ) {
      while( $j<$#m_den  and  $m_num[$i] < $m_den[$j] )  {  $j++;  }
        
      if ( $m_num[$i] == $m_den[$j] ) {
        splice @m_num, $i, 1;
        splice @m_den, $j, 1;
        $i--;
        $j--;
      }
    }
    
    # So, this is all the simplification done:
    
    # Step four is to multiply out groups of integers into
    # larger integers. (assume integers up to INTEGER_TOP)
    
    for( $i=0; $i<$#m_num; $i++ ) {
      while( $i<$#m_num  and  INTEGER_TOP > ($dum1 = $m_num[$i]*$m_num[-1]) ) {
        $m_num[$i] = $dum1;
        pop @m_num;
      }
    }
    for( $i=0; $i<$#m_den; $i++ ) {
      while( $i<$#m_den  and  INTEGER_TOP > ($dum1 = $m_den[$i]*$m_den[-1]) ) {
        $m_den[$i] = $dum1;
        pop @m_den;
      }
    }
    
    if( 0 ) {
      print "denomenator: "; foreach $dum1 (@m_num) { print " $dum1"; }  print "\n";
      print "numerator:   "; foreach $dum2 (@m_den) { print " $dum2"; }  print "\n";
    }
    
    # Now need to calculate the final (float) answer
    
    # Step five is to move to floats and to multiply and divide
    # out to arive at the final answer.
    
    $ans = 1;
#    print "final @m_num / @m_den\n";
    while( $dum1 = @m_num  or  $dum2 = @m_den ) {
      if ($dum1 = @m_num) { $ans *= pop @m_num; }
      if ($dum2 = @m_den) { $ans /= pop @m_den; }
    }    
#    print "$ans\n";
    return( $ans );
  }
