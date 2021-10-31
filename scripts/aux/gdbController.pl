#!/usr/bin/perl -w

use strict;
use IPC::Open2;

my $init = "run";
my $command = "bt";
my $delay = 1;

my $need_int=0;

$init = shift @ARGV;
$delay = shift @ARGV;
$command = shift @ARGV;

die("Usage: gdbpriver.pl '' 0.1 'bt' gdb -q /path/to/proc 33344\n\tgdbdriver.pl init_command period_seconds backtrace_command startup_arguments\n") unless $ARGV[0];

my $pid = open2(\*OUT, \*IN, @ARGV);

print "pid=$pid\n";

print IN "set pagination off\n";
print IN "$init\n";

while(<OUT>) {
    if (/Starting program:/) {
    $need_int=1;
    last;
    }
    last if /\(gdb\)/;
}

sub intr() {
    kill 9, $pid;
    exit(0);
}
$SIG{'INT'} = \&intr;
sub spipe() {
    print "PIPE!\n";
}
$SIG{'PIPE'} = \&spipe;

if($need_int) {
# Allow time for start up.
    select undef, undef, undef, 5;
    kill 2, $pid;
}

for(;;) {
    print IN "$command\n"; # backtrace
    print IN "c\n"; # continue the program
    while(<OUT>) {
    last if /Continuing./;
    print;
    }
    select undef, undef, undef, $delay; # sorry, nanosleep fails
    print "INT\n";
    kill 2, $pid; # SIGINT to gdb to make it interrupt the program
}
