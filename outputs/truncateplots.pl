#!/usr/bin/env perl

use File::Basename qw(fileparse);
use File::Copy qw(move);

opendir my $dir, "." or die "Cannot open directory: $!";
my @files = readdir $dir;
close dir;

foreach $fname (@files)
{
	($name, $path, $suffix) = fileparse($fname, "\.[^.]*"); 
	if ($suffix eq ".pdf")	{
		print "Truncating file $fname\n";
		system("pdfcrop $name.pdf out.pdf > /dev/null 2>&1") == 0 or die "call to pdfcrop failed on $fname";
		move "out.pdf", "$name.pdf";		
	}
}
