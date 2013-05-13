#!/usr/bin/env perl

# Crops (in place) all the pdfs in the current directory
# Requires pdfcrop to be in the path
# Uses bash style output redirection also (there's probably a more portable way of doing this ...)

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
