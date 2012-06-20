#!/usr/bin/perl

my $template = `cat temp_Merge.py`;

while (<>) {
 my @elements = split(/\s+/, $_);
 my $template_copy = $template;
 $template_copy =~ s/AAAA/$elements[0]/g;
 $template_copy =~ s/BBBB/$elements[1]/g;
 $template_copy =~ s/CCCC/$elements[2]/g;
 $template_copy =~ s/DDDD/$elements[3]/g;
 $filename = "tp_PMerge_MC_$elements[0]_$elements[2].py";
 open PYTHON, ">$filename" or die;
 print PYTHON "$template_copy";
 close PYTHON;
# $tmp = `cmsRun $filename`;
}
