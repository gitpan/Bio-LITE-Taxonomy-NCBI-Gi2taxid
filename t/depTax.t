#!perl -T

use strict;
use warnings;
use Test::More;

eval { require Bio::LITE::Taxonomy };
plan skip_all => "Bio::LITE::Taxonomy not installed" if $@;
open my $in, '<:raw', "t/data/dict.bin" or die $!;
is((ref $in), 'GLOB', "Dictionary open for reading");
my $dict = Bio::LITE::Taxonomy::NCBI::Gi2taxid->new(dict => $in);
isa_ok($dict,"Bio::LITE::Taxonomy::NCBI::Gi2taxid","as filehandle");

my $dict2 = Bio::LITE::Taxonomy::NCBI::Gi2taxid->new(dict => 'dict.bin');
isa_ok($dict2,"Bio::LITE::Taxonomy::NCBI::Gi2taxid","as filename");

is($dict->get_taxid(0),0,"Uninitilized values");
is($dict->get_taxid(5),5,"Initilized values");

done_testing();

