package Bio::LITE::Taxonomy::NCBI::Gi2taxid;

=head1 NAME

Bio::LITE::Taxonomy::NCBI::Gi2taxid - Mappings of NCBI GIs to Taxids fast and with very low memory footprint.
=head1 SYNOPSIS

Creation of a new Taxid to GI dictionary (binary lookup file):

  use Bio::LITE::Taxonomy::NCBI::Gi2taxid qw/new_dict/;

  new_dict (in => "gi_taxid_prot.dmp",
            out => "gi_taxid_prot.bin");

Usage of the dictionary:

  use Bio::LITE::Taxonomy::NCBI::Gi2taxid;

  my $dict = Bio::LITE::Taxonomy::NCBI::Gi2taxid->new(dict=>"dict.in");
  my $taxid = $dict->get_taxid(12553);

=head1 DESCRIPTION

The NCBI site offers a file to map gene and protein sequences (GIs) with their corresponding taxon of origin (Taxids). If you want to use this information inside a Perl script you will find that (given the high amount of sequences available) it is fairly inefficient to store this information in a regular hash. Only for creating such a hash you will need more than 10 GBs of system memory.

This is a very simple module that has been designed to efficiently map NCBI GIs to Taxids with speed as the primary goal. It is designed to be able to process a high number of GIs very fast and with low memory usage. It is even faster than using a SQL database to retrieve the mappings or using a local DBHash.

To achieve this, it uses a binary index that can be created with the function C<new_dict>. This index has to be created one time for each mapping file.

The original mapping files can be downloaded from the NCBI site at the following address: L<ftp://ftp.ncbi.nih.gov/pub/taxonomy/>.

=head1 FUNCTIONS

=head2 C<new_dict>

This function creates a new binary dictionary from the NCBI mapping file. The file should be uncompressed before being passed to the script. The function accepts the following parameters:

=over 4

=item in

This is the uncompressed mapping file from the NCBI. The function accepts a filename or a filehandle

=item out

Optional. Where the binary dictionary is going to be printed. The function accepts a filename or a filehandle (that should be opened with writing permissions). If absent STDOUT will be assumed.

=back

=head1 CONSTRUCTOR

=head2 C<new>

Once the binary dictionary is created it can be used as an object using this constructor. It accepts the following parameters

=over 4

=item dict

This is the binary dictionary obtained with the C<new_dict> function. The name of the file or a filehandle is accepted.

=item save_mem

Optional. Use this option to avoid to load the binary dictionary into memory. This will save almost 1GB of system memory but looking up for Taxids will be ~20% slower. This option of I<off> by default.

=back

=head1 METHODS

=head2 C<get_taxid>

This method receives a GI and returns the corresponding Taxid.

=head1 SEE ALSO

L<DBHash>
L<Bio::DB::Taxonomy>

=head1 AUTHOR

Miguel Pignatelli

Any comments should be addressed to emepyc@gmail.com

=head1 LICENSE

Copyright 2009 Miguel Pignatelli, all rights reserved.

This library is free software; you may redistribute it and/or modify it under the same terms as Perl itself.


=cut

use strict;
use warnings;
use Carp qw/croak/;
use File::Tail;

use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
@ISA = qw(Exporter);

@EXPORT = ();   # Only qualified exports are allowed
@EXPORT_OK = qw(new_dict);
$VERSION = '0.01';

sub new
  {
    my ($class, %args) = @_;
    my %opts;

    $args{'dict'} or croak "Need the dictionary file";
    croak "\nERROR\n$args{dict}: File not found\n" unless -e $args{dict};

    my $save_mem = $args{'save_mem'} || 0;
    my $dictfh;
    if ((UNIVERSAL::isa($args{dict}, 'GLOB')) or (ref \$args{dict} eq 'GLOB')) {
      $dictfh = $args{dict}; # TO DO: Test flags, the filehandle must be readable
    } else {
      croak "\nERROR\nThe file containing the gi <-> taxid correspondences must be converted to binary format. See the documentation of this module for more details." unless (-B $args{dict});
      open $dictfh, "<:raw", $args{dict} or croak "$!: $args{dict}";
    }

    $opts{fh} = $dictfh;
    $opts{save_mem} = $save_mem;
    my $self = bless \%opts;

    $self -> _build_dict();
    return $self;
  }

sub _build_dict
  {
    my ($self) = @_;
    my $data;
    $self->{save_mem} && return;
    my $n = sysread( $self->{fh}, $data, -s $self->{fh}  );
    croak "Can't read input file" unless ($n);
    $self->{dict} = $data;
  }

sub get_taxid
  {
    my ($self, $gi) = @_;
    return $self->_direct_lookup ($gi);
  }

sub _direct_lookup {
  my ($self,$gi) = @_;
  if ($self->{save_mem}){
    my $taxid;
    sysseek ($self->{fh},$gi*4,0);
    sysread($self->{fh},$taxid,4,);
    return (unpack "N",$taxid);
  } else {
    return (unpack "N",substr($self->{dict},$gi*4,4));
  }
}

sub new_dict {
  my (%args) = @_;

  # TO DO -> Multiple inputs should be allowed
  defined $args{in} or croak "No input dictionary file provided";

  my $outfh;
  if (! defined $args{out}) {
    $outfh = \*STDOUT
  } elsif ((UNIVERSAL::isa($args{out}, 'GLOB')) or (ref \$args{out} eq 'GLOB')) {
    $outfh = $args{out};  ## TO DO: Test flags, the filehandle must be writeable
  } else {
    open $outfh, ">", $args{out} or croak $!;
  }

  my $infh;
  if ((UNIVERSAL::isa($args{in}, 'GLOB')) or (ref \$args{in} eq 'GLOB')) {
    $infh = $args{in}; ## TO DO: Test flags, the filehandle must be readable
  } else {
    open $infh, "<", $args{in} or croak $!;
  }

  my $ftail = File::Tail->new(name=>$args{in},tail=>1) or croak "$!: $args{in}";
  my $last_line = $ftail->read();
  croak "$args{in} is empty" unless (defined $last_line);
  my ($last_val) = split /\t/, $last_line;

  my $bin = 0;
  substr($bin,$_*4,4,pack ("N",0)) for (0..$last_val);

  while (<$infh>) {
    my ($key,$val) = split /\t/;
    substr($bin,$key*4,4,pack("N",$val));
  }
  close ($infh);
  print {$outfh} $bin;
  if (defined $args{out} && ref \$args{out} eq "SCALAR") {
    close($outfh)}
}


1;
