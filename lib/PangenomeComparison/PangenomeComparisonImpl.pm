package PangenomeComparison::PangenomeComparisonImpl;
use strict;
use Bio::KBase::Exceptions;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org 
our $VERSION = "0.1.0";

=head1 NAME

PangenomeComparison

=head1 DESCRIPTION

A KBase module: PangenomeComparison
This sample module contains one small method - filter_contigs.

=cut

#BEGIN_HEADER
use Bio::KBase::AuthToken;
use Bio::KBase::workspace::Client;
use Config::IniFiles;
use Data::Dumper;
#END_HEADER

sub new
{
    my($class, @args) = @_;
    my $self = {
    };
    bless $self, $class;
    #BEGIN_CONSTRUCTOR
    
    my $config_file = $ENV{ KB_DEPLOYMENT_CONFIG };
    my $cfg = Config::IniFiles->new(-file=>$config_file);
    my $wsInstance = $cfg->val('PangenomeComparison','workspace-url');
    die "no workspace-url defined" unless $wsInstance;
    
    $self->{'workspace-url'} = $wsInstance;
    
    #END_CONSTRUCTOR

    if ($self->can('_init_instance'))
    {
	$self->_init_instance();
    }
    return $self;
}

=head1 METHODS



=head2 build_pangenome

  $return = $obj->build_pangenome($input)

=over 4

=item Parameter and return types

=begin html

<pre>
$input is a PangenomeComparison.BuildPangenomeParams
$return is a PangenomeComparison.BuildPangenomeResult
BuildPangenomeParams is a reference to a hash where the following keys are defined:
	genome_refs has a value which is a reference to a list where each element is a string
	genomeset_ref has a value which is a string
	workspace has a value which is a string
	output_id has a value which is a string
BuildPangenomeResult is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string
	pg_ref has a value which is a string

</pre>

=end html

=begin text

$input is a PangenomeComparison.BuildPangenomeParams
$return is a PangenomeComparison.BuildPangenomeResult
BuildPangenomeParams is a reference to a hash where the following keys are defined:
	genome_refs has a value which is a reference to a list where each element is a string
	genomeset_ref has a value which is a string
	workspace has a value which is a string
	output_id has a value which is a string
BuildPangenomeResult is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string
	pg_ref has a value which is a string


=end text



=item Description



=back

=cut

sub build_pangenome
{
    my $self = shift;
    my($input) = @_;

    my @_bad_arguments;
    (ref($input) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"input\" (value was \"$input\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to build_pangenome:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'build_pangenome');
    }

    my $ctx = $PangenomeComparison::PangenomeComparisonServer::CallContext;
    my($return);
    #BEGIN build_pangenome
    my $token=$ctx->token;
    my $wsClient=Bio::KBase::workspace::Client->new($self->{'workspace-url'},token=>$token);
    my $provenance = [{}];
    $provenance = $ctx->provenance if defined $ctx->provenance;

    if (!exists $input->{'output_id'}) {
        die "Parameter output_id is not set in input arguments";
    }
    my $id = $input->{'output_id'};
    if (!exists $input->{'workspace'}) {
        die "Parameter workspace is not set in input arguments";
    }
    my $workspace_name=$input->{'workspace'};

    my @genomes;
    if (defined $input->{genomeset_ref}) {
	eval {
	    my $genomeset=$wsClient->get_objects([{ref=>$input->{genomeset_ref}}])->[0]{data};
	    push @{$provenance->[0]->{'input_ws_objects'}}, $input->{genomeset_ref};
	    map { push @genomes, $_->{ref} } values %{$genomeset->{elements}};
	};
    }
    if ($@) {
	die "Error loading genomeset from workspace:\n".$@;
    }
    if (defined $input->{genome_refs}) {
	    foreach my $ref (@{$input->{genome_refs}}) {
		next if ! defined $ref;
		push @genomes, $ref;
		push @{$provenance->[0]->{'input_ws_objects'}}, $ref;
	    }
    }

    my $orthlist = [];
    my $okdb;
    my $pangenome = {
    	id => $id,
		type => "kmer",
		genome_refs => [],
		orthologs => [],
    };
    my $proteins = {};

    print STDERR "Processing ", scalar @genomes, " genomes\n";
    my $i = 0;
    foreach my $currgenome_ref (@genomes) {
    	my $gkdb = {};
    	my $genepairs;
    	my $bestorthos = [];
	my $currgenome = undef;
	eval {
	    print STDERR "Getting object from workspace with ref $currgenome_ref\n";
	    my $obj = $wsClient->get_objects([{ref=>$currgenome_ref}])->[0];
	    $currgenome_ref = $obj->{info}->[6]."/".$obj->{info}->[0]."/".$obj->{info}->[4]; # widget needs this kind of ref
	    $currgenome=$obj->{data};
	    push @{$provenance->[0]->{'input_ws_objects'}}, $currgenome_ref;
	};
	if ($@) {
	    die "Error loading genome from workspace:\n".$@;
	}
	
    	push(@{$pangenome->{genome_refs}},$currgenome_ref);
    	if ($i == 1) {
    		my $array = [split(/\s/,$currgenome->{scientific_name})];
    		$pangenome->{name} = $array->[0]." pangenome";
    	}
    	my $ftrs = $currgenome->{features};
    	for (my $j=0; $j < @{$ftrs}; $j++) {
    		my $feature = $ftrs->[$j];
    		if (defined($feature->{protein_translation})) {
    			$proteins->{$feature->{id}} = $feature->{protein_translation};
    			my $matchortho;
    			my $bestortho;
    			my $bestscore = 0;
    			my $seq = $feature->{protein_translation};
    			for (my $k=	0; $k < (length($seq)-8); $k++) {
    				my $kmer = substr($seq,$k,8);
    				if ($i > 0) {
	    				if (defined($okdb->{$kmer})) {
	    					if (!defined($matchortho->{$okdb->{$kmer}})) {
	    						$matchortho->{$okdb->{$kmer}} = 0;
	    					}
	    					$matchortho->{$okdb->{$kmer}}++;
	    					if ($matchortho->{$okdb->{$kmer}} > $bestscore) {
	    						$bestscore = $matchortho->{$okdb->{$kmer}};
	    						$bestortho = $okdb->{$kmer};
	    					}
	    				}
    				}
    				if (defined($gkdb->{$kmer}) && !defined($gkdb->{$kmer}->{-1})) {
    					if (keys(%{$gkdb->{$kmer}}) >= 5) {
    						my $keylist = [keys(%{$gkdb->{$kmer}})];
    						for (my $m=0; $m < 4; $m++) {
    							for (my $n=($m+1); $n < 5; $n++) {
    								$genepairs->{$keylist->[$m]}->{$keylist->[$n]}--;
    								$genepairs->{$keylist->[$n]}->{$keylist->[$m]}--;
    							}
    						}
    						$gkdb->{$kmer} = {-1 => 0};
    					} else {
    						foreach my $key (keys(%{$gkdb->{$kmer}})) {
    							if ($key ne $j) {
    								if (!defined($genepairs->{$key}->{$j})) {
    									$genepairs->{$key}->{$j} = 0;
    									$genepairs->{$j}->{$key} = 0;
    								}
    								$genepairs->{$key}->{$j}++;
    								$genepairs->{$j}->{$key}++;
    							}
    						}
    						$gkdb->{$kmer}->{$j} = 1;
    					}
    				} else {
    					$gkdb->{$kmer}->{$j} = 1;
    				}
    			}
    			if ($bestscore < 10) {
    				$bestorthos->[$j] = -1;
    			} else {
    				$bestorthos->[$j] = $bestortho;
    				push(@{$pangenome->{orthologs}->[$bestortho]->{orthologs}},[$ftrs->[$j]->{id},0,$currgenome_ref]);
    			}
    		}
    	};
    	foreach my $kmer (keys(%{$gkdb})) {
    		if (!defined($gkdb->{$kmer}->{-1})) {
	    		my $keep = 1;
	    		if (keys(%{$gkdb->{$kmer}}) > 1) {
	    			my $keylist = [keys(%{$gkdb->{$kmer}})];
	    			for (my $m=0; $m < (@{$keylist}-1); $m++) {
	    				for (my $n=($m+1); $n < @{$keylist}; $n++) {
	    					if ($genepairs->{$keylist->[$m]}->{$keylist->[$n]} < 10 && $bestorthos->[$keylist->[$m]] == $bestorthos->[$keylist->[$n]]) {
	    						$keep = 0;
	    						$m = 1000;
	    						last;
	    					};
	    				}
	    			}
	    		}
	    		if ($keep == 1) {
	    			foreach my $gene (keys(%{$gkdb->{$kmer}})) {
	    				if ($bestorthos->[$gene] == -1) {
	    					$bestorthos->[$gene] = @{$pangenome->{orthologs}};
	    					my $list = [[$ftrs->[$gene]->{id},0,$currgenome_ref]];
	    					foreach my $partner (keys(%{$genepairs->{$gene}})) {
	    						if ($genepairs->{$gene}->{$partner} >= 10 && $bestorthos->[$partner] == -1) {
	    							$bestorthos->[$partner] = @{$pangenome->{orthologs}};
	    							push(@{$list},[$ftrs->[$partner]->{id},0,$currgenome_ref]);
	    						}
	    					}
	    					my $seq = $ftrs->[$gene]->{protein_translation};
	    					my $index = @{$pangenome->{orthologs}};
	    					my $neworthofam = {
						    	id => $ftrs->[$gene]->{id},
						    	type => $ftrs->[$gene]->{type},
						    	function => $ftrs->[$gene]->{function},
								protein_translation => $ftrs->[$gene]->{protein_translation},
								orthologs => $list
						    };
						    if (!defined($neworthofam->{function})) {
						    	$neworthofam->{function} = "unknown";
						    }
	    					push(@{$pangenome->{orthologs}},$neworthofam);
	    				}
	    				$okdb->{$kmer} = $bestorthos->[$gene];
	    			}
	    		}
    		}
    	}
	$i++;
    }
    print STDERR "Final score computing!\n";
    foreach my $kmer (keys(%{$okdb})) {
    	my $index = $okdb->{$kmer};
    	my $list = $pangenome->{orthologs}->[$index]->{orthologs};
    	my $hits = [];
    	for (my $i=0; $i < @{$list}; $i++) {
    		if (index($proteins->{$list->[$i]->[0]},$kmer) >= 0) {
    			push(@{$hits},$i);
    		}
    	}
    	my $numhits = @{$hits};
    	my $numorthos = @{$list};
    	if ((2*$numhits) >= $numorthos) {
    		foreach my $item (@{$hits}) {
    			$list->[$item]->[1]++;
    		}
    	}
    }

    my $pg_metadata = $wsClient->save_objects({
	'workspace' => $workspace_name,
	'objects' => [{
	    type => 'KBaseGenomes.Pangenome',
	    name => $id,
	    data => $pangenome
		      }]});

    my $report = "Pangenome saved to $workspace_name/$id\n";
    my $reportObj = { "objects_created"=>[{'ref'=>"$workspace_name/$id", "description"=>"Pangenome"}],
		      "text_message"=>$report };
    my $reportName = "pangenome_report_${id}";

    my $metadata = $wsClient->save_objects({
	'id' => $pg_metadata->[0]->[6],
	'objects' => [{
	    type => 'KBaseReport.Report',
	    data => $reportObj,
	    name => $reportName,
	    'meta' => {},
	    'hidden' => 1,
	    'provenance' => $provenance
		      }]});

    $return = { 'report_name'=>$reportName, 'report_ref', $metadata->[0]->[6]."/".$metadata->[0]->[0]."/".$metadata->[0]->[4], 'pg_ref' => $workspace_name."/".$id};
    #END build_pangenome
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to build_pangenome:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'build_pangenome');
    }
    return($return);
}




=head2 compare_genomes

  $return = $obj->compare_genomes($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a PangenomeComparison.CompareGenomesParams
$return is a PangenomeComparison.CompareGenomesResult
CompareGenomesParams is a reference to a hash where the following keys are defined:
	pangenome_id has a value which is a string
	pangenome_ws has a value which is a string
	protcomp_id has a value which is a string
	protcomp_ws has a value which is a string
	output_id has a value which is a string
	workspace has a value which is a string
CompareGenomesResult is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string
	cg_ref has a value which is a string

</pre>

=end html

=begin text

$params is a PangenomeComparison.CompareGenomesParams
$return is a PangenomeComparison.CompareGenomesResult
CompareGenomesParams is a reference to a hash where the following keys are defined:
	pangenome_id has a value which is a string
	pangenome_ws has a value which is a string
	protcomp_id has a value which is a string
	protcomp_ws has a value which is a string
	output_id has a value which is a string
	workspace has a value which is a string
CompareGenomesResult is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string
	cg_ref has a value which is a string


=end text



=item Description

Compares the specified genomes and computes unique features and core features

=back

=cut

sub compare_genomes
{
    my $self = shift;
    my($params) = @_;

    my @_bad_arguments;
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to compare_genomes:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'compare_genomes');
    }

    my $ctx = $PangenomeComparison::PangenomeComparisonServer::CallContext;
    my($return);
    #BEGIN compare_genomes
    #END compare_genomes
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to compare_genomes:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'compare_genomes');
    }
    return($return);
}




=head2 version 

  $return = $obj->version()

=over 4

=item Parameter and return types

=begin html

<pre>
$return is a string
</pre>

=end html

=begin text

$return is a string

=end text

=item Description

Return the module version. This is a Semantic Versioning number.

=back

=cut

sub version {
    return $VERSION;
}

=head1 TYPES



=head2 BuildPangenomeParams

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
genome_refs has a value which is a reference to a list where each element is a string
genomeset_ref has a value which is a string
workspace has a value which is a string
output_id has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
genome_refs has a value which is a reference to a list where each element is a string
genomeset_ref has a value which is a string
workspace has a value which is a string
output_id has a value which is a string


=end text

=back



=head2 BuildPangenomeResult

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string
pg_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string
pg_ref has a value which is a string


=end text

=back



=head2 CompareGenomesParams

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
pangenome_id has a value which is a string
pangenome_ws has a value which is a string
protcomp_id has a value which is a string
protcomp_ws has a value which is a string
output_id has a value which is a string
workspace has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
pangenome_id has a value which is a string
pangenome_ws has a value which is a string
protcomp_id has a value which is a string
protcomp_ws has a value which is a string
output_id has a value which is a string
workspace has a value which is a string


=end text

=back



=head2 CompareGenomesResult

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string
cg_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string
cg_ref has a value which is a string


=end text

=back



=cut

1;
