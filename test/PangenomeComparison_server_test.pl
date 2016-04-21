use strict;
use Data::Dumper;
use Test::More;
use Config::Simple;
use Time::HiRes qw(time);
use Bio::KBase::AuthToken;
use Bio::KBase::workspace::Client;
use PangenomeComparison::PangenomeComparisonImpl;

local $| = 1;
my $token = $ENV{'KB_AUTH_TOKEN'};
my $config_file = $ENV{'KB_DEPLOYMENT_CONFIG'};
my $config = new Config::Simple($config_file)->get_block('PangenomeComparison');
my $ws_url = $config->{"workspace-url"};
my $ws_name = undef;
my $ws_client = new Bio::KBase::workspace::Client($ws_url,token => $token);
my $auth_token = Bio::KBase::AuthToken->new(token => $token, ignore_authrc => 1);
my $ctx = LocalCallContext->new($token, $auth_token->user_id);
$PangenomeComparison::PangenomeComparisonServer::CallContext = $ctx;
my $impl = new PangenomeComparison::PangenomeComparisonImpl();

sub get_ws_name {
    if (!defined($ws_name)) {
        my $suffix = int(time * 1000);
        $ws_name = 'test_PangenomeComparison_' . $suffix;
        $ws_client->create_workspace({workspace => $ws_name});
    }
    return $ws_name;
}

eval {
    my $obj_name = "genome.1";
    my $obj_name2 = "genome.2";
    my $contig1 = {id => '1', length => 10, md5 => 'md5', sequence => 'agcttttcat'};
    my $contig2 = {id => '2', length => 5, md5 => 'md5', sequence => 'agctt'};
    my $contig3 = {id => '3', length => 12, md5 => 'md5', sequence => 'agcttttcatgg'};
    my $obj = {contigs => [$contig1,$contig2,$contig3], id => 'id', md5 => 'md5',
            name => 'name', source => 'source', source_id => 'source_id', type => 'type', "scientific_name" => 'Scientific name', "domain" => "Bacteria", "genetic_code" => 1, "features" => [{"id"=>"ftr1","protein_translation"=>"abcdefghijklmnop","type"=>"CDS"}]};
    $ws_client->save_objects({workspace => get_ws_name(), objects =>
				  [{type => 'KBaseGenomes.Genome', name => $obj_name, data => $obj}, {type => 'KBaseGenomes.Genome', name => $obj_name2, data => $obj}]});
    eval { 
	my $ret = $impl->build_pangenome({workspace=>get_ws_name(), output_id=>"pg.1", genome_refs=>[get_ws_name()."/".$obj_name,get_ws_name()."/".$obj_name2]});
    };
    if ($@) {
	print("Error while running build_pangenome on two genomes: $@\n");
    }
    my $mset = {description => "mset", elements => {param0 => {ref => get_ws_name()."/".$obj_name}, param1 => {ref => get_ws_name()."/".$obj_name2}}};
    $ws_client->save_objects({workspace => get_ws_name(), objects =>
				  [{type => 'KBaseSearch.GenomeSet', name => "mset", data => $mset}]});
    eval { 
	my $ret = $impl->build_pangenome({workspace=>get_ws_name(), output_id=>"pg.1", genomeset_ref=>get_ws_name()."/mset"});
    };
    if ($@) {
	print("Error while running build_pangenome on genomeset: $@\n");
    }

    done_testing(0);
};
my $err = undef;
if ($@) {
    $err = $@;
}
eval {
    if (defined($ws_name)) {
        $ws_client->delete_workspace({workspace => $ws_name});
        print("Test workspace was deleted\n");
    }
};
if (defined($err)) {
    if(ref($err) eq "Bio::KBase::Exceptions::KBaseException") {
        die("Error while running tests: " . $err->trace->as_string);
    } else {
        die $err;
    }
}

{
    package LocalCallContext;
    use strict;
    sub new {
        my($class,$token,$user) = @_;
        my $self = {
            token => $token,
            user_id => $user
        };
        return bless $self, $class;
    }
    sub user_id {
        my($self) = @_;
        return $self->{user_id};
    }
    sub token {
        my($self) = @_;
        return $self->{token};
    }
    sub provenance {
        my($self) = @_;
        return [{'service' => 'PangenomeComparison', 'method' => 'please_never_use_it_in_production', 'method_params' => []}];
    }
    sub authenticated {
        return 1;
    }
    sub log_debug {
        my($self,$msg) = @_;
        print STDERR $msg."\n";
    }
    sub log_info {
        my($self,$msg) = @_;
        print STDERR $msg."\n";
    }
}
