use PangenomeComparison::PangenomeComparisonImpl;

use PangenomeComparison::PangenomeComparisonServer;
use Plack::Middleware::CrossOrigin;



my @dispatch;

{
    my $obj = PangenomeComparison::PangenomeComparisonImpl->new;
    push(@dispatch, 'PangenomeComparison' => $obj);
}


my $server = PangenomeComparison::PangenomeComparisonServer->new(instance_dispatch => { @dispatch },
				allow_get => 0,
			       );

my $handler = sub { $server->handle_input(@_) };

$handler = Plack::Middleware::CrossOrigin->wrap( $handler, origins => "*", headers => "*");
