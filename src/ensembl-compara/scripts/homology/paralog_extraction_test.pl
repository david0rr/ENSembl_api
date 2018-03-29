use strict;
use warnings;

use Bio::EnsEMBL::Registry;

# Load the registry automatically
my $reg = "Bio::EnsEMBL::Registry";

#Use mysql server as input
$reg->load_registry_from_url('mysql://anonymous@ensembldb.ensembl.org/91');

# Get the compara gene member adaptor
my $gene_member_adaptor = $reg->get_adaptor( "Multi", "compara", "GeneMember" ) || die "Could not fetch GeneMemberAdaptor";

# Get
my $gene_tree_adaptor = $reg->get_adaptor( "Multi", "compara", "GeneTree" ) || die "Could not fetch GeneTreeAdaptor";

my %species_list = ( "oryctolagus_cuniculus"    => 108,
                    "loxodonta_africana"   => 98 );

my @species_list_array = keys(%species_list);

foreach my $species ( sort keys %species_list ) {
   print "Fetching: $species\n";

   my $all_genes = $gene_member_adaptor->fetch_all_by_GenomeDB( $species_list{$species} );

   my $gene_count = 0;
   foreach my $gene_member ( @{$all_genes} ) {
       my $tree = $gene_tree_adaptor->fetch_default_for_Member($gene_member);
       if ($tree) {
           my $aln = $tree->get_alignment_of_homologues( $gene_member, 'ENSEMBL_PARALOGUES' , \@species_list_array ) || die "Could not get alignments";
           my $simple_align = $aln->get_SimpleAlign( -SEQ_TYPE => ($tree->member_type eq
'protein' ? 'cds' : undef), -APPEND_SP_SHORT_NAME => 1, -ID_TYPE => 'STABLE' );
           my $all_members = $aln->get_all_Members;

           my $number_of_taxa = scalar( @{$all_members} );
           my $dir            = "/Users/bsdo/Desktop/homologs/$species/paralogs/taxa_$number_of_taxa";
           system("mkdir -p $dir");

           my $alignIO = Bio::AlignIO->newFh( -file => ">$dir/$gene_count.fasta", -format => "fasta" );
           print $alignIO $simple_align;
           $gene_count++;
       }
   }
}