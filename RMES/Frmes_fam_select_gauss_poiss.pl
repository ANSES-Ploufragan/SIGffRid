#!/usr/bin/perl -w

use strict;
use Env qw(SIGFFRID_DIR);
use lib $SIGFFRID_DIR.'libperl/';
use get_dir;

## PERMET D'OBTENIR UN FICHIER DE SORTIE UTILISABLE A PARTIR D'RMES(RESULTATS CLASSES PAR ORDRE DE MANIERE DECROISSANTE), SANS AVOIR A LANCER RMES NI RMES.FORMAT
#UTILISE L'OPTION -FAM!

## Parametres ##
my $prog_tag = '[Frmes_fam_select_gauss_poiss.pl]';

# Type de fichier rmes utilise
my @rmes = (' --gauss',' --compoundpoisson');

# Test de la validite des parametres
# charu's changes.. 
# the change is that the program will get 3 inputs 1.the name of created r'mes files (and this will be bac-id)
#2. path/name of the formatted input sequence file
# 3 input will be the working directory input by the user in the interface

if (($#ARGV < 2))
{
    die "To use this programme : $0\n<suffix for created r'mes files: can give bac id>\n<name of the sequence file>\n< directory>\n";
}

#---------------- Charu's changes start--------------------#
# Contient le nom du fichier cree par 
my $fic = $ARGV[0]; ## we are using bac id as $fic - the suffix used for naming RMES result files


# Variable indiquant le dossier de destination des resultats final, $dossier= working directory were we the input data and the results will be saved for all SIGffRid programs
my $dossier = $ARGV[2]; 
$dossier =~ s/([^\/]+)\/$/$1/;


# Test si le chemin du dossier entre existe, sinon renvoie d'un message d'erreur 
if(! -d $dossier)
{
    die "$prog_tag [Error] $dossier/ directory does not exist, line ".__LINE__."\n";
}

# Contient l'adresse de la sequence a analyser
 # $sequence= $working_directory."/RMES_".$fic.".txt" == # output after using the extrait sequence script run
my $sequence = $ARGV[1]; 


#~~~~~
# Contiennent les valeurs de hmin et hmax (ici 2 et 8)
my $hmin = 3;
my $hmax = 8;


# Variable indiquant quel fichier de famille de mots est traite (8 fichiers au total, de famC_2 a famC_8)
my $indice;
# Variable contenant l'adresse des fichiers de famille de mots
# my $fam = '/users/adage/touzain/familles_de_mots/famC_';


## i have made a directory called RMES that will be a sub folder in the folder where SIGffRid program is placed
# I did this to make it easier the use of the various files that is stored by R'MES program for future use if there is and to make the files #arrangement neater now before i launch this script i change the cwd to the directory where the RMES program + this script is present. So the familles_de_mot is present in  the cwd ( and cwd is ..../RMES/)
my $dir= get_dir(__FILE__); # Cwd::cwd().'/';
print "current dir:$dir, line ".__LINE__."\n";
my $fam = $dir."familles_de_mots/famC_";
$fam =~ /(.*?)[^\/]+$/; 
(-d $1)or die "$prog_tag [Error] $1 direstory does not exist line ".__LINE__."!\n";

###############################################################################
############################  R'MES : do not need this test   ##################
# ??
# open(RMES,"which $rmes[0] |")or die "We can not lacalize $rmes[0] line ".__LINE__."!\n"; 
# ## Variable indiquant l'adresse du logiciel R'MES
# my $commande_rmes = <RMES>;
# $commande_rmes =~ s/$rmes[0]$//;
# chomp $commande_rmes;
# close RMES;
#  #$commande_rmes = $dossier."/"; # we will give the rmes.guassein and other needed programs
my $commande_rmes = 'rmes';
my $nice;
if($^O eq 'darwin')
{
    # Mac OS X
    $nice = 'nice -n 19 ';
}
else
{
    # PC
    $nice = 'nice +19 '; 
}
######################

my $length_seq = 0;
open(SEQ,"< $sequence")or die "Can not open $sequence file line ".__LINE__."\n";
while(<SEQ>){
    /^>/ and next;
    $length_seq += length($_);
    # if we have several sequences in the same file previous seq have Z letter(s) at the end to show rmes
    # we hev several seq
    $_ =~ /([Zz+])$/;
    defined($1) and $length_seq -= length($1);
}
close SEQ;

my $rmes; my @suffix;
foreach $indice($hmin..$hmax)
{ 
    my $eval_gauss_poiss = $length_seq >> ($indice << 1);

    if($eval_gauss_poiss > 100){ $rmes = $rmes[0]; &run_rmes($indice,'gaussien');  }
    elsif($eval_gauss_poiss < 10){ $rmes = $rmes[1]; &run_rmes($indice,'compoundpoisson')or die "PB with $rmes[1] for eval_gauss_poiss ($eval_gauss_poiss ) < 10 line ".__LINE__."\n";  }
    else{
		# print "we here in else 2\n";
		$rmes = $rmes[1];
		&run_rmes($indice,'compoundpoisson')or do{ $rmes = $rmes[0]; pop @suffix; &run_rmes($indice,'gaussien'); };
    }
}



sub run_rmes(@){
    my ($indice) = shift @_;
    push @suffix, shift @_;
# #   system("nice -n 19  /usr/bin/rmes.gaussien -o rmestest.txt6 -seq ../donnees_brutes/130306/RMES_input_files/RMES_SCO.txt -hmin 3 -hmax 8 -max -fam /Volumes/applications/Applications/rmes-2.1.6/familles_de_mots/famC_6")== 0 or print "Bad!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ;

    ## Lancement du programme $rmes ##
	my $fic_indice_f = "${dossier}/$fic$indice";

#     system("$commande_rmes $rmes -o $fic$indice -s $sequence --fam $fam$indice --max --dna") == 0 #--lmin $hmin --lmax $hmax can only be used with guassian when we dont use famille files
	system("$commande_rmes $rmes -s $sequence -o $fic_indice_f  -l $indice --max --dna  ") == 0 #--lmin $hmin --lmax $hmax can only be used with guassian when we dont use famille files	
		or do{
			print "!!! Le programme rmes n'a pas pu etre lance : $?\ncommande: $nice $commande_rmes$rmes -o $fic_indice_f -s $sequence --max --fam $fam$indice\n";
			return 0;
			};
    
    print "- le fichier $fic_indice_f.0 a ete genere\n";
    
    ## Lancement du programme rmes.format qui met les resultats de $rmes sous forme de tableaux exploitables ##
    my $tableau_f = "${dossier}/$fic$indice\_rmes.$suffix[$#suffix]\_tableau";
    system("rmes.format < $fic_indice_f.0 > $tableau_f") == 0
	or die "!!! Le programme rmes.format n'a pas pu etre lance : $?\n";
    
    print "- le fichier $tableau_f a ete genere\n";

    return 1;
}

## Concatenation des 8 fichiers $fic\_$rmes\_tableau obtenus
my $sortie_rmesfam_f = "${dossier}/sortie_$fic\_rmes-fam";
open (CONCAT2,'>', $sortie_rmesfam_f)
    or die "$prog_tag [Error] $sortie_rmesfam_f file cannot be created:$!, line ".__LINE__."\n";
    
my $cptfile = 0;
foreach $indice ($hmin..$hmax)
{
	my $tableau_f = "${dossier}/$fic$indice\_rmes.$suffix[$cptfile]\_tableau";
    open (CONCAT1,'<', $tableau_f)
		or die "$prog_tag [Error] $tableau_f file cannot be opened:$!, line ".__LINE__."\n";
    # Saut des 8 premieres lignes des fichiers
    foreach my $z(0..8)
    {
		<CONCAT1>;
    }
    # Recopiage des fichiers a la suite dans le fichier sortie_$fic_rmes-fam
    while (my $concat = <CONCAT1>)
    {
        print CONCAT2  "$concat";
    }
    $cptfile++;
}

print "- $sortie_rmesfam_f file created\n";

###############################################################################
############################ PROG #############################################;

## OUVERTURE ##

# Test si la longueur mini entree est bien < a la longueur maxi entree
if ($hmin > $hmax)
{
    die "$prog_tag [Error] La longueur minimal (hmin:hmin) doit etre inferieure a la longueur maximale (hmax:$hmax), line ".__LINE__."\n";
}

# Ouverture du fichier $ARGV[0] contenant les resultats a traiter
open (FIC,'<', $sortie_rmesfam_f) or die "$prog_tag [Error] Cannot open $sortie_rmesfam_f file$!, line ".__LINE__."\n";


## VARIABLES ##

my $ligne;

# Tableau contenant les mots (colonne "mot" du fichier de resultats)
my @mots;
# Tableau contenant les scores (colonne "stat" du fichier de resultats)
my @score;
# Tableau contenant la longueur des mots du tableau @mot
my @longueur;

my $p;
my $i=0;

# Variable donnant la valeur de l'ordre ($longueurmot-1)
my $ordre;
# Taille des mots
my $longueurmot;



#### PROGRAMME ####

# On parcours le fichier ligne par ligne
while ($ligne = <FIC>)
{
# Recherche des motifs "mots" formes de 'atcg'->$1 et de nombre dans la colonne stats (apres le 4eme '|')->$3
    if ($ligne =~ /\s([acgt]{$hmin,$hmax})(\s\|[^\|]+){3}\s\|([^\|]+)/)
    {
	# Ajouts des differents elements dans les tableaux correspondants
		push @longueur, length ($1);
		push @mots, $1;
		push @score, $3;
    }
}

# Ouverture du fichier de sortie
my $sortie_f = "$dossier/sortie_$fic\_rmes\_$hmax\_$hmin";
open (FSOR, '>', $sortie_f)or die "$prog_tag [Error] Cannot open $sortie_f file, line ".__LINE__."\n";

print "- Le fichier $dossier/sortie_$fic\_rmes\_$hmax\_$hmin a ete genere\n";

# Test pour les differentes longueurs possible (entre $max et $min -> ordre decroissant)
for ($longueurmot=$hmax; $longueurmot>=$hmin; $longueurmot--)
{
	# Calcul et ecriture de l'ordre (longueur mot -1)
    $ordre = $longueurmot - 1 ;
    print FSOR "#ordre $ordre\n";
    
	# On parcours les tableau en testant si la longueur du mot contenu dans le tableau @longueur, correspond a la variable $longueur mot. Si c'est le cas, le mot et le score s'affiche.
    foreach $p(0..$#mots)
    {
		if($longueur[$p] == $longueurmot)
		{
		    print FSOR "$mots[$p] $score[$p]\n";
		}
    }
}
close FSOR;

##############################################################################################################
# PERMET DE GENERER LES SEQUENCES COMPLEM ET D'ELIMINER LES SEQ PALINDROMIQUES A PARTIR DU FICHIER DE SORTIE # 

my $motif;
my $motif_inv;
my $motif2;
my $reste_ligne;
my $l;

my $mots_f   = "$dossier/mots_plus_comp_sans_palind_$fic\_rmes\_$hmax\_$hmin.txt";
open (FIC, '<', $sortie_f)or die "$prog_tag [Error] Cannot open $sortie_f file, line ".__LINE__."\n";
open (FSOR,'>', $mots_f)  or die "$prog_tag [Error] Cannot create $mots_f file, line ".__LINE__."\n";

while (my $line = <FIC>)
{
    if ($line =~ /([atcg]{2,8})(.+)/)
    {
	$motif = $1;
	$reste_ligne = $2;

	########## seq complem;
	$l= length($motif);
	for ($p = $l; $p>=0; $p--)
	{
	    $motif_inv .= substr($motif,$p,1);
	}
	foreach $p(0..$l)
	{
	    if (substr($motif_inv,$p,1) eq 'a') {$motif2 .= 't';next}; 
	    if (substr($motif_inv,$p,1) eq 't') {$motif2 .= 'a';next};
	    if (substr($motif_inv,$p,1) eq 'g') {$motif2 .= 'c';next};
	    if (substr($motif_inv,$p,1) eq 'c') {$motif2 .= 'g';next};
	}
	####

	print FSOR "$line";
	
	if ($motif2 ne $motif)
	{
	    print FSOR "$motif2 $reste_ligne\n";
	}
    }

    else 
    {
	print FSOR $line;
    }
	
    $motif ='';
    $motif_inv='';
    $motif2='';
}
close FIC;
close FSOR;
print "- le fichier $mots_f a ete genere\n";

