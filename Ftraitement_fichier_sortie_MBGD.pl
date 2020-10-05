#!/usr/bin/perl -w

use strict; 
$|++;

## TRAITE LES FICHIERS DE SORTIE DE LA BDD MBGD ET DONNE UN FICHIER DE SORTIE POUR CHAQUE PAIRES DE BACTERIES ORTHOLOGUES
## FONCTIONNE CORRECTEMENT

# VERIF TRI_RAPIDE ********************************************
# my @toto1 = qw(premier deuxieme troisieme trois cinquieme septieme);
# my @toto2 = qw(t o a t q s);
# my @toto3 = qw(1 2 3 4 3 4);
# my @toto4 = qw(p d t q c s);

# &TRI_RAPIDE_LETTRE(0,$#toto1,\@toto1,\@toto2,\@toto3,\@toto4);

# print "@toto1\n@toto2\n@toto3\n";   
# die;
# ******************************************************************************

my $prog_tag = '[Ftraitement_fichier_sortie_MBGD.pl]';
my $not_verbose = 1; # flag for debug
my $treat_two_bact_only = 0; # we treat only two sequence (otherwise create every combinations possibles between bacterial IDs)

#### OUVERTURE ####
# the program takes in 2 inputs: the mbgd file (a txt file)  and 2.  results directory

if(($#ARGV < 3) or (!-e $ARGV[0]))
{
      die "to use this programme: $0 <file of orthologues(MBDG file)>, <result directory>, <first bacterium id>, <second bacterium id>\n";
}

my $bac1_user= $ARGV[2];
my $bac2_user= $ARGV[3];
my $dossier = $ARGV[1];
$dossier =~ s/([^\/])$/$1\//;

# Test si le chemin du dossier entre existe, sinon renvoie d'un message d'erreur
if(! -d $dossier)
{
    die "le dossier $dossier n'existe pas";
}

open (FIC,"< $ARGV[0]") or die "$prog_tag [Error] impossible d\'ouvrir ce fichier $ARGV[0] :$!";



#### DATE ####

#   0      1     2        3      4       5       6       7       8
# my($sec, $min, $heure, $mjour, $mois, $annee, $sjour, $ajour,
my $est_dst = localtime;
# die localtime( time());
my ($ce_jour_lettre, $ce_jour_num, $ce_mois, $ce_mois_num,
$cette_annee, $date);
$ce_jour_lettre = ( qw(sunday monday tuesday wednesday thursday friday
saturday))[(localtime)[6]];
$ce_jour_num    = (localtime)[3];
$ce_mois_num    = sprintf "%02u",(localtime)[4]+1;
$ce_mois        = ( qw(january february march april may june july august
september october november december))[(localtime)[4]];
$cette_annee    = 1900 + (localtime)[5];

my $infos .= "DATE: $ce_jour_lettre, the $ce_jour_num of $ce_mois $cette_annee\n";
$date = $ce_jour_num.$ce_mois_num.substr($cette_annee,2,2);
print $infos;

# Charu's changes
# I am not putting the time stamp in the saved file so the above date part ($date) is not being used

#### VARIABLES ####


# Incremente le nombre de ligne
my $i = 0;

# Variables utilisees dans differentes boucles...
my $j = 0;
my $p = 0;
my $q = 0;

# Tableau contenant les ID des bacteries
my @ID_bact;

# Tableau contenant le motif correspondant a la fonction des genes
my @signature;

# Tableau contenant le chiffre correspondant a la motif des genes
my @fonction;

# Tableau contenant le nom des genes
my @name = ();

# Tableau de tableau contenant les ID de genes pour toutes les bacteries
my @bacterie = ();


# Contient les ID des gene, qui sont ajoute au tableau @bacterie (en fait, contient pour 1 ligne, l'ID de bacterie puis ID de gene correspondant puis ID de bacterie... etc...
my @ID_gene = ();

# Variable permettant d'afficher sign ou match dans le nom du fichier de sortie, selon le cas (elle ecrira "match" par defaut)
my $sign_match = "match";

# Variable contenant le separateur entre les ID de genes paralogues
my $ID_separateur = ",";

# Variable contenant l'ID precedent dans le tableau (sert a eviter les repetitions
my $doublon = '';

# Tableau contenant toutes les donnees de sortie pour un fichier, ce qui permet de trier ces donnees et d'eviter les repetitions
my @sortie = ();

####### PROGRAMME ########
# not needed my $tag; # will tell if motif is there or not
my $prem_line = <FIC>;

@ID_bact = split/[\s+]/, $prem_line;         # separation de tous les elements de cette 1ere ligne dans @ID

# define which bacteria id we are interested in to decrease amount of data stored
my $index_bac1_user    = undef;
my $index_bac2_user    = undef;
my $index_FuncCat_mbgd = undef;
my $index_Gene         = undef;
my $index_Description  = undef;
my $ind_first_col = 3; # 2 previously, changed the 2017 06 13
my $ind_last_col;
foreach my $i(0..$#ID_bact){ 
    if($ID_bact[$i] eq 'Gene'){
	$ind_last_col = $i;
	# added 20160619
	$ind_last_col--;
	$index_Gene = $i;
    }
    # for following indexes, we add +ind_first_col because the ind_first_col first columns were removed in
    # ID_bact and these indexes will be used on raw lines to get info
    elsif($ID_bact[$i] eq 'FuncCat(mbgd)')
    {
	$index_FuncCat_mbgd = $i;
    }
   
    elsif($ID_bact[$i] eq 'Description')
    {
	$index_Description = $i;
    }    
}

@ID_bact = @ID_bact[$ind_first_col..$ind_last_col];


foreach(0..$#ID_bact)
{
    if($ID_bact[$_] eq $bac1_user)
    {
	$index_bac1_user = $_;
    }
    elsif($ID_bact[$_] eq $bac2_user)
    {
	$index_bac2_user = $_;
    }
}
defined $index_bac1_user or die "$prog_tag [Error] bac1_user ($bac1_user) not found in file $ARGV[0], line ".__LINE__."\n";
defined $index_bac2_user or die "$prog_tag [Error] bac2_user ($bac2_user) not found in file $ARGV[0], line ".__LINE__."\n";
defined $index_FuncCat_mbgd or die "$prog_tag [Error] 'FuncCat(mbgd)' header not found in file $ARGV[0], line ".__LINE__."\n";
defined $index_Gene or die "$prog_tag [Error] 'Gene' not found in file $ARGV[0], line ".__LINE__."\n";
defined $index_Description or die "$prog_tag [Error] 'Description' not found in file $ARGV[0], line ".__LINE__."\n";

print "$prog_tag index FuncCat_mbgd:$index_FuncCat_mbgd\n";
print "$prog_tag index Gene:$index_Gene\n";
print "$prog_tag index Description:$index_Description\n";

## Charu's changes
## we cannot be sure if the user is inputting the MBGD file withthe motifs or without motifs so we have to do some additional checking
# if ( ($ID_bact[$#ID_bact-2] eq "FuncCat") &&	($ID_bact[$#ID_bact-1] eq "Motifs") &&	($ID_bact[$#ID_bact] eq "Gene") )
#   {@ID_bact = @ID_bact[$ind_first_col..($#ID_bact-3)];}       # extraction des ID des bacteries seulement que l'on place dans @ID_bact

# elsif ( ($ID_bact[$#ID_bact-1] eq "FuncCat")  && ($ID_bact[$#ID_bact] eq "Gene") )
#   {@ID_bact = @ID_bact[$ind_first_col..($#ID_bact-2)];}       # extraction des ID des bacteries seulement que l'on place dans @ID_bact

## we check if the bacterial ids the user has input are contained in the mbgd file or not...
# die "ID_bact @ID_bact\n";

my $cpt_gene_without_name = 1;
my $deb_noname = "noName";
$i++;

my $b_bac1_found = undef;
my $b_bac2_found = undef;
my $nb_lines_with_both_ids = 0;

# Defilement du fichier ligne par ligne
while (my $line = <FIC>)
{


    # code to retain only lines with IDs we are interested in
    $b_bac1_found = 0;
    $b_bac2_found = 0;
    if($line =~ /$bac1_user\:/){ $b_bac1_found = 1; }
    if($line =~ /$bac2_user\:/){ $b_bac2_found = 1; }
    ($b_bac1_found and $b_bac2_found)and $nb_lines_with_both_ids++;
    ($b_bac1_found and $b_bac2_found)or next;
    
    # A la 1ere ligne
    $i++;
    
    # Recherche du motif xxx:XX(X)0000(0)
    if ((@ID_gene) = $line =~ /([a-z]{2,3}\d*):([^\s]+)/g)
    {
	# print $prog_tag.' '.join(',', "ID_gene:", @ID_gene, "\nfor line $line, line ".__LINE__."\n");
	
	
	# et XXX0000 =>$2
	# $ID_gene=$2;
	
	# Elimination du (chiffre) a la fin de certains ID de genes
	# in ID_gene you have:
	# bacterium_id1 gene_id1 bacterium_id2 gene_id2 bacterium_id3 gene_id3
	for(my $nID_g = 0; $nID_g < $#ID_gene; $nID_g+=2)
	{ 
	    $ID_gene[$nID_g+1] =~ s/\(\d+\)$//;

	    # print "$prog_tag ID_gene of $nID_g ($ID_gene[$nID_g] with gene $ID_gene[$nID_g+1]), line ".__LINE__."\n";
	    
	    foreach $j (0..$#ID_bact)
	    {               
		# Recherche de l'ID de bacterie correspondant a xxx
		if ($ID_gene[$nID_g] eq $ID_bact[$j])                   
		{
		    # print "$prog_tag ID_gene of $nID_g ($ID_gene[$nID_g]) eq ID_bact of $j ($ID_bact[$j]), we are about to record $ID_gene[$nID_g+1], line ".__LINE__."\n";
		    # Si la case correspondante est deja remplie, on ajoute le nouvel ID gene (paralogue) a l'ID gene deja present
		    if (exists $bacterie[$j][$i-2])
		    {
			$bacterie[$j][$i-2] .= "$ID_separateur$ID_gene[$nID_g+1]";
			# print "$ID_separateur$ID_gene[$nID_g+1] recorded in  bacterie of $j $i-2 (".$bacterie[$j][$i-2]."), line ".__LINE__."\n";
		    }
		    else 
		    { 
			# Ajout dans la colonne appropriee de l'ID du gene s'il correspond a la bacterie de cette colonne
			push @{$bacterie[$j]},$ID_gene[$nID_g+1];
			# print "$prog_tag we record ".$ID_gene[$nID_g+1]." in bacterie of $j array:".join(',',  @{$bacterie[$j]}).", line ".__LINE__."\n";
		    }
		} 
	    }
	} 
    } 

    # print "$prog_tag bacterie array sizes for first bacterium: $#bacterie $#{$bacterie[0]}, line ".__LINE__."\n";
    
    # Recherche pour chacune des colonnes
    foreach $j (0..$#ID_bact)                      
    {          
	# s'il existe des cases vides apres le remplissage prealable avec les ID des genes                
	if (!exists $bacterie[$j][$i-2])              
	{         
	    # et si c'est le cas, on place un vide ("") dans cette case afin que le push suivant n'ait
	    # pas une ligne de retard      
	    push @{$bacterie[$j]},"";                  
	}                                             
    }
    # print "$prog_tag recorded gene ids:\n";
    # foreach $j (0..$#ID_bact)                      
    # {
    # 	print "$prog_tag size for bacterium $j: $#{$bacterie[$j]}, line ".__LINE__."\n";
    # 	foreach my $z(0..$#{$bacterie[$j]})
    # 	{
    # 	    print $bacterie[$j][$z].',';
    # 	}
    # }
    # print "\nfor line\n$line\nline ".__LINE__."\n";
    # die;

    # Recherche du motif "1chiffre" correspondant aux fonctions des groupes d'orthologues    
##    if ($line=~ /\s+(\d{1,3})\s*[^\s]*\s*$/)      OLD: before 200508
#    if ($line=~ /\s+([\d\.]+)\s*$/)            
#    {                                         
#	# Quand le motif est trouve, il est ajoute a @fonction 
#	push (@fonction,$1);                      ## function stores the "motif" given in the MBGD file
#    }
#    # mais si le motif n'est pas present                                      
#    else                                    
#    {           
#	die "For the line: $line, we do not find the function number\n";
#	# dans ce cas on ajoute "" a @function pour qu'il n'y est pas de decalage              
#	push (@fonction, " ");                        
#    }
    # replaced by ... the 2017 06 13 
    push @fonction, (split /\t/, $line)[$index_FuncCat_mbgd];
    
    # Recherche du motif "xxx(A)(A)" correspondant aux noms des genes
# #    if ($line=~ /\s+\d+\s*(\w{1,3}[^\s]*)?\s*$/)      OLD: before 200508     
#    if ($line=~ /(\s\w{1,3}\d*|)\s+([\d\.]+\s+){3}[\d\.]+\s*$/)      
#    {                          
#	if(not defined $1){ 
#	    print "name not defined for line $line";
#	    push @name, "";
#	}
#	else{
#	    # Quand le motif est trouve, il est ajoute a @name
#	    push (@name, $1);                        
#	    $name[$#name] =~ tr/ //;
#	}
    #    }
#       # mais le nom n'est pas toujours present                                      
#    else                                    
#    {           
#	# dans ce cas on ajoute "" a @name pour qu'il n'y ait pas de decalage             
#	die "EXPRESSION REGULIERE POUR LE NOM PAS RECONNUE pour ligne $line";
#	push (@name, " ");                        
#    }
    # replaced by ... the 2017 06 13
    my $tmp_name = (split /\t/, $line)[$index_Gene];
    my $tmp_description = (split /\t/, $line)[$index_Description];
    if(defined $tmp_name)
    {
	push @name, $tmp_name;
	$name[$#name] =~ tr/ //;
    }
    else
    {
	push @name, '';
    }
    if(defined $tmp_description)
    {
	chomp($tmp_description);
	$tmp_description =~ s/ /_/g;
	if($name[$#name] ne ''){ $name[$#name] .= ';'; }
	$name[$#name] .= $tmp_description;
    }
    
    # Recherche des motifs "xxxx:xxxx:xxxx" correspondants aux signatures des genes
    while ($line=~ /\s*[^\s]+:([\w]+:[\w]+)\s*/g)
    {
	# Si la case correspondante a la ligne en cours est vide
	if (!exists $signature[$i-2])	#@signature stores the motif given in the file
	{
	    # la signature est ajoutee
	    push (@signature, $1);
	    # $sign_match affichera sign dans le nom de fichier
	    $sign_match = "sign";
	}
	# Par contre si elle est deja pleine et qu'il y a plusieurs signatures sur la ligne 
	else
	{
	    # la nouvelle signature est rajoutee au bout de la chaîne en separant par un ";"
	    $signature[$i-2] .= ";$1";
	}
    }
    # Si apres passage sur la ligne, la recherche n'a rien donnee
    if (!exists $signature[$i-2])
    {
	# Un vide "" est place dans le tableau
	push (@signature, "");
    }
}
print "$prog_tag number of lines with $bac1_user and $bac2_user IDs:$nb_lines_with_both_ids\n";


if(! $not_verbose)
{
    print "verif taille tableaux\n";
    foreach my $j (0..$#ID_bact)                      
    {          
		print "taille bacterie $j $#{$bacterie[$j]}\n";
    }
    print "end verif taille tableaux\n";
    
    print "verif contenu tableaux\n";
    for my $i(0..$#{$bacterie[0]})
    { 
		for my $j (0..$#ID_bact)                      
		{          
		    print "$ID_bact[$j]: $bacterie[$j][$i] ";
		}
		print "fct: $fonction[$i]";
		if(defined $signature[$i]){ print " sign: $signature[$i]"; }
		if(defined $name[$i])     { print " name: $name[$i]\n";    }
		else{ print "\n"; }
    }
    print "end verif contenu tableaux\n";
    # exit;
}


 
#### SORTIES ####

# Variable contenant le nombre de combinaisons de paires de bacteries possible
my $nb;

# Variable contenant la valeur de l'ID de la 1ere bacterie de la paire
my $bact1=0;

# Variable contenant la valeur de l'ID de la 2nde bacterie de la paire
my $bact2=0;

# Variable permettant le calul des combinaison de paires de bacterie possible
my $combin=0;



# Calcul du nombre de combinaisons de paires possible (si nb bact = n, le nb de comb de paires est de (n-1)+(n-2)+...+1 ) 
$p=0;


foreach $q(1..$#ID_bact)
{
    $p+=$q;
}


# Pour les $p combinaisons possible
foreach $nb(1..$p)
{

    # les ID des bacteries prennent toutes les valeurs possibles correspondant aux differentes combinaisons possibles
    if ($bact2 == $#ID_bact) 
    {
		$bact1++;
		$combin++;
		$bact2=0;
		$bact2+=$combin;
    }
    $bact2++;

    if(($ID_bact[$bact1] eq $bac2_user)and($ID_bact[$bact2] eq $bac1_user))
    {
		my $tmpid = $bact1;
		$bact1 = $bact2;
		$bact2 = $tmpid;
#      foreach my $id(0..$#ID_bact){
#	if($ID_bact[$id] eq $bac1_user){
#	  $bact1 = $id;
#	}
#	if($ID_bact[$id] eq $bac2_user){
#	  $bact2 = $id;
#	}
#      }
      
    }
    elsif($treat_two_bact_only and (($ID_bact[$bact1] ne $bac1_user) or ($ID_bact[$bact2] ne $bac2_user)))
    {
		# added 2017 06 13
		die "ID_bact of $bact1 ($ID_bact[$bact1]) ne bac1_user ($bac1_user) OR ID_bact of $bact2 ($ID_bact[$bact2]) ne bac2_user ($bac2_user) line ".__LINE__."\n";
		next;
    }
    elsif(($ID_bact[$bact1] ne $bac1_user) or ($ID_bact[$bact2] ne $bac2_user))
    {
	# added 2019 02 26
	# avoid treatment of all possible pairs, only analyse the requested
	next;
    }

    # ???? # path will be given (we do not need this variable if stored in SIGffRid directory
    # Creation des differents fichiers de sortie (un pour chaque paire), dont le nom sera IDbact1_IDbact2  
		#first $sign_match=sign and then $sign_match= match;
    

    open (FSOR,"> $dossier/ortho_$ID_bact[$bact1]_$ID_bact[$bact2]_fct_$sign_match.txt") or die "impossible de creer le fichier :$!";

    # Affichage du nom de fichier dans la fenetre du terminal
    print "le fichier ortho_$ID_bact[$bact1]_$ID_bact[$bact2]_fct_$sign_match.txt a ete genere\n";

    # Affichage en premiere ligne de IDbact1 IDbact2
    print FSOR "$ID_bact[$bact1] $ID_bact[$bact2] fct name\n";
    print "ID bact $ID_bact[$bact1] bact2 $ID_bact[$bact2]\n";
    &TRI_RAPIDE_LETTRE(0,$#{$bacterie[$bact1]},\@bacterie,$bact1,\@signature,\@name,\@fonction);

    if(! $not_verbose)
    {
       # if($bact1 == 2)
	# {
	    print "verification TRI\n";
	    for my $l (0..$#{$bacterie[$bact1]})
	    {
	       print "$bacterie[0][$l] $bacterie[1][$l]\n";
	    }
		print "verification TRI\n";
	   
	   # exit;
	# }
    }

    # $doublon contient le premier ID des l'ID precedents de la bacterie
    $doublon = '';
    # Affichage des resultats correspondant, pour chaque fichier de sortie cree
    foreach $p(0..$i-2)
    {
	# Seuls les couples d'orthologues entre les deux bacteries comparees sont affiches

	if (($bacterie[$bact1][$p] =~ /\w+/) and ($bacterie[$bact2][$p] =~ /\w+/))
	{
	    # on recupere le premier indice
	    $bacterie[$bact1][$p] =~ /^([^\, ]+).*/;
	    my $firstID = $1;
	    print "firstID vaut $firstID, ID $bacterie[$bact1][$p]\n";
	    
	    
	    #$doublon = $bacterie[$bact1][$p-1];
	    #$doublon =~ s/([^,]+).*/$1/;
	    
	    # si le 1er ID de la ligne precedente (concernee par cette bacterie) n'est pas identique (donc si y a pas de repetitions),on ecrit la ligne
	    if ($firstID ne $doublon)
	    {
		# Affichage sous la forme: IDgene0A1_IDgeneA2_IDgeneAn IDgeneB1_IDgeneB2_IDgeneBn  sign1;sign2 fonction nom
		# avec geneA orthologue geneB   et   geneA1 paralogue gene A2 et An
		# print FSOR "$bacterie[$bact1][$p] $bacterie[$bact2][$p] $fonction[$p] $name[$p] $signature[$p]\n"; 
		print FSOR "$bacterie[$bact1][$p] $bacterie[$bact2][$p] $fonction[$p] $name[$p]\n"; 
		# die "b1:$bacterie[$bact1][$p] b2:$bacterie[$bact2][$p] fct:$fonction[$p] name:$name[$p] sign:$signature[$p]\n"; 				
		
		if(defined $firstID)
		{
		    $doublon = $firstID;
		}
		else
		{
		    die "ID $bacterie[$bact1][$p] not recognized line 327";
		}
	    }
	    else
	    {
		print "firstID $firstID already recorded\n";
	    }
	}
    }
}




#### FERMETURE ####



close FIC;
close FSOR;

exit;

{
    no strict;
	sub TRI_RAPIDE_LETTRE($$\@$\@\@\@)
	{
	    local($min, $max, *letter_tab, $b, *arrayVal0, *arrayVal1, *arrayVal2) = @_;
	    local($pos);
	    
	    # Est-ce que la position actuelle du tri est plus petite que la taille du tableau a trier ?
	    # Oui alors on commence/continue le tri
	    if ( $min < $max )
	    {
			# we get the position where last sort stop
			$pos = &PARTITIONNER_LETTRE($min, $max, *letter_tab, $b, *arrayVal0, *arrayVal1, *arrayVal2);
			($not_verbose)or print "min $min pos $pos\n";
			# we sort the left section of the remaining section in the table
			&TRI_RAPIDE_LETTRE($min, $pos, *letter_tab, $b, *arrayVal0, *arrayVal1, *arrayVal2);
			($not_verbose)or print "pos+1 ".($pos+1)." max $max\n";
			# we sort the right section of the remaining section in the table
			&TRI_RAPIDE_LETTRE($pos+1, $max, *letter_tab, $b, *arrayVal0, *arrayVal1, *arrayVal2);
	    }
	    # Non, on peut donc s'arreter, le tri est terminé !
	    else
	    {
			return( 1 );
	    }
	}

	sub PARTITIONNER_LETTRE($$\@$\@\@\@)
	{
	    local($min, $max, *letter_tab, $b, *arrayVal0, *arrayVal1, *arrayVal2) = @_;
	    local($i_sub, $j_sub, $centre, $tmp);
	    
	    # we record as temp. value for the center this of the left value
	    $centre = $letter_tab[$b][$min];
	    $i_sub = $min - 1;
	    $j_sub = $max + 1;
	    
	    
	    while ( 1 )
	    {
			$j_sub--;
			$i_sub++;
			# until (( $letter_tab[$j_sub] le $centre )or( $j_sub == 0 )) { $j_sub--; }
			# until (( $letter_tab[$i_sub] ge $centre )or( $i_sub == $#letter_tab )) { $i_sub++; }
		
			until (( $letter_tab[$b][$j_sub] le $centre )or( $j_sub == 0 )) { $j_sub--; }
			until (( $letter_tab[$b][$i_sub] ge $centre )or( $i_sub == $#{$letter_tab[$b]} )) { $i_sub++; }
			 
			if ( $i_sub < $j_sub )
			{
		
			    for my $nbb(0..$#bacterie)
			    {
					$tmp = $letter_tab[$nbb][$j_sub];
					$letter_tab[$nbb][$j_sub] = $letter_tab[$nbb][$i_sub];
					$letter_tab[$nbb][$i_sub] = $tmp;
			    }
			     
			    # $tmp = $letter_tab[$j_sub];
			    # $letter_tab[$j_sub] = $letter_tab[$i_sub];
			    # $letter_tab[$i_sub] = $tmp;
			     
			    # $tmp = $arrayTri[$j_sub];
			    # $arrayTri[$j_sub] = $arrayTri[$i_sub];
			    # $arrayTri[$i_sub] = $tmp;
			    if($#arrayVal0 != -1)
			    {
					$tmp = $arrayVal0[$j_sub];
					$arrayVal0[$j_sub] = $arrayVal0[$i_sub];
					$arrayVal0[$i_sub] = $tmp;
			    }
			    if($#arrayVal1 != -1)
			    {
				 $tmp = $arrayVal1[$j_sub];
				 $arrayVal1[$j_sub] = $arrayVal1[$i_sub];
				 $arrayVal1[$i_sub] = $tmp;
			    }
			    if($#arrayVal2 != -1)
			    {
				 $tmp = $arrayVal2[$j_sub];
				 $arrayVal2[$j_sub] = $arrayVal2[$i_sub];
				 $arrayVal2[$i_sub] = $tmp;
			    }
			}
			else
			{
		
			    return( $j_sub );
			}
	    }
	}
}
