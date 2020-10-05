#!/usr/bin/perl -w

# permet l'extraction des sequences promotrices (un ficher par groupe d'orthologues)
# traite les deux bacteries, regroupe les genes de differentes fontions dans des repertoires dediees si la fonction est precisee

# CAUTION: want for first line: 'first_bact_ID sec_bact_ID'

# CELUI A UTILISER POUR OBTENIR LES FICHIERS POUR CHAQUE RELATION D'ORTHOLOGIE A PARTIR DU FICHIER ISSU DE MGBD ET DES FICHIERS DE PROMOTEURS
use strict;
use Env qw(SIGFFRID_DIR);
use lib $SIGFFRID_DIR.'libperl';
use SeqUtilities;          # qw(print50);

my $prog_tag      = '[Fprend_seq_prom_ortho_pairs_fct.pl]';
my $bool_fasta    = 0;
my $chemin_dest   = '';  # chemin pour atteindre le repertoire ou l'on veut les resultats
my @nom_fic       = '';  # noms des fichier de sequences amont (incluant le chemin pour y acceder)
my $IDgbk         = '';  # identifiant de l'espece dont on va recuperer les sequences promotrices qui nous interessent
my $IDgbk2        = ''; # identifiant de la deuxieme espece dont on va recuperer les sequences promotrices qui nous interessent
my $length_min    = 30;   # longueur minimale des sequences promotrices que l'on souhaite
my $length_max    = 350;  # longueur maximale des sequences promotrices que l'on souhaite
my $dist_deb_trad = 0; # distance par rapport au debut de la traduction, et qui ne sera pas prise en compte pour l'extraction de promoteur
my $not_verbose   = 1;   # option facon linux
my $separateur_ID_regexp = '\,_';  # separe les identifiants de paralogues d'une meme bacterie ou les identifiants diff d'un meme gene
my $separateur_ID = ',';  # separe les identifiants de paralogues d'une meme bacterie
my $bool_fct_used = 0; # booleen permettant de savoir si on genere un repertoire par fonction ou si on met toutes les sequences dans le meme repertoire
my $MAX_FIC_NAME_LENGTH = 30;

# verification et enregistrement des parametres
# 
# FABRICE ORIGINAL CODE
if(($#ARGV == 6)or($#ARGV == 5)){
	  # fichier d'orthologues
	if(-e $ARGV[0]){
		open(F_ORTHO,"< $ARGV[0]") or die "$prog_tag [Error] $ARGV[0] file for genes ID can't be opened (first parameter), check name (and pathway maybe)!, line ".__LINE__."\n";
	}
	else{
		die "$prog_tag [Error] $ARGV[0] does not exist!, check name (and pathway maybe)!, line ".__LINE__."\n";
	}
	
	# identifiant d'espece
	$IDgbk = quotemeta($ARGV[1]);
	$IDgbk2 = quotemeta($ARGV[2]);
	
	# fichier comportant toutes les sequences promotrices de l'espece qui nous interesse
	if(-e $ARGV[3]){
		$ARGV[3] =~ /\/?([\w\.]+)$/;
		$nom_fic[0] = $ARGV[3];
	}
	else{
		die "$prog_tag [Error] $ARGV[3] does not exist!, check name (and pathway maybe)!, line ".__LINE__."\n";
	}
	
	# fichier comportant toutes les sequences promotrices de l'espece qui nous interesse
	if(-e $ARGV[4]){
		$ARGV[4] =~ /\/([\w\.]+)$/;
		$nom_fic[1] = $ARGV[4];
	}
	else{
		die "$prog_tag [Error] $ARGV[4] does not exist!, check name (and pathway maybe)!\n";
	}
	
	#($ARGV[5] =~ /^[01]$/)or die "$prog_tag [Error] $ARGV[5], fifth parameter has to be a boolean";
	$bool_fct_used = $ARGV[5];
	
	# chemin pour atteindre le repertoire ou l'utilisateur veut ses resultats
	if($#ARGV == 6)
	{
		$chemin_dest = $ARGV[6];
		if (-e "$chemin_dest/SIGffRid_orthologs")
		{ 
			print " the directory SIGffRid_orthologs already exists:  \n\t\t\t We delete it and are creating a new one. All earlier data is erased \n";
			system ("rm -fr $chemin_dest/SIGffRid_orthologs");
		# 		system ("mkdir $chemin_dest/SIGffRid_orthologs ");
		}
		
		if ($chemin_dest !~ /\/$/)
		{
		  $chemin_dest .= '/';
		}
		
		system ("mkdir $chemin_dest"."SIGffRid_orthologs");
		$chemin_dest .= 'SIGffRid_orthologs/';
		
		if ($bool_fct_used) {	$chemin_dest .= 'SIGffRid_ortholog_files'; }	
	}
}

else{
	die "$prog_tag Give good parameters: this program (for bacteria) must be used as following\n$0\n<file with gene ID (blast_donne_orth... or ortho...)>\n<first species ID>\n<second species ID>\n<file containing upstream sequences of the first given organism>\n<file containing upstream sequences of the second given organism>\n<boolean: we create directories for functions (1) or not (0)>\n[output directory], line ".__LINE__."\n";
}



# nous chargeons le fichier des identifiants de genes orthologues
my @f_ortho = ();
while(<F_ORTHO>)
{
    if(/^(?:\s*$|>)/){ next; }
    push @f_ortho, $_;
}
close F_ORTHO;

# print "verif f_ortho\n";
# foreach(@f_ortho)
# {
  #  print $_;
# }
# print "END verif f_ortho\n";
# die;
# tri par rapport a la colonne contenant l'identifiant de l'espece qui nous interesse
my @colonne_gene_a_trier = (-1,-1);

my @colonnes = split /\s+/, $f_ortho[0];
pop @colonnes; # we remove "name" chain
pop @colonnes; # we remove "fct" chain

# die "$#colonnes\n";
colsearch:for my $colonnes_cherche(0..$#colonnes)
{
	if($colonnes[$colonnes_cherche]=~ /(\|REFGENOMES\|)?$IDgbk/){
		$colonne_gene_a_trier[0] = $colonnes_cherche;
	}
	if($colonnes[$colonnes_cherche]=~ /(\|REFGENOMES\|)?$IDgbk2/){
		$colonne_gene_a_trier[1] = $colonnes_cherche;
	}
	
	if(($colonne_gene_a_trier[0] != -1)and($colonne_gene_a_trier[1] != -1))
	{ 
	  last colsearch;
	}
}
if(($colonne_gene_a_trier[0] == -1)or($colonne_gene_a_trier[1] == -1)){
	die "$prog_tag [Error] One column of your orthologues file has to own the ID given $IDgbk and another $IDgbk2 ID!, line ".__LINE__."\n";
}

shift @f_ortho; # to eliminate 'SAV SCO' first line (only species ID)
# my $first_valid_line = 0;
# suppression des lignes inutiles en tete du fichier d'orthologues
# SUPP:for my $empty_line(@f_ortho){
  # if($empty_line =~ /$IDgbk\d+/){ last SUPP; }
  # else{ $first_valid_line++;}
# }
# (!$first_valid_line)or splice @f_ortho, 0, $first_valid_line++;

# die "colonnes @colonnes @colonne_gene_a_trier\n";

# sous programme qui enregistre dans les fichiers (uniquement pour la premiere relation d'orthologie concernee par le gene traite)
# ATTENTION, il utilise des variables globales declarees ci-dessus


open(NOEXISTING,'> '.$chemin_dest."No_existing_promot_seq.txt")or die "$prog_tag [Error] Impossible to create $chemin_dest"."No_existing_promot_seq.txt file, line ".__LINE__."\n";

my @fct = (); # store fucntion (integers) corresponding to created directories

#  pour chaque bacterie
for my $col_a_trier(0..$#colonne_gene_a_trier)
{
    print "bacterie correspondant au fichier $nom_fic[$col_a_trier] $colonnes[$col_a_trier]\n";
    system('sleep 3');   
    open(F_SEQ_PROM,"< $nom_fic[$col_a_trier]") or die "$prog_tag [Error] $nom_fic[$col_a_trier] file with promoting sequences can't be opened (parameter), line ".__LINE__."\n";
    # nous chargeons le fichier de toutes les sequences promotrices de l'espece qui nous interesse
    my @f_seq_prom = <F_SEQ_PROM>;
    close F_SEQ_PROM;

    # ADAPTATION REQUISE POUR PRISE EN COMPTE DE LA TOTALITE DE L'IDENTIFIANT (QUOTEMETA NECESSAIRE)
    # tri suivant la colonne correspondant a l'identifiant de l'espece
    my(@compare1,@compare2);
    @f_ortho = sort{
		@compare1 = split /\s+/,$a;
		@compare2 = split /\s+/,$b;
		# ce deuxieme "split" s'explique sur la necessite de faire le tri a partir du numero final identifiant le gene (sinon le tri n'est pas vraiment bon (un cas ne fonctionne pas)(a part l'identifiant, la partie suivant REFGENOME) est identique pour tous les genes)
		if($compare1[$colonne_gene_a_trier[$col_a_trier]]=~ /\|REFGENOMES\|/){
		    @compare1 = split /\|REFGENOMES\|/,$compare1[$colonne_gene_a_trier[$col_a_trier]];
		    @compare2 = split /\|REFGENOMES\|/,$compare2[$colonne_gene_a_trier[$col_a_trier]];
		    $compare1[1] cmp $compare2[1]	 
		}
		else{
		    $compare1[$colonne_gene_a_trier[$col_a_trier]] cmp $compare2[$colonne_gene_a_trier[$col_a_trier]]
		}
    }@f_ortho;
    # die "@f_ortho\n";

    my $deb_boucle      = 0; # point d'initiation de la boucle B2, cela evite de reparcourir des seq promotrices qui ne nous serons plus utiles
    my $prev_tempo      = 0; # memorise l'identifiant du dernier gene concerne par la sequence sauvegardee (evite de faire n fois une meme extraction)
    my $memoire_seq     = ''; # memorise le dernier identifiant de gene et la sequence promotrice correspondante
    my $good_length     = 0; # indique si le dernier identiifiant concerne avait une sequence promotrice de taille adequate (la encore, evite des traitements inutiles si on retrouve cet ID dans les orthologues)
    my $intitule_sortie = ''; # memorise l'intitule
    my $nom_fic_sor     = ''; # recueillera le nom du fichier de sortie (seq amonts des orthologues)

    my $reg_expr_rech       = '';
    my $reg_expr_rech_apres = '';
    my $alphabet_ID_strict  = '\w\d\-_\|\.';
    my @intitule_plus_seq   = ();
    my %H_ID_ln             = ();
    my @ID_redondant        = ();# liste qui va permettre de supprimer les identifiants de genes pouvant correspondre a plusieurs genes (en raison de reannotation par exemple)
    
    # recuperation de tous les intitules de sequences promotrices precedes du numero de ligne ou il a ete trouve dans le fichier 
    open(INTITULEPROM,"grep -n '>' $nom_fic[$col_a_trier]|");

    # creation d'une table de hachage associant a chaque identifiants (plusieurs possibles pour un meme gene), le numero de ligne de l'intitule du promoteur correspondant (dans le fichier de sequences promotrices)
    for my $intitule(<INTITULEPROM>)
    {
		$intitule =~ /^(\d+):>[^:]+:\s*([$alphabet_ID_strict$separateur_ID_regexp]+)/;
		my $num_line = $1;
		if(! defined $2){ die "$prog_tag [Error] dol2 non defini ligne $intitule, line ".__LINE__."\n";	}
		my @ID = split /[$separateur_ID_regexp]/, $2; # =~ /([^$separateur]+)/g);
		if(! defined $num_line){ die "$prog_tag [Error] dol1 non defini ligne $intitule, line ".__LINE__."\n";	}
		if($#ID == -1){ die "$prog_tag [Error] ID non defini ligne $intitule, line ".__LINE__."\n";}
		# les identifiants sont tous precedes de (donc separes par) _, d'ou la necessite du shift pour supprimer la premiere case (vide)
		# shift @ID;
		# die "@ID, separateur_ID_regexp $separateur_ID_regexp\n";
		foreach my $ID(@ID)
		{
		    my $uc_ID = uc($ID);
		    if(! exists $H_ID_ln{$uc_ID})
		    {
				$H_ID_ln{$uc_ID} = $num_line;      
				print "on enregistre l'ID $uc_ID de la ligne $num_line\n";
		    }
		    # si le meme ID est present sur la meme ligne, cela ne pose pas de probleme puisqu'il s'agit d'un seul gene. Dans le cas contraire il y a redondance d'identifiant: nous devons supprimer cet identifiant
		    elsif( $num_line !=  $H_ID_ln{$uc_ID}){  push @ID_redondant, $uc_ID; }
		}
    }
    close INTITULEPROM;
    
    # suppression dans la table de hachage correspondant a des identifiants n'etant pas propres a un seul gene
    for my $ID(@ID_redondant)
    {
		if(exists $H_ID_ln{$ID}){ delete $H_ID_ln{$ID}; }
    }

    # generation d'un expression reguliere qui permettra de chercher l'identifiant en fonction de la colonne concernee par l'espece qui nous interesse 
    # si la colonne qui nous interesse est la premiere ($colonne_gene_a_trier == 0), nous n'avons pas besoin de modifier $reg_expr_rech etant donne que ce sera toujours le premier identifiant trouve qui nous interessera.
    
    # if more than 2 species IDs on a line
    for my $cpt_col(1..$colonne_gene_a_trier[$col_a_trier])
    {
		$reg_expr_rech .= '['.$alphabet_ID_strict.$separateur_ID_regexp.']+\s+';
    }

    # for IDs after interesting one
    for my $cpt_col($colonne_gene_a_trier[$col_a_trier]+1..$#colonnes)
    {
		$reg_expr_rech_apres .= '\s+['.$alphabet_ID_strict.$separateur_ID_regexp.']+';
    }
	
	# die "^$reg_expr_rech([$alphabet_ID_strict$separateur_ID_regexp]+)$reg_expr_rech_apres\s*((\d+)\s*(.+?))?$\n";
	B1:for my $line(@f_ortho)
	{
		# print "EXP_REG line $line being treated with ^$reg_expr_rech([$alphabet_ID_strict$separateur_ID_regexp]+)$reg_expr_rech_apres(\d+)?\s+?(.+?)\s*$\n";
		#  $reg_expr_rech correspond a une expression reguliere generee dynamiquement on fonction de la colonne dans laquelle nous devons trouver l'identifiant qui nous interesse
		
		if($line =~ /^$reg_expr_rech([$alphabet_ID_strict$separateur_ID_regexp]+)$reg_expr_rech_apres\s*((\d+)\s*(.+?))?$/)
		{
			# if($line =~ /^$reg_expr_rech([$alphabet_ID_strict$separateur_ID_regexp]+)$reg_expr_rech_apres\s+?/){
			
			# debug eventuel ICI ??
			# print "for line $line, dol1 $1 dol2 $2 dol3 $3 trouves reg_expr_rech_apres $reg_expr_rech_apres\n";
			# exit;
			# recuperation de l'identifiant present dans la colonne correspondant a l'espece qui nous preoccupe
			my $tempoID = uc($1);
			#	  if($tempoID =~ /SCO_/){ die "$prog_tag [Error] Probleme for $line, ID found $tempoID:$!\n";}
			# print "tempoID $tempoID dol1 $1\n";
			my $fct     = '';
			my $autre   = '';
			
			# if we do not have only IDs
			if(defined $3)
			{
				# case we have fct, name, signatures
				($bool_fct_used)and $fct = $3.'/';
			}
			# case we have comments
			if((defined $4)and($4 !~ /^\s*$/))
			{
				$autre = $4;
				# print "autre $autre\n";
			}
			# print "dol1 $1 dol3 $3 dol4 $4\n";
	  
	
			# creation d'un nom de fichier sans espaces ni pipe
			$line =~ /^\s*([$alphabet_ID_strict$separateur_ID_regexp]+\s+[$alphabet_ID_strict$separateur_ID_regexp]+)/;
			$nom_fic_sor = $1;
			chomp($nom_fic_sor);
			$nom_fic_sor =~ s/\|/_/g;
			$nom_fic_sor =~ s/\s+/__/g;
			
				
			# if(defined $1){ $nom_fic_sor = $1; }
			# else{ die "$prog_tag [Error] Problem with regular expression for file name, dol1 not defined nom_fic_sor $nom_fic_sor!\n"; }
			
			
			
			# si l'identifiant du gene orthologue correspond au precedent...
			if($tempoID eq $prev_tempo)
			{
				# si la sequence amont  est suffisamment longue, on ecrit dans le fichier orthologue correspondant
				if($good_length)
				{
					open(FSOR,'>> '.$chemin_dest.$fct.$nom_fic_sor);
					print FSOR $intitule_sortie.&print50($bool_fasta,$memoire_seq)."\n";
					close FSOR;
					($not_verbose)or print "\nIn ".$chemin_dest.$fct.$nom_fic_sor." 1 we write\n$intitule_sortie".&print50($bool_fasta,$memoire_seq)."\n";
				
					print "$tempoID DEJA TROUVE, on le met dans ".$chemin_dest.$fct.$nom_fic_sor."\n";
				}
				# on passe a la relation d'orthologie suivante
				next B1;
			}
			
			$memoire_seq = '';

			# create subdirectory, one by function (numerically defined) and record existing dir
			if(! -d $chemin_dest.$fct)
			{
				mkdir $chemin_dest.$fct, 0755; 
				push @fct, $fct;
			}
			my @ID2 = split /[$separateur_ID_regexp]/, $tempoID;
			my $bool_1ID_trouve = 0;
  
			# if($ID2[0] =~ /SCO/){	  die "@ID2"; }
			# recherche dans le fichier de seq promotrices de chaque identifiant possible
			BOUCLEID2:for my $ID2(@ID2)
			{
			    my $uc_ID2 = uc($ID2);
			    if(exists $H_ID_ln{$uc_ID2})
			    {
					$intitule_sortie = ">Promot_for:$tempoID\n";
					# print "promoteur $tempoID line $H_ID_ln{$uc_ID2} prem-ligne correspond a \n$f_seq_prom[$H_ID_ln{$uc_ID2}]\n";
				    B2:for my $ln_fic_prom($H_ID_ln{$uc_ID2}..$#f_seq_prom)
				    {
					# die "$f_seq_prom[$ln_fic_prom]\n";
						if( $f_seq_prom[$ln_fic_prom] =~ /^>/ )
						{ 
						    # print "fin de memorisation ligne $ln_fic_prom\n";
						    $memoire_seq =~ s/\n//g;
						    print "$tempoID\n$memoire_seq\n";
						    &impression_f(\$memoire_seq,\$fct,\$nom_fic_sor,\$intitule_sortie,\$good_length, \$alphabet_ID_strict);
						    $prev_tempo = $tempoID;
						    $bool_1ID_trouve = 1;
						    last BOUCLEID2; 
						} 
						else
						{  
						    $memoire_seq .= $f_seq_prom[$ln_fic_prom];
						    # print "on memorise $f_seq_prom[$ln_fic_prom] line $ln_fic_prom\n";
						    if($memoire_seq !~ /^[acgtACGT\n]+$/){ die "$prog_tag [Error] Problem $memoire_seq contains other characters than ACGTacgt\\n, line ".__LINE__."\n";}
						}
				    }
			    }
			    else
			    {
					print "$uc_ID2 not found in genome IDs\n";
			    }
			} # fin de BOUCLEID2
			
			if(! $bool_1ID_trouve)
			{  
		      ($not_verbose)or print "No existing promoting sequence for $tempoID gene @ID2 fct $fct\n"; 
		      print NOEXISTING "No existing promoting sequence for $tempoID gene @ID2 fct $fct\n"; 
		      $prev_tempo = $tempoID;
			}
	  
		}
		# fin de if($line...
		else
		{
			die "$prog_tag [Error] $line\nTerm not recognized by our first regular expression in B1 loop\nreg_expr_rech vaut $reg_expr_rech, reg_expr_rech_apres $reg_expr_rech_apres for $line dol1 $1 dol2 $2, line ".__LINE__."\n";
		}
	} # fin de B1

}

close NOEXISTING;

# checking to verify we have two seq in a file and remove it otherelse
if($#fct == -1)
{
    print "Aucune fonction, nous verifions le dossier global\n";
    push @fct, '';
}
print "Directories created: @fct\n";
for my $dir(@fct)
{
    opendir(DIR,$chemin_dest.$dir)or die "$prog_tag [Error] Impossible to open ".$chemin_dest.$dir.":$!, line ".__LINE__."\n";
    my @f = grep !/^\./, readdir DIR;
    closedir(DIR);
    print "on verifie dans le repertoire $dir si on a deux seq par fichier\n";

    for my $f(@f)
    {

		open(F,'< '.$chemin_dest.$dir.$f)or die "$prog_tag [Error] Impossible to open $chemin_dest".$dir.$f.":$!, line ".__LINE__."\n";
		my $nb_seq = 0;
		while(my $l = <F> )
		{
		    if($l =~ /^>/g)
		    {
			$nb_seq++;
			print "on trouve $l, nb_seq vaut $nb_seq\n";
		    }
		}
		close F;
		if($nb_seq < 2)
		{
		    print "on lance rm -f $chemin_dest$dir$f car $nb_seq < 2\n";
		    system("rm -f $chemin_dest$dir$f");	
		}
		elsif($nb_seq == 2)
		{
		    print "$chemin_dest$dir$f comporte 2 seq, OK\n";
		}
		else
		{
		    die "$prog_tag [Error] $chemin_dest$dir$f comporte plus de 2 seq ($nb_seq), PB, line ".__LINE__."\n";
		}
    }
}

### SUB starts


sub impression_f(\$\$\$\$\$\$) 
{
    my ($r_memoire_seq, $r_fct, $r_nom_fic_sor, $r_intitule_sortie, $r_good_length, $r_alphabet_ID_strict) = @_;

	# nous verifions que la sequence promotrice est suffisamment longue pour nous interesser...
	if(length($$r_memoire_seq)>=$length_min + $dist_deb_trad) 
	{
		# si elle est trop longue par rapport a la taille max definie, nous la tronquons
		if (length($$r_memoire_seq)> $length_max + $dist_deb_trad)
		{
			$$r_memoire_seq = substr($$r_memoire_seq,length($$r_memoire_seq)-$length_max -$dist_deb_trad,$length_max);
		}
		else
		{
			$$r_memoire_seq = substr($$r_memoire_seq,0, length($$r_memoire_seq) - $dist_deb_trad);
		}
		
		# impression a la suite dans le fichier de sortie de la sequence promotrice
		if(length($$r_nom_fic_sor) > $MAX_FIC_NAME_LENGTH)
		{
			($not_verbose)or print "length ".length($$r_nom_fic_sor)."\n";
			$$r_nom_fic_sor =~ /__/; 
			my $nb_IDr = 0;
			my $r_str = $'; # '
			my $l_str = $`; # '
			while($r_str =~  /[$$r_alphabet_ID_strict]+/g ) # '
			{
				$nb_IDr++;
			}
			($not_verbose)or print "nb_IDr $nb_IDr dans $r_str\n";
		
			my $nb_IDl = 0;
			while($l_str =~ /[$$r_alphabet_ID_strict]+/g )
			{
				$nb_IDl++;
			}
			($not_verbose)or print "nb_IDl $nb_IDl dans $l_str\n";
			$$r_nom_fic_sor =~ /^([$$r_alphabet_ID_strict]+)/;
			my $new_name = $1.$separateur_ID.'et_al_'.$nb_IDl.'__';
			$$r_nom_fic_sor =~ /__([$$r_alphabet_ID_strict]+)/;
			$new_name .= $1.$separateur_ID.'et_al_'.$nb_IDr;
			$$r_nom_fic_sor = $new_name;
		}
		
		open(FSOR,'>> '.$chemin_dest.$$r_fct.$$r_nom_fic_sor)or die "$prog_tag [Error] Impossible to open $chemin_dest".$$r_fct.$$r_nom_fic_sor."\nMessage:$!, line ".__LINE__."\n";
		print FSOR $$r_intitule_sortie.&print50($bool_fasta,$$r_memoire_seq)."\n";
		close FSOR;
		if(length($$r_memoire_seq)<30)
		{
			die "$prog_tag [Error] trop court, ".length($$r_memoire_seq)." est theoriquement inf a  $length_min + $dist_deb_trad = ".( $length_min + $dist_deb_trad).":$!, line ".__LINE__."\n";
		}
		($not_verbose)or print "\nIn ".$chemin_dest.$$r_nom_fic_sor." 2 we write\n$$r_intitule_sortie".&print50($bool_fasta,$$r_memoire_seq)."\n";
		# nous memorisons que la longueur de cette sequence convient pour enregistrement
		$$r_good_length = 1;
	}
	else
	{
		($not_verbose)or print "\n".(length($$r_memoire_seq))." letters (<".($length_min + $dist_deb_trad)."), promoting sequence to short $$r_memoire_seq\n";
		# die "\n".(length($$r_memoire_seq))." letters (<".($length_min + $dist_deb_trad)."), promoting sequence to short $$r_memoire_seq\n";
		# nous memorisons que la sequence promotrice est trop courte pour enregistrement (evitera de refaire le test plus tard si d'autres relations d'orthologie concernent cette seq promotrice)
		$$r_good_length = 0;
	}
}## end of sub impression

1;
