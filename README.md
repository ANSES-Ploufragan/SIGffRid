# SIGffRid

### Meaning of program name

SIGffRid means:
$$SIG$$ma $$f$$actor (binding sites) $$f$$inder using
$$R$$'MES to select $$i$$nput $$d$$ata.

> SIGffRid name is mainly a recall to __SigR__ _Streptomyces coelicolor_
> Sigma factor binding sites, well known and easily detectable
> even using counts instead of statistics in SIGffRid for
> motif extension. 
> This first "confirmation" helped me to argue to develop
> the statistical extension criteria needed for SIGffRid
> to find other kinds of sigma factors.

## PRESENTATION

SIGffRid try to __deduce Sigma Factor Binding Sites__ using 
two phylogenetically close bacterial genomes (typically, 
showing 16S RNA identity of 97%, with the highest possible number 
of orthologous sigma factors).

SIGffRid runs on UNIX/Linux systems (tested on
Ubuntu 16.04 LTS and [Solus 4.1](https://getsol.us/home/)). 

## PRINCIPLE

SIGffRid performs a simultaneous analysis of pairs of promoter regions of orthologous genes in Bacterial species. SIGffRid uses a prior identification of slighly over-represented 3 to 7 letters long words in the whole genomes of the two bacterial genomes, as selection criteria for potential -35 and -10 boxes using R'MES program.

From each pair of intergenic upstream sequences of orthologs (obtained from MBGD database), it get pairs of conserved (individually over-represented) words with close spacers on the two sequences: patterns. These patterns are then grouped independantly for each bacterium using pairs of short seeds (of which one is possibly gapped), allowing a variable-length spacer between them (one seed matching to one single word of the pattern). Conserved nucleotides in seeds concerned 6 positions (from 2 to 4 for an individual seed composing a pair).

Next, the conserved area corresponding to seed pair is extended guided by statistical considerations, a feature that ensures a selection of motifs with statistically relevant properties.The extension process is based on a Markov models of third order obtained from each whole genome and concerns the position flanking conserved area related to the letter whose the observed number of occurrences is the less expected in the n sequences concerned. If there are several statistically significant letters at this position, the extension is tried using several letters simultaneously.
	In this case, we obtain a motif with a position corresponding to several letters.If this motif is statistically significant, we stop the extension process: we have the consensus and all the motifs and positions in upstream sequences related to this motif. Otherwise, extension is tried on each statisitcally significant letter at this position. Sequences who are not concerned are grouped so as to run a distinct extension process.
	Entension process stop when the motif was extended by 10 letters without interesting results or when less than 8 sequences are concerned.

After each extension by one letter (at one of the 4 positions flanking seeds), the conserved motif is evaluated.
We compute a ratio R of the number of occurrences on upstream sequences on the number of occurrences in whole genome (taking into account the two strands).The higher is R, the more specific is the motif for upstream sequences.

To know if this specificity is statistically significant, we use  Likelyhood Ratio Test (LRT) which measure the difference of concentration of our motif in the two sets of sequences (test made on data corresponding to the bacteria the motif is got from).
If R is greater than a threshold computed according to the proportion of upstream sequences in the studied bacteria, and LRT test gives a score greater than the quartile at 5 percent of the Xhi square law, the motif is considered as an interesting one.

Important Note:
      We have only used a minimum value threshold for LRT test, which decrease greatly the number of results given.
Using our program on some bacteria, we suggest not to consider the motifs with higher LRT score as the best candidates. Known sigma factor binding sites have quite low LRT scores (from 20 to 50 in most cases) but with the highest R.


## REQUIREMENT

SIGffRid 1.0 (command line version) needs:
- Perl (tested with version >= 5.8.8)

install YAML

using CPAN (cpan install <module>)
- XML::DOM::XPath
(you probably will need to run ```sudo cpan -f -i XML:DOM::XPath``` (to force install) to avoid some problem with BioPerl installation (due to testing failures))
- BioPerl
- SeqIO                (to extract data from GenBank files)
- IPC::System::Simple  (to run system command)
- PerlIO::gzip         (to read write zip files)
- Proc::Daemon
- Getopt::Long         (for options)
- Cwd                  (for current directory)
- Statistics::Basic

install rmes-3.1:
- git clone https://forgemia.inra.fr/sophie.schbath/rmes.git
(with rmes.fomat)
- R'MES documentation is provided here: [rmes_doc](http://migale.jouy.inra.fr/?q=rmes)


add to your .bashrc
```
export SIGFFRID_DIR=path_to_your_sigffrid_global_directory
export ANSES_SCRIPTS_DIR=path_to_your_sigffrid_global_directory
export ANSES_USER_MAIL=fabrice.touzain@anses.fr
```
and if using **SGE** (replace long by the name of the SGE queue in which you want to submit jobs)
```
export ANSES_Q_NAME='long.q'
```

or if using **SLURM** (replace long by the name of the slurm queue in which you want to submit jobs)
```
export ANSES_SLURM_Q_NAME='long'
```

then type
```
source ~/.bashrc
```

## WARNINGS

### Memory usage

For a pair of bacterial genomes, __SIGffRid requires at 
least 16 Gbytes (RAM)__ (for genomes of 4 Mbases). The memory 
usage probably can go from 6 to 16 Gbytes for bacterial genomes
if you treat simultaneously 
all the pairs of orthologues found in genomes as long as 8 
Mbases (typically *Streptomyces coelicolor* and 
*Streptomyces avermitilis*). In such a case, use 
"functional category" option (-use_fct) which groups sequences of 
pairs of orthologues by functional categories as defined 
in MBGD files.
By using this option for big genomes, number of results 
(and false positives) will decrease, and memory usage too.

This option ("grouping by functional category") is NOT 
recommanded to treat genomes smaller than 5 Gbases: you 
will obtain more relevant results by using all orthologue 
relationships (in all cases).

Time consuming:

__In all cases__, __SIGffRid will give results after 1 to 4 
days__ (depending on "grouping by functional category" option 
, on number of seeds used) (100% CPU on a Intel core i7 8700K at 
3.7 GHz, ssd NMVE).

**Cluster mod 'SLURM' NOT TESTED YET** (tested only with modes 
'SGE' (Ubuntu) and 'no' (Solus).

### Option -use_fct

This option groups pairs of orthologs by functional categories
as defined in MBGD:
- Amino acid biosynthesis
- Purines, pyrimidines, nucleosides, and nucleotides
- Fatty acid, phospholipid and sterol metabolism
- Biosynthesis of cofactors, prosthetic groups, and carriers
- Central intermediary metabolism
- Energy metabolism
- Transport and binding proteins
- DNA replication, restriction, modification, recombination, and repair
- Transcription
- Translation
- Regulatory functions
- Cell envelope
- Cellular processes
- Other categories
- Hypothetical

It decreases resources needs but restrict comparison to each category. 


## INPUTS

- Ids of the two bacterial species:
VERY IMPORTANT:
__Bacterial IDs MUST be__ the __IDs used to identify the genome 
in MBGD file!!!!!__ (file called cluster when you download 
it, see 3 items below). We call __bacterial ID__ the __group of 
letters at the beginning of a gene ID__ (for example __sco__ for 
__Streptomyces coelicolor__ whose first gene locus is __SCO__0001).
CAUTION: In some cases, genes IDs provided in MBGD are not 
the same as in the ncbi GenBank file.

- Working directory: the directory where all the outputs 
and other intermediary files generated by SIGffRid program 
will be stored. Each working directory has a log file 
associated with it that stores information about the 
processes executed. It is a good advice to __start SIGffRid__ 
program on __a new dataset__ in __a new working directory__.

- Genomic data: For each input species the user needs to 
input the genome file (GenBank format) for the whole 
bacterial genome (obtained from NCBI website at 
http://www.ncbi.nlm.nih.gov/genomes/static/eub_g.html). 
The input will be the sequence file (.gb/.gbk format).
CAUTION: The .gb file MUST contain the sequence (default, 
not). Therefore, when the gb is display in your browser, 
check the box on the right site of the ncbi interface 
"show sequence" and then "update view". Afterwhat you can 
record the gb file by clicking on
"Send",
check "complet record",
check "file",
select "GenBank full".


- Orthologous relationship file (MBGD file): This file 
is the file that contains information about the orthologous 
pairs of genes on the given set of species. This file can 
be obtained from the MBGD database (http://mbgd.genome.ad.jp/) 
or can be user-defined file. Such file can help to restrict 
the number of genes in analysis by using the user chosen 
genes in finding sigma factor binding sites.
HELP TO GET ORTHOLOGOUS RELATIONSHIP AND FUNCTIONAL 
CATEGORY CLUSTERING FROM MBGD DATABASE:
1. Go to MBGD homepage http://mbgd.genome.ad.jp/.
2. Click on "Create Ortholog Table" button in the left column (send 
you to http://mbgd.genome.ad.jp/htbin/SelectOrganism.pl page).
3. Check boxes related to the two species to compare in the displayed 
list (multiple selection or use the link to "taxonomy 
browser").
4. Click on "Next" button at the top of the page.
5. In the new page, Click on "Create/View Cluster Table" button.
6. On the web page displayed, click on "Save Complete 
Table" button.
The file obtained (default named as:'cluster.tab') is the 
one you have to give to SIGffRid as parameter to define 
orthologous relationships (you can rename it, its name 
has no importance).

__PLEASE, CHECK cluster file after download__ (sometimes, 
database record does not work and cluster file has only 
the header line).


## LAUNCH

Type in terminal
```
perl sigffrid_cmd_line.pl
```
and enter.


## RESULTS

They consist in
- motif consensus, occurrences, and upstream sequences 
where motifs are found,
- annotation of genes concerned by each motif,

Located in:
{Working_directory}/{your_id}/SIGffRid_orthologs/  directory

Listing of motifs sorted by score are in files:
- res_sort_R_for_results_grep_{bacteria2_id}.txt
- res_sort_R_for_results_grep_{bacteria1_id}.txt

You will find all motif occurrences and sequences they 
are extracted from in __SIGffRid_results/MOTIF_UPSTR_treshmulti3_bsu_8/__ directory, in the
__SIGffRid_results/MOTIF_UPSTR_treshmulti3_bsu_8/motif{bacterial ID}_{xxxxx}.txt__ file.

You will find annotations for genes related to a motif xxxxx in
- MOTIF_ANNOT_UPSTR       directory, in the
- MOTIF_ANNOT_UPSTR/annot_motif{bacterial_ID}_{xxxxx}.txt file



Advanced users (using command lines) can know how many 
times a same motif is found by the program (with different 
seeds). This information, not provided by default can be used 
to select some interesting motifs (in addition to scores).
For a complete list of motifs, type:
```
grep 'MOTIF' align* 
```
in SIGffRid_orthologs subdirectory if you do not use "grouping 
by functional category" option, otherwise type
```
grep 'MOTIF' */align*
```

## REFERENCE

If you publish or use SIGffRid results, please cite:

Touzain F, Schbath S, Debled-Rennesson I, Aigle B, Kucherov G and 
Leblond P. **SIGffRid: a tool to search for sigma factor binding 
sites in bacterial genomes using comparative approach and 
biologically driven statistics** (article). *BMC Bioinformatics* 2008, 
9(1):73. [DOI:10.1186/1471-2105-9-73](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-73)                   

Questions and comments should be sent to *fabrice.touzain [@] anses.fr*

## FINANCIAL SUPPORT

Fabrice Touzain was supported by the __Région Lorraine__ and by the __ACI IMPBio__ from the 
__Ministère de l'Education Nationale, de l'Enseignement supérieur et de la Recherche__.

An additional __demi-ATER__ was usefull to finish this program and Fabrice Touzain
PhD. Thank you again to the __Région Lorraine__.

<img src="img/logo_loria_complet.jpg" width="150">
<img src="img/UL_LOGO-300x133.jpg" width="150">
<img src="img/LOGO_CNRS_2019_RVB-150x150.png" width="70">
<img src="img/inr_logo_rouge_150.jpg" width="150">
<img src="img/logo_ANSES.png" width="150">


## ACKNOWLEDGMENT

Enormous *THANKS* to Professor __Isabelle Debled-Rennesson__ ([home](https://members.loria.fr/IDebledRennesson/files/)) for her constant
guidance and encouragement, to __Sophie Schbath__ ([home](http://genome.jouy.inra.fr/~schbath/)) for her statistical explanations and pedagogy to design the method,
to __Gregory Kucherov__ for the three hundred kilometers he drived under the snow to make my PhD defense possible in 2007.

Special thanks to PhD who shared by office for their help and kindness: Laurent Noé ([home](https://sites.google.com/view/laurentnoe/home)) and Laurent Provot ([home](http://lprovot.fr/)). 

## LICENCE

Terms of the CeCILL licence can be found here: [CeCILL](http://gitlab.gvb.anses.fr/touzain/sigffrid_cmd_line/blob/master/LICENSE.txt)
