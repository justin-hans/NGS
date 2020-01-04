#####################################
#####	README SCRIPT NGS	#####
#####################################

### adresse du git : https://github.com/justin-hans/NGS ###

Les 3 fichiers .sh attachés sont des scripts bash voués à l'analyse des données 
de séquençage issues de 3 individus ayant participé au projet 1000 génome.
Ces scripts sont prêts à l'emploi : leur utilisation requiert simplement la 
spécification d'un répertoire de travail. L'utilisateur doit donc saisir
le chemin d'accès au répertoire de travail dans l'espace adéquat (WORKDIR)
dans l'en-tête de chacun des 3 scripts, puis executer successivement les 3 scripts
pour procéder à l'analyse, dans l'ordre suivant :
	-> mapping.sh
	-> variant calling.sh
	-> trio-analysis.sh



###	Détail de l'analyse réalisée par les 3 scripts	###

# Objectif général #

Les données issues du programme 1000 génomes (et téléchargées par les scripts) 
sont les données de séquençage provenant de 3 individus d'une même famille
(père, mère et fille). Le but global de l'analyse est de détecter les variants (SNPs)
dans chacun de ces 3 individus. Ces variants sont des variations d'un seul nucléotide
(substitution, addition, délétion).
En pratique, l'analyse de variants chez de nombreux individus peut s'avérer utile dans le cadre
des études d'association entre maladies et SNps : c'est le cas de GWAS par exemple.
Cependant, pour des raisons techniques (temps de calcul), nous nous limitons ici à l'analyse 
des régions exoniques du chromosome 20.
Les 3 scripts utilisent des fonctions issues de bibliothèques dédiées à la bioinformatique,
ert plus précisément à l'analyse de séquence de séquençage. Il s'agit des suites GATK
et Picard, développées par le Broad Institute, et en accès libre (https://software.broadinstitute.org/gatk/)
(https://broadinstitute.github.io/picard/). Une documentation détaillée est notamment disponible sur ces sites.

# mapping.sh #

Le script débute par le téléchargement d'un chromosome de référence, et des données de séquençage 
pour les individus étudiés (directement sur le FTP du projet 1000 génomes : les données sont en open access).
Les données de séquençage, obtenus grâce à la technique illumina, sont de nombreux courts segments (reads) qui
correspondent à une petite partie du génome de l'individu. Ces reads se chevauchent : ainsi, nous disposons de 
l'ensemble de la séquence de l'individu diploïde. On peut spécifier le nombre de runs de séquençage à télécharger
et analyser (par défaut, NUMBER_RUNS=8). Un plus grand nombre de runs correspond à un plus grand nombre de reads,
et donc une meilleure définition/précision, mais aussi à un plus grand temps de téléchargement et d'analyse.
La première étape consiste donc à aligner les différents reads sur le chromosome de référence :

- Après une étape d'indexation du chromosome, on récupère le nom des reads à télécharger à partir d'un 
fichier index obtenu sur le FTP. La variable NUMBER_RUNS spécifie le nombre de runs à obtenir.
- Les reads sont ensuite mappés, filtrés et ordonnés sur le chromosome (bwa mem, samtools view, samtools sort).
Ces 3 fonctions sont exécutées directement à la suite grâce au pipe (|) : en supprimant les étapes intermédiaires,
on réduit le temps de calcul et on optimise la procédure.
-> Après une étape de merge et d'indexation, on obtient alors un fichier .bam contenant l'ensemble des reads alignés 
sur le chromosome de référence.

# variant-calling.sh #

Il s'agit à proprement parler de l'étape de détection de variants. On commence par télécharger une 
liste des variants et indels connus, disponibles sur le site du projet 1000 génomes eux-aussis.

- Après une étape d'indexation et la  création d'un dictionnaire de séquence, on marque les dupliquats
parmis les reads.
- Les éventuels indels présents dans les reads perturbent leur alignement sur le chromosome :
ils provoquent la détection de variants qui n'ont pas lieu d'être (en effet, d'un côté ou de l'autre de
l'indel, l'alignement du read est impossible). Il est donc nécessaire de localiser ces indels pour
effectuer un alignement de bonne qualité. Pour cela, on utilise les fonctions RealignerTargetCreator
et IndelRealigner (détection des cibles, puis réalignement effectif).
- Vient alors l'étape de Base quality recalibration. Les séquenceurs utilisés (ILLUMINA en l'ocurrence)
fournissent une séquence, ainsi que des indices de confiance associés à chaque base séquencée. Ces indices de 
confiance sont particulièrement importement lors de l'étape de recherche de variants. Cependant, la technique de
séquençage utilisée est source d'une erreur technique non-aléatoire, qui impacte directement l'indice de confiance
pour chaque base, notamment en fonction de la position dans le read, de l'environnement (séquence, bases adjacentes...).
Il est donc nécessaire de corriger cette erreur systématique, au moyen d'un algorithme de machine learning : BaseRecalibrator
- Enfin, les fonctions HaplotypeCaller et GenotypeGVCFs permettent la recherche effective de variants, qui sont consignés
dans un fichier .vcf. Ce fichier vcf peut notamment être ouvert dans IGV pour visualiser les variants au sein de 
la séquence de référence.

# trio-analysis.sh #

Dans l'étape précédente, nous avons mappé les reads séquencés chez 3 individus sur le génome de référence. 
Cependant, nous avons traité les 3 individus de façon indépendant, et n'avons pas tiré profit du lien de
parenté qui unit père, mère et fille. Ce lien peut notamment être eploité pour réaliser une détection plus
fiable des variants, en permettant de discriminer les variants hérités des variants nouvellement apparus.
Pour ce faire, nous utilisons notamment un fichier PEDIGREE, qui contient les informations relatives
aux liens généalogiques qui unissent les sujets du projet 1000 génomes. Nous réalisons ensuite une recherche
de variants "jointe" qui prend en compte les liens généalogiques, en utilisant la fonction PhaseByTransmission.
Afin de déterminer l'impact de la prise en compte du pedigree sur l'analyse de variants, on compare ensuite
les fichier vcf avant et après le phasing par PhaseByTransmission, en utilisant GenotypeConcordance.
On obtient alors un fichier texte qui compile les différences entre les 2 types d'analyse.



