#!/bin/bash
#$ -cwd

#script_genresc.sh <d> <NINDTP> <NINDAL> <GEN> <tR> <MIG> <INT>
#To be used only in a machine with /state/partition1 directory

rm script_genresc.sh.*

##########################################################################################
# Outline                                                                                #
#                                                                                        #
#           __ 'NINDTP'R'NINDTP' (line with same size as threatened population and GR)   #
#          |                                                                             #
#          |__ 'NINDTP'R'NINDAL' (line with different size and GR)                       #
# NINDTP---|                                                                             #
#          |__ 'NINDTP'NR'NINDAL' (line with different size and without GR)              #
#          |                                                                             #
#          |                  __ 'NINDTP'RR'NINDTP'                                      #
#          |                 |                       (line with same size as threatened  #
#          |__ 'NINDTP'NR  __|__ 'NINDTP'RR'NINDAL'  population without GR, and          #
#                            |                       subsequent formation of sublines)   #
#                            |__ 'NINDTP'NR'NINDTP'                                      #
#                            |                                                           #
#                            |__ 'NINDTP'NR'NINDAL'                                      #
#                                                                                        #
##########################################################################################


####################################### ARGUMENTS ######################################

#Check number of arguments
if [ $# -ne 7 ]  
then
	echo "Usage: $0 <d> <NINDTP> <NINDAL> <GEN> <tR> <MIG> <INT>" 
	exit 1
fi

#Set arguments
d=$1
NINDTP=$2
NINDAL=$3
GEN=$4
tR=$5
MIG=$6
INT=$7
#REPS=$8

############################### VARIABLES AND DIRECTORIES ##############################

#Parameters
LAMB=0.2
NCRO=300
AVEs=0.2
AVEh=0.283
Vs=0
Opt=0
NR=0

m_TP=5
m_AL=1
dist=0

#Working directory
WDIR=$PWD 
mkdir -p $WDIR/genresc_results/L$LAMB.k$NCRO.s$AVEs.h$AVEh/TP$NINDTP.AL$NINDAL.M$MIG.INT$INT.g$tR
DIR="genresc_results/L$LAMB.k$NCRO.s$AVEs.h$AVEh/TP$NINDTP.AL$NINDAL.M$MIG.INT$INT.g$tR"

#Scratch directory
mkdir -p /state/partition1/noeliaGR$d/$SLURM_JOBID/
SCDIR="/state/partition1/noeliaGR$d" 

############################# TRANSFER OF FILES TO SCRATCH #############################

#Copy all files in scratch directory
cp seedfile $SCDIR/$SLURM_JOBID/
cp genresc $SCDIR/$SLURM_JOBID/
cp POPFILE_L$LAMB.K$NCRO.s$AVEs.h$AVEh $SCDIR/$SLURM_JOBID/popfile
cp DATAFILE_L$LAMB.K$NCRO.s$AVEs.h$AVEh $SCDIR/$SLURM_JOBID/datafile

#File with information of node and directory
touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`

#Move to scratch directory
cd $SCDIR/$SLURM_JOBID

####################################### GENRESC ########################################

START=$(date +%s)
time ./genresc>>out<<@
0
-99
$NINDTP	NIND Threatened Population
$NINDAL	NIND Alternative Lines
$m_TP	Migrants for lines with NINDTP individuals
$m_AL	Migrants for lines with NINDAL individuals
0	Gender of Migrants (males 0, males&females 1)
99	Lenght genome (99=free)
$NCRO	NCRO (max 2000)(Neu=Ncro)
30	NLOCI (2-30)
$LAMB	Lambda_s
0.0	Lambda_L
0.33	beta_s
$AVEs	ave |s|
1	dom model (0=cnt; 1:variable)
$AVEh	ave h
$Vs	Stabilizing selection (Vs)
1	VE
$Opt	Optimal
$GEN	generations
$tR	Formation of R and NR lines (max 'generations')
$tR	Formation of sublines ('tLINES' - 'generations')
$NR	For R lines, initial NR-generations after formation of line
$MIG	Number of migrations (99: periodic)
$INT	Generation intervals of migration (lines NINDTP)
$INT	Generation intervals of migration (lines NINDAL)
$dist	Generations since migration (distribution -w- among replicates)
0.0	Minimum fitness value for extinction (99: without extinction)
1000	Replicates
@

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "genresc took 		$DIFF seconds" >> timefile

###################### TRANSFER OF FILES TO MAIN DIRECTORY ############################

cp -r $SCDIR/$SLURM_JOBID/seedfile $WDIR
cp -r $SCDIR/$SLURM_JOBID/timefile $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/dfilename.dat $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/genfile_L1.dat $WDIR/$DIR/genfileL1_"$NINDTP"R"$NINDTP".dat
cp -r $SCDIR/$SLURM_JOBID/genfile_L2.dat $WDIR/$DIR/genfileL2_"$NINDTP"R"$NINDAL".dat
cp -r $SCDIR/$SLURM_JOBID/genfile_L3.dat $WDIR/$DIR/genfileL3_"$NINDTP"NR"$NINDAL".dat
cp -r $SCDIR/$SLURM_JOBID/genfile_L4subl1.dat $WDIR/$DIR/genfileL4sub1_"$NINDTP"RR"$NINDTP".dat
cp -r $SCDIR/$SLURM_JOBID/genfile_L4subl2.dat $WDIR/$DIR/genfileL4sub2_"$NINDTP"RR"$NINDAL".dat
cp -r $SCDIR/$SLURM_JOBID/genfile_L4subl3.dat $WDIR/$DIR/genfileL4sub3_"$NINDTP"NR"$NINDTP".dat
cp -r $SCDIR/$SLURM_JOBID/genfile_L4subl4.dat $WDIR/$DIR/genfileL4sub4_"$NINDTP"NR"$NINDAL".dat
cp -r $SCDIR/$SLURM_JOBID/distribution_qsh.dat $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/summary_outline.dat $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/out $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/rescfile.dat $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/extinction.dat $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/distribution_w.dat $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/repfile.dat $WDIR/$DIR/

############################# CLEANING OF SCRATCH ####################################

rm -r $SCDIR/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*


