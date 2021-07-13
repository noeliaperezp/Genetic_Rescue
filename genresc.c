/* genresc.c */

/* ********************************************************************* */

#include "libhdr"
#include "ranlib.h"
#define NN 1501  /* max number of NIND and gen */
#define MM 2001  /* max number of NCRO */
#define MI 201  /* max number of migrations */
#define maxmpt 5001
#define normalthreshold 30  

int NIND, NINDNP, NINDTP, NINDAL, N, NCRO, NLOCI, NSEGLOCNP, TOTLOCI;
int NMIG, NMIGTP, NMIGAL, mGEND, numMIG, migrations, migration_events[10], exp_migevents[10], rescue[NN], events;
int gen, tgen, generations, tLINES, tsubL, preNR, endmig, mINT, mINT_TP, mINT_AL, end, init, dist;
int lin, lines, sublines, subl, i, k, l, rep, replicates, dom, mig;
int lastinmutantspoissontable, lastinmutantspoissontableL, lastinrecombinantpoissontable;
int father[NN], mother[NN], sfath[NN], smoth[NN], mark[10001];
int RM[NN], gNP[10001][MM][2], g[NN][MM][2], sgnt[NN][MM][2], sp[NN][MM][2];
int initialgen[MM][31], initialgenNP[MM][31], sinitgen[MM][31], seg_markNP[MM][31];
int rep_extinct[10][NN], repmark;

double L, Lambda_s, Lambda_L, Vs, VE, VA, sVA, OPT, females[NN], males[NN], wEXT;
double ave_s, addedfixed_sNP, addedfixed_sf, alpha_s, beta_s, k_s, ave_hs, AA, Aa, aa, LEQ_NP, TLEQ_NP, LEQ_L_NP, LEQ_NL_NP;
double q[MM][31], qs[MM][31], s[MM][31], ss[MM][31], sNP[MM][31], hs[MM][31], hss[MM][31], hsNP[MM][31];
double smat[NN][NN], ssmt[NN][NN], mat[NN][NN], F[NN], sF[NN], pm_sf[NN], spm_sf[NN], p_a[NN], sp_a[NN];
double mutantspoissontable[maxmpt], mutantspoissontableL[maxmpt], recombinantpoissontable[maxmpt];
double STR_m[MI], STR_l[MI], summ_dq_l[MI], summ_dq_m[MI], mean_w[MI];

struct acc gmean_sf[10][NN], gvar_sf[10][NN], gmean_a[10][NN], gvar_a[10][NN], SK2[10][NN], Hw[10][NN], fp[10][NN], Fp[10][NN];
struct acc NSEGLOCm[10][NN], NSEGLOCl[10][NN], sm[10][NN], sl[10][NN], hm[10][NN], hl[10][NN], qm[10][NN], ql[10][NN];
struct acc AAm[10][NN], Aam[10][NN], aam[10][NN], AAl[10][NN], Aal[10][NN], aal[10][NN];
struct acc LEQ_Lf[10][NN], LEQ_NLf[10][NN], LEQ_T[10][NN];
struct acc mutNP[10][NN], mutShared[10][NN], mutLine[10][NN], mut[10][NN], mutNP_m[10][NN], mutm[10][NN];
struct acc lostNP[10][NN], lostLine[10][NN], lost[10][NN], lostNP_m[10][NN], lostm[10][NN];
struct acc mate_mm[10][MI], mate_ml[10][MI], mate_ll[10][MI], FEC_mig[10][NN], FEC_line[10][NN];
struct acc summ_w[10][MI], ext_summ_w[10][MI], sur_summ_w[10][MI];
struct acc gload_line[10][MI], gload_mig[10][MI], ext_gload_line[10][MI], ext_gload_mig[10][MI], sur_gload_line[10][MI], sur_gload_mig[10][MI];
struct acc STR_line[10][MI], STR_mig[10][MI], ext_STR_line[10][MI], ext_STR_mig[10][MI], sur_STR_line[10][MI], sur_STR_mig[10][MI];

/* ******************* Distribution of contributions ******************* */

struct acc off_0_02[10][NN], off_02_04[10][NN], off_04_06[10][NN], off_06_08[10][NN], off_08_1[10][NN];

/* ********************** Distribution of fitness ********************** */

struct acc w_00[10][101], w_00_01[10][101], w_01_02[10][101], w_02_03[10][101], w_03_04[10][101], w_04_05[10][101];
struct acc w_05_06[10][101], w_06_07[10][101], w_07_08[10][101], w_08_09[10][101], w_09_10[10][101], w_10[10][101];

struct acc ext_w_00[10][101], ext_w_00_01[10][101], ext_w_01_02[10][101], ext_w_02_03[10][101], ext_w_03_04[10][101], ext_w_04_05[10][101];
struct acc ext_w_05_06[10][101], ext_w_06_07[10][101], ext_w_07_08[10][101], ext_w_08_09[10][101], ext_w_09_10[10][101], ext_w_10[10][101];

struct acc sur_w_00[10][101], sur_w_00_01[10][101], sur_w_01_02[10][101], sur_w_02_03[10][101], sur_w_03_04[10][101], sur_w_04_05[10][101];
struct acc sur_w_05_06[10][101], sur_w_06_07[10][101], sur_w_07_08[10][101], sur_w_08_09[10][101], sur_w_09_10[10][101], sur_w_10[10][101];

/* ************************ Extinct replicates ************************* */

struct acc wext_00[10], wext_00_01[10], wext_01_02[10], wext_02_03[10], wext_03_04[10], wext_04_05[10];
struct acc wext_05_06[10], wext_06_07[10], wext_07_08[10], wext_08_09[10], wext_09_10[10], wext_10[10];

struct acc STR_00[10], STR_00_01[10], STR_01_02[10], STR_02_03[10], STR_03_04[10], STR_04_05[10];
struct acc STR_05_06[10], STR_06_07[10], STR_07_08[10], STR_08_09[10], STR_09_10[10], STR_10[10];

struct acc dq_00[10], dq_00_01[10], dq_01_02[10], dq_02_03[10], dq_03_04[10], dq_04_05[10];
struct acc dq_05_06[10], dq_06_07[10], dq_07_08[10], dq_08_09[10], dq_09_10[10], dq_10[10];

/* *********************** Distribution of genes *********************** */

struct acc q_ng_00[10][101], q_ng_00_01[10][101], q_ng_01_02[10][101], q_ng_02_03[10][101];
struct acc q_ng_03_04[10][101], q_ng_04_05[10][101], q_ng_05_06[10][101], q_ng_06_07[10][101];
struct acc q_ng_07_08[10][101], q_ng_08_09[10][101], q_ng_09_10[10][101], q_ng_10[10][101];

struct acc s_ng_0[10][101], s_ng_0_106[10][101], s_ng_106_104[10][101], s_ng_104_102[10][101], s_ng_102_101[10][101];
struct acc s_ng_01_02[10][101], s_ng_02_04[10][101], s_ng_04_06[10][101], s_ng_06_10[10][101], s_ng_10[10][101];
struct acc s_ng_0_fix[10][101], s_ng_0_106_fix[10][101], s_ng_106_104_fix[10][101], s_ng_104_102_fix[10][101], s_ng_102_101_fix[10][101];
struct acc s_ng_01_02_fix[10][101], s_ng_02_04_fix[10][101], s_ng_04_06_fix[10][101], s_ng_06_10_fix[10][101], s_ng_10_fix[10][101];
struct acc s_ng_0_lost[10][101], s_ng_0_106_lost[10][101], s_ng_106_104_lost[10][101], s_ng_104_102_lost[10][101], s_ng_102_101_lost[10][101];
struct acc s_ng_01_02_lost[10][101], s_ng_02_04_lost[10][101], s_ng_04_06_lost[10][101], s_ng_06_10_lost[10][101], s_ng_10_lost[10][101];

struct acc hs_ng_00_02[10][101], hs_ng_02_04[10][101], hs_ng_04_06[10][101], hs_ng_06_08[10][101], hs_ng_08_10[10][101];
struct acc hs_ng_00_02_fix[10][101], hs_ng_02_04_fix[10][101], hs_ng_04_06_fix[10][101], hs_ng_06_08_fix[10][101], hs_ng_08_10_fix[10][101];
struct acc hs_ng_00_02_lost[10][101], hs_ng_02_04_lost[10][101], hs_ng_04_06_lost[10][101], hs_ng_06_08_lost[10][101], hs_ng_08_10_lost[10][101];


FILE *fptr, *fgenL1, *fgenL2, *fgenL3, *fgensubL1, *fgensubL2, *fgensubL3, *fgensubL4, *fdis, *fdisw, *frep, *fdat, *fpop, *fsum, *fmate, *frepext; 

/* ********************************************************************* */

main()
{
	fgenL1 = fopen ("genfile_L1.dat","w");
	fgenL2 = fopen ("genfile_L2.dat","w");
	fgenL3 = fopen ("genfile_L3.dat","w");
	fgensubL1 = fopen ("genfile_L4subl1.dat","w");
	fgensubL2 = fopen ("genfile_L4subl2.dat","w");
	fgensubL3 = fopen ("genfile_L4subl3.dat","w");
	fgensubL4 = fopen ("genfile_L4subl4.dat","w");

	fdis = fopen ("distribution_qsh.dat","w");
	fdisw = fopen ("distribution_w.dat","w");
	fsum = fopen ("summary_outline.dat","w");
	frepext = fopen ("extinction.dat","w");
	fmate = fopen ("rescfile.dat","w");

	getinputs();
	sum_outline();

	if (tracelevel!=0) 	fptr = fopen ("dfilename.dat","w");

	recombination_masks();
	natural_population();

	tgen = 10;
	lines = 9;

	for (rep=1; rep<=replicates; rep++)
	{
		frep = fopen ("repfile.dat","a");
		fprintf (frep,"replicate %d\n", rep);
		fclose(frep);

		if (tracelevel!=0) fprintf (fptr,"\n***********************\n***** Replicate %d *****\n***********************\n", rep);

		for (lin=0; lin<lines; lin++)
		{
			exp_migevents[lin] = 0;
			migration_events[lin] = 0;
			repmark = 0;

			//THREATENED POPULATION
			if (lin==0)	TP();

			//LINES
			else if (lin==1)	L1(); /* -NINDTP- R -NINDTP- */
			else if (lin==2)	L2(); /* -NINDTP- R -NINDAL- */
			else if (lin==3)	L3(); /* -NINDTP- NR -NINDAL- */
			else if (lin==4)	L4(); /* -NINDTP- NR */

			//SUBLINES
			else if (lin==5)	subL1(); /* -NINDTP- RR -NINDTP- */
			else if (lin==6)	subL2(); /* -NINDTP- RR -NINDAL- */
			else if (lin==7)	subL3(); /* -NINDTP- NR -NINDTP- */
			else if (lin==8)	subL4(); /* -NINDTP- NR -NINDAL- */

			events = 0;
			for (mig=1; mig<=exp_migevents[lin]; mig++)	mean_w[mig] = 0.0;

			for (gen=init; gen<end; gen++)
			{
				parents();
				if (rescue[gen]==1)
				{
					if (wEXT!=99.0)	
					{				
						//males as migrants
						if ((mGEND==0) && (females[gen]<=wEXT))
						{
							if (tracelevel!=0)	fprintf(fptr,"\n** Females no available **");
							rep_extinct[lin][gen] += 1;
							repmark = 1;

							if (((lin>=1)&&(lin<=3)) && (gen!=tLINES)) extinction_w();
							else if (((lin>=5)&&(lin<=8)) && (gen!=tsubL))  extinction_w();

							goto labelLINE;
						}
						//males and females as migrants (NMIG > 1)
						else if ((mGEND==1) && ((females[gen]<=wEXT)&&(males[gen]<=wEXT)))
						{
							if (tracelevel!=0)	fprintf(fptr,"\n** Parents no available **");
							rep_extinct[lin][gen] += 1;
							repmark = 1;

							if (((lin>=1)&&(lin<=3)) && (gen!=tLINES)) extinction_w();
							else if (((lin>=5)&&(lin<=8)) && (gen!=tsubL))  extinction_w();

							goto labelLINE;
						}
						//males and females as migrants (NMIG = 1)
						else if ((mGEND==1) && (NMIG==1) && (males[gen]<=wEXT))
						{
							if (tracelevel!=0)	fprintf(fptr,"\n** Males no available **");
							rep_extinct[lin][gen] += 1;
							repmark = 1;

							if (((lin>=1)&&(lin<=3)) && (gen!=tLINES)) extinction_w();
							else if (((lin>=5)&&(lin<=8)) && (gen!=tsubL))  extinction_w();

							goto labelLINE;
						}
					}
					migration_events[lin] += 1;
					if (tracelevel!=0)	fprintf (fptr,"\n***** Generation %d - migration%d - *****\n", gen, migration_events[lin]);
					NIND = NIND + NMIG;
					migration();
				}
				else
				{
					if(wEXT!=99.0)
					{
						if ((females[gen]==0.0)||(males[gen]==0.0) || (females[gen]<wEXT)||(males[gen]<wEXT))
						{
							if (((females[gen]==0.0)||(females[gen]<wEXT)) && (tracelevel!=0))	fprintf(fptr,"\n** Females no available **");
							if (((males[gen]==0.0)||(males[gen]<wEXT)) && (tracelevel!=0))	fprintf(fptr,"\n** Males no available **");
							rep_extinct[lin][gen] += 1;
							repmark = 1;

							if (((lin>=1)&&(lin<=3)) && (gen!=tLINES)) extinction_w();
							else if (((lin>=5)&&(lin<=8)) && (gen!=tsubL))  extinction_w();

							if ((lin==0)||(lin==4))	goto labelREP;
							else	goto labelLINE;
						}
					}
					if (rescue[gen]==2)	migration_events[lin] += 1; /* generation control */
					if (tracelevel!=0)	fprintf (fptr,"\n***** Generation %d *****\n", gen);
				}

				poisson_tables();
				mutation_neutral();
				mutation_selected();
				neutral_genes();
				selected_genes();

				if (rescue[gen]==1)	summary_genresc();
				if (rescue[gen]==2)	summary_NR();
//				if (tracelevel!=0)	dumpoffspringaftermutation();

				genotypic_values();

				if (gen%tgen == 0)	distribution_qsh();
				if (tracelevel!=0)	dumpphenotypes();
				if (Vs!=0.0)		additivevariance();

				coancestry_matrix();
				mating();

//				if (tracelevel!=0)	dumpoffspring();
			}

			if ( (lin==0) || (lin==4) )	save();

			labelLINE: /* end of line */;
			if ((lin!=0)&&(events>0))			distribution_w();
			if ((lin!=0)&&(migration_events[lin]>0))	genetic_load();
		}

		labelREP: /* end of replicate */;
	}

	printout();
	replicates_out();
	distribution_qsh_out();
	distribution_w_out();
	genresc_out();
	writeseed();
}

/* ********************************************************************* */

getinputs()
{
	tracestart();
	getseed();
	getintandskip("NIND Threatened Population 'NINDTP' (max 1000):",&NINDTP,1,1500);
	getintandskip("NIND Alternative Lines 'NINDAL' (max 1000):",&NINDAL,1,1000);
	getintandskip("Number of Migrants for lines with NINDTP individuals (max 1000):",&NMIGTP,0,1000);
	getintandskip("Number of Migrants for lines with NINDAL individuals (max 1000):",&NMIGAL,0,1000);
	getintandskip("Gender of Migrants (males 0, males&females 1):",&mGEND,0,1);
	getrealandskip("Length of genome in Morgans(99:FreeRecom) :",&L,0.0,99.0);
	getintandskip("NCRO (min 1, max 2000):",&NCRO,1,2000);
	getintandskip("NLOCI (first is neutral)(min 2, max 30):",&NLOCI,2,30);
	TOTLOCI = NCRO * NLOCI;
	getrealandskip("Lambda_s :",&Lambda_s,0.0,(double)infinity);
	getrealandskip("Lambda_L (s=1,h=0.02):",&Lambda_L,0.0,(double)infinity);
	getrealandskip("Beta_s :",&beta_s,0.0,(double)infinity);
	getrealandskip("Average |s| :",&ave_s,0.0,1.0);
	alpha_s = beta_s / ave_s;
	getintandskip("dom (constant 0, variable 1):",&dom,0,1);
	getrealandskip("Average h:",&ave_hs,0.0,(double)infinity);
	k_s = alpha_s * (pow((2.0*ave_hs),((-1.0)/beta_s))-1.0);
	getrealandskip("Stabilizing selection (Vs) :",&Vs,0.0,(double)infinity);
	getrealandskip("VE :",&VE,0.0,(double)infinity);
	getrealandskip("If Vs!=0.0, optimal fitness :",&OPT,0.0,(double)infinity);
	getintandskip("Number of generations :",&generations,2,1000);
	getintandskip("Formation of R and NR lines (t) :",&tLINES,1,generations);
	getintandskip("Formation of sublines (t) :",&tsubL,tLINES,generations);
	getintandskip("For R lines, initial NR period after formation of line :",&preNR,0,generations);
	getintandskip("Number of migrations (99: periodic):",&numMIG,0,generations);
	getintandskip("Generation intervals of migration (Lines 'NINDTP'):",&mINT_TP,0,generations);
	getintandskip("Generation intervals of migration (Lines 'NINDAL'):",&mINT_AL,0,generations);
	getintandskip("Generations since migration (distribution -w- among replicates) :",&dist,0,mINT_AL);
	getrealandskip("Minimum fitness value for extinction (99: without extinction) :",&wEXT,0.0,99.0);
	getintandskip("Number of replicates :",&replicates,1,infinity);
}

/* ********************************************************************* */

recombination_masks ()
{
	for (l=0; l<NLOCI; l++)   RM[l]=pow(2.0,(double)l);
}

/* ********************************************************************* */

natural_population ()
{
	int g0, g1, dinitialgen, dN;
	double qNP[MM][31], dadds, dadda, ds, da, dhs, dha;

	/* ***** take genotypic values of natural population ***** */

	fpop=fopen("popfile","r");

	fscanf(fpop,"%d", &dN);
	NINDNP = dN;

	for (i=0; i<dN; i++)
	for (k=0; k<NCRO; k++)
	{
		fscanf(fpop,"%d%d", &g0, &g1);
		gNP[i][k][0] = g0;
		gNP[i][k][1] = g1;
	}

	fclose(fpop);

	/* **************** take effects of genes **************** */

	fdat=fopen("datafile","r");

	fscanf(fdat,"%lf%lf", &dadds, &dadda);
	addedfixed_sNP = dadds;

	for (k=0; k<NCRO; k++) 
	for (l=0; l<NLOCI; l++)
	{
		fscanf(fdat,"%lf%lf%lf%lf%d", &ds, &da, &dhs, &dha, &dinitialgen);
		sNP[k][l] = ds;
		hsNP[k][l] = dhs;
		
		initialgenNP[k][l] = dinitialgen;
	}

	/* ******* estimate LEQ in the natural population ******* */

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<NINDNP; i++)
		{
			if (((gNP[i][k][0] & RM[l])==RM[l])&&((gNP[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
	    		else if (((gNP[i][k][0] & RM[l])!=RM[l])&&((gNP[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		qNP[k][l] = (aa/(double)NINDNP)+(Aa/(2.0*(double)NINDNP));
		if (qNP[k][l] != 0.0) 	seg_markNP[k][l] = 1;
		else			seg_markNP[k][l] = 0;

		LEQ_NP += 2.0 * ((sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0)) * qNP[k][l] * (1.0 - qNP[k][l]);

		if (sNP[k][l]==(-1.0))	LEQ_L_NP += 2.0 * ((sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0)) * qNP[k][l] * (1.0 - qNP[k][l]);
		else			LEQ_NL_NP += 2.0 * ((sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0)) * qNP[k][l] * (1.0 - qNP[k][l]);

		TLEQ_NP += -(sNP[k][l] * qNP[k][l]);
		if (qNP[k][l] > 0.0)	NSEGLOCNP ++;
	} 

	fclose(fdat);

//COMMENT
//	if (tracelevel!=0)	fprintf(fptr,"\nCOMMENT - natural_population\n");
//	if (tracelevel!=0)	for (i=0; i<NINDNP; i++)	fprintf(fptr,"gNP[%d][0][0]=%d	gNP[%d][0][1]=%d\n", i, gNP[i][0][0], i, gNP[i][0][1]);
//	if (tracelevel!=0)	for (k=0; k<NCRO; k++)	for (l=0; l<NLOCI; l++)	fprintf(fptr,"k=%d l=%d sNP=%f  hNP=%f  qNP=%f  inigenNP=%d  seg_markNP=%d\n", k, l, sNP[k][l], hsNP[k][l], qNP[k][l], initialgenNP[k][l], seg_markNP[k][l]);

}

/* ********************************************************************* */

TP()
{
	if (tracelevel!=0)	fprintf (fptr,"\n********** THREATENED POPULATION (NINDTP = %d) **********\n", NINDTP);

	NIND = NINDTP;
	init = 0;
	end = tLINES;

	for (gen=init; gen<end; gen++)	rescue[gen] = 0;

	sample_NP();

	/* initial max fecundity */
	for (i=0; i<NIND; i++) pm_sf[i] = 1.0;
}

/* ********************************************************************* */

L1()
{
	if (tracelevel!=0)	fprintf (fptr,"\n********** FORMATION OF LINES **********\n");
	if (tracelevel!=0)	fprintf (fptr,"\n******** %dR%d (L%d) ********\n", NINDTP, NINDTP, lin);

	NIND = NINDTP;
	NMIG = NMIGTP;
	mINT = mINT_TP;
	init = tLINES;
	end = generations + 1;

	/* genetic rescue */
	if (numMIG==99)	endmig = generations;
	else 		endmig = init + preNR + (mINT * (numMIG - 1));

	for (gen=init; gen<end; gen++)			rescue[gen] = 0;
	for (gen=init+preNR; gen<=endmig; gen+=mINT)  { rescue[gen] = 1; exp_migevents[lin] += 1; }

	sample();
}

/* ********************************************************************* */

L2()
{
	if (tracelevel!=0)	fprintf (fptr,"\n******** %dR%d (L%d) ********\n", NINDTP, NINDAL, lin);

	NIND = NINDAL;
	NMIG = NMIGAL;
	mINT = mINT_AL;
	init = tLINES;
	end = generations + 1;

	/* genetic rescue */
	if (numMIG==99)	endmig = generations;
	else 		endmig = init + preNR + (mINT * (numMIG - 1));

	for (gen=init; gen<end; gen++)			rescue[gen] = 0;
	for (gen=init+preNR; gen<=endmig; gen+=mINT)  { rescue[gen] = 1; exp_migevents[lin] += 1; }

	sample();
}

/* ********************************************************************* */

L3()
{
	if (tracelevel!=0)	fprintf (fptr,"\n******** %dNR%d (L%d) ********\n", NINDTP, NINDAL, lin);

	NIND = NINDAL;
	mINT = mINT_AL;
	init = tLINES;
	end = generations + 1;

	/* control of generations with GR */
	if (numMIG==99)	endmig = generations;
	else 		endmig = init + preNR + (mINT * (numMIG - 1));

	for (gen=init; gen<end; gen++)			rescue[gen] = 0;
	for (gen=init+preNR; gen<=endmig; gen+=mINT)  { rescue[gen] = 2; exp_migevents[lin] += 1; }

	sample();
}

/* ********************************************************************* */

L4()
{
	if (tracelevel!=0)	fprintf (fptr,"\n******** %dNR (L%d) ********\n", NINDTP, lin);

	NIND = NINDTP;
	mINT = mINT_TP;
	init = tLINES;
	end = tsubL;

	/* control of generations with GR */
	if (numMIG==99)	endmig = generations;
	else 		endmig = init + preNR + (mINT * (numMIG - 1));

	for (gen=init; gen<end; gen++)			rescue[gen] = 0;
	for (gen=init+preNR; gen<=endmig; gen+=mINT)  { rescue[gen] = 2; exp_migevents[lin] += 1; }

	sample();
}

/* ********************************************************************* */

subL1()
{
	if (tracelevel!=0)	fprintf (fptr,"\n********** FORMATION OF SUBLINES **********\n");
	if (tracelevel!=0)	fprintf (fptr,"\n******** %dRR%d (subL%d) ********\n", NINDTP, NINDTP, lin-4);

	NIND = NINDTP;
	NMIG = NMIGTP;
	mINT = mINT_TP;
	init = tsubL;
	end = generations + 1;

	/* genetic rescue */
	if (numMIG==99)	endmig = generations;
	else 		endmig = init + preNR + (mINT*(numMIG - 1));

	for (gen=init; gen<end; gen++)			rescue[gen] = 0;
	for (gen=init+preNR; gen<=endmig; gen+=mINT)  { rescue[gen] = 1; exp_migevents[lin] += 1; }

	sample();
}

/* ********************************************************************* */

subL2()
{
	if (tracelevel!=0)	fprintf (fptr,"\n******** %dRR%d (subL%d) ********\n", NINDTP, NINDAL, lin-4);

	NIND = NINDAL;
	NMIG = NMIGAL;
	mINT = mINT_AL;
	init = tsubL;
	end = generations + 1;

	/* genetic rescue */
	if (numMIG==99)	endmig = generations;
	else 		endmig = init + preNR + (mINT*(numMIG - 1));

	for (gen=init; gen<end; gen++)			rescue[gen] = 0;
	for (gen=init+preNR; gen<=endmig; gen+=mINT)  { rescue[gen] = 1; exp_migevents[lin] += 1; }

	sample();
}

/* ********************************************************************* */

subL3()
{
	if (tracelevel!=0)	fprintf (fptr,"\n******** %dNR%d (subL%d) ********\n", NINDTP, NINDTP, lin-4);

	NIND = NINDTP;
	mINT = mINT_TP;
	init = tsubL;
	end = generations + 1;

	/* control of generations with GR */
	if (numMIG==99)	endmig = generations;
	else 		endmig = init + preNR + (mINT * (numMIG - 1));

	for (gen=init; gen<end; gen++)			rescue[gen] = 0;
	for (gen=init+preNR; gen<=endmig; gen+=mINT)  { rescue[gen] = 2; exp_migevents[lin] += 1; }

	sample();
}

/* ********************************************************************* */

subL4()
{
	if (tracelevel!=0)	fprintf (fptr,"\n******** %dNR%d (subL%d) ********\n", NINDTP, NINDAL, lin-4);

	NIND = NINDAL;
	mINT = mINT_AL;
	init = tsubL;
	end = generations + 1;

	/* control of generations with GR */
	if (numMIG==99)	endmig = generations;
	else 		endmig = init + preNR + (mINT * (numMIG - 1));

	for (gen=init; gen<end; gen++)			rescue[gen] = 0;
	for (gen=init+preNR; gen<=endmig; gen+=mINT)  { rescue[gen] = 2; exp_migevents[lin] += 1; }

	sample();
}

/* ********************************************************************* */

sample_NP ()
{
	int np, ran_i;

	/* ****** sample individuals from Base Population ******* */

	for (i=0; i<NINDTP; i++)
	{
		ran_i = (int)(uniform() * NINDNP);

		for (k=0; k<NCRO; k++)
		{
			np=gNP[i][k][0]; gNP[i][k][0]=gNP[ran_i][k][0]; gNP[ran_i][k][0]=np;
			np=gNP[i][k][1]; gNP[i][k][1]=gNP[ran_i][k][1]; gNP[ran_i][k][1]=np;
		}
	}

	for (i=0; i<NINDTP; i++)
	{
		for (k=0; k<NCRO; k++)
		{
			g[i][k][0]=gNP[i][k][0];
			g[i][k][1]=gNP[i][k][1];
		}
		mark[i]=1;
	}

	/* ********************* frequences ********************* */

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<NINDTP; i++)
		{
			if (((g[i][k][0] & RM[l])==RM[l])&&((g[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
			else if (((g[i][k][0] & RM[l])!=RM[l])&&((g[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		q[k][l] = (aa/(double)NINDTP)+(Aa/(2.0*(double)NINDTP));
	}

	/* ***** take effects of genes from Base Population ***** */

	addedfixed_sf = addedfixed_sNP;

	for (k=0; k<NCRO; k++) 
	for (l=0; l<NLOCI; l++)
	{
		s[k][l] = sNP[k][l];
		hs[k][l] = hsNP[k][l];
		
		initialgen[k][l] = initialgenNP[k][l];
	}

//COMMENT
//	if (tracelevel!=0)    fprintf(fptr,"\nCOMMENT - sampleNP\n");
//	if (tracelevel!=0)	for (i=0; i<NINDTP; i++)	fprintf(fptr,"gNP[%d][0][0]=%d	gNP[%d][0][1]=%d\n", i, gNP[i][0][0], i, gNP[i][0][1]);
//	if (tracelevel!=0)	for (k=0; k<NCRO; k++)	for (l=0; l<NLOCI; l++)	fprintf(fptr,"k=%d l=%d s=%f  h=%f  q=%f  inigen=%d  seg_markNP=%d\n", k, l, s[k][l], hs[k][l], q[k][l], initialgen[k][l], seg_markNP[k][l]);

}

/* ********************************************************************* */

save()
{
	int j;

	/* ************ genotypic values ************* */

	if (NINDTP<NINDAL)
	{
		NIND = NINDAL;
		if (tracelevel!=0)	fprintf(fptr,"NINDTP%d<NINDAL%d\n", NINDTP, NINDAL);
	}	
	else	NIND = NINDTP;

	for (i=0; i<NIND; i++)
	for (k=0; k<NCRO; k++)
	{
		sgnt[i][k][0]=g[i][k][0];
		sgnt[i][k][1]=g[i][k][1];
	}

	for (i=0; i<NIND; i++)	spm_sf[i] = pm_sf[i];

	/* ***************** parents ***************** */

	for (i=0; i<NIND; i++)
	{
		smoth[i]=mother[i];
		sfath[i]=father[i];
	}
	for (i=0; i<NIND; i++)
	for (j=0; j<NIND; j++)
	{
		ssmt[sfath[i]][sfath[j]]=smat[father[i]][father[j]];
		ssmt[sfath[i]][smoth[j]]=smat[father[i]][mother[j]];
		ssmt[smoth[i]][sfath[j]]=smat[mother[i]][father[j]];
		ssmt[smoth[i]][smoth[j]]=smat[mother[i]][mother[j]];
	}

	/* ***** effects of genes and frequences ***** */

	for (k=0; k<NCRO; k++) 
	for (l=0; l<NLOCI; l++)
	{
		ss[k][l] = s[k][l];
		hss[k][l] = hs[k][l];
		qs[k][l] = q[k][l];
		
		sinitgen[k][l] = initialgen[k][l];
	}

	/* ********** stabilizing selection ********** */

	if (Vs!=0.0)
	{
		for (i=0; i<NIND; i++)
		{
			sp_a[mother[i]] = p_a[mother[i]];
			sp_a[father[i]] = p_a[father[i]];

			sF[mother[i]] = F[mother[i]];
			sF[father[i]] = F[father[i]];
		}
		sVA = VA;
	}

//COMMENT
//	if (tracelevel!=0)	fprintf(fptr,"\n\nCOMMENT - saveind\n");
//	if (tracelevel!=0)	for (i=0; i<NIND; i++)	fprintf(fptr,"i=%d mother=%d father=%d  sgnt[%d][0][0]=%d  sgnt[%d][0][1]=%d\n", i, smoth[i], sfath[i],  i, sgnt[i][0][0], i, sgnt[i][0][1]); 
//	if (tracelevel!=0)	for (k=0; k<NCRO; k++)	for (l=0; l<NLOCI; l++)	fprintf(fptr,"s[%d][%d]=%f hs[%d][%d]=%f q[%d][%d]=%f  inigen=%d  seg_markNP=%d\n", k, l, ss[k][l], k, l, hss[k][l], k, l, qs[k][l], sinitgen[k][l], seg_markNP[k][l]);

//	if ((tracelevel!=0)&&(Vs!=0.0))	fprintf(fptr,"\n*** stabilizing selection\n");
//	if ((tracelevel!=0)&&(Vs!=0.0))	for (i=0; i<NIND; i++)	fprintf(fptr,"p_a[mother[%d]]=%f\tp_a[father[%d]]=%f\tF[mother[%d]]=%f\tF[father[%d]]=%f\n", i, sp_a[mother[i]], i, sp_a[father[i]], i, sF[mother[i]], i, sF[father[i]]);
//	if ((tracelevel!=0)&&(Vs!=0.0))	fprintf(fptr,"VA=%f\n", sVA);

}

/* ********************************************************************* */

sample()
{
	int j;

	/* ************* genotypic values ************* */

	for (i=0; i<NIND; i++)
	for (k=0; k<NCRO; k++)
	{
		g[i][k][0]=sgnt[i][k][0];
		g[i][k][1]=sgnt[i][k][1];
	}

	for (i=0; i<NIND; i++)	pm_sf[i] = spm_sf[i];

	/* ******** parents of each individual ******** */

	for (i=0; i<NIND; i++)
	{
		mother[i]=smoth[i];
		father[i]=sfath[i];
	}
	for (i=0; i<NIND; i++)
	for (j=0; j<NIND; j++)
	{
		smat[father[i]][father[j]]=ssmt[sfath[i]][sfath[j]];
		smat[father[i]][mother[j]]=ssmt[sfath[i]][smoth[j]];
		smat[mother[i]][father[j]]=ssmt[smoth[i]][sfath[j]];
		smat[mother[i]][mother[j]]=ssmt[smoth[i]][smoth[j]];
	}

	/* ***** effects of genes and frequences ***** */

	for (k=0; k<NCRO; k++) 
	for (l=0; l<NLOCI; l++)
	{
		s[k][l] = ss[k][l];
		hs[k][l] = hss[k][l];
		q[k][l] = qs[k][l];
		
		initialgen[k][l] = sinitgen[k][l];
	}

	/* ********** stabilizing selection ********** */

	if (Vs!=0.0)
	{
		for (i=0; i<NIND; i++)
		{
			p_a[mother[i]] = sp_a[mother[i]];
			p_a[father[i]] = sp_a[father[i]];

			F[mother[i]] = sF[mother[i]];
			F[father[i]] = sF[father[i]];
		}
		VA = sVA;
	}

	/* ***** return to mark 0 those individuals caught from NP in other lines during migration events ***** */

	for (i=NINDTP; i<NINDNP; i++)	mark[i]=0;

//COMMENT
//	if (tracelevel!=0)	fprintf(fptr,"\n\nCOMMENT - sample\n");
//	if (tracelevel!=0)	for (i=0; i<NIND; i++)	fprintf(fptr,"g[%d][0][0]=%d	g[%d][0][1]=%d\n", i, g[i][0][0], i, g[i][0][1]);
//	if (tracelevel!=0)	for (k=0; k<NCRO; k++) 	for (l=0; l<NLOCI; l++)	fprintf(fptr,"s[%d][%d]=%f hs[%d][%d]=%f q[%d][%d]=%f  inigen=%d  seg_markNP=%d\n", k, l, s[k][l], k, l, hs[k][l], k, l, q[k][l], initialgen[k][l], seg_markNP[k][l]);

//	if ((tracelevel!=0)&&(Vs!=0.0))	fprintf(fptr,"\n*** stabilizing selection\n");
//	if ((tracelevel!=0)&&(Vs!=0.0))	for (i=0; i<NIND; i++)	fprintf(fptr,"p_a[mother[%d]]=%f\tp_a[father[%d]]=%f\tF[mother[%d]]=%f\tF[father[%d]]=%f\n", i, p_a[mother[i]], i, p_a[father[i]], i, F[mother[i]], i, F[father[i]]);
//	if ((tracelevel!=0)&&(Vs!=0.0))	fprintf(fptr,"VA=%f\n", VA);

}

/* ********************************************************************* */

parents()
{
	double sum_fem=0.0, sum_mal=0.0;

	for (i=0; i<NIND; i++)
	{
		if (i%2==0)	sum_fem += pm_sf[i];
		else		sum_mal += pm_sf[i];
	}

	females[gen] = sum_fem/(NIND/2);
	males[gen] = sum_mal/(NIND/2);

//COMMENT
	if (tracelevel!=0)	fprintf(fptr,"\nCOMMENT - parents (before new mut)\nfemales = %f  males = %f\n", females[gen], males[gen]);
}

/* ********************************************************************* */

extinction_w()
{
	/* exclusion of migrants and generations tLINES or tsubL */
	double qext[MM][31], parentals, STR=0.0, summ_dq=0.0, d;

	//lethals on parental generation
	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;
		for (i=0; i<NIND; i++)
		{
			if (((sp[i][k][0] & RM[l])==RM[l])&&((sp[i][k][1] & RM[l])==RM[l]))	{aa+=1.0;}
			else if (((sp[i][k][0] & RM[l])!=RM[l])&&((sp[i][k][1] & RM[l])!=RM[l]))	{AA+=1.0;}
			else	{Aa+=1.0;}
		}
		qext[k][l] = (aa/(double)NIND)+(Aa/(2.0*(double)NIND));

		if (s[k][l]!=(-1.0))	/* dq excluding lethals */
		{
			d = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
			summ_dq += d * qext[k][l];
		}
		if (s[k][l]==(-1.0))	STR += qext[k][l];
	}

	/* individuals gen-1 */
	parentals = (females[gen-1] + males[gen-1])/2.0;

	if (parentals == 0.0)					{accum (&wext_00[lin],  1.0); accum (&STR_00[lin], STR); accum (&dq_00[lin], summ_dq);}
	else if ( (parentals > 0.0) && (parentals < 0.1) )	{accum (&wext_00_01[lin],  1.0); accum (&STR_00_01[lin], STR); accum (&dq_00_01[lin], summ_dq);}
	else if ( (parentals >= 0.1) && (parentals < 0.2) )	{accum (&wext_01_02[lin],  1.0); accum (&STR_01_02[lin], STR); accum (&dq_01_02[lin], summ_dq);}
	else if ( (parentals >= 0.2) && (parentals < 0.3) )	{accum (&wext_02_03[lin],  1.0); accum (&STR_02_03[lin], STR); accum (&dq_02_03[lin], summ_dq);}
	else if ( (parentals >= 0.3) && (parentals < 0.4) )	{accum (&wext_03_04[lin],  1.0); accum (&STR_03_04[lin], STR); accum (&dq_03_04[lin], summ_dq);}
	else if ( (parentals >= 0.4) && (parentals < 0.5) )	{accum (&wext_04_05[lin],  1.0); accum (&STR_04_05[lin], STR); accum (&dq_04_05[lin], summ_dq);}
	else if ( (parentals >= 0.5) && (parentals < 0.6) )	{accum (&wext_05_06[lin],  1.0); accum (&STR_05_06[lin], STR); accum (&dq_05_06[lin], summ_dq);}
	else if ( (parentals >= 0.6) && (parentals < 0.7) )	{accum (&wext_06_07[lin],  1.0); accum (&STR_06_07[lin], STR); accum (&dq_06_07[lin], summ_dq);}
	else if ( (parentals >= 0.7) && (parentals < 0.8) )	{accum (&wext_07_08[lin],  1.0); accum (&STR_07_08[lin], STR); accum (&dq_07_08[lin], summ_dq);}
	else if ( (parentals >= 0.8) && (parentals < 0.9) )	{accum (&wext_08_09[lin],  1.0); accum (&STR_08_09[lin], STR); accum (&dq_08_09[lin], summ_dq);}
	else if ( (parentals >= 0.9) && (parentals < 1.0) )	{accum (&wext_09_10[lin],  1.0); accum (&STR_09_10[lin], STR); accum (&dq_09_10[lin], summ_dq);}
	else if (parentals == 1.0)				{accum (&wext_10[lin],  1.0); accum (&STR_10[lin], STR); accum (&dq_10[lin], summ_dq);}
}

/* ********************************************************************* */

migration()
{
	int mi, ran_i;
	N = NIND-NMIG;

	/* ***** take migrants from Base Population ***** */

	if (tracelevel!=0)	fprintf(fptr,"\nMigrants from Base Population (number of migrants: %d)\n", NMIG);

	for (i=N; i<NIND; i++)
	{
		do {ran_i = (int)(uniform() * NINDNP);}		while (mark[ran_i]==1);
		mark[ran_i] = mark[i];

		for (k=0; k<NCRO; k++)
		{
			mi=gNP[i][k][0]; gNP[i][k][0]=gNP[ran_i][k][0]; gNP[ran_i][k][0]=mi;
			mi=gNP[i][k][1]; gNP[i][k][1]=gNP[ran_i][k][1]; gNP[ran_i][k][1]=mi;

			g[i][k][0]=gNP[i][k][0];
			g[i][k][1]=gNP[i][k][1];
		}

		mother[i] = -1;
		father[i] = -1;
		mark[i] = 1;

		if ((tracelevel!=0)&&(mGEND==0)) fprintf(fptr,"(mal) ran_i %d g[%d][0][0]=%d\tg[%d][0][1]=%d\tmark=%d\n", ran_i, i, g[i][0][0], i, g[i][0][1], mark[i]);
		if ((tracelevel!=0)&&(mGEND!=0))
		{
			if(i%2==0)	fprintf(fptr,"(fem) ran_i %d g[%d][0][0]=%d\tg[%d][0][1]=%d\tmark=%d\n", ran_i, i, g[i][0][0], i, g[i][0][1], mark[i]);
			else		fprintf(fptr,"(mal) ran_i %d g[%d][0][0]=%d\tg[%d][0][1]=%d\tmark=%d\n", ran_i, i, g[i][0][0], i, g[i][0][1], mark[i]);
		}
	}

	/* *********** frequences after migration ********** */

	selected_genes();
}

/* ********************************************************************* */

poisson_tables ()
{   
	/* SELECTED LOCI WITH POISSON (2NL) NEW MUTATIONS */

	if ( (exp(-2.0*(double)NIND*Lambda_s) != 0.0)&&(2.0*(double)NIND*Lambda_s < normalthreshold) )
	generatepoissontable(2.0*(double)NIND*Lambda_s, &lastinmutantspoissontable, mutantspoissontable, maxmpt-1);

	if ( (exp(-2.0*(double)NIND*Lambda_L) != 0.0)&&(2.0*(double)NIND*Lambda_L < normalthreshold) )
	generatepoissontable(2.0*(double)NIND*Lambda_L, &lastinmutantspoissontableL, mutantspoissontableL, maxmpt-1);

	/* NUMERO DE RECOMBINACIONES POISSON CON MEDIA L*/

	if ( (exp(-L) != 0.0) && (L < normalthreshold) )
	generatepoissontable(L, &lastinrecombinantpoissontable, recombinantpoissontable, maxmpt-1);
}

/* ********************************************************************* */

mutation_neutral ()
{
	int ran_k, ran_h, ran_i, m, muts;

	/* NEUTRAL GENES: (POISSON) 2N(Lambda_a) NEW MUTATIONS (TWO DIRECTIONS) */

	muts = mutationnumber();

	if (tracelevel!=0)    fprintf(fptr,"\n New neutral mutants = %d\n", muts);

	for (m=0; m<muts; m++)
	{
		ran_i = (int)(uniform()*NIND);
		ran_k = (int)(uniform()*NCRO);
		ran_h = (int)(uniform()*2.0);

		if ( (g[ran_i][ran_k][ran_h] & RM[0])==RM[0] ) 	g[ran_i][ran_k][ran_h]=(g[ran_i][ran_k][ran_h] & (~RM[0]));	
		else						g[ran_i][ran_k][ran_h]=(g[ran_i][ran_k][ran_h] | RM[0]);
	}
}

/* ********************************************************************* */

void mutation_selected ()
{
	int countNoSS, NoSS_k[60001], NoSS_l[60001], mutants_ocurred;
	int ran_k, ran_h, ran_i, ran_l, m, muts, mutsL;
	double zs;

	/* SELECTED GENES: (POISSON) 2N(Lambda_s + Lambda_L) NEW MUTATIONS */

	muts = mutationnumber();
	mutsL = mutationnumberL();

	if ( (muts + mutsL) == 0 ) goto label;

	if (tracelevel!=0)    fprintf(fptr,"\n New selected mutants  muts = %d  mutsL = %d\n", muts, mutsL);

	countNoSS = 0;

	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	{
		if (q[k][l] == 0.0)
		{
			countNoSS += 1;
			NoSS_k[countNoSS-1] = k;
			NoSS_l[countNoSS-1] = l;
		}
	}

	if (tracelevel!=0)    fprintf(fptr," countNoSS=%d\n", countNoSS);

	mutants_ocurred = 0;

	if (countNoSS != 0)
	{
		disorder_NoSS (NoSS_k,NoSS_l,countNoSS);

		for (m=0; m<countNoSS; m++)
		{
			if (mutants_ocurred==(muts+mutsL))    goto label;

			ran_i = (int)(uniform()*NIND);
   			ran_h = (int)(uniform()*2.0);

			g[ran_i][NoSS_k[m]][ran_h]=(g[ran_i][NoSS_k[m]][ran_h] | RM[NoSS_l[m]]); 

		   	mutants_ocurred += 1;
			initialgen[NoSS_k[m]][NoSS_l[m]] = gen;

			/* **** Mutations segregating in NP **** */

			if (seg_markNP[NoSS_k[m]][NoSS_l[m]] == 1)
			{
				s[NoSS_k[m]][NoSS_l[m]] = sNP[NoSS_k[m]][NoSS_l[m]];
				hs[NoSS_k[m]][NoSS_l[m]] = hsNP[NoSS_k[m]][NoSS_l[m]];
			}
			else
			{
				/* ****** Lethal mutations ****** */

				if(mutants_ocurred <= mutsL)
				{
					s[NoSS_k[m]][NoSS_l[m]] = (-1.0);
					hs[NoSS_k[m]][NoSS_l[m]] = 0.02;
				}
				else
				{
	     			/* ****** values of s, hs ****** */

	    				zs = gengam (alpha_s, beta_s);
					if (zs > 1.0)	zs=1.0;
					s[NoSS_k[m]][NoSS_l[m]] = (-zs);
					if (dom==0)	hs[NoSS_k[m]][NoSS_l[m]] = ave_hs;
					else		hs[NoSS_k[m]][NoSS_l[m]] = uniform() * exp(k_s*s[NoSS_k[m]][NoSS_l[m]]);
				}
			}
//COMMENT
			if (tracelevel!=0)    fprintf(fptr,"alpha_s=%f  beta_s=%f ", alpha_s, beta_s);
			if (tracelevel!=0)    fprintf(fptr,"s=%f  h=%f  ", s[NoSS_k[m]][NoSS_l[m]], hs[NoSS_k[m]][NoSS_l[m]]);
			if (tracelevel!=0)    fprintf(fptr,"ran_i=%d  k=%d  l=%d   ran_h=%d  seg_markNP=%d\n", ran_i, NoSS_k[m],  NoSS_l[m], ran_h, seg_markNP[NoSS_k[m]][NoSS_l[m]]);
		}
	}

	for (m=mutants_ocurred; m<muts; m++)
	{
		ran_i = (int)(uniform()*NIND);
		ran_k = (int)(uniform()*NCRO);
		do {ran_l = (int)(uniform()*NLOCI);}   while (ran_l==0);
		ran_h = (int)(uniform()*2.0);

		g[ran_i][ran_k][ran_h]=(g[ran_i][ran_k][ran_h] | RM[ran_l]); 

//COMMENT
//		if (tracelevel!=0)    fprintf(fptr,"(rec. mut.) ran_i=%d  ran_k=%d  ran_l=%d  ran_h=%d\n", ran_i, ran_k, ran_l, ran_h);
	}

	label: /* end of mutations */;

	/* without extinctions */
	if(wEXT==99.0)
	{
		for (k=0; k<NCRO; k++)
		for (l=1; l<NLOCI; l++)
		{
			if (s[k][l] < (-0.99))
			{
//				if (tracelevel!=0)    fprintf(fptr,"(lethals) k=%d l=%d  s=%f  ",k,l,s[k][l]);
				s[k][l] = (-0.99);
//				if (tracelevel!=0)    fprintf(fptr,"s=%f\n",s[k][l]);
			}
		}
	}
}

/* ********************************************************************* */

int mutationnumber()
{
	int r;
	if ((2.0*(double)NIND*Lambda_s < normalthreshold) && (exp(-2.0*(double)NIND*Lambda_s) != 0.0) )
	{
		r = poisson(lastinmutantspoissontable, mutantspoissontable);
	}
	else r = (int)( normal(2.0*(double)NIND*Lambda_s, sqrt(2.0*(double)NIND*Lambda_s)) );
	return(r);
}

/* ********************************************************************* */

int mutationnumberL()
{
	int r;
	if ((2.0*(double)NIND*Lambda_L < normalthreshold) && (exp(-2.0*(double)NIND*Lambda_L) != 0.0) )
	{
		r = poisson(lastinmutantspoissontableL, mutantspoissontableL);
	}
	else r = (int)( normal(2.0*(double)NIND*Lambda_L, sqrt(2.0*(double)NIND*Lambda_L)) );
	return(r);
}

/* ********************************************************************* */

void disorder_NoSS (NoSS_k,NoSS_l,countNoSS)
int NoSS_k[], NoSS_l[];
{
	int a, b, rnd;
	
	for (i=0; i<countNoSS-1; i++)
	{
	   rnd=(int)(uniform()*(countNoSS-i));
	   a=NoSS_k[countNoSS-1-i]; NoSS_k[countNoSS-1-i]=NoSS_k[rnd]; NoSS_k[rnd]=a;
	   b=NoSS_l[countNoSS-1-i]; NoSS_l[countNoSS-1-i]=NoSS_l[rnd]; NoSS_l[rnd]=b;
	}
}

/* ********************************************************************* */

neutral_genes()
{
	double H=0.0;

	//exclude migrants from average values
	if (rescue[gen]==1)	N = NIND - NMIG;
	else 			N = NIND;

	for (k=0; k<NCRO; k++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<N; i++)
	    	{
			if (((g[i][k][0] & RM[0])==RM[0])&&((g[i][k][1] & RM[0])==RM[0]))	aa+=1.0;
			else if (((g[i][k][0] & RM[0])!=RM[0])&&((g[i][k][1] & RM[0])!=RM[0]))	AA+=1.0;
	 		else	Aa+=1.0;
	    	}

		q[k][0] = (aa/(double)N)+(Aa/(2.0*(double)N));
		H += 2.0 * q[k][0] * (1.0 - q[k][0]);
	}

	accum(&Hw[lin][gen], H/(double)NCRO);

//COMMENT
//	if (tracelevel!=0)	fprintf(fptr,"\nCOMMENT - Neutral genes\n");
//	if (tracelevel!=0)	for (k=0; k<NCRO; k++)	fprintf(fptr,"k=%d\t0\tAA=%1.0f\tAa=%1.0f\taa=%1.0f\tq[k][0]=%f\tinitialgen[k][0]=%d\n",k,AA,Aa,aa,q[k][0],initialgen[k][0]);

}

/* ********************************************************************* */

selected_genes()
{
	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<NIND; i++)
		{
			if (((g[i][k][0] & RM[l])==RM[l])&&((g[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
			else if (((g[i][k][0] & RM[l])!=RM[l])&&((g[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
			else	Aa+=1.0;
		}

		q[k][l] = (aa/(double)NIND)+(Aa/(2.0*(double)NIND));
	}
}

/* ********************************************************************* */

summary_genresc()
{
	int mut[MM][31], NSEGLOC=0, mut_NP=0, mut_Line=0, lost_NP=0, lost_Line=0, mut_shared=0;
	double d, qml[MM][31], sum_s=0.0, sum_h=0.0, sum_q=0.0, AAt=0.0, Aat=0.0, aat=0.0;

	N = NIND - NMIG;
	summ_dq_l[migration_events[lin]] = 0.0;
	summ_dq_m[migration_events[lin]] = 0.0;
	STR_l[migration_events[lin]] = 0.0;
	STR_m[migration_events[lin]] = 0.0;

	for (k=0; k<NCRO; k++)	for (l=1; l<NLOCI; l++)	mut[k][l] = 0;

	/* **** Migrants (after new mutation) **** */

	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;
		for (i=N; i<NIND; i++)
		{
			if (((g[i][k][0] & RM[l])==RM[l])&&((g[i][k][1] & RM[l])==RM[l]))	{aa+=1.0; aat+=1.0;}
			else if (((g[i][k][0] & RM[l])!=RM[l])&&((g[i][k][1] & RM[l])!=RM[l]))	{AA+=1.0; AAt+=1.0;}
			else	{Aa+=1.0; Aat+=1.0;}
		}
		qml[k][l] = (aa/(double)NMIG)+(Aa/(2.0*(double)NMIG));

		if (qml[k][l] > 0.0)
		{
			NSEGLOC ++;
			mut[k][l] = 1;
			sum_q += qml[k][l];
			sum_s += s[k][l];
			sum_h += hs[k][l];	
		}

		d = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
		summ_dq_m[migration_events[lin]] += d * qml[k][l];

		if ( (s[k][l]==(-1.0)) && (qml[k][l]!=0.0) )	STR_m[migration_events[lin]] += qml[k][l];
	}

	accum (&NSEGLOCm[lin][gen], (double)NSEGLOC);
	accum (&qm[lin][gen], sum_q/(double)NSEGLOC);
	accum (&sm[lin][gen], sum_s/(double)NSEGLOC);
	accum (&hm[lin][gen], sum_h/(double)NSEGLOC);
	accum (&AAm[lin][gen], AAt/(double)NMIG);
	accum (&Aam[lin][gen], Aat/(double)NMIG);
	accum (&aam[lin][gen], aat/(double)NMIG);

	/* ****** Line (after new mutation) ****** */

	NSEGLOC=0; sum_s=0.0; sum_h=0.0; sum_q=0.0; AAt=0.0; Aat=0.0, aat=0.0;

	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;
		for (i=0; i<N; i++)
		{
			if (((g[i][k][0] & RM[l])==RM[l])&&((g[i][k][1] & RM[l])==RM[l]))	{aa+=1.0; aat+=1.0;}
			else if (((g[i][k][0] & RM[l])!=RM[l])&&((g[i][k][1] & RM[l])!=RM[l]))	{AA+=1.0; AAt+=1.0;}
			else	{Aa+=1.0; Aat+=1.0;}
		}
		qml[k][l] = (aa/(double)N)+(Aa/(2.0*(double)N));

		if (qml[k][l] > 0.0)
		{
			NSEGLOC ++;

			if (seg_markNP[k][l]==1)
			{
				/* mutations NP segregating in line */
				mut_NP ++;
				/* mutations NP segregating in line and migrants */			
				if (mut[k][l] == 1)	mut_shared ++;
			}
			/* new mutations segregating in line */
			else	mut_Line ++;	

			sum_q += qml[k][l];
			sum_s += s[k][l];
			sum_h += hs[k][l];	
		}
		else if ((qml[k][l] == 0.0)&&(initialgen[k][l] != (-99)))
		{
			if (seg_markNP[k][l]==1)	lost_NP ++;
			else				lost_Line ++;
		}

		d = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
		summ_dq_l[migration_events[lin]] += d * qml[k][l];

		if (s[k][l]==(-1.0))	STR_l[migration_events[lin]] += qml[k][l];
	}

	accum (&NSEGLOCl[lin][gen], (double)NSEGLOC);
	accum (&ql[lin][gen], sum_q/(double)NSEGLOC);
	accum (&sl[lin][gen], sum_s/(double)NSEGLOC);
	accum (&hl[lin][gen], sum_h/(double)NSEGLOC);
	accum (&AAl[lin][gen], AAt/(double)N);
	accum (&Aal[lin][gen], Aat/(double)N);
	accum (&aal[lin][gen], aat/(double)N);

	accum (&mutNP[lin][gen], (double)mut_NP);
	accum (&mutShared[lin][gen], (double)mut_shared);
	accum (&mutLine[lin][gen], (double)mut_Line);
	accum (&lostNP[lin][gen], (double)lost_NP);
	accum (&lostLine[lin][gen], (double)lost_Line);
}

/* ********************************************************************* */

summary_NR()
{
	int NSEGLOC=0, mut_NP=0, mut_Line=0, lost_NP=0, lost_Line=0;
	double d, qml[k][l], sum_s=0.0, sum_h=0.0, sum_q=0.0, AAt=0.0, Aat=0.0, aat=0.0;

	summ_dq_l[migration_events[lin]] = 0.0;
	STR_l[migration_events[lin]] = 0;

	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;
		for (i=0; i<NIND; i++)
		{
			if (((g[i][k][0] & RM[l])==RM[l])&&((g[i][k][1] & RM[l])==RM[l]))	{aa+=1.0; aat+=1.0;}
			else if (((g[i][k][0] & RM[l])!=RM[l])&&((g[i][k][1] & RM[l])!=RM[l]))	{AA+=1.0; AAt+=1.0;}
			else	{Aa+=1.0; Aat+=1.0;}
		}
		qml[k][l] = (aa/(double)NIND)+(Aa/(2.0*(double)NIND));

		if (qml[k][l] > 0.0)
		{
			NSEGLOC ++;
			if (seg_markNP[k][l]==1)	mut_NP ++;
			else				mut_Line ++;

			sum_q += qml[k][l];
			sum_s += s[k][l];
			sum_h += hs[k][l];	
		}
		else if ((qml[k][l] == 0.0)&&(initialgen[k][l] != (-99)))
		{
			if (seg_markNP[k][l]==1)	lost_NP ++;
			else				lost_Line ++;
		}

		d = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
		summ_dq_l[migration_events[lin]] += d * qml[k][l];

		if (s[k][l]==(-1.0))	STR_l[migration_events[lin]] += qml[k][l];
	}

	accum (&NSEGLOCl[lin][gen], (double)NSEGLOC);
	accum (&ql[lin][gen], sum_q/(double)NSEGLOC);
	accum (&sl[lin][gen], sum_s/(double)NSEGLOC);
	accum (&hl[lin][gen], sum_h/(double)NSEGLOC);
	accum (&AAl[lin][gen], AAt/(double)NIND);
	accum (&Aal[lin][gen], Aat/(double)NIND);
	accum (&aal[lin][gen], aat/(double)NIND);

	accum (&mutNP[lin][gen], (double)mut_NP);
	accum (&mutLine[lin][gen], (double)mut_Line);
	accum (&lostNP[lin][gen], (double)lost_NP);
	accum (&lostLine[lin][gen], (double)lost_Line);
}

/* ********************************************************************* */

genetic_load()
{
	for (mig=1; mig<=migration_events[lin]; mig++)
	{
		/* All annotated replicates */
		accum (&gload_line[lin][mig], summ_dq_l[mig]);
		accum (&gload_mig[lin][mig], summ_dq_m[mig]);
		accum (&STR_line[lin][mig], STR_l[mig]);
		accum (&STR_mig[lin][mig], STR_m[mig]);

		/* Extinct replicates after annotation */
		if (repmark == 1)
		{
			accum (&ext_gload_line[lin][mig], summ_dq_l[mig]);
			accum (&ext_gload_mig[lin][mig], summ_dq_m[mig]);
			accum (&ext_STR_line[lin][mig], STR_l[mig]);
			accum (&ext_STR_mig[lin][mig], STR_m[mig]);
		}
		/* Surviving replicates */
		else
		{
			accum (&sur_gload_line[lin][mig], summ_dq_l[mig]);
			accum (&sur_gload_mig[lin][mig], summ_dq_m[mig]);
			accum (&sur_STR_line[lin][mig], STR_l[mig]);
			accum (&sur_STR_mig[lin][mig], STR_m[mig]);
		}
	}
}

/* ********************************************************************* */

dumpoffspringaftermutation()
{
	if (tracelevel==0)   return (0);
	fprintf(fptr,"\n Offspring after mutation \n");	

	if ((rescue[gen]==1)&&(mGEND==0))
	{
		for (i=0; i<NIND; i++)
		{
			if ((i%2==0)&&(i<NIND-NMIG))	fprintf(fptr,"(gf0 gf1)	%d  %d\n",g[i][0][0],g[i][0][1]); /* females */
			else				fprintf(fptr,"(gm0 gm1)	%d  %d\n",g[i][0][0],g[i][0][1]); /* males and migrants */
		}
	}
	else
	{
		for (i=0; i<NIND; i++)
		{
			if (i%2==0)			fprintf(fptr,"(gf0 gf1)	%d  %d\n",g[i][0][0],g[i][0][1]); /* females */
			else				fprintf(fptr,"(gm0 gm1)	%d  %d\n",g[i][0][0],g[i][0][1]); /* males */
		}
	}
}

/* ********************************************************************* */

genotypic_values()
{
	double leqLf[NN], leqNLf[NN], leqT[NN], genval_sf[NN], q_l[MM][31];
	double gsum_sf=0.0, gsum2_sf=0.0, sum_leqNLf=0.0, sum_leqLf=0.0, sum_leqT=0.0;
	double pmo_a[NN], pfa_a[NN], ave_pp_a, sq_p_a, gsum_a=0.0, VAW, gsum2_a=0.0;

//	if (tracelevel!=0)	fprintf(fptr,"\nGenotypic values \n");

	//exclude migrants from average values
	if (rescue[gen]==1)
	{
		N = NIND - NMIG;

		for (k=0; k<NCRO; k++)
		for (l=1; l<NLOCI; l++)
		{
			AA=0.0; Aa=0.0; aa=0.0;
			for (i=0; i<N; i++)
			{
				if (((g[i][k][0] & RM[l])==RM[l])&&((g[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
				else if (((g[i][k][0] & RM[l])!=RM[l])&&((g[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
				else	Aa+=1.0;
			}
			q_l[k][l] = (aa/(double)N)+(Aa/(2.0*(double)N));
		}
	}
	else
	{
	 	N = NIND;
		for (k=0; k<NCRO; k++)	for (l=1; l<NLOCI; l++)	 q_l[k][l] = q[k][l];
	}

	for (i=0; i<NIND; i++)
	{
		leqLf[i] = 0.0;
		leqNLf[i] = 0.0;
		leqT[i] = 0.0;

		genval_sf[i] = addedfixed_sf;

		for (k=0; k<NCRO; k++)
		for (l=1; l<NLOCI; l++)
		if (initialgen[k][l] != (-99))
		{

			/* ************* Genotypic values  ************ */

			if (((g[i][k][0] & RM[l])==RM[l])&&((g[i][k][1] & RM[l])==RM[l]))	genval_sf[i] *= (1.0 + s[k][l]);
			else if (((g[i][k][0] & RM[l])!=RM[l])&&((g[i][k][1] & RM[l])!=RM[l])) /* AA */;
			else	genval_sf[i] *= (1.0 + (s[k][l]*hs[k][l]));

			/* *************** estimate leq ************** */

			if (s[k][l]==(-1.0))	/* Lethal mutations */
			{
				leqT[i] += -(s[k][l] * q_l[k][l]);
				leqLf[i] += 2.0 * ((s[k][l]*hs[k][l]) - (s[k][l]/2.0)) * q_l[k][l] * (1.0 - q_l[k][l]);
			}
			else			/* Non Lethal mutations */
			{
				leqT[i] += -(s[k][l] * q_l[k][l]);
				leqNLf[i] += 2.0 * ((s[k][l]*hs[k][l]) - (s[k][l]/2.0)) * q_l[k][l] * (1.0 - q_l[k][l]);
			}

//COMMENT
//			if (tracelevel!=0)
//			{
//				if ( ((g[i][k][0] & RM[l])==RM[l]) || ((g[i][k][1] & RM[l])==RM[l]) )	fprintf(fptr,"gen=%d i=%d k=%d l=%d s=%f h=%f q=%f ini=%d\n", gen, i, k, l, s[k][l], hs[k][l], q_l[k][l], initialgen[k][l]);
//			}
		}
//COMMENT
//		if (tracelevel!=0)	fprintf(fptr,"%d    genval_sf = %f\n", i, genval_sf[i]);

	}

	/* ************************************** */
	/* ************* phenotypeB ************* */

	if (Vs!=0.0)
	{
		for (i=0; i<NIND; i++)
		{
			pmo_a[mother[i]] = p_a[mother[i]];
			pfa_a[father[i]] = p_a[father[i]];
		}
//COMMENT
//		if (tracelevel!=0)	fprintf(fptr,"\nCOMMENT - phenotypeB\n*** stabilizing selection\n");
//		if (tracelevel!=0)	for (i=0; i<NIND; i++)	fprintf(fptr,"i=%d\tpmo_a(%d)=%f\tpfa_a(%d)=%f\n", i, mother[i], pmo_a[i], father[i], pfa_a[i]);
	}

	for (i=0; i<NIND; i++)
	{
		if (Vs==0.0)	pm_sf[i] = genval_sf[i];
		else
		{
			if ((gen==0) || ((rescue[gen]==1)&&(i>=NIND-NMIG)))	p_a[i] = normal(0.0, VE);		
			else
			{
				ave_pp_a = (pmo_a[mother[i]]+pfa_a[father[i]])/2.0;		
				VAW = VA * (1 - ((F[mother[i]]+F[father[i]])/2.0));
				p_a[i] = normal(ave_pp_a, VAW);
			}

			sq_p_a = pow ((p_a[i] - OPT), 2.0);
			pm_sf[i] = genval_sf[i] * exp(-0.5*(sq_p_a)/Vs);
		}
	}

	for (i=0; i<N; i++)
	{
		sum_leqNLf += leqNLf[i];
		sum_leqLf += leqLf[i];
		sum_leqT += leqT[i];

		gsum_sf += pm_sf[i];
		gsum2_sf += (pm_sf[i]*pm_sf[i]);
		if (Vs!=0.0)
		{
			gsum_a += p_a[i];
			gsum2_a += (p_a[i]*p_a[i]);
		}
	}

	if (rescue[gen-dist]!=0)
	{
		events = migration_events[lin];
		mean_w[events] = gsum_sf/(double)N;
	}

	accum (&LEQ_NLf[lin][gen], sum_leqNLf/(double)N);
	accum (&LEQ_Lf[lin][gen], sum_leqLf/(double)N);
	accum (&LEQ_T[lin][gen], sum_leqT/(double)N);

	accum (&gmean_sf[lin][gen], gsum_sf/(double)N);
	accum (&gvar_sf[lin][gen], (gsum2_sf - (gsum_sf*gsum_sf / (double)N)) / ((double)N - 1.0));

	if (Vs!=0.0)
	{
		accum (&gmean_a[lin][gen], gsum_a/(double)N);
		accum (&gvar_a[lin][gen], (gsum2_a - (gsum_a*gsum_a / (double)N)) / ((double)N - 1.0));	
	}

	if (tracelevel!=0)   fprintf(fptr,"\nLEQ_NLf = %f  LEQ_Lf = %f LEQ_T = %f  mean fec = %f  var_fec = %f\n", sum_leqNLf/(double)N, sum_leqLf/(double)N, sum_leqT/(double)N, gsum_sf/(double)N, (gsum2_sf - (gsum_sf*gsum_sf / (double)N)) / ((double)N - 1.0));

//COMMENT
	if ((tracelevel!=0)&&(Vs!=0.0))
	{
		if (gen==0)	fprintf(fptr,"\nVE = %f\n", VE);
		else
		{
			if (rescue[gen]==1)	fprintf(fptr,"\nVA (parents TP)  = %f\nVE (migrants)    = %f\n", VA, VE);
			else			fprintf(fptr,"\nVA (parents) = %f\n", VA);
		}
	}
	if ((tracelevel!=0)&&(Vs!=0.0))	for (i=0; i<NIND; i++)	fprintf(fptr,"i=%d	p_a=%f	pm_sf=%f\n", i, p_a[i], pm_sf[i]);
}

/* ********************************************************************* */

distribution_qsh()
{
	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	if (initialgen[k][l] != (-99)) 
	{
		/* distribution q */
		if (q[k][l] == 0.0)				accum (&q_ng_00[lin][gen/tgen],  1.0);
		else if ( (q[k][l] > 0.0) && (q[k][l] < 0.1) )	accum (&q_ng_00_01[lin][gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.1) && (q[k][l] < 0.2) )	accum (&q_ng_01_02[lin][gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.2) && (q[k][l] < 0.3) )	accum (&q_ng_02_03[lin][gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.3) && (q[k][l] < 0.4) )	accum (&q_ng_03_04[lin][gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.4) && (q[k][l] < 0.5) )	accum (&q_ng_04_05[lin][gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.5) && (q[k][l] < 0.6) )	accum (&q_ng_05_06[lin][gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.6) && (q[k][l] < 0.7) )	accum (&q_ng_06_07[lin][gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.7) && (q[k][l] < 0.8) )	accum (&q_ng_07_08[lin][gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.8) && (q[k][l] < 0.9) )	accum (&q_ng_08_09[lin][gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.9) && (q[k][l] < 1.0) )	accum (&q_ng_09_10[lin][gen/tgen],  1.0);
		else if (q[k][l] == 1.0)			accum (&q_ng_10[lin][gen/tgen],  1.0);

		/* distribution s */
		if (q[k][l] == 0.0) 
		{
			if (s[k][l] == 0.0)						accum (&s_ng_0_lost[lin][gen/tgen],  1.0);		
			else if ( (s[k][l] < 0.0) && (s[k][l] > (-0.000001)) )		accum (&s_ng_0_106_lost[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.000001)) && (s[k][l] > (-0.0001)) )	accum (&s_ng_106_104_lost[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.0001)) && (s[k][l] > (-0.01)) )	accum (&s_ng_104_102_lost[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.01)) && (s[k][l] > (-0.1)) )		accum (&s_ng_102_101_lost[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.1)) && (s[k][l] > (-0.2)) )		accum (&s_ng_01_02_lost[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.2)) && (s[k][l] > (-0.4)) )		accum (&s_ng_02_04_lost[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.4)) && (s[k][l] > (-0.6)) )		accum (&s_ng_04_06_lost[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.6)) && (s[k][l] > (-1.0)) )		accum (&s_ng_06_10_lost[lin][gen/tgen],  1.0);
			else if (s[k][l] == (-1.0))					accum (&s_ng_10_lost[lin][gen/tgen],  1.0);
		}
		else if (q[k][l] == 1.0) 
		{
			if (s[k][l] == 0.0)						accum (&s_ng_0_fix[lin][gen/tgen],  1.0);		
			else if ( (s[k][l] < 0.0) && (s[k][l] > (-0.000001)) )		accum (&s_ng_0_106_fix[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.000001)) && (s[k][l] > (-0.0001)) )	accum (&s_ng_106_104_fix[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.0001)) && (s[k][l] > (-0.01)) )	accum (&s_ng_104_102_fix[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.01)) && (s[k][l] > (-0.1)) )		accum (&s_ng_102_101_fix[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.1)) && (s[k][l] > (-0.2)) )		accum (&s_ng_01_02_fix[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.2)) && (s[k][l] > (-0.4)) )		accum (&s_ng_02_04_fix[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.4)) && (s[k][l] > (-0.6)) )		accum (&s_ng_04_06_fix[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.6)) && (s[k][l] > (-1.0)) )		accum (&s_ng_06_10_fix[lin][gen/tgen],  1.0);
			else if (s[k][l] == (-1.0))					accum (&s_ng_10_fix[lin][gen/tgen],  1.0);
		}
		else
		{
			if (s[k][l] == 0.0)						accum (&s_ng_0[lin][gen/tgen],  1.0);		
			else if ( (s[k][l] < 0.0) && (s[k][l] > (-0.000001)) )		accum (&s_ng_0_106[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.000001)) && (s[k][l] > (-0.0001)) )	accum (&s_ng_106_104[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.0001)) && (s[k][l] > (-0.01)) )	accum (&s_ng_104_102[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.01)) && (s[k][l] > (-0.1)) )		accum (&s_ng_102_101[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.1)) && (s[k][l] > (-0.2)) )		accum (&s_ng_01_02[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.2)) && (s[k][l] > (-0.4)) )		accum (&s_ng_02_04[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.4)) && (s[k][l] > (-0.6)) )		accum (&s_ng_04_06[lin][gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.6)) && (s[k][l] > (-1.0)) )		accum (&s_ng_06_10[lin][gen/tgen],  1.0);
			else if (s[k][l] == (-1.0))					accum (&s_ng_10[lin][gen/tgen],  1.0);
		}

		/* distribution h */
		if (q[k][l] == 0.0)
		{
			if ( (hs[k][l] > 0.0) && (hs[k][l] <= 0.2) )		accum (&hs_ng_00_02_lost[lin][gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.2) && (hs[k][l] <= 0.4) )	accum (&hs_ng_02_04_lost[lin][gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.4) && (hs[k][l] <= 0.6) )	accum (&hs_ng_04_06_lost[lin][gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.6) && (hs[k][l] <= 0.8) )	accum (&hs_ng_06_08_lost[lin][gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.8) && (hs[k][l] <= 1.0) )	accum (&hs_ng_08_10_lost[lin][gen/tgen],  1.0);
		}
		else if (q[k][l] == 1.0)
		{
			if ( (hs[k][l] > 0.0) && (hs[k][l] <= 0.2) )		accum (&hs_ng_00_02_fix[lin][gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.2) && (hs[k][l] <= 0.4) )	accum (&hs_ng_02_04_fix[lin][gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.4) && (hs[k][l] <= 0.6) )	accum (&hs_ng_04_06_fix[lin][gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.6) && (hs[k][l] <= 0.8) )	accum (&hs_ng_06_08_fix[lin][gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.8) && (hs[k][l] <= 1.0) )	accum (&hs_ng_08_10_fix[lin][gen/tgen],  1.0);
		}
		else
		{
			if ( (hs[k][l] > 0.0) && (hs[k][l] <= 0.2) )		accum (&hs_ng_00_02[lin][gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.2) && (hs[k][l] <= 0.4) )	accum (&hs_ng_02_04[lin][gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.4) && (hs[k][l] <= 0.6) )	accum (&hs_ng_04_06[lin][gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.6) && (hs[k][l] <= 0.8) )	accum (&hs_ng_06_08[lin][gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.8) && (hs[k][l] <= 1.0) )	accum (&hs_ng_08_10[lin][gen/tgen],  1.0);
		}
	}
}

/* ********************************************************************* */

distribution_w()
{
	for (mig=1; mig<=events; mig++)
	{
		/* All annotated replicates */
		accum (&summ_w[lin][mig], mean_w[mig]);

		if (mean_w[mig] == 0.0)						accum (&w_00[lin][mig],  1.0);
		else if ( (mean_w[mig] > 0.0) && (mean_w[mig] < 0.1) )		accum (&w_00_01[lin][mig],  1.0);
		else if ( (mean_w[mig] >= 0.1) && (mean_w[mig] < 0.2) )		accum (&w_01_02[lin][mig],  1.0);
		else if ( (mean_w[mig] >= 0.2) && (mean_w[mig] < 0.3) )		accum (&w_02_03[lin][mig],  1.0);
		else if ( (mean_w[mig] >= 0.3) && (mean_w[mig] < 0.4) )		accum (&w_03_04[lin][mig],  1.0);
		else if ( (mean_w[mig] >= 0.4) && (mean_w[mig] < 0.5) )		accum (&w_04_05[lin][mig],  1.0);
		else if ( (mean_w[mig] >= 0.5) && (mean_w[mig] < 0.6) )		accum (&w_05_06[lin][mig],  1.0);
		else if ( (mean_w[mig] >= 0.6) && (mean_w[mig] < 0.7) )		accum (&w_06_07[lin][mig],  1.0);
		else if ( (mean_w[mig] >= 0.7) && (mean_w[mig] < 0.8) )		accum (&w_07_08[lin][mig],  1.0);
		else if ( (mean_w[mig] >= 0.8) && (mean_w[mig] < 0.9) )		accum (&w_08_09[lin][mig],  1.0);
		else if ( (mean_w[mig] >= 0.9) && (mean_w[mig] < 1.0) )		accum (&w_09_10[lin][mig],  1.0);
		else if (mean_w[mig] == 1.0)					accum (&w_10[lin][mig],  1.0);

		/* Extinct replicates after annotation */
		if (repmark == 1)
		{
			accum (&ext_summ_w[lin][mig], mean_w[mig]);

			if (mean_w[mig] == 0.0)						accum (&ext_w_00[lin][mig],  1.0);
			else if ( (mean_w[mig] > 0.0) && (mean_w[mig] < 0.1) )		accum (&ext_w_00_01[lin][mig],  1.0);
			else if ( (mean_w[mig] >= 0.1) && (mean_w[mig] < 0.2) )		accum (&ext_w_01_02[lin][mig],  1.0);
			else if ( (mean_w[mig] >= 0.2) && (mean_w[mig] < 0.3) )		accum (&ext_w_02_03[lin][mig],  1.0);
			else if ( (mean_w[mig] >= 0.3) && (mean_w[mig] < 0.4) )		accum (&ext_w_03_04[lin][mig],  1.0);
			else if ( (mean_w[mig] >= 0.4) && (mean_w[mig] < 0.5) )		accum (&ext_w_04_05[lin][mig],  1.0);
			else if ( (mean_w[mig] >= 0.5) && (mean_w[mig] < 0.6) )		accum (&ext_w_05_06[lin][mig],  1.0);
			else if ( (mean_w[mig] >= 0.6) && (mean_w[mig] < 0.7) )		accum (&ext_w_06_07[lin][mig],  1.0);
			else if ( (mean_w[mig] >= 0.7) && (mean_w[mig] < 0.8) )		accum (&ext_w_07_08[lin][mig],  1.0);
			else if ( (mean_w[mig] >= 0.8) && (mean_w[mig] < 0.9) )		accum (&ext_w_08_09[lin][mig],  1.0);
			else if ( (mean_w[mig] >= 0.9) && (mean_w[mig] < 1.0) )		accum (&ext_w_09_10[lin][mig],  1.0);
			else if (mean_w[mig] == 1.0)					accum (&ext_w_10[lin][mig],  1.0);
		}
		/* Surviving replicates */
		else
		{
			accum (&sur_summ_w[lin][mig], mean_w[mig]);

			if (mean_w[mig] == 0.0)						accum (&sur_w_00[lin][mig],  1.0);
			else if ( (mean_w[mig] > 0.0) && (mean_w[mig] < 0.1) )		accum (&sur_w_00_01[lin][mig],  1.0);
			else if ( (mean_w[mig] >= 0.1) && (mean_w[mig] < 0.2) )		accum (&sur_w_01_02[lin][mig],  1.0);
			else if ( (mean_w[mig] >= 0.2) && (mean_w[mig] < 0.3) )		accum (&sur_w_02_03[lin][mig],  1.0);
			else if ( (mean_w[mig] >= 0.3) && (mean_w[mig] < 0.4) )		accum (&sur_w_03_04[lin][mig],  1.0);
			else if ( (mean_w[mig] >= 0.4) && (mean_w[mig] < 0.5) )		accum (&sur_w_04_05[lin][mig],  1.0);
			else if ( (mean_w[mig] >= 0.5) && (mean_w[mig] < 0.6) )		accum (&sur_w_05_06[lin][mig],  1.0);
			else if ( (mean_w[mig] >= 0.6) && (mean_w[mig] < 0.7) )		accum (&sur_w_06_07[lin][mig],  1.0);
			else if ( (mean_w[mig] >= 0.7) && (mean_w[mig] < 0.8) )		accum (&sur_w_07_08[lin][mig],  1.0);
			else if ( (mean_w[mig] >= 0.8) && (mean_w[mig] < 0.9) )		accum (&sur_w_08_09[lin][mig],  1.0);
			else if ( (mean_w[mig] >= 0.9) && (mean_w[mig] < 1.0) )		accum (&sur_w_09_10[lin][mig],  1.0);
			else if (mean_w[mig] == 1.0)					accum (&sur_w_10[lin][mig],  1.0);
		}
	}
}

/* ********************************************************************* */

dumpphenotypes()
{
	if (tracelevel==0)   return (0);
	fprintf(fptr,"\nFecundity of parents after mutation\n\n");

	for (i=0; i<NIND; i++)		
	{
		if ((rescue[gen]==1)&&(mGEND==0))
		{
			if ((i%2==0)&&(i<NIND-NMIG))	fprintf(fptr,"%d fem %f\n", i, pm_sf[i]); /* females */
			else 				fprintf(fptr,"%d mal %f\n", i, pm_sf[i]); /* males and migrants */
		}
		else
		{
			if (i%2==0)			fprintf(fptr,"%d fem %f\n", i, pm_sf[i]); /* females */
			else 				fprintf(fptr,"%d mal %f\n", i, pm_sf[i]); /* males */
		}
	}
}

/* ********************************************************************* */

additivevariance()
{
	double d_s, alfa_s;

	VA=0.0;

	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	{
		d_s = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);

		if (s[k][l]>=0.0)	alfa_s = (s[k][l]/2.0) + ( d_s * (1.0 - 2.0*q[k][l]) );
		else			alfa_s = (-s[k][l]/2.0) + ( d_s * (2.0*q[k][l] - 1.0) );

		VA += 2.0 * alfa_s * alfa_s * q[k][l] * (1.0 - q[k][l]);
	} 

//COMMENT
	if (tracelevel!=0)	fprintf(fptr, "\nCOMMENT - additive variance (parents)\nVA = %f\n", VA);
}

/* ********************************************************************* */

coancestry_matrix()
{
	int j, N;

	if (gen==0)
	{
		for (i=0; i<NIND; i++)
		for (j=0; j<NIND; j++)
		{
			if (i == j)	mat[i][j] = 0.5;
			else		mat[i][j] = 0.0;
		}
	}
	else
	{
		for (i=0; i<NIND; i++)
		for (j=0; j<NIND; j++)
		{
			if (i == j)	mat[i][j] = 0.5 * (1.0 + smat[father[i]][mother[i]]);
			else		mat[i][j] = 0.25 * (smat[father[i]][father[j]] + smat[father[i]][mother[j]] + smat[mother[i]][father[j]] + smat[mother[i]][mother[j]]);
		}
	}

	for (i=0; i<NIND; i++)	F[i] = 2.0*mat[i][i]-1.0;

	for (i=0; i<NIND; i++)
	for (j=0; j<NIND; j++)
	{
		smat[i][j] = mat[i][j];
	}

	//exclude migrants from average values
	if (rescue[gen]==1)	N = NIND - NMIG;
	else 			N = NIND;

	for (i=0; i<N; i++)
	for (j=0; j<N; j++)
	{
		accum(&fp[lin][gen], mat[i][j]); 
		accum(&Fp[lin][gen], (2.0*mat[i][i]-1.0)); 
	}

//COMMENT
//	if (tracelevel!=0)
//	{
//		fprintf(fptr, "\nCOMMENT - coancestry_matrix\n");   
//		for (i=0; i<NIND; i++)
//		{
//			fprintf(fptr, "gen=%d   i=%d   fat=%d   mot=%d   Fp=%f\n", gen, i, father[i], mother[i], 2.0*mat[i][i]-1.0);
//		}
//	}

	if (tracelevel!=0)	dumpparents(mat);
}

/* ********************************************************************* */

dumpparents(double mat[][NN])
{
	if (tracelevel==0)   return (0);

	fprintf(fptr,"\nParents\n\n");
	for (i=0; i<NIND; i++)   fprintf(fptr,"i=%d\tmother=%d\tfather=%d\tpm_sf=%f\tFp=%f\n", i, mother[i], father[i], pm_sf[i], 2.0*mat[i][i]-1.0);
}

/* ********************************************************************* */

int fecunditysearch_mo(double *array, int size)
{
	int i, p;
	double r;

	r = uniform(); 
	if ( (r >= 0.0) && (r <= array[0]) )	p = 0; 
	else for (i=2; i<size; i+=2)		if ( (r > array[i-2]) && (r <= array[i]) )	p = i; 

	return(p);
}

/* ********************************************************************* */

int fecunditysearch_fa(double *array, int size)
{
	int i, p;
	double r;

	r = uniform(); 

	if ((rescue[gen]==1)&&(mGEND==0))
	{
		if ( (r >= 0.0) && (r <= array[1]) )	p = 1; 
		else
		{
			for (i=3; i<size-NMIG; i+=2)	if ( (r > array[i-2]) && (r <= array[i]) )	p = i; /* males - no migrants */
			for (i=size-NMIG; i<size; i++)	if ( (r > array[i-1]) && (r <= array[i]) )	p = i; /* males - migrants */
		}
	}
	else
	{
		if ( (r >= 0.0) && (r <= array[1]) )	p = 1; 
		else for (i=3; i<size; i+=2)	if ( (r > array[i-2]) && (r <= array[i]) )	p = i; 
	}

	return(p);
}

/* ********************************************************************* */

mating()
{
	int family[NN], mo, fa, NOFF, numberrecs, nr, rndk, rndl, marker;
	int EE[MM], FF[MM], ncrorec[MM], pointrec[MM][31];

	double rnd, sk2, cum_pm_sf_fem[NN], cum_pm_sf_mal[NN], pm_sfo[NN];
	double family_sum=0.0, family_sum2=0.0, cum_fem=0.0, cum_mal=0.0;
	double matmm=0.0, matml=0.0, matll=0.0;

	/* ************** KEEP GENOTYPES OF PARENTS ************** */ 

	for (i=0; i<NIND; i++)
	for (k=0; k<NCRO; k++)
	{
		sp[i][k][0]=g[i][k][0];
		sp[i][k][1]=g[i][k][1];
	}

	/* *********** ACCUMULATE psf FECUNDITY VALUES ********** */

	if ((rescue[gen]==1)&&(mGEND==0))
	{
		/* if there is migration, all migrants are males */
		for (i=0; i<NIND; i++)
		{
			if ((i%2==0)&&(i<NIND-NMIG))	/* females */
			{
				cum_fem += pm_sf[i];
				cum_pm_sf_fem[i] = cum_fem;
			}
			else				/* males and migrants */
			{
				cum_mal += pm_sf[i];
				cum_pm_sf_mal[i] = cum_mal;
			}
		}
	}
	else
	{
		/* if there is migration, migrants can be males or females */
		for (i=0; i<NIND; i++)
		{
			if (i%2==0)	/* females */
			{
				cum_fem += pm_sf[i];
				cum_pm_sf_fem[i] = cum_fem;
			}
			else		/* males */
			{
				cum_mal += pm_sf[i];
				cum_pm_sf_mal[i] = cum_mal;
			}
		}
	}

	if ((rescue[gen]==1)&&(mGEND==0))
	{
		for (i=0; i<NIND; i++)
		{
			if ((i%2==0)&&(i<NIND-NMIG))	cum_pm_sf_fem[i] = cum_pm_sf_fem[i] / cum_fem;	/* females */
			else				cum_pm_sf_mal[i] = cum_pm_sf_mal[i] / cum_mal;	/* males and migrants */
		}
	}
	else
	{
		for (i=0; i<NIND; i++)
		{
			if (i%2==0)	cum_pm_sf_fem[i] = cum_pm_sf_fem[i] / cum_fem;	/* females */
			else		cum_pm_sf_mal[i] = cum_pm_sf_mal[i] / cum_mal;	/* males */
		}
	}

//COMMENT
	if (tracelevel!=0)
	{
		fprintf (fptr,"\ncum_fem = %f   cum_mal = %f\n", cum_fem, cum_mal); 
		if ((rescue[gen]==1)&&(mGEND==0)) for (i=0; i<NIND-NMIG; i+=2)  fprintf (fptr,"i=%d cum_pm_sf_fem[%d] = %f\n",i,i,cum_pm_sf_fem[i]);
		else			   	  for (i=0; i<NIND; i+=2)	fprintf (fptr,"i=%d cum_pm_sf_fem[%d] = %f\n",i,i,cum_pm_sf_fem[i]);	

		fprintf (fptr,"\n"); 
		if ((rescue[gen]==1)&&(mGEND==0))
		{
			for (i=1; i<NIND-NMIG; i+=2)	fprintf (fptr,"i=%d cum_pm_sf_mal[%d] = %f\n",  i, i, cum_pm_sf_mal[i]);
			for (i=NIND-NMIG; i<NIND; i++)	fprintf (fptr,"i=%d cum_pm_sf_mal[%d] = %f\n",  i, i, cum_pm_sf_mal[i]);
		}
		else	for (i=1; i<NIND; i+=2)		fprintf (fptr,"i=%d cum_pm_sf_mal[%d] = %f\n",  i, i, cum_pm_sf_mal[i]);
	}

	/* ****************** ROUNDS OF MATING ***************** */

	if (tracelevel!=0)	fprintf (fptr,"\n\nRounds of mating and genotypic values of progeny before mutation\n\n");

	/* ************ offspring ************ */

	if ( ((lin==0)&&(gen==(tLINES-1))) || ((lin==4)&&(gen==(tsubL-1))) )
	{
		//generate enought offspring to cover each starting line
		if (NINDTP<NINDAL)
		{
			NOFF = NINDAL;
			if (tracelevel!=0)	fprintf(fptr,"NINDTP%d<NINDAL%d\n", NINDTP, NINDAL);
		}
		else	NOFF = NIND;	
	}
	else
	{
		if (rescue[gen]==1)	NOFF = NIND - NMIG;
		else			NOFF = NIND;
	}
	
	for (i=0; i<NIND; i++)   family[i] = 0;
	for (i=0; i<NOFF; i++)
	{
		/* ************ parents ************ */

		if ((rescue[gen]==1)&&(mGEND==0))
		{
			mo = fecunditysearch_mo(cum_pm_sf_fem, NIND-NMIG);
			fa = fecunditysearch_fa(cum_pm_sf_mal, NIND);		
		}
		else
		{
			mo = fecunditysearch_mo(cum_pm_sf_fem, NIND);
			fa = fecunditysearch_fa(cum_pm_sf_mal, NIND);		
		}

		mother[i] = mo;
		father[i] = fa;
		family[mo] += 1;
		family[fa] += 1;

		if (tracelevel!=0)	fprintf (fptr,"Progeny %d mo = %d    fa = %d", i, mo, fa);

		/* ******* OUTFILE - mating ******* */

		if (rescue[gen]==1)
		{
			if ((mo>=NIND-NMIG)&&(fa>=NIND-NMIG))	matmm += 1.0;
			else if ((mo<NIND-NMIG)&&(fa<NIND-NMIG))	matll += 1.0;
			else	matml += 1.0;
		}

		/* ******* Free recombination ******* */

		if(L==99.0)
		{
			for (k=0; k<NCRO; k++)
			{
			   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
				FF[k] = ~EE[k];
			   	g[i][k][0]=((EE[k]&sp[mo][k][0])|(FF[k]&sp[mo][k][1]));
			}
// COMMENT
//			if (tracelevel!=0)   fprintf (fptr,"  i=%d EE[0]=%d EE[1]=%d EE[2]=%d sf00=%d sf01=%d sf10=%d sf11=%d sf20=%d sf21=%d g00=%d g10=%d g20=%d \n", i, EE[0], EE[1], EE[2], sp[mo][0][0], sp[mo][0][1], sp[mo][1][0], sp[mo][1][1], sp[mo][2][0], sp[mo][2][1], g[i][0][0], g[i][1][0], g[i][2][0]);

			for (k=0; k<NCRO; k++)
			{
			   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
			   	FF[k] = ~EE[k];
			   	g[i][k][1]=((EE[k]&sp[fa][k][0])|(FF[k]&sp[fa][k][1]));
			}
// COMMENT
//			if (tracelevel!=0)   fprintf (fptr,"i=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", i, EE[0], EE[1], EE[2], sp[fa][0][0], sp[fa][0][1], sp[fa][1][0], sp[fa][1][1], sp[fa][2][0], sp[fa][2][1], g[i][0][1], g[i][1][1], g[i][2][1]);

		}

		/* **** Restricted recombination **** */ 

		else
		{
			/* Chromosome from mother */

			numberrecs = recombinationnumber();
			for (k=0; k<NCRO; k++)
			{
				ncrorec[k] = 0;
				for (l=0; l<NLOCI; l++)  pointrec[k][l] = 0;
			}
			for (nr=0; nr<numberrecs; nr++)
			{
				rndk = (int)(uniform()*NCRO);
				rndl = (int)(uniform()*NLOCI);
				ncrorec[rndk] = 1;
				pointrec[rndk][rndl] = 1;
			}

			marker = 1;

			for (k=0; k<NCRO; k++)
			{
				EE[k]=0;
				if (ncrorec[k] == 0)
				{
					if (marker==(-1))
					{
						EE[k] = ~EE[k];
					}
				}
				else
				{
					for (l=0; l<NLOCI; l++)
			      		{
						if (pointrec[k][l] == 0)
						{
							if (marker==(-1))  EE[k] = EE[k] | RM[l];
						}
						else
						{
							if (marker==1)
							{
								EE[k] = EE[k] | RM[l];
								marker = marker * (-1);
							}
							else
							{
								marker = marker * (-1);
							}
						}
					}
				}
			}

			rnd = uniform();
			for (k=0; k<NCRO; k++)
			{
				if (rnd < 0.5) 	EE[k] = ~EE[k];
				FF[k] = ~EE[k];
				g[i][k][0]=((EE[k]&sp[mo][k][0])|(FF[k]&sp[mo][k][1]));
			}

// COMMENT
//			if (tracelevel!=0)   fprintf (fptr,"i=%d EE[0]=%d EE[1]=%d EE[2]=%d sf00=%d sf01=%d sf10=%d sf11=%d sf20=%d sf21=%d g00=%d g10=%d g20=%d \n", i, EE[0], EE[1], EE[2], sp[mo][0][0], sp[mo][0][1], sp[mo][1][0], sp[mo][1][1], sp[mo][2][0], sp[mo][2][1], g[i][0][0], g[i][1][0], g[i][2][0]);
/*
			if (tracelevel!=0)
			{
				fprintf (fptr,"EE\n");
				for (k=0; k<50; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((EE[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
				fprintf (fptr,"FF\n");
				for (k=0; k<50; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((FF[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
			}
*/

			/* Chromosome from father */

			numberrecs = recombinationnumber();

			for (k=0; k<NCRO; k++)
			{
				ncrorec[k] = 0;
				for (l=0; l<NLOCI; l++)  pointrec[k][l] = 0;
			}
			for (nr=0; nr<numberrecs; nr++)
			{
				rndk = (int)(uniform()*NCRO);
				rndl = (int)(uniform()*NLOCI);
				ncrorec[rndk] = 1;
				pointrec[rndk][rndl] = 1;
			}

			marker = 1;

			for (k=0; k<NCRO; k++)
			{
				EE[k]=0;
				if (ncrorec[k] == 0)
				{
					if (marker==(-1))
					{
						EE[k] = ~EE[k];
					}
				}
				else
				{
					for (l=0; l<NLOCI; l++)
			      		{
						if (pointrec[k][l] == 0)
						{
							if (marker==(-1))  EE[k] = EE[k] | RM[l];
						}
						else
						{
							if (marker==1)
							{
								EE[k] = EE[k] | RM[l];
								marker = marker * (-1);
							}
							else
							{
								marker = marker * (-1);
							}
						}
					}
				}
			}

			rnd = uniform();
			for (k=0; k<NCRO; k++)
			{
				if (rnd < 0.5) 	EE[k] = ~EE[k];
				FF[k] = ~EE[k];
				g[i][k][1]=((EE[k]&sp[fa][k][0])|(FF[k]&sp[fa][k][1]));
			}
// COMMENT
//			if (tracelevel!=0)   fprintf (fptr,"i=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", i, EE[0], EE[1], EE[2], sp[fa][0][0], sp[fa][0][1], sp[fa][1][0], sp[fa][1][1], sp[fa][2][0], sp[fa][2][1], g[i][0][1], g[i][1][1], g[i][2][1]);

/*			if (tracelevel!=0)
			{
				fprintf (fptr,"EE\n");
				for (k=0; k<50; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((EE[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
				fprintf (fptr,"FF\n");
				for (k=0; k<50; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((FF[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
			}
*/
		}

		/* **** Genotypic value of offspring **** */

		pm_sfo[i]=addedfixed_sf;

		for (k=0; k<NCRO; k++)
		for (l=1; l<NLOCI; l++)
		if (initialgen[k][l] != (-99))
		{
	    		if (((g[i][k][0] & RM[l])==RM[l])&&((g[i][k][1] & RM[l])==RM[l]))	pm_sfo[i] *= (1.0 + s[k][l]);
			else    if (((g[i][k][0] & RM[l])!=RM[l])&&((g[i][k][1] & RM[l])!=RM[l])) 	/* AA */;
			else	pm_sfo[i] *= (1.0 + (s[k][l]*hs[k][l]));
		}

		if (tracelevel!=0)	fprintf (fptr,"  pm_sf = %f\n", pm_sfo[i]);
	}

	/* ************** VARIANCE OF FAMILY SIZE ************** */

	for (i=0; i<NIND; i++)
	{
		family_sum += family[i];
		family_sum2 += (family[i] * family[i]);
	}

	sk2 = ((family_sum2-(family_sum*family_sum/NIND))/(NIND-1.0)); 

	accum(&SK2[lin][gen], sk2);

	if (tracelevel!=0)	fprintf (fptr,"\nContributions from parents\n\n");
	if (tracelevel!=0)	for (i=0; i<NIND; i++)	fprintf (fptr,"family(%d) = %d\n", i, family[i]);
	if (tracelevel!=0)	fprintf (fptr,"\nSK2 = %5.3f\n", sk2);
	
	/* *** OUTFILE - mating *** */

	if ((lin==1)||(lin==2)||(lin==5)||(lin==6))	summ_matings(family, matmm, matml, matll);

	/* ************ */
	NIND = NOFF;
	for (i=0; i<NIND; i++)	pm_sf[i] = pm_sfo[i];
}

/* ********************************************************************* */

summ_matings(int family[], double matmm, double matml, double matll)
{
	double sum_sf_L=0.0, sum_sf_mig=0.0;

	for (i=0; i<NIND; i++)
	{
		if ( (pm_sf[i] > 0.0) && (pm_sf[i] <= 0.2) )		accum (&off_0_02[lin][gen],  family[i]);
		if ( (pm_sf[i] > 0.2) && (pm_sf[i] <= 0.4) )		accum (&off_02_04[lin][gen],  family[i]);
		if ( (pm_sf[i] > 0.4) && (pm_sf[i] <= 0.6) )		accum (&off_04_06[lin][gen],  family[i]);
		if ( (pm_sf[i] > 0.6) && (pm_sf[i] <= 0.8) )		accum (&off_06_08[lin][gen],  family[i]);
		if ( (pm_sf[i] > 0.8) && (pm_sf[i] <= 1.0) )		accum (&off_08_1[lin][gen],  family[i]);

		if (rescue[gen]==1)
		{
			if (i<NIND-NMIG)	sum_sf_L += pm_sf[i];	
			else			sum_sf_mig += pm_sf[i];
		}
	}

	if (rescue[gen]==1)
	{
		accum (&FEC_line[lin][gen],  sum_sf_L/(NIND-NMIG));
		accum (&FEC_mig[lin][gen],  sum_sf_mig/(NMIG));

		accum (&mate_mm[lin][gen],  matmm);
		accum (&mate_ml[lin][gen],  matml);
		accum (&mate_ll[lin][gen],  matll);
	}
}

/* ********************************************************************* */

int recombinationnumber ()
{
	int r;
	if ((L < normalthreshold) && (exp(-L) != 0.0) )
	{
		r = poisson(lastinrecombinantpoissontable, recombinantpoissontable);
	}
	else r = (int)normal(L, sqrt(L));
	return(r);
}

/* ********************************************************************* */

dumpoffspring()
{
	if (tracelevel==0)   return (0);
	fprintf(fptr,"\n Offspring before mutation \n");

	if ( ((lin==0)&&(gen==(tLINES-1))) || ((lin==4)&&(gen==(tsubL-1))) )
	{
		if (NINDTP<NINDAL)
		{
			NIND = NINDAL;
			if (tracelevel!=0)	fprintf(fptr,"NINDTP%d<NINDAL%d\n", NINDTP, NINDAL);
		}
	}
	
	for (i=0; i<NIND; i++)
	{
		if (i%2==0)	fprintf(fptr,"(gf0 gf1)	%d  %d\n",g[i][0][0],g[i][0][1]); /* females */
		else		fprintf(fptr,"(gm0 gm1)	%d  %d\n",g[i][0][0],g[i][0][1]); /* males */		
	}
}

/* ********************************************************************* */

printout()
{
	/* DESCRIPTION OF LINES */

	fprintf(fgenL1,"******************************************\n************* %dR%d (line 1) *************\nLine mantained with %d individuals throughout the period analyzed.\nGenetic rescue from generation %d to end of simulation.\n******************************************", NINDTP, NINDTP, NINDTP, tLINES, NMIGTP, mINT_TP);
	fprintf(fgenL2,"******************************************\n************ %dR%d (line 2) ************\nLine mantained with %d individuals until generation %d, and %d individuals from generation %d to end of simulation.\nGenetic rescue from generation %d to end of simulation.\n******************************************", NINDTP, NINDAL, NINDTP, tLINES-1, NINDAL, tLINES, tLINES, NMIGAL, mINT_AL);
	fprintf(fgenL3,"******************************************\n************ %dNR%d (line 3) ************\nLine mantained with %d individuals until generation %d, and %d individuals from generation %d to end of simulation.\nLine without genetic rescue.\n******************************************", NINDTP, NINDAL, NINDTP, tLINES-1, NINDAL, tLINES);
	fprintf(fgensubL1,"******************************************\n************ %dRR%d (line 4, subline1) ************\nLine mantained with %d individuals throughout the period analyzed.\nNo rescue until generation %d. Retarded rescue from generation %d to end of simulation.\n******************************************", NINDTP, NINDTP, NINDTP, tsubL-1, tsubL, NMIGTP, mINT_TP);
	fprintf(fgensubL2,"******************************************\n************ %dRR%d (line 4, subline2) ************\nLine mantained with %d individuals until generation %d, and %d individuals from generation %d to end of simulation.\nNo rescue until generation %d. Retarded rescue from generation %d to end of simulation.\n******************************************", NINDTP, NINDAL, NINDTP, tsubL-1, NINDAL, tsubL, tsubL-1, tsubL, NMIGAL, mINT_AL);
	fprintf(fgensubL3,"******************************************\n************ %dNR%d (line 4, subline3) ************\nLine mantained with %d individuals throughout the period analyzed.\nLine without genetic rescue.\n******************************************", NINDTP, NINDTP, NINDTP);
	fprintf(fgensubL4,"******************************************\n************ %dNR%d (line 4, subline4) ************\nLine mantained with %d individuals until generation %d, and %d individuals from generation %d to end of simulation.\nLine without genetic rescue.\n******************************************", NINDTP, NINDAL, NINDTP, tsubL-1, NINDAL, tsubL);

	/* PARAMETERS */

	fprintf(fgenL1,"\n\nNINDTP=%d  NINDAL=%d  NMIGTP=%d  NMIGAL=%d migGENDER=%d  L=%4.2f  N.S-LOCI=%d  N.N-LOCI=%d  extinction(w)=%f  reps=%d\ngens=%d  tlines=%d  tsubL=%d tmigration=%d mig_intervals_TP=%d mig_intervals_AL=%d migrations=%d\n\n    Vs=%f   VE=%f   Opt=%f   Lambda_s=%f   Lambda_L=%f    ave_s=%f\n    beta_s=%f   alpha_s=%f  dom=%d   k_s=%f   ave_hs=%f    addedfixed_sNP=%f\n", NINDTP, NINDAL, NMIGTP, NMIGAL, mGEND, L, TOTLOCI-NCRO, NCRO, wEXT, replicates, generations, tLINES, tsubL, tLINES+preNR, mINT_TP, mINT_AL, numMIG, Vs, VE, OPT, Lambda_s, Lambda_L, ave_s, beta_s, alpha_s, dom, k_s, ave_hs, addedfixed_sNP);
	fprintf(fgenL2,"\n\nNINDTP=%d  NINDAL=%d  NMIGTP=%d  NMIGAL=%d migGENDER=%d  L=%4.2f  N.S-LOCI=%d  N.N-LOCI=%d  extinction(w)=%f  reps=%d\ngens=%d  tlines=%d  tsubL=%d tmigration=%d mig_intervals_TP=%d mig_intervals_AL=%d migrations=%d\n\n    Vs=%f   VE=%f   Opt=%f   Lambda_s=%f   Lambda_L=%f    ave_s=%f\n    beta_s=%f   alpha_s=%f  dom=%d   k_s=%f   ave_hs=%f    addedfixed_sNP=%f\n", NINDTP, NINDAL, NMIGTP, NMIGAL, mGEND, L, TOTLOCI-NCRO, NCRO, wEXT, replicates, generations, tLINES, tsubL, tLINES+preNR, mINT_TP, mINT_AL, numMIG, Vs, VE, OPT, Lambda_s, Lambda_L, ave_s, beta_s, alpha_s, dom, k_s, ave_hs, addedfixed_sNP);
	fprintf(fgenL3,"\n\nNINDTP=%d  NINDAL=%d  NMIGTP=%d  NMIGAL=%d migGENDER=%d  L=%4.2f  N.S-LOCI=%d  N.N-LOCI=%d  extinction(w)=%f  reps=%d\ngens=%d  tlines=%d  tsubL=%d tmigration=%d  mig_intervals_TP=%d mig_intervals_AL=%d migrations=%d\n\n    Vs=%f   VE=%f   Opt=%f   Lambda_s=%f   Lambda_L=%f    ave_s=%f\n    beta_s=%f   alpha_s=%f  dom=%d   k_s=%f   ave_hs=%f    addedfixed_sNP=%f\n", NINDTP, NINDAL, NMIGTP, NMIGAL, mGEND, L, TOTLOCI-NCRO, NCRO, wEXT, replicates, generations, tLINES, tsubL, tLINES+preNR, mINT_TP, mINT_AL, numMIG, Vs, VE, OPT, Lambda_s, Lambda_L, ave_s, beta_s, alpha_s, dom, k_s, ave_hs, addedfixed_sNP);
	fprintf(fgensubL1,"\n\nNINDTP=%d  NINDAL=%d  NMIGTP=%d  NMIGAL=%d migGENDER=%d  L=%4.2f  N.S-LOCI=%d  N.N-LOCI=%d  extinction(w)=%f  reps=%d\ngens=%d  tlines=%d  tsubL=%d tmigration=%d  mig_intervals_TP=%d mig_intervals_AL=%d migrations=%d\n\n    Vs=%f   VE=%f   Opt=%f   Lambda_s=%f   Lambda_L=%f    ave_s=%f\n    beta_s=%f   alpha_s=%f  dom=%d   k_s=%f   ave_hs=%f    addedfixed_sNP=%f\n", NINDTP, NINDAL, NMIGTP, NMIGAL, mGEND, L, TOTLOCI-NCRO, NCRO, wEXT, replicates, generations, tLINES, tsubL, tsubL+preNR, mINT_TP, mINT_AL, numMIG, Vs, VE, OPT, Lambda_s, Lambda_L, ave_s, beta_s, alpha_s, dom, k_s, ave_hs, addedfixed_sNP);
	fprintf(fgensubL2,"\n\nNINDTP=%d  NINDAL=%d  NMIGTP=%d  NMIGAL=%d migGENDER=%d  L=%4.2f  N.S-LOCI=%d  N.N-LOCI=%d  extinction(w)=%f  reps=%d\ngens=%d  tlines=%d  tsubL=%d tmigration=%d  mig_intervals_TP=%d mig_intervals_AL=%d migrations=%d\n\n    Vs=%f   VE=%f   Opt=%f   Lambda_s=%f   Lambda_L=%f    ave_s=%f\n    beta_s=%f   alpha_s=%f  dom=%d   k_s=%f   ave_hs=%f    addedfixed_sNP=%f\n", NINDTP, NINDAL, NMIGTP, NMIGAL, mGEND, L, TOTLOCI-NCRO, NCRO, wEXT, replicates, generations, tLINES, tsubL, tsubL+preNR, mINT_TP, mINT_AL, numMIG, Vs, VE, OPT, Lambda_s, Lambda_L, ave_s, beta_s, alpha_s, dom, k_s, ave_hs, addedfixed_sNP);
	fprintf(fgensubL3,"\n\nNINDTP=%d  NINDAL=%d  NMIGTP=%d  NMIGAL=%d migGENDER=%d  L=%4.2f  N.S-LOCI=%d  N.N-LOCI=%d  extinction(w)=%f  reps=%d\ngens=%d  tlines=%d  tsubL=%d tmigration=%d  mig_intervals_TP=%d mig_intervals_AL=%d migrations=%d\n\n    Vs=%f   VE=%f   Opt=%f   Lambda_s=%f   Lambda_L=%f    ave_s=%f\n    beta_s=%f   alpha_s=%f  dom=%d   k_s=%f   ave_hs=%f    addedfixed_sNP=%f\n", NINDTP, NINDAL, NMIGTP, NMIGAL, mGEND, L, TOTLOCI-NCRO, NCRO, wEXT, replicates, generations, tLINES, tsubL, tsubL+preNR, mINT_TP, mINT_AL, numMIG, Vs, VE, OPT, Lambda_s, Lambda_L, ave_s, beta_s, alpha_s, dom, k_s, ave_hs, addedfixed_sNP);
	fprintf(fgensubL4,"\n\nNINDTP=%d  NINDAL=%d  NMIGTP=%d  NMIGAL=%d migGENDER=%d  L=%4.2f  N.S-LOCI=%d  N.N-LOCI=%d  extinction(w)=%f  reps=%d\ngens=%d  tlines=%d  tsubL=%d tmigration=%d  mig_intervals_TP=%d mig_intervals_AL=%d migrations=%d\n\n    Vs=%f   VE=%f   Opt=%f   Lambda_s=%f   Lambda_L=%f    ave_s=%f\n    beta_s=%f   alpha_s=%f  dom=%d   k_s=%f   ave_hs=%f    addedfixed_sNP=%f\n", NINDTP, NINDAL, NMIGTP, NMIGAL, mGEND, L, TOTLOCI-NCRO, NCRO, wEXT, replicates, generations, tLINES, tsubL, tsubL+preNR, mINT_TP, mINT_AL, numMIG, Vs, VE, OPT, Lambda_s, Lambda_L, ave_s, beta_s, alpha_s, dom, k_s, ave_hs, addedfixed_sNP);

	fprintf(fgenL1,"\nLEQ(2dpq) = %f   LEQ_L(2dpq) = %f   LEQ_NL(2dpq) = %f   TLEQ(sq) = %f  NSEGLOCNP=%d\n", LEQ_NP, LEQ_L_NP, LEQ_NL_NP, TLEQ_NP, NSEGLOCNP);
	fprintf(fgenL2,"\nLEQ(2dpq) = %f   LEQ_L(2dpq) = %f   LEQ_NL(2dpq) = %f   TLEQ(sq) = %f  NSEGLOCNP=%d\n", LEQ_NP, LEQ_L_NP, LEQ_NL_NP, TLEQ_NP, NSEGLOCNP);
	fprintf(fgenL3,"\nLEQ(2dpq) = %f   LEQ_L(2dpq) = %f   LEQ_NL(2dpq) = %f   TLEQ(sq) = %f  NSEGLOCNP=%d\n", LEQ_NP, LEQ_L_NP, LEQ_NL_NP, TLEQ_NP, NSEGLOCNP);
	fprintf(fgensubL1,"\nLEQ(2dpq) = %f   LEQ_L(2dpq) = %f   LEQ_NL(2dpq) = %f   TLEQ(sq) = %f  NSEGLOCNP=%d\n", LEQ_NP, LEQ_L_NP, LEQ_NL_NP, TLEQ_NP, NSEGLOCNP);
	fprintf(fgensubL2,"\nLEQ(2dpq) = %f   LEQ_L(2dpq) = %f   LEQ_NL(2dpq) = %f   TLEQ(sq) = %f  NSEGLOCNP=%d\n", LEQ_NP, LEQ_L_NP, LEQ_NL_NP, TLEQ_NP, NSEGLOCNP);
	fprintf(fgensubL3,"\nLEQ(2dpq) = %f   LEQ_L(2dpq) = %f   LEQ_NL(2dpq) = %f   TLEQ(sq) = %f  NSEGLOCNP=%d\n", LEQ_NP, LEQ_L_NP, LEQ_NL_NP, TLEQ_NP, NSEGLOCNP);
	fprintf(fgensubL4,"\nLEQ(2dpq) = %f   LEQ_L(2dpq) = %f   LEQ_NL(2dpq) = %f   TLEQ(sq) = %f  NSEGLOCNP=%d\n", LEQ_NP, LEQ_L_NP, LEQ_NL_NP, TLEQ_NP, NSEGLOCNP);


	/* ******* NINDTP_R_NINDTP (L1) ******* */

	fprintf(fgenL1,"\ngen     Wf     repvarWf     p_a     var(p_a)    fp       Fp      LEQ_NL     LEQ_L     LEQ_T      SK2       Hw      repvarHw\n");

	for (gen=0; gen<tLINES; gen++)	/* TP */
	{
		fprintf(fgenL1, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", gen, accmean(&gmean_sf[0][gen]), variance(&gmean_sf[0][gen]), accmean(&gmean_a[0][gen]), accmean(&gvar_a[0][gen]), accmean(&fp[0][gen]), accmean(&Fp[0][gen]), accmean(&LEQ_NLf[0][gen]), accmean(&LEQ_Lf[0][gen]), accmean(&LEQ_T[0][gen]), accmean(&SK2[0][gen]), accmean(&Hw[0][gen]), variance(&Hw[0][gen]));
	}

	for (gen=tLINES; gen<=generations; gen++)  /* Formation of lines */
	{
		fprintf(fgenL1, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", gen, accmean(&gmean_sf[1][gen]), variance(&gmean_sf[1][gen]), accmean(&gmean_a[1][gen]), accmean(&gvar_a[1][gen]), accmean(&fp[1][gen]), accmean(&Fp[1][gen]), accmean(&LEQ_NLf[1][gen]), accmean(&LEQ_Lf[1][gen]), accmean(&LEQ_T[1][gen]), accmean(&SK2[1][gen]), accmean(&Hw[1][gen]), variance(&Hw[1][gen]));
	}

	fprintf(fgenL1, "\n");

	/* ******* NINDTP_R_NINDAL (L2) ******* */

	fprintf(fgenL2,"\ngen     Wf     repvarWf    p_a     var(p_a)    fp       Fp      LEQ_NL     LEQ_L     LEQ_T      SK2       Hw      repvarHw\n");

	for (gen=0; gen<tLINES; gen++)	/* TP */
	{
		fprintf(fgenL2, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", gen, accmean(&gmean_sf[0][gen]), variance(&gmean_sf[0][gen]), accmean(&gmean_a[0][gen]), accmean(&gvar_a[0][gen]), accmean(&fp[0][gen]), accmean(&Fp[0][gen]), accmean(&LEQ_NLf[0][gen]), accmean(&LEQ_Lf[0][gen]), accmean(&LEQ_T[0][gen]), accmean(&SK2[0][gen]), accmean(&Hw[0][gen]), variance(&Hw[0][gen]));
	}

	for (gen=tLINES; gen<=generations; gen++)  /* Formation of lines */
	{
		fprintf(fgenL2, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", gen, accmean(&gmean_sf[2][gen]), variance(&gmean_sf[2][gen]), accmean(&gmean_a[2][gen]), accmean(&gvar_a[2][gen]), accmean(&fp[2][gen]), accmean(&Fp[2][gen]), accmean(&LEQ_NLf[2][gen]), accmean(&LEQ_Lf[2][gen]), accmean(&LEQ_T[2][gen]), accmean(&SK2[2][gen]), accmean(&Hw[2][gen]), variance(&Hw[2][gen]));
	}

	fprintf(fgenL2, "\n");

	/* ******* NINDTP_NR_NINDAL (L3) ******* */

	fprintf(fgenL3,"\ngen     Wf     repvarWf    p_a     var(p_a)    fp       Fp      LEQ_NL     LEQ_L     LEQ_T      SK2       Hw      repvarHw\n");

	for (gen=0; gen<tLINES; gen++)	/* TP */
	{
		fprintf(fgenL3, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", gen, accmean(&gmean_sf[0][gen]), variance(&gmean_sf[0][gen]), accmean(&gmean_a[0][gen]), accmean(&gvar_a[0][gen]), accmean(&fp[0][gen]), accmean(&Fp[0][gen]), accmean(&LEQ_NLf[0][gen]), accmean(&LEQ_Lf[0][gen]), accmean(&LEQ_T[0][gen]), accmean(&SK2[0][gen]), accmean(&Hw[0][gen]), variance(&Hw[0][gen]));
	}

	for (gen=tLINES; gen<=generations; gen++)  /* Formation of lines */
	{
		fprintf(fgenL3, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", gen, accmean(&gmean_sf[3][gen]), variance(&gmean_sf[3][gen]), accmean(&gmean_a[3][gen]), accmean(&gvar_a[3][gen]), accmean(&fp[3][gen]), accmean(&Fp[3][gen]), accmean(&LEQ_NLf[3][gen]), accmean(&LEQ_Lf[3][gen]), accmean(&LEQ_T[3][gen]), accmean(&SK2[3][gen]), accmean(&Hw[3][gen]), variance(&Hw[3][gen]));
	}

	fprintf(fgenL3, "\n");


	/* ******* NINDTP_RR_NINDTP (subL1) ******* */

	fprintf(fgensubL1,"\ngen     Wf     repvarWf    p_a     var(p_a)    fp       Fp      LEQ_NL     LEQ_L     LEQ_T      SK2       Hw      repvarHw\n");

	for (gen=0; gen<tLINES; gen++)	/* TP */
	{
		fprintf(fgensubL1, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", gen, accmean(&gmean_sf[0][gen]), variance(&gmean_sf[0][gen]), accmean(&gmean_a[0][gen]), accmean(&gvar_a[0][gen]), accmean(&fp[0][gen]), accmean(&Fp[0][gen]), accmean(&LEQ_NLf[0][gen]), accmean(&LEQ_Lf[0][gen]), accmean(&LEQ_T[0][gen]), accmean(&SK2[0][gen]), accmean(&Hw[0][gen]), variance(&Hw[0][gen]));
	}

	for (gen=tLINES; gen<tsubL; gen++)  /* Formation of lines - NR */
	{
		fprintf(fgensubL1, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", gen, accmean(&gmean_sf[4][gen]), variance(&gmean_sf[4][gen]), accmean(&gmean_a[4][gen]), accmean(&gvar_a[4][gen]), accmean(&fp[4][gen]), accmean(&Fp[4][gen]), accmean(&LEQ_NLf[4][gen]), accmean(&LEQ_Lf[4][gen]), accmean(&LEQ_T[4][gen]), accmean(&SK2[4][gen]), accmean(&Hw[4][gen]), variance(&Hw[4][gen]));
	}

	for (gen=tsubL; gen<=generations; gen++)	/* Formation of sublines - Retarded Rescue */
	{
		fprintf(fgensubL1, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", gen, accmean(&gmean_sf[5][gen]), variance(&gmean_sf[5][gen]), accmean(&gmean_a[5][gen]), accmean(&gvar_a[5][gen]), accmean(&fp[5][gen]), accmean(&Fp[5][gen]), accmean(&LEQ_NLf[5][gen]), accmean(&LEQ_Lf[5][gen]), accmean(&LEQ_T[5][gen]), accmean(&SK2[5][gen]), accmean(&Hw[5][gen]), variance(&Hw[5][gen]));
	}

	fprintf(fgensubL1, "\n");

	/* ******* NINDTP_RR_NINDAL (subL2) ******* */

	fprintf(fgensubL2,"\ngen     Wf     repvarWf    p_a     var(p_a)    fp       Fp      LEQ_NL     LEQ_L     LEQ_T      SK2       Hw      repvarHw\n");

	for (gen=0; gen<tLINES; gen++)	/* TP */
	{
		fprintf(fgensubL2, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", gen, accmean(&gmean_sf[0][gen]), variance(&gmean_sf[0][gen]), accmean(&gmean_a[0][gen]), accmean(&gvar_a[0][gen]), accmean(&fp[0][gen]), accmean(&Fp[0][gen]), accmean(&LEQ_NLf[0][gen]), accmean(&LEQ_Lf[0][gen]), accmean(&LEQ_T[0][gen]), accmean(&SK2[0][gen]), accmean(&Hw[0][gen]), variance(&Hw[0][gen]));
	}

	for (gen=tLINES; gen<tsubL; gen++)  /* Formation of lines - NR */
	{
		fprintf(fgensubL2, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", gen, accmean(&gmean_sf[4][gen]), variance(&gmean_sf[4][gen]), accmean(&gmean_a[4][gen]), accmean(&gvar_a[4][gen]), accmean(&fp[4][gen]), accmean(&Fp[4][gen]), accmean(&LEQ_NLf[4][gen]), accmean(&LEQ_Lf[4][gen]), accmean(&LEQ_T[4][gen]), accmean(&SK2[4][gen]), accmean(&Hw[4][gen]), variance(&Hw[4][gen]));
	}

	for (gen=tsubL; gen<=generations; gen++) /* Fromation of sublines - Retarded Rescue */
	{
		fprintf(fgensubL2, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", gen, accmean(&gmean_sf[6][gen]), variance(&gmean_sf[6][gen]), accmean(&gmean_a[6][gen]), accmean(&gvar_a[6][gen]), accmean(&fp[6][gen]), accmean(&Fp[6][gen]), accmean(&LEQ_NLf[6][gen]), accmean(&LEQ_Lf[6][gen]), accmean(&LEQ_T[6][gen]), accmean(&SK2[6][gen]), accmean(&Hw[6][gen]), variance(&Hw[6][gen]));
	}

	fprintf(fgensubL2, "\n");

	/* ******* NINDTP_NR_NINDTP (subL3) ******* */

	fprintf(fgensubL3,"\ngen     Wf     repvarWf    p_a     var(p_a)    fp       Fp      LEQ_NL     LEQ_L     LEQ_T      SK2       Hw      repvarHw\n");

	for (gen=0; gen<tLINES; gen++)	/* TP */
	{
		fprintf(fgensubL3, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", gen, accmean(&gmean_sf[0][gen]), variance(&gmean_sf[0][gen]), accmean(&gmean_a[0][gen]), accmean(&gvar_a[0][gen]), accmean(&fp[0][gen]), accmean(&Fp[0][gen]), accmean(&LEQ_NLf[0][gen]), accmean(&LEQ_Lf[0][gen]), accmean(&LEQ_T[0][gen]), accmean(&SK2[0][gen]), accmean(&Hw[0][gen]), variance(&Hw[0][gen]));
	}

	for (gen=tLINES; gen<tsubL; gen++)  /* Formation of lines - NR */
	{
		fprintf(fgensubL3, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", gen, accmean(&gmean_sf[4][gen]), variance(&gmean_sf[4][gen]), accmean(&gmean_a[4][gen]), accmean(&gvar_a[4][gen]), accmean(&fp[4][gen]), accmean(&Fp[4][gen]), accmean(&LEQ_NLf[4][gen]), accmean(&LEQ_Lf[4][gen]), accmean(&LEQ_T[4][gen]), accmean(&SK2[4][gen]), accmean(&Hw[4][gen]), variance(&Hw[4][gen]));
	}

	for (gen=tsubL; gen<=generations; gen++) /* Formation of sublines - NR */
	{
		fprintf(fgensubL3, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", gen, accmean(&gmean_sf[7][gen]), variance(&gmean_sf[7][gen]), accmean(&gmean_a[7][gen]), accmean(&gvar_a[7][gen]), accmean(&fp[7][gen]), accmean(&Fp[7][gen]), accmean(&LEQ_NLf[7][gen]), accmean(&LEQ_Lf[7][gen]), accmean(&LEQ_T[7][gen]), accmean(&SK2[7][gen]), accmean(&Hw[7][gen]), variance(&Hw[7][gen]));
	}

	fprintf(fgensubL3, "\n");

	/* ******* NINDTP_NR_NINDAL (subL4) ******* */

	fprintf(fgensubL4,"\ngen     Wf     repvarWf    p_a     var(p_a)    fp       Fp      LEQ_NL     LEQ_L     LEQ_T      SK2       Hw      repvarHw\n");

	for (gen=0; gen<tLINES; gen++)	/* TP */
	{
		fprintf(fgensubL4, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", gen, accmean(&gmean_sf[0][gen]), variance(&gmean_sf[0][gen]), accmean(&gmean_a[0][gen]), accmean(&gvar_a[0][gen]), accmean(&fp[0][gen]), accmean(&Fp[0][gen]), accmean(&LEQ_NLf[0][gen]), accmean(&LEQ_Lf[0][gen]), accmean(&LEQ_T[0][gen]), accmean(&SK2[0][gen]), accmean(&Hw[0][gen]), variance(&Hw[0][gen]));
	}

	for (gen=tLINES; gen<tsubL; gen++)  /* Formation of lines - NR */
	{
		fprintf(fgensubL4, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", gen, accmean(&gmean_sf[4][gen]), variance(&gmean_sf[4][gen]), accmean(&gmean_a[4][gen]), accmean(&gvar_a[4][gen]), accmean(&fp[4][gen]), accmean(&Fp[4][gen]), accmean(&LEQ_NLf[4][gen]), accmean(&LEQ_Lf[4][gen]), accmean(&LEQ_T[4][gen]), accmean(&SK2[4][gen]), accmean(&Hw[4][gen]), variance(&Hw[4][gen]));
	}

	for (gen=tsubL; gen<=generations; gen++) /* Formation of sublines - NR */
	{
		fprintf(fgensubL4, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", gen, accmean(&gmean_sf[8][gen]), variance(&gmean_sf[8][gen]), accmean(&gmean_a[8][gen]), accmean(&gvar_a[8][gen]), accmean(&fp[8][gen]), accmean(&Fp[8][gen]), accmean(&LEQ_NLf[8][gen]), accmean(&LEQ_Lf[8][gen]), accmean(&LEQ_T[8][gen]), accmean(&SK2[8][gen]), accmean(&Hw[8][gen]), variance(&Hw[8][gen]));
	}

	fprintf(fgensubL4, "\n");
}

/* ********************************************************************* */

replicates_out()
{
	fprintf(frepext,"Number of extinct replicates per generation\n");
	fprintf(frepext,"\ngen  %dR%d(L1)  %dR%d(L2)  %dNR%d(L3_g%d)  %dRR%d(subL1)  %dRR%d(subL2)  %dNR%d(subL3)  %dNR%d(subL4_g%d)\n",NINDTP,NINDTP,NINDTP,NINDAL,NINDTP,NINDAL,tLINES,NINDTP,NINDTP,NINDTP,NINDAL,NINDTP,NINDTP,NINDTP,NINDAL,tsubL);

	for (gen=0; gen<tLINES; gen++)	/* TP */
	{
		fprintf(frepext, "%d      %d          %d            %d              %d              %d              %d              %d\n",gen,rep_extinct[0][gen],rep_extinct[0][gen],rep_extinct[0][gen],rep_extinct[0][gen],rep_extinct[0][gen],rep_extinct[0][gen],rep_extinct[0][gen]);
	}
	for (gen=tLINES; gen<tsubL; gen++)
	{
		fprintf(frepext, "%d      %d          %d            %d              %d              %d              %d              %d\n",gen,rep_extinct[1][gen],rep_extinct[2][gen],rep_extinct[3][gen],rep_extinct[4][gen],rep_extinct[4][gen],rep_extinct[4][gen],rep_extinct[4][gen]);
	}
	for (gen=tsubL; gen<=generations; gen++)
	{
		fprintf(frepext, "%d      %d          %d            %d              %d              %d              %d              %d\n",gen,rep_extinct[1][gen],rep_extinct[2][gen],rep_extinct[3][gen],rep_extinct[5][gen],rep_extinct[6][gen],rep_extinct[7][gen],rep_extinct[8][gen]);
	}

	fprintf(frepext,"\nDistribution w (gen before extinction) - total number of extinct replicates\n");

	fprintf(frepext,"\ndistribution	%dR%d(L1)  %dR%d(L2)  %dNR%d(L3_g%d)  %dRR%d(subL1)  %dRR%d(subL2)  %dNR%d(subL3)  %dNR%d(subL4_g%d)\n",NINDTP,NINDTP,NINDTP,NINDAL,NINDTP,NINDAL,tLINES,NINDTP,NINDTP,NINDTP,NINDAL,NINDTP,NINDTP,NINDTP,NINDAL,tsubL);
	fprintf (frepext, "w=0.0	%f    %f    %f    %f    %f    %f    %f\n", accsum(&wext_00[1]), accsum(&wext_00[2]), accsum(&wext_00[3]), accsum(&wext_00[5]), accsum(&wext_00[6]), accsum(&wext_00[7]), accsum(&wext_00[8]));
	fprintf (frepext, "w=0.0-0.1	%f    %f    %f    %f    %f    %f    %f\n", accsum(&wext_00_01[1]), accsum(&wext_00_01[2]), accsum(&wext_00_01[3]), accsum(&wext_00_01[5]), accsum(&wext_00_01[6]), accsum(&wext_00_01[7]), accsum(&wext_00_01[8]));
	fprintf (frepext, "w=0.1-0.2	%f    %f    %f    %f    %f    %f    %f\n", accsum(&wext_01_02[1]), accsum(&wext_01_02[2]), accsum(&wext_01_02[3]), accsum(&wext_01_02[5]), accsum(&wext_01_02[6]), accsum(&wext_01_02[7]), accsum(&wext_01_02[8]));
	fprintf (frepext, "w=0.2-0.3	%f    %f    %f    %f    %f    %f    %f\n", accsum(&wext_02_03[1]), accsum(&wext_02_03[2]), accsum(&wext_02_03[3]), accsum(&wext_02_03[5]), accsum(&wext_02_03[6]), accsum(&wext_02_03[7]), accsum(&wext_02_03[8]));
	fprintf (frepext, "w=0.3-0.4	%f    %f    %f    %f    %f    %f    %f\n", accsum(&wext_03_04[1]), accsum(&wext_03_04[2]), accsum(&wext_03_04[3]), accsum(&wext_03_04[5]), accsum(&wext_03_04[6]), accsum(&wext_03_04[7]), accsum(&wext_03_04[8]));
	fprintf (frepext, "w=0.4-0.5	%f    %f    %f    %f    %f    %f    %f\n", accsum(&wext_04_05[1]), accsum(&wext_04_05[2]), accsum(&wext_04_05[3]), accsum(&wext_04_05[5]), accsum(&wext_04_05[6]), accsum(&wext_04_05[7]), accsum(&wext_04_05[8]));
	fprintf (frepext, "w=0.5-0.6	%f    %f    %f    %f    %f    %f    %f\n", accsum(&wext_05_06[1]), accsum(&wext_05_06[2]), accsum(&wext_05_06[3]), accsum(&wext_05_06[5]), accsum(&wext_05_06[6]), accsum(&wext_05_06[7]), accsum(&wext_05_06[8]));
	fprintf (frepext, "w=0.6-0.7	%f    %f    %f    %f    %f    %f    %f\n", accsum(&wext_06_07[1]), accsum(&wext_06_07[2]), accsum(&wext_06_07[3]), accsum(&wext_06_07[5]), accsum(&wext_06_07[6]), accsum(&wext_06_07[7]), accsum(&wext_06_07[8]));
	fprintf (frepext, "w=0.7-0.8	%f    %f    %f    %f    %f    %f    %f\n", accsum(&wext_07_08[1]), accsum(&wext_07_08[2]), accsum(&wext_07_08[3]), accsum(&wext_07_08[5]), accsum(&wext_07_08[6]), accsum(&wext_07_08[7]), accsum(&wext_07_08[8]));
	fprintf (frepext, "w=0.8-0.9	%f    %f    %f    %f    %f    %f    %f\n", accsum(&wext_08_09[1]), accsum(&wext_08_09[2]), accsum(&wext_08_09[3]), accsum(&wext_08_09[5]), accsum(&wext_08_09[6]), accsum(&wext_08_09[7]), accsum(&wext_08_09[8]));
	fprintf (frepext, "w=0.9-1.0	%f    %f    %f    %f    %f    %f    %f\n", accsum(&wext_09_10[1]), accsum(&wext_09_10[2]), accsum(&wext_09_10[3]), accsum(&wext_09_10[5]), accsum(&wext_09_10[6]), accsum(&wext_09_10[7]), accsum(&wext_09_10[8]));
	fprintf (frepext, "w=1.0	%f    %f    %f    %f    %f    %f    %f\n", accsum(&wext_10[1]), accsum(&wext_10[2]), accsum(&wext_10[3]), accsum(&wext_10[5]), accsum(&wext_10[6]), accsum(&wext_10[7]), accsum(&wext_10[8]));

	fprintf(frepext,"\nDistribution w (gen before extinction) - average number of lethals per individual\n");

	fprintf (frepext,"\ndistribution	%dR%d(L1)  %dR%d(L2)  %dNR%d(L3_g%d)  %dRR%d(subL1)  %dRR%d(subL2)  %dNR%d(subL3)  %dNR%d(subL4_g%d)\n",NINDTP,NINDTP,NINDTP,NINDAL,NINDTP,NINDAL,tLINES,NINDTP,NINDTP,NINDTP,NINDAL,NINDTP,NINDTP,NINDTP,NINDAL,tsubL);
	fprintf (frepext, "w=0.0	%f    %f    %f    %f    %f    %f    %f\n", accmean(&STR_00[1]), accmean(&STR_00[2]), accmean(&STR_00[3]), accmean(&STR_00[5]), accmean(&STR_00[6]), accmean(&STR_00[7]), accmean(&STR_00[8]));
	fprintf (frepext, "w=0.0-0.1	%f    %f    %f    %f    %f    %f    %f\n", accmean(&STR_00_01[1]), accmean(&STR_00_01[2]), accmean(&STR_00_01[3]), accmean(&STR_00_01[5]), accmean(&STR_00_01[6]), accmean(&STR_00_01[7]), accmean(&STR_00_01[8]));
	fprintf (frepext, "w=0.1-0.2	%f    %f    %f    %f    %f    %f    %f\n", accmean(&STR_01_02[1]), accmean(&STR_01_02[2]), accmean(&STR_01_02[3]), accmean(&STR_01_02[5]), accmean(&STR_01_02[6]), accmean(&STR_01_02[7]), accmean(&STR_01_02[8]));
	fprintf (frepext, "w=0.2-0.3	%f    %f    %f    %f    %f    %f    %f\n", accmean(&STR_02_03[1]), accmean(&STR_02_03[2]), accmean(&STR_02_03[3]), accmean(&STR_02_03[5]), accmean(&STR_02_03[6]), accmean(&STR_02_03[7]), accmean(&STR_02_03[8]));
	fprintf (frepext, "w=0.3-0.4	%f    %f    %f    %f    %f    %f    %f\n", accmean(&STR_03_04[1]), accmean(&STR_03_04[2]), accmean(&STR_03_04[3]), accmean(&STR_03_04[5]), accmean(&STR_03_04[6]), accmean(&STR_03_04[7]), accmean(&STR_03_04[8]));
	fprintf (frepext, "w=0.4-0.5	%f    %f    %f    %f    %f    %f    %f\n", accmean(&STR_04_05[1]), accmean(&STR_04_05[2]), accmean(&STR_04_05[3]), accmean(&STR_04_05[5]), accmean(&STR_04_05[6]), accmean(&STR_04_05[7]), accmean(&STR_04_05[8]));
	fprintf (frepext, "w=0.5-0.6	%f    %f    %f    %f    %f    %f    %f\n", accmean(&STR_05_06[1]), accmean(&STR_05_06[2]), accmean(&STR_05_06[3]), accmean(&STR_05_06[5]), accmean(&STR_05_06[6]), accmean(&STR_05_06[7]), accmean(&STR_05_06[8]));
	fprintf (frepext, "w=0.6-0.7	%f    %f    %f    %f    %f    %f    %f\n", accmean(&STR_06_07[1]), accmean(&STR_06_07[2]), accmean(&STR_06_07[3]), accmean(&STR_06_07[5]), accmean(&STR_06_07[6]), accmean(&STR_06_07[7]), accmean(&STR_06_07[8]));
	fprintf (frepext, "w=0.7-0.8	%f    %f    %f    %f    %f    %f    %f\n", accmean(&STR_07_08[1]), accmean(&STR_07_08[2]), accmean(&STR_07_08[3]), accmean(&STR_07_08[5]), accmean(&STR_07_08[6]), accmean(&STR_07_08[7]), accmean(&STR_07_08[8]));
	fprintf (frepext, "w=0.8-0.9	%f    %f    %f    %f    %f    %f    %f\n", accmean(&STR_08_09[1]), accmean(&STR_08_09[2]), accmean(&STR_08_09[3]), accmean(&STR_08_09[5]), accmean(&STR_08_09[6]), accmean(&STR_08_09[7]), accmean(&STR_08_09[8]));
	fprintf (frepext, "w=0.9-1.0	%f    %f    %f    %f    %f    %f    %f\n", accmean(&STR_09_10[1]), accmean(&STR_09_10[2]), accmean(&STR_09_10[3]), accmean(&STR_09_10[5]), accmean(&STR_09_10[6]), accmean(&STR_09_10[7]), accmean(&STR_09_10[8]));
	fprintf (frepext, "w=1.0	%f    %f    %f    %f    %f    %f    %f\n", accmean(&STR_10[1]), accmean(&STR_10[2]), accmean(&STR_10[3]), accmean(&STR_10[5]), accmean(&STR_10[6]), accmean(&STR_10[7]), accmean(&STR_10[8]));

	fprintf(frepext,"\nDistribution w (gen before extinction) - dq\n");

	fprintf (frepext,"\ndistribution	%dR%d(L1)  %dR%d(L2)  %dNR%d(L3_g%d)  %dRR%d(subL1)  %dRR%d(subL2)  %dNR%d(subL3)  %dNR%d(subL4_g%d)\n",NINDTP,NINDTP,NINDTP,NINDAL,NINDTP,NINDAL,tLINES,NINDTP,NINDTP,NINDTP,NINDAL,NINDTP,NINDTP,NINDTP,NINDAL,tsubL);
	fprintf (frepext, "w=0.0	%f    %f    %f    %f    %f    %f    %f\n", accmean(&dq_00[1]), accmean(&dq_00[2]), accmean(&dq_00[3]), accmean(&dq_00[5]), accmean(&dq_00[6]), accmean(&dq_00[7]), accmean(&dq_00[8]));
	fprintf (frepext, "w=0.0-0.1	%f    %f    %f    %f    %f    %f    %f\n", accmean(&dq_00_01[1]), accmean(&dq_00_01[2]), accmean(&dq_00_01[3]), accmean(&dq_00_01[5]), accmean(&dq_00_01[6]), accmean(&dq_00_01[7]), accmean(&dq_00_01[8]));
	fprintf (frepext, "w=0.1-0.2	%f    %f    %f    %f    %f    %f    %f\n", accmean(&dq_01_02[1]), accmean(&dq_01_02[2]), accmean(&dq_01_02[3]), accmean(&dq_01_02[5]), accmean(&dq_01_02[6]), accmean(&dq_01_02[7]), accmean(&dq_01_02[8]));
	fprintf (frepext, "w=0.2-0.3	%f    %f    %f    %f    %f    %f    %f\n", accmean(&dq_02_03[1]), accmean(&dq_02_03[2]), accmean(&dq_02_03[3]), accmean(&dq_02_03[5]), accmean(&dq_02_03[6]), accmean(&dq_02_03[7]), accmean(&dq_02_03[8]));
	fprintf (frepext, "w=0.3-0.4	%f    %f    %f    %f    %f    %f    %f\n", accmean(&dq_03_04[1]), accmean(&dq_03_04[2]), accmean(&dq_03_04[3]), accmean(&dq_03_04[5]), accmean(&dq_03_04[6]), accmean(&dq_03_04[7]), accmean(&dq_03_04[8]));
	fprintf (frepext, "w=0.4-0.5	%f    %f    %f    %f    %f    %f    %f\n", accmean(&dq_04_05[1]), accmean(&dq_04_05[2]), accmean(&dq_04_05[3]), accmean(&dq_04_05[5]), accmean(&dq_04_05[6]), accmean(&dq_04_05[7]), accmean(&dq_04_05[8]));
	fprintf (frepext, "w=0.5-0.6	%f    %f    %f    %f    %f    %f    %f\n", accmean(&dq_05_06[1]), accmean(&dq_05_06[2]), accmean(&dq_05_06[3]), accmean(&dq_05_06[5]), accmean(&dq_05_06[6]), accmean(&dq_05_06[7]), accmean(&dq_05_06[8]));
	fprintf (frepext, "w=0.6-0.7	%f    %f    %f    %f    %f    %f    %f\n", accmean(&dq_06_07[1]), accmean(&dq_06_07[2]), accmean(&dq_06_07[3]), accmean(&dq_06_07[5]), accmean(&dq_06_07[6]), accmean(&dq_06_07[7]), accmean(&dq_06_07[8]));
	fprintf (frepext, "w=0.7-0.8	%f    %f    %f    %f    %f    %f    %f\n", accmean(&dq_07_08[1]), accmean(&dq_07_08[2]), accmean(&dq_07_08[3]), accmean(&dq_07_08[5]), accmean(&dq_07_08[6]), accmean(&dq_07_08[7]), accmean(&dq_07_08[8]));
	fprintf (frepext, "w=0.8-0.9	%f    %f    %f    %f    %f    %f    %f\n", accmean(&dq_08_09[1]), accmean(&dq_08_09[2]), accmean(&dq_08_09[3]), accmean(&dq_08_09[5]), accmean(&dq_08_09[6]), accmean(&dq_08_09[7]), accmean(&dq_08_09[8]));
	fprintf (frepext, "w=0.9-1.0	%f    %f    %f    %f    %f    %f    %f\n", accmean(&dq_09_10[1]), accmean(&dq_09_10[2]), accmean(&dq_09_10[3]), accmean(&dq_09_10[5]), accmean(&dq_09_10[6]), accmean(&dq_09_10[7]), accmean(&dq_09_10[8]));
	fprintf (frepext, "w=1.0	%f    %f    %f    %f    %f    %f    %f\n", accmean(&dq_10[1]), accmean(&dq_10[2]), accmean(&dq_10[3]), accmean(&dq_10[5]), accmean(&dq_10[6]), accmean(&dq_10[7]), accmean(&dq_10[8]));
}

/* ********************************************************************* */

sum_outline()
{

	fprintf(fsum,"\n\n                  %dR%d (N=%d +mig)\n                ----------------------------------------> %d generations\n               |\n               |  %dR%d (N=%d +mig)\n  Threat.Pop   |----------------------------------------> %d generations\n -------------*|\n    N=%d       |  %dNR%d (N=%d)\n               |----------------------------------------> %d generations (CONTROL)\n               |\n               |                     %dRR%d (N=%d +mig)\n               |                  ----------------------> %d generations\n               |                 |\n               |                 |   %dRR%d (N=%d +mig)\n               |                 |----------------------> %d generations\n               |  %dNR%d (N=%d)  |\n                ---------------**|   %dNR%d (N=%d)\n                                 |----------------------> %d generations (CONTROL)\n                                 |\n                                 |   %dNR%d (N=%d)\n                                  ----------------------> %d generations (CONTROL)\n", NINDTP, NINDTP, NINDTP, generations, NINDTP, NINDAL, NINDAL, generations, NINDTP, NINDTP, NINDAL, NINDAL, generations, NINDTP, NINDTP, NINDTP, generations, NINDTP, NINDAL, NINDAL, generations, NINDTP, NINDTP, NINDTP, NINDTP, NINDTP, NINDTP, generations, NINDTP, NINDAL, NINDAL, generations); 

	fprintf(fsum,"\n\n              * Formation of lines: generation %d\n             ** Formation of sublines: generation %d\n              R Genetic Rescue\n             RR Retarded Rescue\n             NR Non-Rescue\n", tLINES, tsubL);

	fprintf(fsum,"\n\nSeed: %d\nNIND Threatened Population 'NINDTP' (max 1000): %d\nNIND Alternative Lines 'NINDAL' (max 1000): %d\nNumber of Migrants for lines with NINDTP individuals (max 1000): %d\nNumber of Migrants for lines with NINDAL individuals (max 1000): %d\nGender of Migrants (males 0, males&females 1): %d\nLength of genome in Morgans(99:FreeRecom): %4.2f\nNCRO (min 1, max 2000): %d\nNLOCI (first is neutral)(min 2, max 30): %d\nLambda_s: %f\nLambda_L (s=1,h=0.02): %f\nBeta_s: %f\nAverage |s|: %f\ndom (constant 0, variable 1): %d\nAverage h: %f\nStabilizing selection (Vs): %f\nVE: %f\nIf Vs!=0.0, optimal fitness: %f\nNumber of generations: %d\nFormation of R and NR lines (t): %d\nFormation of sublines (t): %d\nFor R lines, initial NR period after formation of line : %d\nNumber of migrations (99: periodic): %d\nGeneration intervals of migration (Lines 'NINDTP'): %d\nGeneration intervals of migration (Lines 'NINDAL'): %d\nGenerations since migration (distribution -w- among replicates): %d\nMinimum fitness value for extinction (99: without extinction): %f\nNumber of replicates: %d\n\n", seed, NINDTP, NINDAL, NMIGTP, NMIGAL, mGEND, L, NCRO, NLOCI, Lambda_s, Lambda_L, beta_s, ave_s, dom, ave_hs, Vs, VE, OPT, generations, tLINES, tsubL, preNR, numMIG, mINT_TP, mINT_AL, dist, wEXT, replicates);

}

/* ********************************************************************* */

distribution_qsh_out()
{
	fprintf(fdis,"\nDISTRIBUTION qsh\n");

	/* LINES */
	for (lin=1; lin<=3; lin++)
	{
		if (lin==1)		fprintf(fdis,"\n *********** %dR%d (line 1)***********\n", NINDTP, NINDTP);
		else if (lin==2)	fprintf(fdis,"\n *********** %dR%d (line 2) ***********\n", NINDTP, NINDAL);
		else			fprintf(fdis,"\n *********** %dNR%d (line 3) ***********\n", NINDTP, NINDAL);

		fprintf (fdis, "\nGen		   "); 
		for (i=0; i<=(generations/tgen); i++)	fprintf (fdis, "%d        ", i*tgen);
		fprintf (fdis, "\n");

		// distribution of q values in the population

		fprintf (fdis, "q=0.0          ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_00[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_00[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.0-0.1      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_00_01[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_00_01[lin][i])/replicates);		
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.1-0.2      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_01_02[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_01_02[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.2-0.3      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_02_03[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_02_03[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.3-0.4      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_03_04[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_03_04[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.4-0.5      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_04_05[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_04_05[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.5-0.6      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_05_06[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_05_06[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.6-0.7      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_06_07[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_06_07[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.7-0.8      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_07_08[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_07_08[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.8-0.9      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_08_09[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_08_09[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.9-1.0      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_09_10[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_09_10[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=1.0          ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_10[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_10[lin][i])/replicates);
		fprintf (fdis, "\n");

		// distribution of s values in the population

		fprintf (fdis, "LOST\n");

		fprintf (fdis, "s=0.0          ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=0.0-10-6     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_106_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_106_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-10-6-10-4   ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_106_104_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_106_104_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-10-4-0.01   ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_104_102_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_104_102_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.01-0.1    ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_102_101_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_102_101_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.1-0.2     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_01_02_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_01_02_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.2-0.4     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_02_04_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_02_04_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.4-0.6     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_04_06_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_04_06_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.6-1.0     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_06_10_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_06_10_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-1.0         ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_10_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_10_lost[lin][i])/replicates);
		fprintf (fdis, "\n");

		fprintf (fdis, "SEGREGATING\n");

		fprintf (fdis, "s=0.0          ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=0.0-10-6     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_106[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_106[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-10-6-10-4   ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_106_104[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_106_104[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-10-4-0.01  ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_104_102[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_104_102[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.01-0.1   ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_102_101[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_102_101[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.1-0.2     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_01_02[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_01_02[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.2-0.4     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_02_04[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_02_04[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.4-0.6     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_04_06[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_04_06[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.6-1.0     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_06_10[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_06_10[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-1.0         ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_10[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_10[lin][i])/replicates);
		fprintf (fdis, "\n");

		fprintf (fdis, "FIXED\n");

		fprintf (fdis, "s=0.0          ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=0.0-10-6     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_106_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_106_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-10-6-10-4   ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_106_104_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_106_104_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-10-4-0.01  ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_104_102_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_104_102_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.01-0.1   ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_102_101_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_102_101_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.1-0.2     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_01_02_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_01_02_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.2-0.4     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_02_04_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_02_04_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.4-0.6     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_04_06_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_04_06_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.6-1.0     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_06_10_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_06_10_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-1.0         ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_10_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_10_fix[lin][i])/replicates);
		fprintf (fdis, "\n\n");

		// distribution of hs values in the population

		fprintf (fdis, "LOST\n");

		fprintf (fdis, "hs=0.0-0.2     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_00_02_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_00_02_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.2-0.4     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_02_04_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_02_04_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.4-0.6     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_04_06_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_04_06_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.6-0.8     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_06_08_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_06_08_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.8-1.0     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_08_10_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_08_10_lost[lin][i])/replicates);
		fprintf (fdis, "\n");

		fprintf (fdis, "SEGREGATING\n");

		fprintf (fdis, "hs=0.0-0.2     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_00_02[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_00_02[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.2-0.4     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_02_04[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_02_04[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.4-0.6     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_04_06[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_04_06[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.6-0.8     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_06_08[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_06_08[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.8-1.0     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_08_10[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_08_10[lin][i])/replicates);
		fprintf (fdis, "\n");

		fprintf (fdis, "FIXED\n");

		fprintf (fdis, "hs=0.0-0.2     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_00_02_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_00_02_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.2-0.4     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_02_04_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_02_04_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.4-0.6     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_04_06_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_04_06_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.6-0.8     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_06_08_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_06_08_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.8-1.0     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_08_10_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_08_10_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
	}

	/* SUBLINES */
	for (lin=5; lin<=8; lin++)
	{
		if (lin==5)		fprintf(fdis,"\n *********** %dRR%d (line 4 - subline 1) ***********\n", NINDTP, NINDTP);	
		else if (lin==6)	fprintf(fdis,"\n *********** %dRR%d (line 4 - subline 2) ***********\n", NINDTP, NINDAL);
		else if (lin==7)	fprintf(fdis,"\n *********** %dNR%d (line 4 - subline 3) ***********\n", NINDTP, NINDTP);
		else			fprintf(fdis,"\n *********** %dNR%d (line 4 - subline 4) ***********\n", NINDTP, NINDAL);

		fprintf (fdis, "\nGen		   "); 
		for (i=0; i<=(generations/tgen); i++)	fprintf (fdis, "%d        ", i*tgen);
		fprintf (fdis, "\n");

		// distribution of q values in the population

		fprintf (fdis, "q=0.0          ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_00[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_00[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_00[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.0-0.1      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_00_01[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_00_01[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_00_01[lin][i])/replicates);		
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.1-0.2      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_01_02[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_01_02[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_01_02[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.2-0.3      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_02_03[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_02_03[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_02_03[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.3-0.4      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_03_04[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_03_04[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_03_04[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.4-0.5      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_04_05[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_04_05[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_04_05[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.5-0.6      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_05_06[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_05_06[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_05_06[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.6-0.7      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_06_07[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_06_07[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_06_07[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.7-0.8      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_07_08[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_07_08[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_07_08[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.8-0.9      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_08_09[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_08_09[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_08_09[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=0.9-1.0      ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_09_10[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_09_10[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_09_10[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "q=1.0          ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_10[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_10[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&q_ng_10[lin][i])/replicates);
		fprintf (fdis, "\n");

	// distribution of s values in the population

		fprintf (fdis, "LOST\n");

		fprintf (fdis, "s=0.0          ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_lost[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=0.0-10-6     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_106_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_106_lost[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_106_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-10-6-10-4   ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_106_104_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_106_104_lost[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_106_104_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-10-4-0.01   ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_104_102_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_104_102_lost[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_104_102_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.01-0.1    ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_102_101_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_102_101_lost[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_102_101_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.1-0.2     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_01_02_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_01_02_lost[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_01_02_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.2-0.4     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_02_04_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_02_04_lost[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_02_04_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.4-0.6     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_04_06_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_04_06_lost[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_04_06_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.6-1.0     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_06_10_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_06_10_lost[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_06_10_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-1.0         ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_10_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_10_lost[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_10_lost[lin][i])/replicates);
		fprintf (fdis, "\n");

		fprintf (fdis, "SEGREGATING\n");

		fprintf (fdis, "s=0.0          ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=0.0-10-6     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_106[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_106[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_106[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-10-6-10-4   ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_106_104[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_106_104[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_106_104[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-10-4-0.01  ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_104_102[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_104_102[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_104_102[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.01-0.1   ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_102_101[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_102_101[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_102_101[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.1-0.2     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_01_02[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_01_02[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_01_02[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.2-0.4     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_02_04[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_02_04[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_02_04[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.4-0.6     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_04_06[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_04_06[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_04_06[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.6-1.0     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_06_10[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_06_10[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_06_10[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-1.0         ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_10[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_10[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_10[lin][i])/replicates);
		fprintf (fdis, "\n");

		fprintf (fdis, "FIXED\n");

		fprintf (fdis, "s=0.0          ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_fix[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=0.0-10-6     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_106_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_106_fix[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_0_106_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-10-6-10-4   ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_106_104_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_106_104_fix[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_106_104_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-10-4-0.01  ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_104_102_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_104_102_fix[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_104_102_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.01-0.1   ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_102_101_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_102_101_fix[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_102_101_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.1-0.2     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_01_02_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_01_02_fix[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_01_02_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.2-0.4     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_02_04_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_02_04_fix[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_02_04_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.4-0.6     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_04_06_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_04_06_fix[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_04_06_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-0.6-1.0     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_06_10_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_06_10_fix[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_06_10_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "s=-1.0         ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_10_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_10_fix[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&s_ng_10_fix[lin][i])/replicates);
		fprintf (fdis, "\n\n");

		// distribution of hs values in the population

		fprintf (fdis, "LOST\n");

		fprintf (fdis, "hs=0.0-0.2     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_00_02_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_00_02_lost[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_00_02_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.2-0.4     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_02_04_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_02_04_lost[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_02_04_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.4-0.6     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_04_06_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_04_06_lost[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_04_06_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.6-0.8     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_06_08_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_06_08_lost[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_06_08_lost[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.8-1.0     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_08_10_lost[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_08_10_lost[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_08_10_lost[lin][i])/replicates);
		fprintf (fdis, "\n");

		fprintf (fdis, "SEGREGATING\n");

		fprintf (fdis, "hs=0.0-0.2     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_00_02[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_00_02[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_00_02[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.2-0.4     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_02_04[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_02_04[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_02_04[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.4-0.6     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_04_06[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_04_06[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_04_06[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.6-0.8     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_06_08[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_06_08[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_06_08[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.8-1.0     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_08_10[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_08_10[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_08_10[lin][i])/replicates);
		fprintf (fdis, "\n");

		fprintf (fdis, "FIXED\n");

		fprintf (fdis, "hs=0.0-0.2     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_00_02_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_00_02_fix[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_00_02_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.2-0.4     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_02_04_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_02_04_fix[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_02_04_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.4-0.6     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_04_06_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_04_06_fix[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_04_06_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.6-0.8     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_06_08_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_06_08_fix[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_06_08_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
		fprintf (fdis, "hs=0.8-1.0     ");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_08_10_fix[0][i])/replicates);
		for (i=(tLINES/tgen); i<(tsubL/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_08_10_fix[4][i])/replicates);
		for (i=(tsubL/tgen); i<=(generations/tgen); i++) fprintf (fdis, "%8.2f  ", accsum(&hs_ng_08_10_fix[lin][i])/replicates);
		fprintf (fdis, "\n");
	}
}

/* ********************************************************************* */

distribution_w_out()
{
	fprintf(fdisw, "\n*** DISTRIBUTION w ***\n");

	for (lin=1; lin<lines; lin++)
	{
		if (lin==1)	 { fprintf(fdisw,"\n *********** %dR%d (line 1)***********\n", NINDTP, NINDTP); init = tLINES; mINT = mINT_TP; }
		else if (lin==2) { fprintf(fdisw,"\n *********** %dR%d (line 2) ***********\n", NINDTP, NINDAL); init = tLINES; mINT = mINT_AL; }
		else if (lin==3) { fprintf(fdisw,"\n *********** %dNR%d (line 3) ***********\n", NINDTP, NINDAL); init = tLINES; mINT = mINT_AL; }
		else if (lin==4) { fprintf(fdisw,"\n *********** %dNR (line 4) ***********\n", NINDTP); init = tLINES; mINT = mINT_TP; }
		else if (lin==5) { fprintf(fdisw,"\n *********** %dRR%d (subL1) ***********\n", NINDTP, NINDTP); init = tsubL; mINT = mINT_TP; }
		else if (lin==6) { fprintf(fdisw,"\n *********** %dRR%d (subL2) ***********\n", NINDTP, NINDAL); init = tsubL; mINT = mINT_AL; }
		else if (lin==7) { fprintf(fdisw,"\n *********** %dNR%d (subL3) ***********\n", NINDTP, NINDTP); init = tsubL; mINT = mINT_TP; }
		else if (lin==8) { fprintf(fdisw,"\n *********** %dNR%d (subL4) ***********\n", NINDTP, NINDAL); init = tsubL; mINT = mINT_AL; }

		fprintf(fdisw, "\nAverage w:\n");
		for (i=1; i<=exp_migevents[lin]; i++)	fprintf (fdisw, "Gen %d	%f\n", init+preNR+mINT*(i-1)+dist, accmean(&summ_w[lin][i]));

		fprintf (fdisw,"\nDistribution	    "); 
		for (i=1; i<=exp_migevents[lin]; i++)	fprintf (fdisw, "gen%d        ", init+preNR+mINT*(i-1)+dist);
		fprintf (fdisw,"\n");

		fprintf (fdisw, "w=0.0          ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&w_00[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.0-0.1      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&w_00_01[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.1-0.2      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&w_01_02[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.2-0.3      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&w_02_03[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.3-0.4      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&w_03_04[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.4-0.5      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&w_04_05[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.5-0.6      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&w_05_06[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.6-0.7      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&w_06_07[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.7-0.8      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&w_07_08[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.8-0.9      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&w_08_09[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.9-1.0      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&w_09_10[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=1.0          ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&w_10[lin][i]));
		fprintf (fdisw, "\n");

		//EXTINCT REPLICATES
		fprintf (fdisw,"\nEXTINCT REPLICATES\n");

		fprintf(fdisw, "\nAverage w:\n");
		for (i=1; i<=exp_migevents[lin]; i++)	fprintf (fdisw, "Gen %d	%f\n", init+preNR+mINT*(i-1)+dist, accmean(&ext_summ_w[lin][i]));

		fprintf (fdisw,"\nDistribution	    "); 
		for (i=1; i<=exp_migevents[lin]; i++)	fprintf (fdisw, "gen%d        ", init+preNR+mINT*(i-1)+dist);
		fprintf (fdisw,"\n");

		fprintf (fdisw, "w=0.0          ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&ext_w_00[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.0-0.1      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&ext_w_00_01[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.1-0.2      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&ext_w_01_02[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.2-0.3      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&ext_w_02_03[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.3-0.4      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&ext_w_03_04[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.4-0.5      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&ext_w_04_05[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.5-0.6      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&ext_w_05_06[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.6-0.7      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&ext_w_06_07[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.7-0.8      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&ext_w_07_08[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.8-0.9      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&ext_w_08_09[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.9-1.0      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&ext_w_09_10[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=1.0          ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&ext_w_10[lin][i]));
		fprintf (fdisw, "\n");

		//SURVIVING REPLICATES
		fprintf (fdisw,"\nSURVIVING REPLICATES\n");

		fprintf(fdisw, "\nAverage w:\n");
		for (i=1; i<=exp_migevents[lin]; i++)	fprintf (fdisw, "Gen %d	%f\n", init+preNR+mINT*(i-1)+dist, accmean(&sur_summ_w[lin][i]));

		fprintf (fdisw,"\nDistribution	    "); 
		for (i=1; i<=exp_migevents[lin]; i++)	fprintf (fdisw, "gen%d        ", init+preNR+mINT*(i-1)+dist);
		fprintf (fdisw,"\n");

		fprintf (fdisw, "w=0.0          ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&sur_w_00[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.0-0.1      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&sur_w_00_01[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.1-0.2      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&sur_w_01_02[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.2-0.3      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&sur_w_02_03[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.3-0.4      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&sur_w_03_04[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.4-0.5      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&sur_w_04_05[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.5-0.6      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&sur_w_05_06[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.6-0.7      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&sur_w_06_07[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.7-0.8      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&sur_w_07_08[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.8-0.9      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&sur_w_08_09[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=0.9-1.0      ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&sur_w_09_10[lin][i]));
		fprintf (fdisw, "\n");
		fprintf (fdisw, "w=1.0          ");
		for (i=1; i<=exp_migevents[lin]; i++) fprintf (fdisw, "%8.2f  ", accsum(&sur_w_10[lin][i]));
		fprintf (fdisw, "\n");
	}
}

/* ********************************************************************* */

genresc_out() 
{
	int nlin, mig;

	/* LINES AND SUBLINES */
	for (nlin=1; nlin<=4; nlin++)
	{
		fprintf(fmate,"\n**********************************\n");
		if (nlin==1)	  {lin=1; fprintf(fmate,"\n********* %dR%d (line 1) *********\n", NINDTP, NINDTP); init=tLINES; mINT=mINT_TP;}
		else if (nlin==2) {lin=2; fprintf(fmate,"\n********* %dR%d (line 2) *********\n", NINDTP, NINDAL); init=tLINES; mINT=mINT_AL;}
		else if (nlin==3) {lin=5; fprintf(fmate,"\n********* %dRR%d (subline 1) *********\n", NINDTP, NINDTP); init=tsubL; mINT=mINT_TP;}
		else 		  {lin=6; fprintf(fmate,"\n********* %dRR%d (subline 2) *********\n", NINDTP, NINDAL); init=tsubL; mINT=mINT_AL;}

		if (numMIG==99)	end = generations;
		else		end = init+preNR+mINT*(numMIG-1);	

		fprintf(fmate,"\ngen  aveWf_mig   aveWf_line   NSEGLOC_mig   NSEGLOC_line   mut_NP(line)   shared(line&mig)   lost_NP(line)   mut_line   lost_line     AA_mig      AA_line       Aa_mig     Aa_line     aa_mig     aa_line    ave_s_mig   ave_s_line   ave_h_mig   ave_h_line     ave_q_mig    ave_q_line     couples_mig-mig     couples_mig-line     couples_line-line\n");
		for (gen=init; gen<=end; gen+=mINT)
		{
			fprintf(fmate, "%d    %f    %f     %f    %f    %f       %f     %f     %f    %f    %f  %f  %f  %f  %f   %f    %f   %f    %f     %f      %f      %f          %f             %f           %f\n", gen, accmean(&FEC_mig[lin][gen]), accmean(&FEC_line[lin][gen]), accmean(&NSEGLOCm[lin][gen]), accmean(&NSEGLOCl[lin][gen]), accmean(&mutNP[lin][gen]), accmean(&mutShared[lin][gen]), accmean(&lostNP[lin][gen]), accmean(&mutLine[lin][gen]), accmean(&lostLine[lin][gen]), accmean(&AAm[lin][gen]), accmean(&AAl[lin][gen]), accmean(&Aam[lin][gen]), accmean(&Aal[lin][gen]), accmean(&aam[lin][gen]), accmean(&aal[lin][gen]), accmean(&sm[lin][gen]), accmean(&sl[lin][gen]), accmean(&hm[lin][gen]), accmean(&hl[lin][gen]), accmean(&qm[lin][gen]), accmean(&ql[lin][gen]), accmean(&mate_mm[lin][gen]), accmean(&mate_ml[lin][gen]), accmean(&mate_ll[lin][gen]));
		}

		mig = 0;
		fprintf (fmate,"\nGENETIC LOAD\n");
		fprintf(fmate,"\nGen  dq_line(ALL)   dq_mig(ALL)  dq_line(EXTINCT)   dq_mig(EXTINCT)  dq_line(SURVIVORS)   dq_mig(SURVIVORS)\n");
		for (gen=init; gen<=end; gen+=mINT)
		{
			mig += 1;
			fprintf(fmate, "%d    %f       %f       %f           %f          %f            %f\n", gen, accmean(&gload_line[lin][mig]), accmean(&gload_mig[lin][mig]), accmean(&ext_gload_line[lin][mig]), accmean(&ext_gload_mig[lin][mig]), accmean(&sur_gload_line[lin][mig]), accmean(&sur_gload_mig[lin][mig]));
		}

		mig = 0;
		fprintf (fmate,"\nSTERILITY(average number of segregating loci with s=-1.0)\n");
		fprintf(fmate,"\nGen  line(ALL)   mig(ALL)  line(EXTINCT)   mig(EXTINCT)  line(SURVIVORS)   mig(SURVIVORS)\n");
		for (gen=init; gen<=end; gen+=mINT)
		{
			mig += 1;
			fprintf(fmate, "%d    %f   %f    %f        %f       %f          %f\n", gen, accmean(&STR_line[lin][mig]), accmean(&STR_mig[lin][mig]), accmean(&ext_STR_line[lin][mig]), accmean(&ext_STR_mig[lin][mig]), accmean(&sur_STR_line[lin][mig]), accmean(&sur_STR_mig[lin][mig]));
		}

		fprintf (fmate,"\nCONTRIBUTIONS\n");
		fprintf (fmate, "\nGen		   "); 
		for (i=init; i<=generations; i++)	fprintf (fmate, "%d        ", i);
		fprintf (fmate, "\n");
		fprintf (fmate, "wf=0.0-0.2     ");
		for (i=init; i<=generations; i++) fprintf (fmate, "%8.2f  ", accmean(&off_0_02[lin][i]));
		fprintf (fmate, "\n");
		fprintf (fmate, "wf=0.2-0.4     ");
		for (i=init; i<=generations; i++) fprintf (fmate, "%8.2f  ", accmean(&off_02_04[lin][i]));
		fprintf (fmate, "\n");
		fprintf (fmate, "wf=0.4-0.6     ");
		for (i=init; i<=generations; i++) fprintf (fmate, "%8.2f  ", accmean(&off_04_06[lin][i]));
		fprintf (fmate, "\n");
		fprintf (fmate, "wf=0.6-0.8     ");
		for (i=init; i<=generations; i++) fprintf (fmate, "%8.2f  ", accmean(&off_06_08[lin][i]));
		fprintf (fmate, "\n");
		fprintf (fmate, "wf=0.8-1.0     ");
		for (i=init; i<=generations; i++) fprintf (fmate, "%8.2f  ", accmean(&off_08_1[lin][i]));
		fprintf (fmate, "\n");
	}

	/* CONTROLS */
	for (nlin=1; nlin<=4; nlin++)
	{
		if (nlin==1)	  {lin=3; fprintf(fmate,"\n********* %dNR%d (line 3) *********\n", NINDTP, NINDAL); init=tLINES; mINT=mINT_AL;}
		else if (nlin==2) {lin=4; fprintf(fmate,"\n********* %dNR (line 4) *********\n", NINDTP); init=tLINES; mINT=mINT_TP;}
		else if (nlin==3) {lin=7; fprintf(fmate,"\n********* %dNR%d (subline 3) *********\n", NINDTP, NINDTP); init=tsubL; mINT=mINT_TP;}
		else		  {lin=8; fprintf(fmate,"\n********* %dNR%d (subline 4) *********\n", NINDTP, NINDAL); init=tsubL; mINT=mINT_AL; }

		if (numMIG==99)
		{
			if (lin==4)	end = tsubL-1;	
			else		end = generations;
		}
		else	end = init+preNR+mINT*(numMIG-1);	

		fprintf(fmate,"\ngen   NSEGLOC      mut_NP(line)   lost_NP(line)   mut_line   lost_line        AA             Aa           aa          ave_s      ave_h     ave_q\n");
		for (gen=init; gen<=end; gen+=mINT)
		{
			fprintf(fmate, "%d    %f    %f    %f    %f    %f    %f     %f    %f    %f  %f  %f\n", gen, accmean(&NSEGLOCl[lin][gen]), accmean(&mutNP[lin][gen]), accmean(&lostNP[lin][gen]), accmean(&mutLine[lin][gen]), accmean(&lostLine[lin][gen]), accmean(&AAl[lin][gen]), accmean(&Aal[lin][gen]), accmean(&aal[lin][gen]), accmean(&sl[lin][gen]), accmean(&hl[lin][gen]), accmean(&ql[lin][gen]));
		}

		mig = 0;
		fprintf (fmate,"\nGENETIC LOAD\n");
		fprintf(fmate,"\nGen  dq_line(ALL)   dq_line(EXTINCT)   dq_line(SURVIVORS)\n");
		for (gen=init; gen<=end; gen+=mINT)
		{
			mig += 1;
			fprintf(fmate, "%d    %f         %f            %f\n", gen, accmean(&gload_line[lin][mig]), accmean(&ext_gload_line[lin][mig]), accmean(&sur_gload_line[lin][mig]));
		}

		mig = 0;
		fprintf (fmate,"\nSTERILITY(average number of segregating loci with s=-1.0)\n");
		fprintf(fmate,"\nGen  line(ALL)   line(EXTINCT)   line(SURVIVORS)\n");
		for (gen=init; gen<=end; gen+=mINT)
		{
			mig += 1;
			fprintf(fmate, "%d    %f       %f       %f\n", gen, accmean(&STR_line[lin][mig]), accmean(&ext_STR_line[lin][mig]), accmean(&sur_STR_line[lin][mig]));
		}
	}
}

/* ********************************************************************* */
