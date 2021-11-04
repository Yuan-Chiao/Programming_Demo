/* Biomechanical Considerations for Abdominal Loading by Seat Belt Pretensioners*/
/* Author: Yuan-Chiao Lu */
/* Data source: University of Virginia TRW Project*/

DATA Abdomen_Raw;
   INPUT ID     Source $ Age Gender $ Height Weight MAIS Depth PBF Disp Comp Vel VCMAX VC FCMAX;
   DATALINES;
1 Trosseille      76      M 1.75 78      2      230      6.11      68      29.56521739      8.2      2.42      1.13      1.8
2 Trosseille      81      M 1.68 70      2      272      7.04      68      25      11.3      2.83      1.69      1.76
3 Trosseille      85      M 1.65 51      4      235      10.3      76      32.34042553      11.5      3.72      2.1      3.33
4 Trosseille      64      F 1.49 45      0      206      7.47      58      28.15533981      11.4      3.21      1.15      2.1
5 Trosseille      86      F 1.5 50      2      207      7.59      62      29.95169082      11.7      3.5      0.85      2.27
6 Foster      24      M 1.8 95.7    2      302      9.48      132.2      43.77483444      9.42      4.12      1.8      4.15
7 Foster      58      M 1.88 111.2   0      356      10      126.1      35.42134831      6.89      2.44      1.11      3.55
8 Foster      80      M 1.668 57.6    3      252      7.94      138.4      54.92063492      13.3      7.3      4.13      4.36
9 Foster      83      M 1.73 81.6    3      259      9.64      130.3      50.30888031      8.51      4.28      2.03      4.85
10 Foster      85      M 1.685 80.7    0      360      5.72      98.7      27.41666667      6.31      1.73      0.58      1.57
11 Foster      45      M 1.735 74.8    0      258      4.98      92.3      35.7751938      6.13      2.19      1.1      1.78
12 Foster      59      M 1.69 62.1    0      261      5.77      99.5      38.12260536      7.54      2.87      1.17      2.2
13 Foster      86      F 1.591 37.9    0      191      3.1      49.4      25.86387435      5.35      1.38      0.57      0.8
14 Foster      86      F .   .  0      227      2.78      55.7      24.53744493      3.95      0.97      0.44      0.68
15 Hardy      77      F 1.68 53      0      239      4      79      33.05439331      3.2      1.06      0.7      1.32
16 Hardy      77      F 1.68 53      0      310      6.1      90      29.03225806      3.4      0.99      1.1      1.77
17 Hardy      78      M 1.7 52      0      208      3.1      75      36.05769231      2.1      0.76      0.8      1.12
18 Hardy      78      M 1.7 52      0      215      4.1      56      26.04651163      2.9      0.75      0.6      1.07
19 Hardy      88      M 1.56 72      0      288      4.3      95      32.98611111      3.7      1.22      0.9      1.42
20 Hardy      88      M 1.56 72      0      308      4.5      114      37.01298701      5.6      2.07      1.2      1.67
21 Steffan      47      M 1.9 75      0      229      .      76.9      33.58078603      .      .      .      .
22 Steffan      49      M 1.82 85      0      270      6      95.8      35.48148148      .      .      .      2.13
23 Steffan      73      F 1.62 71      0      262      4.9      109.4      41.75572519      .      .      .      2.05
24 Steffan      42.5      F 1.66 60      0      239      7.4      75      31.38075314      .      .      .      2.32
25 Steffan      58      M 1.74 100     0      349      5.5      165.6      47.44985673      .      .      .      2.61
26 Steffan      59      M 1.8 95.5    0      294      10.6      103.7      35.27210884      .      .      .      3.74
27 Steffan      50      M 1.75 75      0      286      11.3      113.6      39.72027972      .      .      .      4.49
28 Steffan      87      F 1.71 70.1    0      266      .      150.2      56.46616541      .      .      .      .
29 Steffan      66      M 1.69 52      3      225      12.6      122.7      54.53333333      .      .      .      6.88
30 Steffan      54      M 1.72 79.3    0      277      10.9      157.3      56.78700361      .      .      .      6.2
31 Steffan      95      F 1.51 43.5    3      187      .      107.4      57.43315508      .      .      .      .
32 Steffan      69      M 1.83 72.2    3      278      9.2      153.4      55.17985612      .      .      .      5.08
33 Steffan      84      F 1.48 52.2    0      210      12.4      127      60.47619048      .      .      .      7.5
34 UVATRW      79      M         1.81  71      2      236      8.22      84        35.59      .        .      .          .
35 UVATRW      64      M         1.96  95      2      252      5.75      98        38.89      .        .      .          .
36 UVATRW      75      M         1.78  86      3      321      5.20      125       38.94      .        .      .          .
37 UVATRW      61      M         1.85  105     3      363      6.16      144.2     39.72      .        .      .          .
38 UVATRW      71      F 1.57 77      2      367      6.49      89.7      24.44      5.08     .      0.374      .
39 UVATRW      70      M 1.78 77      2      328      5.12      91.3      27.84      6.57     .      0.557      .
40 UVATRW      64      M 1.91 57      0      226      4.75      47.3      20.95      4.47     .      0.478      .
41 UVATRW      70      M 1.8 68      0      226      4.62      53.6      23.75      4.4      .      0.667      .
;
run;

Data Abdomen_Raw;
	set Abdomen_Raw;
    if Source='UVATRW' then
		VCMAX=(Comp/100)*Vel;
		FCMAX=(Comp/100)*PBF;
run;

/* MAIS2+ MAIS3+ */
Data Abdomen_MAIS;
	set Abdomen_Raw;
	if MAIS>=2 then
		MAIS2=1;
	else
		MAIS2=0;

	if MAIS>=3 then
		MAIS3=1;
	else
		MAIS3=0;
	BMI=Weight/(Height*Height);
run;

/*Pearson Correlation Coefficients*/
PROC CORR data=Abdomen_MAIS;
      VAR MAIS Age Height Weight BMI Depth PBF Disp Comp Vel VCMAX VC FCMAX;
run;

/*Summary statistics*/
proc means data=Abdomen_MAIS N mean median min max;
	var Age;
run;
proc freq data=Abdomen_MAIS;
	tables Gender MAIS MAIS2 MAIS3;
run;

/*********************************************/
/*    Association between BMI and injury     */
/*********************************************/
/*Logistic Regression*/
/* https://support.sas.com/documentation/cdl/en/statug/63962/HTML/default/viewer.htm#statug_logistic_sect016.htm#statug.logistic.logisticlackfit
AGGREGATE: Determines subpopulations for Pearson chi-square and deviance
LACKFIT: Requests the Hosmer and Lemeshow goodness-of-fit test
SCALE: Specifies the method to correct overdispersion
*/
proc logistic data=Abdomen_MAIS descending; /* To model 1s rather than 0s, we use the descending option. */
	model MAIS2=BMI / SCALE=NONE  aggregate lackfit;
run;

proc logistic data=Abdomen_MAIS descending;
	model MAIS3=BMI / SCALE=NONE  aggregate lackfit;
run;

/****************************************************/
/* Associations between biomech measures and injury */
/****************************************************/
/***********************/
/* Logistic Regression */
/***********************/
/*Logistic Regression for MAIS2*/
goptions reset=all;
PROC GPLOT DATA=Abdomen_MAIS;
        PLOT MAIS2*VCMAX;
run;

proc logistic data=Abdomen_MAIS descending;
		class Gender; /*Variable Gender should be either numeric or in CLASSES list.*/
        model MAIS2=Gender PBF / SCALE=NONE  aggregate lackfit;
run;
proc logistic data=Abdomen_MAIS descending;
		class Gender;
        model MAIS2=Gender Disp / SCALE=NONE  aggregate lackfit;
run;
proc logistic data=Abdomen_MAIS descending;
		class Gender;
        model MAIS2=Gender Comp / SCALE=NONE  aggregate lackfit;
run;
proc logistic data=Abdomen_MAIS descending;
		class Gender;
        model MAIS2=Gender Vel / SCALE=NONE  aggregate lackfit;
run;
proc logistic data=Abdomen_MAIS descending;
		class Gender;
        model MAIS2=Gender VCMAX / SCALE=NONE  aggregate lackfit;
run;
proc logistic data=Abdomen_MAIS descending;
		class Gender;
        model MAIS2=Gender VC / SCALE=NONE  aggregate lackfit;
run;
proc logistic data=Abdomen_MAIS descending;
		class Gender;
        model MAIS2=Gender FCMAX / SCALE=NONE  aggregate lackfit;
run;

/* Best model is fitted by Vel */
/*Plot of Estimated Prob vs. Vel*/
proc logistic data=Abdomen_MAIS descending;
        model MAIS2=Vel / SCALE=NONE  aggregate lackfit;
        Output OUT=Abdomen_MAIS2_Vel P=PRED_MAIS2_Vel L=LOWER_MAIS2_Vel U=UPPER_MAIS2_Vel;
run;
PROC SORT DATA=Abdomen_MAIS2_Vel;
        BY Vel;
run;
symbol1 i = join v=circle l=32  c = black;
symbol2 i = join v=star l=32  c = black;
symbol3 i = join v=star l=32  c = black;
PROC GPLOT DATA=Abdomen_MAIS2_Vel;
        PLOT PRED_MAIS2_Vel*Vel      LOWER_MAIS2_Vel*Vel      UPPER_MAIS2_Vel*Vel / OVERLAY;
RUN;

/* EXPORT DATA */
PROC EXPORT DATA= Abdomen_MAIS2_Vel
            OUTFILE= "/home/luyuanchiao/Abdomen_MAIS2_Vel.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

/*Logistic Regression for MAIS3*/
proc logistic data=Abdomen_MAIS descending;
		class Gender;
        model MAIS3=Gender PBF / SCALE=NONE  aggregate lackfit;
run;
proc logistic data=Abdomen_MAIS descending;
		class Gender;
        model MAIS3=Gender Disp / SCALE=NONE  aggregate lackfit;
run;
proc logistic data=Abdomen_MAIS descending;
		class Gender;
        model MAIS3=Gender Comp / SCALE=NONE  aggregate lackfit;
run;
proc logistic data=Abdomen_MAIS descending;
		class Gender;
        model MAIS3=Gender Vel / SCALE=NONE  aggregate lackfit;
run;
proc logistic data=Abdomen_MAIS descending;
		class Gender;
        model MAIS3=Gender VCMAX / SCALE=NONE  aggregate lackfit;
run;
proc logistic data=Abdomen_MAIS descending;
		class Gender;
        model MAIS3=Gender VC / SCALE=NONE  aggregate lackfit;
run;
proc logistic data=Abdomen_MAIS descending;
		class Gender;
        model MAIS3=Gender FCMAX / SCALE=NONE  aggregate lackfit;
run;

/* Best model is fitted by FCMAX */
/*Plot of Estimated Prob vs. FCMAX*/
proc logistic data=Abdomen_MAIS descending;
        model MAIS3=FCMAX / SCALE=NONE  aggregate lackfit;
        Output OUT=Abdomen_MAIS3_FCMAX P=PRED_MAIS3_FCMAX L=LOWER_MAIS3_FCMAX U=UPPER_MAIS3_FCMAX;
run;
PROC SORT DATA=Abdomen_MAIS3_FCMAX;
        BY FCMAX;
run;
symbol1 i = join v=circle l=32  c = black;
symbol2 i = join v=star l=32  c = black;
symbol3 i = join v=star l=32  c = black;
PROC GPLOT DATA=Abdomen_MAIS3_FCMAX;
        PLOT PRED_MAIS3_FCMAX*FCMAX      LOWER_MAIS3_FCMAX*FCMAX      UPPER_MAIS3_FCMAX*FCMAX / OVERLAY;
RUN;
PROC SORT DATA=Abdomen_MAIS3_FCMAX;
        BY FCMAX;
run;

/* EXPORT DATA */
PROC EXPORT DATA= Abdomen_MAIS3_FCMAX
            OUTFILE= "/home/luyuanchiao/Abdomen_MAIS3_FCMAX.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

/***********************/
/*  Survival Analysis  */
/***********************/
PROC SORT DATA=Abdomen_MAIS;
        BY MAIS2 Vel;
run;
PROC SORT DATA=Abdomen_MAIS;
        BY MAIS3 FCMAX;
run;

/* Look up the Lower Bound and Upper Bound of Vel and FCMAX for Suvival Interval Censor */
Data Abdomen_Survival;
        set Abdomen_MAIS;
        Vel_MAIS2_Lower=Vel;
        Vel_MAIS2_Upper=Vel;
        FCMAX_MAIS3_Lower=FCMAX;
        FCMAX_MAIS3_Upper=FCMAX;

        if MAIS>=2 then
                Vel_MAIS2_Lower=5.08;
        else
                Vel_MAIS2_Upper=11.4;

        if MAIS>=3 then
                FCMAX_MAIS3_Lower=3.33;
        else
                FCMAX_MAIS3_Upper=7.5;

        if Vel=. then
        do
                Vel_MAIS2_Lower=.;
                Vel_MAIS2_Upper=.;
        end;
        if FCMAX=. then
        do
                FCMAX_MAIS3_Lower=.;
                FCMAX_MAIS3_Upper=.;
        end;
run;

/**********************/
/*  interval-censored */
/**********************/
proc lifereg data=Abdomen_Survival;
		model (Vel_MAIS2_Lower,Vel_MAIS2_Upper)= / dist=LOGNORMAL CORRB;
/* expecting one of the following: EXP, EXPONENTIAL, GAMMA, LLOGISTIC, LNORMAL, LOGISTIC, LOGNORMAL, NORMAL, WEIBULL. */
        Output OUT=Abdomen_MAIS2_Survival CDF=Pred_MAIS2_Survival SRES=Res_MAIS2_Survival;
/* expecting one of the following: ;, CDF, CENSORED, CONTROL, CRES, CRESIDUAL, OUT, P, PREDICTED, Q, QUANTILE, QUANTILES, SRES, SRESIDUAL, STD, STD_ERR, XBETA. */
run;

proc lifereg data=Abdomen_Survival;
        model (FCMAX_MAIS3_Lower,FCMAX_MAIS3_Upper)= / dist=LOGNORMAL CORRB;
		Output OUT=Abdomen_MAIS3_Survival CDF=Pred_MAIS3_Survival SRES=Res_MAIS3_Survival;
run;

proc lifetest data=Abdomen_Survival(where=(MAIS3=1)) plots=survival(atrisk);
		time FCMAX*MAIS3(0);
run; 

/**********************/
/*   Test Normality   */
/**********************/
/* Anderson-Darling Statistic */
/* The Anderson-Darling test (Stephens, 1974) is used to test if a sample of data came from a population with a specific distribution. */
proc univariate data=Abdomen_MAIS2_Survival;
        var Res_MAIS2_Survival;
        histogram Res_MAIS2_Survival / normal;
        /* expecting one of the following: ;, ANNOKEY, ANNOTATE, BARLABEL, */
		/*               BARWIDTH, BETA, CAXIS, CBARLINE, CFILL, CFRAME, CFRAMESIDE, CFRAMETOP, CGRID, */
		/*               CHREF, CLIPREF, CONTENTS, CPROP, CTEXT, CV, CVREF, DESCRIPTION, ENDPOINTS, */
		/*               EXPONENTIAL, FONT, FRONTREF, GAMMA, GRID, HANGING, HAXIS, HEIGHT, HMINOR, */
		/*               HOFFSET, HREF, HREFLABELS, HREFLABPOS, INFONT, INHEIGHT, INTERBAR, INTERTILE, */
		/*               KERNEL, LGRID, LHREF, LOGNORMAL, LVREF, MAXNBIN, MAXSIGMAS, MIDPERCENTS, */
		/*               MIDPOINTS, NAME, NCOL, NCOLS, NENDPOINTS, NMIDPOINTS, NOBARS, NOCHART, NOFRAME, */
		/*               NOHLABEL, NOPLOT, NORMAL, NOTABCONTENTS, NOVLABEL, NOVTICK, NROW, NROWS, */
		/*               OUTHISTOGRAM, OUTKERNEL, OVERLAY, PFILL, RTINCLUDE, SB, SU, TILELEGLABEL, */
		/*               TURNVLABELS, VAXIS, VAXISLABEL, VMINOR, VOFFSET, VREF, VREFLABELS, VREFLABPOS, */
		/*               VSCALE, WAXIS, WBARLINE, WEIBULL, WGRID. */
run;

