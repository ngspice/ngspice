/**** BSIM4v7.6.2 Released by Wenwei Yang 07/31/2008 ****/

/**********
Copyright 2006 Regents of the University of California.  All rights reserved.
File: bsim4v7def.h
Author: 2000 Weidong Liu.
Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
Authors: 2008- Wenwei Yang, Ali Niknejad, Chenming Hu
Modified by Xuemei Xi, 11/15/2002.
Modified by Xuemei Xi, 05/09/2003.
Modified by Xuemei Xi, 03/04/2004.
Modified by Xuemei Xi, Mohan Dunga, 09/24/2004.
Modified by Xuemei Xi, 07/29/2005.
Modified by Mohan Dunga, 12/13/2006
Modified by Mohan Dunga, Wenwei Yang, 05/18/2007.
Modified by Wenwei Yang, 07/31/2008.
 * Modified by Tanvir Morshed, Darsen Lu 03/27/2011
**********/

#ifndef BSIM4v7
#define BSIM4v7

#include "ngspice/ifsim.h"
#include "ngspice/gendefs.h"
#include "ngspice/cktdefs.h"
#include "ngspice/complex.h"
#include "ngspice/noisedef.h"

typedef struct sBSIM4v7instance
{
    struct sBSIM4v7model *BSIM4v7modPtr;
    struct sBSIM4v7instance *BSIM4v7nextInstance;
    IFuid BSIM4v7name;
    int BSIM4v7states;     /* index into state table for this device */
    int BSIM4v7dNode;
    int BSIM4v7gNodeExt;
    int BSIM4v7sNode;
    int BSIM4v7bNode;
    int BSIM4v7dNodePrime;
    int BSIM4v7gNodePrime;
    int BSIM4v7gNodeMid;
    int BSIM4v7sNodePrime;
    int BSIM4v7bNodePrime;
    int BSIM4v7dbNode;
    int BSIM4v7sbNode;
    int BSIM4v7qNode;

    double BSIM4v7ueff;
    double BSIM4v7thetavth;
    double BSIM4v7von;
    double BSIM4v7vdsat;
    double BSIM4v7cgdo;
    double BSIM4v7qgdo;
    double BSIM4v7cgso;
    double BSIM4v7qgso;
    double BSIM4v7grbsb;
    double BSIM4v7grbdb;
    double BSIM4v7grbpb;
    double BSIM4v7grbps;
    double BSIM4v7grbpd;

    double BSIM4v7vjsmFwd;
    double BSIM4v7vjsmRev;
    double BSIM4v7vjdmFwd;
    double BSIM4v7vjdmRev;
    double BSIM4v7XExpBVS;
    double BSIM4v7XExpBVD;
    double BSIM4v7SslpFwd;
    double BSIM4v7SslpRev;
    double BSIM4v7DslpFwd;
    double BSIM4v7DslpRev;
    double BSIM4v7IVjsmFwd;
    double BSIM4v7IVjsmRev;
    double BSIM4v7IVjdmFwd;
    double BSIM4v7IVjdmRev;

    double BSIM4v7grgeltd;
    double BSIM4v7Pseff;
    double BSIM4v7Pdeff;
    double BSIM4v7Aseff;
    double BSIM4v7Adeff;

    double BSIM4v7l;
    double BSIM4v7w;
    double BSIM4v7drainArea;
    double BSIM4v7sourceArea;
    double BSIM4v7drainSquares;
    double BSIM4v7sourceSquares;
    double BSIM4v7drainPerimeter;
    double BSIM4v7sourcePerimeter;
    double BSIM4v7sourceConductance;
    double BSIM4v7drainConductance;
     /* stress effect instance param */
    double BSIM4v7sa;
    double BSIM4v7sb;
    double BSIM4v7sd;
    double BSIM4v7sca;
    double BSIM4v7scb;
    double BSIM4v7scc;
    double BSIM4v7sc;

    double BSIM4v7rbdb;
    double BSIM4v7rbsb;
    double BSIM4v7rbpb;
    double BSIM4v7rbps;
    double BSIM4v7rbpd;

    double BSIM4v7delvto;
    double BSIM4v7xgw;
    double BSIM4v7ngcon;

     /* added here to account stress effect instance dependence */
    double BSIM4v7u0temp;
    double BSIM4v7vsattemp;
    double BSIM4v7vth0;
    double BSIM4v7vfb;
    double BSIM4v7vfbzb;
    double BSIM4v7vtfbphi1;
    double BSIM4v7vtfbphi2;
    double BSIM4v7k2;
    double BSIM4v7vbsc;
    double BSIM4v7k2ox;
    double BSIM4v7eta0;

    double BSIM4v7icVDS;
    double BSIM4v7icVGS;
    double BSIM4v7icVBS;
    double BSIM4v7m;
    double BSIM4v7nf;
    int BSIM4v7off;
    int BSIM4v7mode;
    int BSIM4v7trnqsMod;
    int BSIM4v7acnqsMod;
    int BSIM4v7rbodyMod;
    int BSIM4v7rgateMod;
    int BSIM4v7geoMod;
    int BSIM4v7rgeoMod;
    int BSIM4v7min;


    /* OP point */
    double BSIM4v7Vgsteff;
    double BSIM4v7vgs_eff;
    double BSIM4v7vgd_eff;
    double BSIM4v7dvgs_eff_dvg;
    double BSIM4v7dvgd_eff_dvg;
    double BSIM4v7Vdseff;
    double BSIM4v7nstar;
    double BSIM4v7Abulk;
    double BSIM4v7EsatL;
    double BSIM4v7AbovVgst2Vtm;
    double BSIM4v7qinv;
    double BSIM4v7cd;
    double BSIM4v7cbs;
    double BSIM4v7cbd;
    double BSIM4v7csub;
    double BSIM4v7Igidl;
    double BSIM4v7Igisl;
    double BSIM4v7gm;
    double BSIM4v7gds;
    double BSIM4v7gmbs;
    double BSIM4v7gbd;
    double BSIM4v7gbs;
    double BSIM4v7noiGd0;   /* tnoiMod=2 (v4.7) */
    double BSIM4v7Coxeff;

    double BSIM4v7gbbs;
    double BSIM4v7gbgs;
    double BSIM4v7gbds;
    double BSIM4v7ggidld;
    double BSIM4v7ggidlg;
    double BSIM4v7ggidls;
    double BSIM4v7ggidlb;
    double BSIM4v7ggisld;
    double BSIM4v7ggislg;
    double BSIM4v7ggisls;
    double BSIM4v7ggislb;

    double BSIM4v7Igcs;
    double BSIM4v7gIgcsg;
    double BSIM4v7gIgcsd;
    double BSIM4v7gIgcss;
    double BSIM4v7gIgcsb;
    double BSIM4v7Igcd;
    double BSIM4v7gIgcdg;
    double BSIM4v7gIgcdd;
    double BSIM4v7gIgcds;
    double BSIM4v7gIgcdb;

    double BSIM4v7Igs;
    double BSIM4v7gIgsg;
    double BSIM4v7gIgss;
    double BSIM4v7Igd;
    double BSIM4v7gIgdg;
    double BSIM4v7gIgdd;

    double BSIM4v7Igb;
    double BSIM4v7gIgbg;
    double BSIM4v7gIgbd;
    double BSIM4v7gIgbs;
    double BSIM4v7gIgbb;

    double BSIM4v7grdsw;
    double BSIM4v7IdovVds;
    double BSIM4v7gcrg;
    double BSIM4v7gcrgd;
    double BSIM4v7gcrgg;
    double BSIM4v7gcrgs;
    double BSIM4v7gcrgb;

    double BSIM4v7gstot;
    double BSIM4v7gstotd;
    double BSIM4v7gstotg;
    double BSIM4v7gstots;
    double BSIM4v7gstotb;

    double BSIM4v7gdtot;
    double BSIM4v7gdtotd;
    double BSIM4v7gdtotg;
    double BSIM4v7gdtots;
    double BSIM4v7gdtotb;

    double BSIM4v7cggb;
    double BSIM4v7cgdb;
    double BSIM4v7cgsb;
    double BSIM4v7cbgb;
    double BSIM4v7cbdb;
    double BSIM4v7cbsb;
    double BSIM4v7cdgb;
    double BSIM4v7cddb;
    double BSIM4v7cdsb;
    double BSIM4v7csgb;
    double BSIM4v7csdb;
    double BSIM4v7cssb;
    double BSIM4v7cgbb;
    double BSIM4v7cdbb;
    double BSIM4v7csbb;
    double BSIM4v7cbbb;
    double BSIM4v7capbd;
    double BSIM4v7capbs;

    double BSIM4v7cqgb;
    double BSIM4v7cqdb;
    double BSIM4v7cqsb;
    double BSIM4v7cqbb;

    double BSIM4v7qgate;
    double BSIM4v7qbulk;
    double BSIM4v7qdrn;
    double BSIM4v7qsrc;
    double BSIM4v7qdef;

    double BSIM4v7qchqs;
    double BSIM4v7taunet;
    double BSIM4v7gtau;
    double BSIM4v7gtg;
    double BSIM4v7gtd;
    double BSIM4v7gts;
    double BSIM4v7gtb;
    double BSIM4v7SjctTempRevSatCur;
    double BSIM4v7DjctTempRevSatCur;
    double BSIM4v7SswTempRevSatCur;
    double BSIM4v7DswTempRevSatCur;
    double BSIM4v7SswgTempRevSatCur;
    double BSIM4v7DswgTempRevSatCur;

    struct bsim4SizeDependParam  *pParam;

    unsigned BSIM4v7lGiven :1;
    unsigned BSIM4v7wGiven :1;
    unsigned BSIM4v7mGiven :1;
    unsigned BSIM4v7nfGiven :1;
    unsigned BSIM4v7minGiven :1;
    unsigned BSIM4v7drainAreaGiven :1;
    unsigned BSIM4v7sourceAreaGiven    :1;
    unsigned BSIM4v7drainSquaresGiven  :1;
    unsigned BSIM4v7sourceSquaresGiven :1;
    unsigned BSIM4v7drainPerimeterGiven    :1;
    unsigned BSIM4v7sourcePerimeterGiven   :1;
    unsigned BSIM4v7saGiven :1;
    unsigned BSIM4v7sbGiven :1;
    unsigned BSIM4v7sdGiven :1;
    unsigned BSIM4v7scaGiven :1;
    unsigned BSIM4v7scbGiven :1;
    unsigned BSIM4v7sccGiven :1;
    unsigned BSIM4v7scGiven :1;
    unsigned BSIM4v7rbdbGiven   :1;
    unsigned BSIM4v7rbsbGiven   :1;
    unsigned BSIM4v7rbpbGiven   :1;
    unsigned BSIM4v7rbpdGiven   :1;
    unsigned BSIM4v7rbpsGiven   :1;
    unsigned BSIM4v7delvtoGiven   :1;
    unsigned BSIM4v7xgwGiven   :1;
    unsigned BSIM4v7ngconGiven   :1;
    unsigned BSIM4v7icVDSGiven :1;
    unsigned BSIM4v7icVGSGiven :1;
    unsigned BSIM4v7icVBSGiven :1;
    unsigned BSIM4v7trnqsModGiven :1;
    unsigned BSIM4v7acnqsModGiven :1;
    unsigned BSIM4v7rbodyModGiven :1;
    unsigned BSIM4v7rgateModGiven :1;
    unsigned BSIM4v7geoModGiven :1;
    unsigned BSIM4v7rgeoModGiven :1;


    double *BSIM4v7DPdPtr;
    double *BSIM4v7DPdpPtr;
    double *BSIM4v7DPgpPtr;
    double *BSIM4v7DPgmPtr;
    double *BSIM4v7DPspPtr;
    double *BSIM4v7DPbpPtr;
    double *BSIM4v7DPdbPtr;

    double *BSIM4v7DdPtr;
    double *BSIM4v7DdpPtr;

    double *BSIM4v7GPdpPtr;
    double *BSIM4v7GPgpPtr;
    double *BSIM4v7GPgmPtr;
    double *BSIM4v7GPgePtr;
    double *BSIM4v7GPspPtr;
    double *BSIM4v7GPbpPtr;

    double *BSIM4v7GMdpPtr;
    double *BSIM4v7GMgpPtr;
    double *BSIM4v7GMgmPtr;
    double *BSIM4v7GMgePtr;
    double *BSIM4v7GMspPtr;
    double *BSIM4v7GMbpPtr;

    double *BSIM4v7GEdpPtr;
    double *BSIM4v7GEgpPtr;
    double *BSIM4v7GEgmPtr;
    double *BSIM4v7GEgePtr;
    double *BSIM4v7GEspPtr;
    double *BSIM4v7GEbpPtr;

    double *BSIM4v7SPdpPtr;
    double *BSIM4v7SPgpPtr;
    double *BSIM4v7SPgmPtr;
    double *BSIM4v7SPsPtr;
    double *BSIM4v7SPspPtr;
    double *BSIM4v7SPbpPtr;
    double *BSIM4v7SPsbPtr;

    double *BSIM4v7SspPtr;
    double *BSIM4v7SsPtr;

    double *BSIM4v7BPdpPtr;
    double *BSIM4v7BPgpPtr;
    double *BSIM4v7BPgmPtr;
    double *BSIM4v7BPspPtr;
    double *BSIM4v7BPdbPtr;
    double *BSIM4v7BPbPtr;
    double *BSIM4v7BPsbPtr;
    double *BSIM4v7BPbpPtr;

    double *BSIM4v7DBdpPtr;
    double *BSIM4v7DBdbPtr;
    double *BSIM4v7DBbpPtr;
    double *BSIM4v7DBbPtr;

    double *BSIM4v7SBspPtr;
    double *BSIM4v7SBbpPtr;
    double *BSIM4v7SBbPtr;
    double *BSIM4v7SBsbPtr;

    double *BSIM4v7BdbPtr;
    double *BSIM4v7BbpPtr;
    double *BSIM4v7BsbPtr;
    double *BSIM4v7BbPtr;

    double *BSIM4v7DgpPtr;
    double *BSIM4v7DspPtr;
    double *BSIM4v7DbpPtr;
    double *BSIM4v7SdpPtr;
    double *BSIM4v7SgpPtr;
    double *BSIM4v7SbpPtr;

    double *BSIM4v7QdpPtr;
    double *BSIM4v7QgpPtr;
    double *BSIM4v7QspPtr;
    double *BSIM4v7QbpPtr;
    double *BSIM4v7QqPtr;
    double *BSIM4v7DPqPtr;
    double *BSIM4v7GPqPtr;
    double *BSIM4v7SPqPtr;

#ifdef USE_OMP
    /* per instance storage of results, to update matrix at a later stge */
    double BSIM4v7rhsdPrime;
    double BSIM4v7rhsgPrime;
    double BSIM4v7rhsgExt;
    double BSIM4v7grhsMid;
    double BSIM4v7rhsbPrime;
    double BSIM4v7rhssPrime;
    double BSIM4v7rhsdb;
    double BSIM4v7rhssb;
    double BSIM4v7rhsd;
    double BSIM4v7rhss;
    double BSIM4v7rhsq;

    double BSIM4v7_1;
    double BSIM4v7_2;
    double BSIM4v7_3;
    double BSIM4v7_4;
    double BSIM4v7_5;
    double BSIM4v7_6;
    double BSIM4v7_7;
    double BSIM4v7_8;
    double BSIM4v7_9;
    double BSIM4v7_10;
    double BSIM4v7_11;
    double BSIM4v7_12;
    double BSIM4v7_13;
    double BSIM4v7_14;
    double BSIM4v7_15;
    double BSIM4v7_16;
    double BSIM4v7_17;
    double BSIM4v7_18;
    double BSIM4v7_19;
    double BSIM4v7_20;
    double BSIM4v7_21;
    double BSIM4v7_22;
    double BSIM4v7_23;
    double BSIM4v7_24;
    double BSIM4v7_25;
    double BSIM4v7_26;
    double BSIM4v7_27;
    double BSIM4v7_28;
    double BSIM4v7_29;
    double BSIM4v7_30;
    double BSIM4v7_31;
    double BSIM4v7_32;
    double BSIM4v7_33;
    double BSIM4v7_34;
    double BSIM4v7_35;
    double BSIM4v7_36;
    double BSIM4v7_37;
    double BSIM4v7_38;
    double BSIM4v7_39;
    double BSIM4v7_40;
    double BSIM4v7_41;
    double BSIM4v7_42;
    double BSIM4v7_43;
    double BSIM4v7_44;
    double BSIM4v7_45;
    double BSIM4v7_46;
    double BSIM4v7_47;
    double BSIM4v7_48;
    double BSIM4v7_49;
    double BSIM4v7_50;
    double BSIM4v7_51;
    double BSIM4v7_52;
    double BSIM4v7_53;
    double BSIM4v7_54;
    double BSIM4v7_55;
    double BSIM4v7_56;
    double BSIM4v7_57;
    double BSIM4v7_58;
    double BSIM4v7_59;
    double BSIM4v7_60;
    double BSIM4v7_61;
    double BSIM4v7_62;
    double BSIM4v7_63;
    double BSIM4v7_64;
    double BSIM4v7_65;
    double BSIM4v7_66;
    double BSIM4v7_67;
    double BSIM4v7_68;
    double BSIM4v7_69;
    double BSIM4v7_70;
    double BSIM4v7_71;
    double BSIM4v7_72;
    double BSIM4v7_73;
    double BSIM4v7_74;
    double BSIM4v7_75;
    double BSIM4v7_76;
    double BSIM4v7_77;
    double BSIM4v7_78;
    double BSIM4v7_79;
    double BSIM4v7_80;
    double BSIM4v7_81;
    double BSIM4v7_82;
    double BSIM4v7_83;
    double BSIM4v7_84;
    double BSIM4v7_85;
    double BSIM4v7_86;
    double BSIM4v7_87;
    double BSIM4v7_88;
    double BSIM4v7_89;
    double BSIM4v7_90;
    double BSIM4v7_91;
    double BSIM4v7_92;
    double BSIM4v7_93;
    double BSIM4v7_94;
    double BSIM4v7_95;
    double BSIM4v7_96;
    double BSIM4v7_97;
    double BSIM4v7_98;
    double BSIM4v7_99;
    double BSIM4v7_100;
    double BSIM4v7_101;
    double BSIM4v7_102;
    double BSIM4v7_103;
#endif

#define BSIM4v7vbd BSIM4v7states+ 0
#define BSIM4v7vbs BSIM4v7states+ 1
#define BSIM4v7vgs BSIM4v7states+ 2
#define BSIM4v7vds BSIM4v7states+ 3
#define BSIM4v7vdbs BSIM4v7states+ 4
#define BSIM4v7vdbd BSIM4v7states+ 5
#define BSIM4v7vsbs BSIM4v7states+ 6
#define BSIM4v7vges BSIM4v7states+ 7
#define BSIM4v7vgms BSIM4v7states+ 8
#define BSIM4v7vses BSIM4v7states+ 9
#define BSIM4v7vdes BSIM4v7states+ 10

#define BSIM4v7qb BSIM4v7states+ 11
#define BSIM4v7cqb BSIM4v7states+ 12
#define BSIM4v7qg BSIM4v7states+ 13
#define BSIM4v7cqg BSIM4v7states+ 14
#define BSIM4v7qd BSIM4v7states+ 15
#define BSIM4v7cqd BSIM4v7states+ 16
#define BSIM4v7qgmid BSIM4v7states+ 17
#define BSIM4v7cqgmid BSIM4v7states+ 18

#define BSIM4v7qbs  BSIM4v7states+ 19
#define BSIM4v7cqbs  BSIM4v7states+ 20
#define BSIM4v7qbd  BSIM4v7states+ 21
#define BSIM4v7cqbd  BSIM4v7states+ 22

#define BSIM4v7qcheq BSIM4v7states+ 23
#define BSIM4v7cqcheq BSIM4v7states+ 24
#define BSIM4v7qcdump BSIM4v7states+ 25
#define BSIM4v7cqcdump BSIM4v7states+ 26
#define BSIM4v7qdef BSIM4v7states+ 27
#define BSIM4v7qs BSIM4v7states+ 28

#define BSIM4v7numStates 29


/* indices to the array of BSIM4v7 NOISE SOURCES */

#define BSIM4v7RDNOIZ       0
#define BSIM4v7RSNOIZ       1
#define BSIM4v7RGNOIZ       2
#define BSIM4v7RBPSNOIZ     3
#define BSIM4v7RBPDNOIZ     4
#define BSIM4v7RBPBNOIZ     5
#define BSIM4v7RBSBNOIZ     6
#define BSIM4v7RBDBNOIZ     7
#define BSIM4v7IDNOIZ       8
#define BSIM4v7FLNOIZ       9
#define BSIM4v7IGSNOIZ      10
#define BSIM4v7IGDNOIZ      11
#define BSIM4v7IGBNOIZ      12
#define BSIM4v7CORLNOIZ     13
#define BSIM4v7TOTNOIZ      14

#define BSIM4v7NSRCS        15  /* Number of BSIM4v7 noise sources */

#ifndef NONOISE
    double BSIM4v7nVar[NSTATVARS][BSIM4v7NSRCS];
#else /* NONOISE */
        double **BSIM4v7nVar;
#endif /* NONOISE */

} BSIM4v7instance ;

struct bsim4SizeDependParam
{
    double Width;
    double Length;
    double NFinger;

    double BSIM4v7cdsc;
    double BSIM4v7cdscb;
    double BSIM4v7cdscd;
    double BSIM4v7cit;
    double BSIM4v7nfactor;
    double BSIM4v7xj;
    double BSIM4v7vsat;
    double BSIM4v7at;
    double BSIM4v7a0;
    double BSIM4v7ags;
    double BSIM4v7a1;
    double BSIM4v7a2;
    double BSIM4v7keta;
    double BSIM4v7nsub;
    double BSIM4v7ndep;
    double BSIM4v7nsd;
    double BSIM4v7phin;
    double BSIM4v7ngate;
    double BSIM4v7gamma1;
    double BSIM4v7gamma2;
    double BSIM4v7vbx;
    double BSIM4v7vbi;
    double BSIM4v7vbm;
    double BSIM4v7xt;
    double BSIM4v7phi;
    double BSIM4v7litl;
    double BSIM4v7k1;
    double BSIM4v7kt1;
    double BSIM4v7kt1l;
    double BSIM4v7kt2;
    double BSIM4v7k2;
    double BSIM4v7k3;
    double BSIM4v7k3b;
    double BSIM4v7w0;
    double BSIM4v7dvtp0;
    double BSIM4v7dvtp1;
    double BSIM4v7dvtp2;	/* New DIBL/Rout */
    double BSIM4v7dvtp3;
    double BSIM4v7dvtp4;
    double BSIM4v7dvtp5;
    double BSIM4v7lpe0;
    double BSIM4v7lpeb;
    double BSIM4v7dvt0;
    double BSIM4v7dvt1;
    double BSIM4v7dvt2;
    double BSIM4v7dvt0w;
    double BSIM4v7dvt1w;
    double BSIM4v7dvt2w;
    double BSIM4v7drout;
    double BSIM4v7dsub;
    double BSIM4v7vth0;
    double BSIM4v7ua;
    double BSIM4v7ua1;
    double BSIM4v7ub;
    double BSIM4v7ub1;
    double BSIM4v7uc;
    double BSIM4v7uc1;
    double BSIM4v7ud;
    double BSIM4v7ud1;
    double BSIM4v7up;
    double BSIM4v7lp;
    double BSIM4v7u0;
    double BSIM4v7eu;
  	double BSIM4v7ucs;
    double BSIM4v7ute;
	  double BSIM4v7ucste;
    double BSIM4v7voff;
    double BSIM4v7tvoff;
    double BSIM4v7tnfactor; 	/* v4.7 Temp dep of leakage current */
    double BSIM4v7teta0;   	/* v4.7 temp dep of leakage current */
    double BSIM4v7tvoffcv;	/* v4.7 temp dep of leakage current */
    double BSIM4v7minv;
    double BSIM4v7minvcv;
    double BSIM4v7vfb;
    double BSIM4v7delta;
    double BSIM4v7rdsw;
    double BSIM4v7rds0;
    double BSIM4v7rs0;
    double BSIM4v7rd0;
    double BSIM4v7rsw;
    double BSIM4v7rdw;
    double BSIM4v7prwg;
    double BSIM4v7prwb;
    double BSIM4v7prt;
    double BSIM4v7eta0;
    double BSIM4v7etab;
    double BSIM4v7pclm;
    double BSIM4v7pdibl1;
    double BSIM4v7pdibl2;
    double BSIM4v7pdiblb;
    double BSIM4v7fprout;
    double BSIM4v7pdits;
    double BSIM4v7pditsd;
    double BSIM4v7pscbe1;
    double BSIM4v7pscbe2;
    double BSIM4v7pvag;
    double BSIM4v7wr;
    double BSIM4v7dwg;
    double BSIM4v7dwb;
    double BSIM4v7b0;
    double BSIM4v7b1;
    double BSIM4v7alpha0;
    double BSIM4v7alpha1;
    double BSIM4v7beta0;
    double BSIM4v7agidl;
    double BSIM4v7bgidl;
    double BSIM4v7cgidl;
    double BSIM4v7egidl;
    double BSIM4v7fgidl; /* v4.7 New GIDL/GISL */
    double BSIM4v7kgidl; /* v4.7 New GIDL/GISL */
    double BSIM4v7rgidl; /* v4.7 New GIDL/GISL */
    double BSIM4v7agisl;
    double BSIM4v7bgisl;
    double BSIM4v7cgisl;
    double BSIM4v7egisl;
    double BSIM4v7fgisl; /* v4.7 New GIDL/GISL */
    double BSIM4v7kgisl; /* v4.7 New GIDL/GISL */
    double BSIM4v7rgisl; /* v4.7 New GIDL/GISL */
    double BSIM4v7aigc;
    double BSIM4v7bigc;
    double BSIM4v7cigc;
    double BSIM4v7aigs;
    double BSIM4v7bigs;
    double BSIM4v7cigs;
    double BSIM4v7aigd;
    double BSIM4v7bigd;
    double BSIM4v7cigd;
    double BSIM4v7aigbacc;
    double BSIM4v7bigbacc;
    double BSIM4v7cigbacc;
    double BSIM4v7aigbinv;
    double BSIM4v7bigbinv;
    double BSIM4v7cigbinv;
    double BSIM4v7nigc;
    double BSIM4v7nigbacc;
    double BSIM4v7nigbinv;
    double BSIM4v7ntox;
    double BSIM4v7eigbinv;
    double BSIM4v7pigcd;
    double BSIM4v7poxedge;
    double BSIM4v7xrcrg1;
    double BSIM4v7xrcrg2;
    double BSIM4v7lambda; /* overshoot */
    double BSIM4v7vtl; /* thermal velocity limit */
    double BSIM4v7xn; /* back scattering parameter */
    double BSIM4v7lc; /* back scattering parameter */
    double BSIM4v7tfactor;  /* ballistic transportation factor  */
    double BSIM4v7vfbsdoff;  /* S/D flatband offset voltage  */
    double BSIM4v7tvfbsdoff;

/* added for stress effect */
    double BSIM4v7ku0;
    double BSIM4v7kvth0;
    double BSIM4v7ku0temp;
    double BSIM4v7rho_ref;
    double BSIM4v7inv_od_ref;
/* added for well proximity effect */
    double BSIM4v7kvth0we;
    double BSIM4v7k2we;
    double BSIM4v7ku0we;

    /* CV model */
    double BSIM4v7cgsl;
    double BSIM4v7cgdl;
    double BSIM4v7ckappas;
    double BSIM4v7ckappad;
    double BSIM4v7cf;
    double BSIM4v7clc;
    double BSIM4v7cle;
    double BSIM4v7vfbcv;
    double BSIM4v7noff;
    double BSIM4v7voffcv;
    double BSIM4v7acde;
    double BSIM4v7moin;

/* Pre-calculated constants */

    double BSIM4v7dw;
    double BSIM4v7dl;
    double BSIM4v7leff;
    double BSIM4v7weff;

    double BSIM4v7dwc;
    double BSIM4v7dlc;
    double BSIM4v7dwj;
    double BSIM4v7leffCV;
    double BSIM4v7weffCV;
    double BSIM4v7weffCJ;
    double BSIM4v7abulkCVfactor;
    double BSIM4v7cgso;
    double BSIM4v7cgdo;
    double BSIM4v7cgbo;

    double BSIM4v7u0temp;
    double BSIM4v7vsattemp;
    double BSIM4v7sqrtPhi;
    double BSIM4v7phis3;
    double BSIM4v7Xdep0;
    double BSIM4v7sqrtXdep0;
    double BSIM4v7theta0vb0;
    double BSIM4v7thetaRout;
    double BSIM4v7mstar;
	  double BSIM4v7VgsteffVth;
    double BSIM4v7mstarcv;
    double BSIM4v7voffcbn;
    double BSIM4v7voffcbncv;
    double BSIM4v7rdswmin;
    double BSIM4v7rdwmin;
    double BSIM4v7rswmin;
    double BSIM4v7vfbsd;

    double BSIM4v7cof1;
    double BSIM4v7cof2;
    double BSIM4v7cof3;
    double BSIM4v7cof4;
    double BSIM4v7cdep0;
    double BSIM4v7ToxRatio;
    double BSIM4v7Aechvb;
    double BSIM4v7Bechvb;
    double BSIM4v7ToxRatioEdge;
    double BSIM4v7AechvbEdgeS;
    double BSIM4v7AechvbEdgeD;
    double BSIM4v7BechvbEdge;
    double BSIM4v7ldeb;
    double BSIM4v7k1ox;
    double BSIM4v7k2ox;
    double BSIM4v7vfbzbfactor;
    double BSIM4v7dvtp2factor; /* v4.7 */
    struct bsim4SizeDependParam  *pNext;
};


typedef struct sBSIM4v7model
{
    int BSIM4v7modType;
    struct sBSIM4v7model *BSIM4v7nextModel;
    BSIM4v7instance *BSIM4v7instances;
    IFuid BSIM4v7modName;

    /* --- end of generic struct GENmodel --- */

    int BSIM4v7type;

    int    BSIM4v7mobMod;
    int    BSIM4v7cvchargeMod;
    int    BSIM4v7capMod;
    int    BSIM4v7dioMod;
    int    BSIM4v7trnqsMod;
    int    BSIM4v7acnqsMod;
    int    BSIM4v7fnoiMod;
    int    BSIM4v7tnoiMod;
    int    BSIM4v7rdsMod;
    int    BSIM4v7rbodyMod;
    int    BSIM4v7rgateMod;
    int    BSIM4v7perMod;
    int    BSIM4v7geoMod;
    int    BSIM4v7rgeoMod;
    int    BSIM4v7mtrlMod;
    int    BSIM4v7mtrlCompatMod; /* v4.7 */
    int    BSIM4v7gidlMod; /* v4.7 New GIDL/GISL */
    int    BSIM4v7igcMod;
    int    BSIM4v7igbMod;
    int    BSIM4v7tempMod;
    int    BSIM4v7binUnit;
    int    BSIM4v7paramChk;
    char   *BSIM4v7version;
    double BSIM4v7eot;
    double BSIM4v7vddeot;
  	double BSIM4v7tempeot;
  	double BSIM4v7leffeot;
  	double BSIM4v7weffeot;
    double BSIM4v7ados;
    double BSIM4v7bdos;
    double BSIM4v7toxe;
    double BSIM4v7toxp;
    double BSIM4v7toxm;
    double BSIM4v7dtox;
    double BSIM4v7epsrox;
    double BSIM4v7cdsc;
    double BSIM4v7cdscb;
    double BSIM4v7cdscd;
    double BSIM4v7cit;
    double BSIM4v7nfactor;
    double BSIM4v7xj;
    double BSIM4v7vsat;
    double BSIM4v7at;
    double BSIM4v7a0;
    double BSIM4v7ags;
    double BSIM4v7a1;
    double BSIM4v7a2;
    double BSIM4v7keta;
    double BSIM4v7nsub;
    double BSIM4v7phig;
    double BSIM4v7epsrgate;
    double BSIM4v7easub;
    double BSIM4v7epsrsub;
    double BSIM4v7ni0sub;
    double BSIM4v7bg0sub;
    double BSIM4v7tbgasub;
    double BSIM4v7tbgbsub;
    double BSIM4v7ndep;
    double BSIM4v7nsd;
    double BSIM4v7phin;
    double BSIM4v7ngate;
    double BSIM4v7gamma1;
    double BSIM4v7gamma2;
    double BSIM4v7vbx;
    double BSIM4v7vbm;
    double BSIM4v7xt;
    double BSIM4v7k1;
    double BSIM4v7kt1;
    double BSIM4v7kt1l;
    double BSIM4v7kt2;
    double BSIM4v7k2;
    double BSIM4v7k3;
    double BSIM4v7k3b;
    double BSIM4v7w0;
    double BSIM4v7dvtp0;
    double BSIM4v7dvtp1;
    double BSIM4v7dvtp2;	/* New DIBL/Rout */
    double BSIM4v7dvtp3;
    double BSIM4v7dvtp4;
    double BSIM4v7dvtp5;
    double BSIM4v7lpe0;
    double BSIM4v7lpeb;
    double BSIM4v7dvt0;
    double BSIM4v7dvt1;
    double BSIM4v7dvt2;
    double BSIM4v7dvt0w;
    double BSIM4v7dvt1w;
    double BSIM4v7dvt2w;
    double BSIM4v7drout;
    double BSIM4v7dsub;
    double BSIM4v7vth0;
    double BSIM4v7eu;
  	double BSIM4v7ucs;
    double BSIM4v7ua;
    double BSIM4v7ua1;
    double BSIM4v7ub;
    double BSIM4v7ub1;
    double BSIM4v7uc;
    double BSIM4v7uc1;
    double BSIM4v7ud;
    double BSIM4v7ud1;
    double BSIM4v7up;
    double BSIM4v7lp;
    double BSIM4v7u0;
    double BSIM4v7ute;
  	double BSIM4v7ucste;
    double BSIM4v7voff;
    double BSIM4v7tvoff;
    double BSIM4v7tnfactor; 	/* v4.7 Temp dep of leakage current */
    double BSIM4v7teta0;   	/* v4.7 temp dep of leakage current */
    double BSIM4v7tvoffcv;	/* v4.7 temp dep of leakage current */
    double BSIM4v7minv;
    double BSIM4v7minvcv;
    double BSIM4v7voffl;
    double BSIM4v7voffcvl;
    double BSIM4v7delta;
    double BSIM4v7rdsw;
    double BSIM4v7rdswmin;
    double BSIM4v7rdwmin;
    double BSIM4v7rswmin;
    double BSIM4v7rsw;
    double BSIM4v7rdw;
    double BSIM4v7prwg;
    double BSIM4v7prwb;
    double BSIM4v7prt;
    double BSIM4v7eta0;
    double BSIM4v7etab;
    double BSIM4v7pclm;
    double BSIM4v7pdibl1;
    double BSIM4v7pdibl2;
    double BSIM4v7pdiblb;
    double BSIM4v7fprout;
    double BSIM4v7pdits;
    double BSIM4v7pditsd;
    double BSIM4v7pditsl;
    double BSIM4v7pscbe1;
    double BSIM4v7pscbe2;
    double BSIM4v7pvag;
    double BSIM4v7wr;
    double BSIM4v7dwg;
    double BSIM4v7dwb;
    double BSIM4v7b0;
    double BSIM4v7b1;
    double BSIM4v7alpha0;
    double BSIM4v7alpha1;
    double BSIM4v7beta0;
    double BSIM4v7agidl;
    double BSIM4v7bgidl;
    double BSIM4v7cgidl;
    double BSIM4v7egidl;
    double BSIM4v7fgidl; /* v4.7 New GIDL/GISL */
    double BSIM4v7kgidl; /* v4.7 New GIDL/GISL */
    double BSIM4v7rgidl; /* v4.7 New GIDL/GISL */
    double BSIM4v7agisl;
    double BSIM4v7bgisl;
    double BSIM4v7cgisl;
    double BSIM4v7egisl;
    double BSIM4v7fgisl; /* v4.7 New GIDL/GISL */
    double BSIM4v7kgisl; /* v4.7 New GIDL/GISL */
    double BSIM4v7rgisl; /* v4.7 New GIDL/GISL */
    double BSIM4v7aigc;
    double BSIM4v7bigc;
    double BSIM4v7cigc;
    double BSIM4v7aigsd;
    double BSIM4v7bigsd;
    double BSIM4v7cigsd;
    double BSIM4v7aigs;
    double BSIM4v7bigs;
    double BSIM4v7cigs;
    double BSIM4v7aigd;
    double BSIM4v7bigd;
    double BSIM4v7cigd;
    double BSIM4v7aigbacc;
    double BSIM4v7bigbacc;
    double BSIM4v7cigbacc;
    double BSIM4v7aigbinv;
    double BSIM4v7bigbinv;
    double BSIM4v7cigbinv;
    double BSIM4v7nigc;
    double BSIM4v7nigbacc;
    double BSIM4v7nigbinv;
    double BSIM4v7ntox;
    double BSIM4v7eigbinv;
    double BSIM4v7pigcd;
    double BSIM4v7poxedge;
    double BSIM4v7toxref;
    double BSIM4v7ijthdfwd;
    double BSIM4v7ijthsfwd;
    double BSIM4v7ijthdrev;
    double BSIM4v7ijthsrev;
    double BSIM4v7xjbvd;
    double BSIM4v7xjbvs;
    double BSIM4v7bvd;
    double BSIM4v7bvs;

    double BSIM4v7jtss;
    double BSIM4v7jtsd;
    double BSIM4v7jtssws;
    double BSIM4v7jtsswd;
    double BSIM4v7jtsswgs;
    double BSIM4v7jtsswgd;
    double BSIM4v7jtweff;
    double BSIM4v7njts;
    double BSIM4v7njtssw;
    double BSIM4v7njtsswg;
    double BSIM4v7njtsd;
    double BSIM4v7njtsswd;
    double BSIM4v7njtsswgd;
    double BSIM4v7xtss;
    double BSIM4v7xtsd;
    double BSIM4v7xtssws;
    double BSIM4v7xtsswd;
    double BSIM4v7xtsswgs;
    double BSIM4v7xtsswgd;
    double BSIM4v7tnjts;
    double BSIM4v7tnjtssw;
    double BSIM4v7tnjtsswg;
    double BSIM4v7tnjtsd;
    double BSIM4v7tnjtsswd;
    double BSIM4v7tnjtsswgd;
    double BSIM4v7vtss;
    double BSIM4v7vtsd;
    double BSIM4v7vtssws;
    double BSIM4v7vtsswd;
    double BSIM4v7vtsswgs;
    double BSIM4v7vtsswgd;

    double BSIM4v7xrcrg1;
    double BSIM4v7xrcrg2;
    double BSIM4v7lambda;
    double BSIM4v7vtl;
    double BSIM4v7lc;
    double BSIM4v7xn;
    double BSIM4v7vfbsdoff;  /* S/D flatband offset voltage  */
    double BSIM4v7lintnoi;  /* lint offset for noise calculation  */
    double BSIM4v7tvfbsdoff;

    double BSIM4v7vfb;
    double BSIM4v7gbmin;
    double BSIM4v7rbdb;
    double BSIM4v7rbsb;
    double BSIM4v7rbpb;
    double BSIM4v7rbps;
    double BSIM4v7rbpd;

    double BSIM4v7rbps0;
    double BSIM4v7rbpsl;
    double BSIM4v7rbpsw;
    double BSIM4v7rbpsnf;

    double BSIM4v7rbpd0;
    double BSIM4v7rbpdl;
    double BSIM4v7rbpdw;
    double BSIM4v7rbpdnf;

    double BSIM4v7rbpbx0;
    double BSIM4v7rbpbxl;
    double BSIM4v7rbpbxw;
    double BSIM4v7rbpbxnf;
    double BSIM4v7rbpby0;
    double BSIM4v7rbpbyl;
    double BSIM4v7rbpbyw;
    double BSIM4v7rbpbynf;

    double BSIM4v7rbsbx0;
    double BSIM4v7rbsby0;
    double BSIM4v7rbdbx0;
    double BSIM4v7rbdby0;

    double BSIM4v7rbsdbxl;
    double BSIM4v7rbsdbxw;
    double BSIM4v7rbsdbxnf;
    double BSIM4v7rbsdbyl;
    double BSIM4v7rbsdbyw;
    double BSIM4v7rbsdbynf;

    double BSIM4v7tnoia;
    double BSIM4v7tnoib;
    double BSIM4v7tnoic;
    double BSIM4v7rnoia;
    double BSIM4v7rnoib;
    double BSIM4v7rnoic;
    double BSIM4v7ntnoi;

    /* CV model and Parasitics */
    double BSIM4v7cgsl;
    double BSIM4v7cgdl;
    double BSIM4v7ckappas;
    double BSIM4v7ckappad;
    double BSIM4v7cf;
    double BSIM4v7vfbcv;
    double BSIM4v7clc;
    double BSIM4v7cle;
    double BSIM4v7dwc;
    double BSIM4v7dlc;
    double BSIM4v7xw;
    double BSIM4v7xl;
    double BSIM4v7dlcig;
    double BSIM4v7dlcigd;
    double BSIM4v7dwj;
    double BSIM4v7noff;
    double BSIM4v7voffcv;
    double BSIM4v7acde;
    double BSIM4v7moin;
    double BSIM4v7tcj;
    double BSIM4v7tcjsw;
    double BSIM4v7tcjswg;
    double BSIM4v7tpb;
    double BSIM4v7tpbsw;
    double BSIM4v7tpbswg;
    double BSIM4v7dmcg;
    double BSIM4v7dmci;
    double BSIM4v7dmdg;
    double BSIM4v7dmcgt;
    double BSIM4v7xgw;
    double BSIM4v7xgl;
    double BSIM4v7rshg;
    double BSIM4v7ngcon;

    /* Length Dependence */
    double BSIM4v7lcdsc;
    double BSIM4v7lcdscb;
    double BSIM4v7lcdscd;
    double BSIM4v7lcit;
    double BSIM4v7lnfactor;
    double BSIM4v7lxj;
    double BSIM4v7lvsat;
    double BSIM4v7lat;
    double BSIM4v7la0;
    double BSIM4v7lags;
    double BSIM4v7la1;
    double BSIM4v7la2;
    double BSIM4v7lketa;
    double BSIM4v7lnsub;
    double BSIM4v7lndep;
    double BSIM4v7lnsd;
    double BSIM4v7lphin;
    double BSIM4v7lngate;
    double BSIM4v7lgamma1;
    double BSIM4v7lgamma2;
    double BSIM4v7lvbx;
    double BSIM4v7lvbm;
    double BSIM4v7lxt;
    double BSIM4v7lk1;
    double BSIM4v7lkt1;
    double BSIM4v7lkt1l;
    double BSIM4v7lkt2;
    double BSIM4v7lk2;
    double BSIM4v7lk3;
    double BSIM4v7lk3b;
    double BSIM4v7lw0;
    double BSIM4v7ldvtp0;
    double BSIM4v7ldvtp1;
    double BSIM4v7ldvtp2;        /* New DIBL/Rout */
    double BSIM4v7ldvtp3;
    double BSIM4v7ldvtp4;
    double BSIM4v7ldvtp5;
    double BSIM4v7llpe0;
    double BSIM4v7llpeb;
    double BSIM4v7ldvt0;
    double BSIM4v7ldvt1;
    double BSIM4v7ldvt2;
    double BSIM4v7ldvt0w;
    double BSIM4v7ldvt1w;
    double BSIM4v7ldvt2w;
    double BSIM4v7ldrout;
    double BSIM4v7ldsub;
    double BSIM4v7lvth0;
    double BSIM4v7lua;
    double BSIM4v7lua1;
    double BSIM4v7lub;
    double BSIM4v7lub1;
    double BSIM4v7luc;
    double BSIM4v7luc1;
    double BSIM4v7lud;
    double BSIM4v7lud1;
    double BSIM4v7lup;
    double BSIM4v7llp;
    double BSIM4v7lu0;
    double BSIM4v7leu;
    double BSIM4v7lucs;
    double BSIM4v7lute;
    double BSIM4v7lucste;
    double BSIM4v7lvoff;
    double BSIM4v7ltvoff;
    double BSIM4v7ltnfactor;         /* v4.7 Temp dep of leakage current */
    double BSIM4v7lteta0;           /* v4.7 temp dep of leakage current */
    double BSIM4v7ltvoffcv;        /* v4.7 temp dep of leakage current */
    double BSIM4v7lminv;
    double BSIM4v7lminvcv;
    double BSIM4v7ldelta;
    double BSIM4v7lrdsw;
    double BSIM4v7lrsw;
    double BSIM4v7lrdw;
    double BSIM4v7lprwg;
    double BSIM4v7lprwb;
    double BSIM4v7lprt;
    double BSIM4v7leta0;
    double BSIM4v7letab;
    double BSIM4v7lpclm;
    double BSIM4v7lpdibl1;
    double BSIM4v7lpdibl2;
    double BSIM4v7lpdiblb;
    double BSIM4v7lfprout;
    double BSIM4v7lpdits;
    double BSIM4v7lpditsd;
    double BSIM4v7lpscbe1;
    double BSIM4v7lpscbe2;
    double BSIM4v7lpvag;
    double BSIM4v7lwr;
    double BSIM4v7ldwg;
    double BSIM4v7ldwb;
    double BSIM4v7lb0;
    double BSIM4v7lb1;
    double BSIM4v7lalpha0;
    double BSIM4v7lalpha1;
    double BSIM4v7lbeta0;
    double BSIM4v7lvfb;
    double BSIM4v7lagidl;
    double BSIM4v7lbgidl;
    double BSIM4v7lcgidl;
    double BSIM4v7legidl;
    double BSIM4v7lfgidl; /* v4.7 New GIDL/GISL */
    double BSIM4v7lkgidl; /* v4.7 New GIDL/GISL */
    double BSIM4v7lrgidl; /* v4.7 New GIDL/GISL */
    double BSIM4v7lagisl;
    double BSIM4v7lbgisl;
    double BSIM4v7lcgisl;
    double BSIM4v7legisl;
    double BSIM4v7lfgisl; /* v4.7 New GIDL/GISL */
    double BSIM4v7lkgisl; /* v4.7 New GIDL/GISL */
    double BSIM4v7lrgisl; /* v4.7 New GIDL/GISL */
    double BSIM4v7laigc;
    double BSIM4v7lbigc;
    double BSIM4v7lcigc;
    double BSIM4v7laigsd;
    double BSIM4v7lbigsd;
    double BSIM4v7lcigsd;
    double BSIM4v7laigs;
    double BSIM4v7lbigs;
    double BSIM4v7lcigs;
    double BSIM4v7laigd;
    double BSIM4v7lbigd;
    double BSIM4v7lcigd;
    double BSIM4v7laigbacc;
    double BSIM4v7lbigbacc;
    double BSIM4v7lcigbacc;
    double BSIM4v7laigbinv;
    double BSIM4v7lbigbinv;
    double BSIM4v7lcigbinv;
    double BSIM4v7lnigc;
    double BSIM4v7lnigbacc;
    double BSIM4v7lnigbinv;
    double BSIM4v7lntox;
    double BSIM4v7leigbinv;
    double BSIM4v7lpigcd;
    double BSIM4v7lpoxedge;
    double BSIM4v7lxrcrg1;
    double BSIM4v7lxrcrg2;
    double BSIM4v7llambda;
    double BSIM4v7lvtl;
    double BSIM4v7lxn;
    double BSIM4v7lvfbsdoff;
    double BSIM4v7ltvfbsdoff;

    /* CV model */
    double BSIM4v7lcgsl;
    double BSIM4v7lcgdl;
    double BSIM4v7lckappas;
    double BSIM4v7lckappad;
    double BSIM4v7lcf;
    double BSIM4v7lclc;
    double BSIM4v7lcle;
    double BSIM4v7lvfbcv;
    double BSIM4v7lnoff;
    double BSIM4v7lvoffcv;
    double BSIM4v7lacde;
    double BSIM4v7lmoin;

    /* Width Dependence */
    double BSIM4v7wcdsc;
    double BSIM4v7wcdscb;
    double BSIM4v7wcdscd;
    double BSIM4v7wcit;
    double BSIM4v7wnfactor;
    double BSIM4v7wxj;
    double BSIM4v7wvsat;
    double BSIM4v7wat;
    double BSIM4v7wa0;
    double BSIM4v7wags;
    double BSIM4v7wa1;
    double BSIM4v7wa2;
    double BSIM4v7wketa;
    double BSIM4v7wnsub;
    double BSIM4v7wndep;
    double BSIM4v7wnsd;
    double BSIM4v7wphin;
    double BSIM4v7wngate;
    double BSIM4v7wgamma1;
    double BSIM4v7wgamma2;
    double BSIM4v7wvbx;
    double BSIM4v7wvbm;
    double BSIM4v7wxt;
    double BSIM4v7wk1;
    double BSIM4v7wkt1;
    double BSIM4v7wkt1l;
    double BSIM4v7wkt2;
    double BSIM4v7wk2;
    double BSIM4v7wk3;
    double BSIM4v7wk3b;
    double BSIM4v7ww0;
    double BSIM4v7wdvtp0;
    double BSIM4v7wdvtp1;
    double BSIM4v7wdvtp2;        /* New DIBL/Rout */
    double BSIM4v7wdvtp3;
    double BSIM4v7wdvtp4;
    double BSIM4v7wdvtp5;
    double BSIM4v7wlpe0;
    double BSIM4v7wlpeb;
    double BSIM4v7wdvt0;
    double BSIM4v7wdvt1;
    double BSIM4v7wdvt2;
    double BSIM4v7wdvt0w;
    double BSIM4v7wdvt1w;
    double BSIM4v7wdvt2w;
    double BSIM4v7wdrout;
    double BSIM4v7wdsub;
    double BSIM4v7wvth0;
    double BSIM4v7wua;
    double BSIM4v7wua1;
    double BSIM4v7wub;
    double BSIM4v7wub1;
    double BSIM4v7wuc;
    double BSIM4v7wuc1;
    double BSIM4v7wud;
    double BSIM4v7wud1;
    double BSIM4v7wup;
    double BSIM4v7wlp;
    double BSIM4v7wu0;
    double BSIM4v7weu;
    double BSIM4v7wucs;
    double BSIM4v7wute;
    double BSIM4v7wucste;
    double BSIM4v7wvoff;
    double BSIM4v7wtvoff;
    double BSIM4v7wtnfactor;         /* v4.7 Temp dep of leakage current */
    double BSIM4v7wteta0;           /* v4.7 temp dep of leakage current */
    double BSIM4v7wtvoffcv;        /* v4.7 temp dep of leakage current */
    double BSIM4v7wminv;
    double BSIM4v7wminvcv;
    double BSIM4v7wdelta;
    double BSIM4v7wrdsw;
    double BSIM4v7wrsw;
    double BSIM4v7wrdw;
    double BSIM4v7wprwg;
    double BSIM4v7wprwb;
    double BSIM4v7wprt;
    double BSIM4v7weta0;
    double BSIM4v7wetab;
    double BSIM4v7wpclm;
    double BSIM4v7wpdibl1;
    double BSIM4v7wpdibl2;
    double BSIM4v7wpdiblb;
    double BSIM4v7wfprout;
    double BSIM4v7wpdits;
    double BSIM4v7wpditsd;
    double BSIM4v7wpscbe1;
    double BSIM4v7wpscbe2;
    double BSIM4v7wpvag;
    double BSIM4v7wwr;
    double BSIM4v7wdwg;
    double BSIM4v7wdwb;
    double BSIM4v7wb0;
    double BSIM4v7wb1;
    double BSIM4v7walpha0;
    double BSIM4v7walpha1;
    double BSIM4v7wbeta0;
    double BSIM4v7wvfb;
    double BSIM4v7wagidl;
    double BSIM4v7wbgidl;
    double BSIM4v7wcgidl;
    double BSIM4v7wegidl;
    double BSIM4v7wfgidl; /* v4.7 New GIDL/GISL */
    double BSIM4v7wkgidl; /* v4.7 New GIDL/GISL */
    double BSIM4v7wrgidl; /* v4.7 New GIDL/GISL */
    double BSIM4v7wagisl;
    double BSIM4v7wbgisl;
    double BSIM4v7wcgisl;
    double BSIM4v7wegisl;
    double BSIM4v7wfgisl; /* v4.7 New GIDL/GISL */
    double BSIM4v7wkgisl; /* v4.7 New GIDL/GISL */
    double BSIM4v7wrgisl; /* v4.7 New GIDL/GISL */
    double BSIM4v7waigc;
    double BSIM4v7wbigc;
    double BSIM4v7wcigc;
    double BSIM4v7waigsd;
    double BSIM4v7wbigsd;
    double BSIM4v7wcigsd;
    double BSIM4v7waigs;
    double BSIM4v7wbigs;
    double BSIM4v7wcigs;
    double BSIM4v7waigd;
    double BSIM4v7wbigd;
    double BSIM4v7wcigd;
    double BSIM4v7waigbacc;
    double BSIM4v7wbigbacc;
    double BSIM4v7wcigbacc;
    double BSIM4v7waigbinv;
    double BSIM4v7wbigbinv;
    double BSIM4v7wcigbinv;
    double BSIM4v7wnigc;
    double BSIM4v7wnigbacc;
    double BSIM4v7wnigbinv;
    double BSIM4v7wntox;
    double BSIM4v7weigbinv;
    double BSIM4v7wpigcd;
    double BSIM4v7wpoxedge;
    double BSIM4v7wxrcrg1;
    double BSIM4v7wxrcrg2;
    double BSIM4v7wlambda;
    double BSIM4v7wvtl;
    double BSIM4v7wxn;
    double BSIM4v7wvfbsdoff;
    double BSIM4v7wtvfbsdoff;

    /* CV model */
    double BSIM4v7wcgsl;
    double BSIM4v7wcgdl;
    double BSIM4v7wckappas;
    double BSIM4v7wckappad;
    double BSIM4v7wcf;
    double BSIM4v7wclc;
    double BSIM4v7wcle;
    double BSIM4v7wvfbcv;
    double BSIM4v7wnoff;
    double BSIM4v7wvoffcv;
    double BSIM4v7wacde;
    double BSIM4v7wmoin;

    /* Cross-term Dependence */
    double BSIM4v7pcdsc;
    double BSIM4v7pcdscb;
    double BSIM4v7pcdscd;
    double BSIM4v7pcit;
    double BSIM4v7pnfactor;
    double BSIM4v7pxj;
    double BSIM4v7pvsat;
    double BSIM4v7pat;
    double BSIM4v7pa0;
    double BSIM4v7pags;
    double BSIM4v7pa1;
    double BSIM4v7pa2;
    double BSIM4v7pketa;
    double BSIM4v7pnsub;
    double BSIM4v7pndep;
    double BSIM4v7pnsd;
    double BSIM4v7pphin;
    double BSIM4v7pngate;
    double BSIM4v7pgamma1;
    double BSIM4v7pgamma2;
    double BSIM4v7pvbx;
    double BSIM4v7pvbm;
    double BSIM4v7pxt;
    double BSIM4v7pk1;
    double BSIM4v7pkt1;
    double BSIM4v7pkt1l;
    double BSIM4v7pkt2;
    double BSIM4v7pk2;
    double BSIM4v7pk3;
    double BSIM4v7pk3b;
    double BSIM4v7pw0;
    double BSIM4v7pdvtp0;
    double BSIM4v7pdvtp1;
    double BSIM4v7pdvtp2;        /* New DIBL/Rout */
    double BSIM4v7pdvtp3;
    double BSIM4v7pdvtp4;
    double BSIM4v7pdvtp5;
    double BSIM4v7plpe0;
    double BSIM4v7plpeb;
    double BSIM4v7pdvt0;
    double BSIM4v7pdvt1;
    double BSIM4v7pdvt2;
    double BSIM4v7pdvt0w;
    double BSIM4v7pdvt1w;
    double BSIM4v7pdvt2w;
    double BSIM4v7pdrout;
    double BSIM4v7pdsub;
    double BSIM4v7pvth0;
    double BSIM4v7pua;
    double BSIM4v7pua1;
    double BSIM4v7pub;
    double BSIM4v7pub1;
    double BSIM4v7puc;
    double BSIM4v7puc1;
    double BSIM4v7pud;
    double BSIM4v7pud1;
    double BSIM4v7pup;
    double BSIM4v7plp;
    double BSIM4v7pu0;
    double BSIM4v7peu;
    double BSIM4v7pucs;
    double BSIM4v7pute;
    double BSIM4v7pucste;
    double BSIM4v7pvoff;
    double BSIM4v7ptvoff;
    double BSIM4v7ptnfactor;         /* v4.7 Temp dep of leakage current */
    double BSIM4v7pteta0;           /* v4.7 temp dep of leakage current */
    double BSIM4v7ptvoffcv;        /* v4.7 temp dep of leakage current */
    double BSIM4v7pminv;
    double BSIM4v7pminvcv;
    double BSIM4v7pdelta;
    double BSIM4v7prdsw;
    double BSIM4v7prsw;
    double BSIM4v7prdw;
    double BSIM4v7pprwg;
    double BSIM4v7pprwb;
    double BSIM4v7pprt;
    double BSIM4v7peta0;
    double BSIM4v7petab;
    double BSIM4v7ppclm;
    double BSIM4v7ppdibl1;
    double BSIM4v7ppdibl2;
    double BSIM4v7ppdiblb;
    double BSIM4v7pfprout;
    double BSIM4v7ppdits;
    double BSIM4v7ppditsd;
    double BSIM4v7ppscbe1;
    double BSIM4v7ppscbe2;
    double BSIM4v7ppvag;
    double BSIM4v7pwr;
    double BSIM4v7pdwg;
    double BSIM4v7pdwb;
    double BSIM4v7pb0;
    double BSIM4v7pb1;
    double BSIM4v7palpha0;
    double BSIM4v7palpha1;
    double BSIM4v7pbeta0;
    double BSIM4v7pvfb;
    double BSIM4v7pagidl;
    double BSIM4v7pbgidl;
    double BSIM4v7pcgidl;
    double BSIM4v7pegidl;
    double BSIM4v7pfgidl; /* v4.7 New GIDL/GISL */
    double BSIM4v7pkgidl; /* v4.7 New GIDL/GISL */
    double BSIM4v7prgidl; /* v4.7 New GIDL/GISL */
    double BSIM4v7pagisl;
    double BSIM4v7pbgisl;
    double BSIM4v7pcgisl;
    double BSIM4v7pegisl;
    double BSIM4v7pfgisl; /* v4.7 New GIDL/GISL */
    double BSIM4v7pkgisl; /* v4.7 New GIDL/GISL */
    double BSIM4v7prgisl; /* v4.7 New GIDL/GISL */
    double BSIM4v7paigc;
    double BSIM4v7pbigc;
    double BSIM4v7pcigc;
    double BSIM4v7paigsd;
    double BSIM4v7pbigsd;
    double BSIM4v7pcigsd;
    double BSIM4v7paigs;
    double BSIM4v7pbigs;
    double BSIM4v7pcigs;
    double BSIM4v7paigd;
    double BSIM4v7pbigd;
    double BSIM4v7pcigd;
    double BSIM4v7paigbacc;
    double BSIM4v7pbigbacc;
    double BSIM4v7pcigbacc;
    double BSIM4v7paigbinv;
    double BSIM4v7pbigbinv;
    double BSIM4v7pcigbinv;
    double BSIM4v7pnigc;
    double BSIM4v7pnigbacc;
    double BSIM4v7pnigbinv;
    double BSIM4v7pntox;
    double BSIM4v7peigbinv;
    double BSIM4v7ppigcd;
    double BSIM4v7ppoxedge;
    double BSIM4v7pxrcrg1;
    double BSIM4v7pxrcrg2;
    double BSIM4v7plambda;
    double BSIM4v7pvtl;
    double BSIM4v7pxn;
    double BSIM4v7pvfbsdoff;
    double BSIM4v7ptvfbsdoff;

    /* CV model */
    double BSIM4v7pcgsl;
    double BSIM4v7pcgdl;
    double BSIM4v7pckappas;
    double BSIM4v7pckappad;
    double BSIM4v7pcf;
    double BSIM4v7pclc;
    double BSIM4v7pcle;
    double BSIM4v7pvfbcv;
    double BSIM4v7pnoff;
    double BSIM4v7pvoffcv;
    double BSIM4v7pacde;
    double BSIM4v7pmoin;

    double BSIM4v7tnom;
    double BSIM4v7cgso;
    double BSIM4v7cgdo;
    double BSIM4v7cgbo;
    double BSIM4v7xpart;
    double BSIM4v7cFringOut;
    double BSIM4v7cFringMax;

    double BSIM4v7sheetResistance;
    double BSIM4v7SjctSatCurDensity;
    double BSIM4v7DjctSatCurDensity;
    double BSIM4v7SjctSidewallSatCurDensity;
    double BSIM4v7DjctSidewallSatCurDensity;
    double BSIM4v7SjctGateSidewallSatCurDensity;
    double BSIM4v7DjctGateSidewallSatCurDensity;
    double BSIM4v7SbulkJctPotential;
    double BSIM4v7DbulkJctPotential;
    double BSIM4v7SbulkJctBotGradingCoeff;
    double BSIM4v7DbulkJctBotGradingCoeff;
    double BSIM4v7SbulkJctSideGradingCoeff;
    double BSIM4v7DbulkJctSideGradingCoeff;
    double BSIM4v7SbulkJctGateSideGradingCoeff;
    double BSIM4v7DbulkJctGateSideGradingCoeff;
    double BSIM4v7SsidewallJctPotential;
    double BSIM4v7DsidewallJctPotential;
    double BSIM4v7SGatesidewallJctPotential;
    double BSIM4v7DGatesidewallJctPotential;
    double BSIM4v7SunitAreaJctCap;
    double BSIM4v7DunitAreaJctCap;
    double BSIM4v7SunitLengthSidewallJctCap;
    double BSIM4v7DunitLengthSidewallJctCap;
    double BSIM4v7SunitLengthGateSidewallJctCap;
    double BSIM4v7DunitLengthGateSidewallJctCap;
    double BSIM4v7SjctEmissionCoeff;
    double BSIM4v7DjctEmissionCoeff;
    double BSIM4v7SjctTempExponent;
    double BSIM4v7DjctTempExponent;
    double BSIM4v7njtsstemp;
    double BSIM4v7njtsswstemp;
    double BSIM4v7njtsswgstemp;
    double BSIM4v7njtsdtemp;
    double BSIM4v7njtsswdtemp;
    double BSIM4v7njtsswgdtemp;

    double BSIM4v7Lint;
    double BSIM4v7Ll;
    double BSIM4v7Llc;
    double BSIM4v7Lln;
    double BSIM4v7Lw;
    double BSIM4v7Lwc;
    double BSIM4v7Lwn;
    double BSIM4v7Lwl;
    double BSIM4v7Lwlc;
    double BSIM4v7Lmin;
    double BSIM4v7Lmax;

    double BSIM4v7Wint;
    double BSIM4v7Wl;
    double BSIM4v7Wlc;
    double BSIM4v7Wln;
    double BSIM4v7Ww;
    double BSIM4v7Wwc;
    double BSIM4v7Wwn;
    double BSIM4v7Wwl;
    double BSIM4v7Wwlc;
    double BSIM4v7Wmin;
    double BSIM4v7Wmax;

    /* added for stress effect */
    double BSIM4v7saref;
    double BSIM4v7sbref;
    double BSIM4v7wlod;
    double BSIM4v7ku0;
    double BSIM4v7kvsat;
    double BSIM4v7kvth0;
    double BSIM4v7tku0;
    double BSIM4v7llodku0;
    double BSIM4v7wlodku0;
    double BSIM4v7llodvth;
    double BSIM4v7wlodvth;
    double BSIM4v7lku0;
    double BSIM4v7wku0;
    double BSIM4v7pku0;
    double BSIM4v7lkvth0;
    double BSIM4v7wkvth0;
    double BSIM4v7pkvth0;
    double BSIM4v7stk2;
    double BSIM4v7lodk2;
    double BSIM4v7steta0;
    double BSIM4v7lodeta0;

    double BSIM4v7web;
    double BSIM4v7wec;
    double BSIM4v7kvth0we;
    double BSIM4v7k2we;
    double BSIM4v7ku0we;
    double BSIM4v7scref;
    double BSIM4v7wpemod;
    double BSIM4v7lkvth0we;
    double BSIM4v7lk2we;
    double BSIM4v7lku0we;
    double BSIM4v7wkvth0we;
    double BSIM4v7wk2we;
    double BSIM4v7wku0we;
    double BSIM4v7pkvth0we;
    double BSIM4v7pk2we;
    double BSIM4v7pku0we;

/* Pre-calculated constants
 * move to size-dependent param */
    double BSIM4v7Eg0;
    double BSIM4v7vtm;
    double BSIM4v7vtm0;
    double BSIM4v7coxe;
    double BSIM4v7coxp;
    double BSIM4v7cof1;
    double BSIM4v7cof2;
    double BSIM4v7cof3;
    double BSIM4v7cof4;
    double BSIM4v7vcrit;
    double BSIM4v7factor1;
    double BSIM4v7PhiBS;
    double BSIM4v7PhiBSWS;
    double BSIM4v7PhiBSWGS;
    double BSIM4v7SjctTempSatCurDensity;
    double BSIM4v7SjctSidewallTempSatCurDensity;
    double BSIM4v7SjctGateSidewallTempSatCurDensity;
    double BSIM4v7PhiBD;
    double BSIM4v7PhiBSWD;
    double BSIM4v7PhiBSWGD;
    double BSIM4v7DjctTempSatCurDensity;
    double BSIM4v7DjctSidewallTempSatCurDensity;
    double BSIM4v7DjctGateSidewallTempSatCurDensity;
    double BSIM4v7SunitAreaTempJctCap;
    double BSIM4v7DunitAreaTempJctCap;
    double BSIM4v7SunitLengthSidewallTempJctCap;
    double BSIM4v7DunitLengthSidewallTempJctCap;
    double BSIM4v7SunitLengthGateSidewallTempJctCap;
    double BSIM4v7DunitLengthGateSidewallTempJctCap;

    double BSIM4v7oxideTrapDensityA;
    double BSIM4v7oxideTrapDensityB;
    double BSIM4v7oxideTrapDensityC;
    double BSIM4v7em;
    double BSIM4v7ef;
    double BSIM4v7af;
    double BSIM4v7kf;

    double BSIM4v7vgsMax;
    double BSIM4v7vgdMax;
    double BSIM4v7vgbMax;
    double BSIM4v7vdsMax;
    double BSIM4v7vbsMax;
    double BSIM4v7vbdMax;

    struct bsim4SizeDependParam *pSizeDependParamKnot;

#ifdef USE_OMP
    int BSIM4v7InstCount;
    struct sBSIM4v7instance **BSIM4v7InstanceArray;
#endif

    /* Flags */
    unsigned  BSIM4v7mobModGiven :1;
    unsigned  BSIM4v7binUnitGiven :1;
    unsigned  BSIM4v7cvchargeModGiven :1;
    unsigned  BSIM4v7capModGiven :1;
    unsigned  BSIM4v7dioModGiven :1;
    unsigned  BSIM4v7rdsModGiven :1;
    unsigned  BSIM4v7rbodyModGiven :1;
    unsigned  BSIM4v7rgateModGiven :1;
    unsigned  BSIM4v7perModGiven :1;
    unsigned  BSIM4v7geoModGiven :1;
    unsigned  BSIM4v7rgeoModGiven :1;
    unsigned  BSIM4v7paramChkGiven :1;
    unsigned  BSIM4v7trnqsModGiven :1;
    unsigned  BSIM4v7acnqsModGiven :1;
    unsigned  BSIM4v7fnoiModGiven :1;
    unsigned  BSIM4v7tnoiModGiven :1;
    unsigned  BSIM4v7mtrlModGiven :1;
    unsigned  BSIM4v7mtrlCompatModGiven :1;
    unsigned  BSIM4v7gidlModGiven :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7igcModGiven :1;
    unsigned  BSIM4v7igbModGiven :1;
    unsigned  BSIM4v7tempModGiven :1;
    unsigned  BSIM4v7typeGiven   :1;
    unsigned  BSIM4v7toxrefGiven   :1;
    unsigned  BSIM4v7eotGiven   :1;
    unsigned  BSIM4v7vddeotGiven   :1;
    unsigned  BSIM4v7tempeotGiven  :1;
    unsigned  BSIM4v7leffeotGiven  :1;
    unsigned  BSIM4v7weffeotGiven  :1;
    unsigned  BSIM4v7adosGiven   :1;
    unsigned  BSIM4v7bdosGiven   :1;
    unsigned  BSIM4v7toxeGiven   :1;
    unsigned  BSIM4v7toxpGiven   :1;
    unsigned  BSIM4v7toxmGiven   :1;
    unsigned  BSIM4v7dtoxGiven   :1;
    unsigned  BSIM4v7epsroxGiven   :1;
    unsigned  BSIM4v7versionGiven   :1;
    unsigned  BSIM4v7cdscGiven   :1;
    unsigned  BSIM4v7cdscbGiven   :1;
    unsigned  BSIM4v7cdscdGiven   :1;
    unsigned  BSIM4v7citGiven   :1;
    unsigned  BSIM4v7nfactorGiven   :1;
    unsigned  BSIM4v7xjGiven   :1;
    unsigned  BSIM4v7vsatGiven   :1;
    unsigned  BSIM4v7atGiven   :1;
    unsigned  BSIM4v7a0Given   :1;
    unsigned  BSIM4v7agsGiven   :1;
    unsigned  BSIM4v7a1Given   :1;
    unsigned  BSIM4v7a2Given   :1;
    unsigned  BSIM4v7ketaGiven   :1;
    unsigned  BSIM4v7nsubGiven   :1;
    unsigned  BSIM4v7phigGiven   :1;
    unsigned  BSIM4v7epsrgateGiven   :1;
    unsigned  BSIM4v7easubGiven   :1;
    unsigned  BSIM4v7epsrsubGiven   :1;
    unsigned  BSIM4v7ni0subGiven   :1;
    unsigned  BSIM4v7bg0subGiven   :1;
    unsigned  BSIM4v7tbgasubGiven   :1;
    unsigned  BSIM4v7tbgbsubGiven   :1;
    unsigned  BSIM4v7ndepGiven   :1;
    unsigned  BSIM4v7nsdGiven    :1;
    unsigned  BSIM4v7phinGiven   :1;
    unsigned  BSIM4v7ngateGiven   :1;
    unsigned  BSIM4v7gamma1Given   :1;
    unsigned  BSIM4v7gamma2Given   :1;
    unsigned  BSIM4v7vbxGiven   :1;
    unsigned  BSIM4v7vbmGiven   :1;
    unsigned  BSIM4v7xtGiven   :1;
    unsigned  BSIM4v7k1Given   :1;
    unsigned  BSIM4v7kt1Given   :1;
    unsigned  BSIM4v7kt1lGiven   :1;
    unsigned  BSIM4v7kt2Given   :1;
    unsigned  BSIM4v7k2Given   :1;
    unsigned  BSIM4v7k3Given   :1;
    unsigned  BSIM4v7k3bGiven   :1;
    unsigned  BSIM4v7w0Given   :1;
    unsigned  BSIM4v7dvtp0Given :1;
    unsigned  BSIM4v7dvtp1Given :1;
    unsigned  BSIM4v7dvtp2Given :1;        /* New DIBL/Rout */
    unsigned  BSIM4v7dvtp3Given :1;
    unsigned  BSIM4v7dvtp4Given :1;
    unsigned  BSIM4v7dvtp5Given :1;
    unsigned  BSIM4v7lpe0Given   :1;
    unsigned  BSIM4v7lpebGiven   :1;
    unsigned  BSIM4v7dvt0Given   :1;
    unsigned  BSIM4v7dvt1Given   :1;
    unsigned  BSIM4v7dvt2Given   :1;
    unsigned  BSIM4v7dvt0wGiven   :1;
    unsigned  BSIM4v7dvt1wGiven   :1;
    unsigned  BSIM4v7dvt2wGiven   :1;
    unsigned  BSIM4v7droutGiven   :1;
    unsigned  BSIM4v7dsubGiven   :1;
    unsigned  BSIM4v7vth0Given   :1;
    unsigned  BSIM4v7euGiven   :1;
    unsigned  BSIM4v7ucsGiven  :1;
    unsigned  BSIM4v7uaGiven   :1;
    unsigned  BSIM4v7ua1Given   :1;
    unsigned  BSIM4v7ubGiven   :1;
    unsigned  BSIM4v7ub1Given   :1;
    unsigned  BSIM4v7ucGiven   :1;
    unsigned  BSIM4v7uc1Given   :1;
    unsigned  BSIM4v7udGiven     :1;
    unsigned  BSIM4v7ud1Given     :1;
    unsigned  BSIM4v7upGiven     :1;
    unsigned  BSIM4v7lpGiven     :1;
    unsigned  BSIM4v7u0Given   :1;
    unsigned  BSIM4v7uteGiven   :1;
    unsigned  BSIM4v7ucsteGiven :1;
    unsigned  BSIM4v7voffGiven   :1;
    unsigned  BSIM4v7tvoffGiven   :1;
    unsigned  BSIM4v7tnfactorGiven  :1;         /* v4.7 Temp dep of leakage current */
    unsigned  BSIM4v7teta0Given   :1;           /* v4.7 temp dep of leakage current */
    unsigned  BSIM4v7tvoffcvGiven   :1;        /* v4.7 temp dep of leakage current */
    unsigned  BSIM4v7vofflGiven  :1;
    unsigned  BSIM4v7voffcvlGiven  :1;
    unsigned  BSIM4v7minvGiven   :1;
    unsigned  BSIM4v7minvcvGiven   :1;
    unsigned  BSIM4v7rdswGiven   :1;
    unsigned  BSIM4v7rdswminGiven :1;
    unsigned  BSIM4v7rdwminGiven :1;
    unsigned  BSIM4v7rswminGiven :1;
    unsigned  BSIM4v7rswGiven   :1;
    unsigned  BSIM4v7rdwGiven   :1;
    unsigned  BSIM4v7prwgGiven   :1;
    unsigned  BSIM4v7prwbGiven   :1;
    unsigned  BSIM4v7prtGiven   :1;
    unsigned  BSIM4v7eta0Given   :1;
    unsigned  BSIM4v7etabGiven   :1;
    unsigned  BSIM4v7pclmGiven   :1;
    unsigned  BSIM4v7pdibl1Given   :1;
    unsigned  BSIM4v7pdibl2Given   :1;
    unsigned  BSIM4v7pdiblbGiven   :1;
    unsigned  BSIM4v7fproutGiven   :1;
    unsigned  BSIM4v7pditsGiven    :1;
    unsigned  BSIM4v7pditsdGiven    :1;
    unsigned  BSIM4v7pditslGiven    :1;
    unsigned  BSIM4v7pscbe1Given   :1;
    unsigned  BSIM4v7pscbe2Given   :1;
    unsigned  BSIM4v7pvagGiven   :1;
    unsigned  BSIM4v7deltaGiven  :1;
    unsigned  BSIM4v7wrGiven   :1;
    unsigned  BSIM4v7dwgGiven   :1;
    unsigned  BSIM4v7dwbGiven   :1;
    unsigned  BSIM4v7b0Given   :1;
    unsigned  BSIM4v7b1Given   :1;
    unsigned  BSIM4v7alpha0Given   :1;
    unsigned  BSIM4v7alpha1Given   :1;
    unsigned  BSIM4v7beta0Given   :1;
    unsigned  BSIM4v7agidlGiven   :1;
    unsigned  BSIM4v7bgidlGiven   :1;
    unsigned  BSIM4v7cgidlGiven   :1;
    unsigned  BSIM4v7egidlGiven   :1;
    unsigned  BSIM4v7fgidlGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7kgidlGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7rgidlGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7agislGiven   :1;
    unsigned  BSIM4v7bgislGiven   :1;
    unsigned  BSIM4v7cgislGiven   :1;
    unsigned  BSIM4v7egislGiven   :1;
    unsigned  BSIM4v7fgislGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7kgislGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7rgislGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7aigcGiven   :1;
    unsigned  BSIM4v7bigcGiven   :1;
    unsigned  BSIM4v7cigcGiven   :1;
    unsigned  BSIM4v7aigsdGiven   :1;
    unsigned  BSIM4v7bigsdGiven   :1;
    unsigned  BSIM4v7cigsdGiven   :1;
    unsigned  BSIM4v7aigsGiven   :1;
    unsigned  BSIM4v7bigsGiven   :1;
    unsigned  BSIM4v7cigsGiven   :1;
    unsigned  BSIM4v7aigdGiven   :1;
    unsigned  BSIM4v7bigdGiven   :1;
    unsigned  BSIM4v7cigdGiven   :1;
    unsigned  BSIM4v7aigbaccGiven   :1;
    unsigned  BSIM4v7bigbaccGiven   :1;
    unsigned  BSIM4v7cigbaccGiven   :1;
    unsigned  BSIM4v7aigbinvGiven   :1;
    unsigned  BSIM4v7bigbinvGiven   :1;
    unsigned  BSIM4v7cigbinvGiven   :1;
    unsigned  BSIM4v7nigcGiven   :1;
    unsigned  BSIM4v7nigbinvGiven   :1;
    unsigned  BSIM4v7nigbaccGiven   :1;
    unsigned  BSIM4v7ntoxGiven   :1;
    unsigned  BSIM4v7eigbinvGiven   :1;
    unsigned  BSIM4v7pigcdGiven   :1;
    unsigned  BSIM4v7poxedgeGiven   :1;
    unsigned  BSIM4v7ijthdfwdGiven  :1;
    unsigned  BSIM4v7ijthsfwdGiven  :1;
    unsigned  BSIM4v7ijthdrevGiven  :1;
    unsigned  BSIM4v7ijthsrevGiven  :1;
    unsigned  BSIM4v7xjbvdGiven   :1;
    unsigned  BSIM4v7xjbvsGiven   :1;
    unsigned  BSIM4v7bvdGiven   :1;
    unsigned  BSIM4v7bvsGiven   :1;

    unsigned  BSIM4v7jtssGiven   :1;
    unsigned  BSIM4v7jtsdGiven   :1;
    unsigned  BSIM4v7jtsswsGiven   :1;
    unsigned  BSIM4v7jtsswdGiven   :1;
    unsigned  BSIM4v7jtsswgsGiven   :1;
    unsigned  BSIM4v7jtsswgdGiven   :1;
        unsigned  BSIM4v7jtweffGiven    :1;
    unsigned  BSIM4v7njtsGiven   :1;
    unsigned  BSIM4v7njtsswGiven   :1;
    unsigned  BSIM4v7njtsswgGiven   :1;
    unsigned  BSIM4v7njtsdGiven   :1;
    unsigned  BSIM4v7njtsswdGiven   :1;
    unsigned  BSIM4v7njtsswgdGiven   :1;
    unsigned  BSIM4v7xtssGiven   :1;
    unsigned  BSIM4v7xtsdGiven   :1;
    unsigned  BSIM4v7xtsswsGiven   :1;
    unsigned  BSIM4v7xtsswdGiven   :1;
    unsigned  BSIM4v7xtsswgsGiven   :1;
    unsigned  BSIM4v7xtsswgdGiven   :1;
    unsigned  BSIM4v7tnjtsGiven   :1;
    unsigned  BSIM4v7tnjtsswGiven   :1;
    unsigned  BSIM4v7tnjtsswgGiven   :1;
    unsigned  BSIM4v7tnjtsdGiven   :1;
    unsigned  BSIM4v7tnjtsswdGiven   :1;
    unsigned  BSIM4v7tnjtsswgdGiven   :1;
    unsigned  BSIM4v7vtssGiven   :1;
    unsigned  BSIM4v7vtsdGiven   :1;
    unsigned  BSIM4v7vtsswsGiven   :1;
    unsigned  BSIM4v7vtsswdGiven   :1;
    unsigned  BSIM4v7vtsswgsGiven   :1;
    unsigned  BSIM4v7vtsswgdGiven   :1;

    unsigned  BSIM4v7vfbGiven   :1;
    unsigned  BSIM4v7gbminGiven :1;
    unsigned  BSIM4v7rbdbGiven :1;
    unsigned  BSIM4v7rbsbGiven :1;
    unsigned  BSIM4v7rbpsGiven :1;
    unsigned  BSIM4v7rbpdGiven :1;
    unsigned  BSIM4v7rbpbGiven :1;

    unsigned BSIM4v7rbps0Given :1;
    unsigned BSIM4v7rbpslGiven :1;
    unsigned BSIM4v7rbpswGiven :1;
    unsigned BSIM4v7rbpsnfGiven :1;

    unsigned BSIM4v7rbpd0Given :1;
    unsigned BSIM4v7rbpdlGiven :1;
    unsigned BSIM4v7rbpdwGiven :1;
    unsigned BSIM4v7rbpdnfGiven :1;

    unsigned BSIM4v7rbpbx0Given :1;
    unsigned BSIM4v7rbpbxlGiven :1;
    unsigned BSIM4v7rbpbxwGiven :1;
    unsigned BSIM4v7rbpbxnfGiven :1;
    unsigned BSIM4v7rbpby0Given :1;
    unsigned BSIM4v7rbpbylGiven :1;
    unsigned BSIM4v7rbpbywGiven :1;
    unsigned BSIM4v7rbpbynfGiven :1;

    unsigned BSIM4v7rbsbx0Given :1;
    unsigned BSIM4v7rbsby0Given :1;
    unsigned BSIM4v7rbdbx0Given :1;
    unsigned BSIM4v7rbdby0Given :1;

    unsigned BSIM4v7rbsdbxlGiven :1;
    unsigned BSIM4v7rbsdbxwGiven :1;
    unsigned BSIM4v7rbsdbxnfGiven :1;
    unsigned BSIM4v7rbsdbylGiven :1;
    unsigned BSIM4v7rbsdbywGiven :1;
    unsigned BSIM4v7rbsdbynfGiven :1;

    unsigned  BSIM4v7xrcrg1Given   :1;
    unsigned  BSIM4v7xrcrg2Given   :1;
    unsigned  BSIM4v7tnoiaGiven    :1;
    unsigned  BSIM4v7tnoibGiven    :1;
    unsigned  BSIM4v7tnoicGiven    :1;
    unsigned  BSIM4v7rnoiaGiven    :1;
    unsigned  BSIM4v7rnoibGiven    :1;
    unsigned  BSIM4v7rnoicGiven    :1;
    unsigned  BSIM4v7ntnoiGiven    :1;

    unsigned  BSIM4v7lambdaGiven    :1;
    unsigned  BSIM4v7vtlGiven    :1;
    unsigned  BSIM4v7lcGiven    :1;
    unsigned  BSIM4v7xnGiven    :1;
    unsigned  BSIM4v7vfbsdoffGiven    :1;
    unsigned  BSIM4v7lintnoiGiven    :1;
    unsigned  BSIM4v7tvfbsdoffGiven    :1;

    /* CV model and parasitics */
    unsigned  BSIM4v7cgslGiven   :1;
    unsigned  BSIM4v7cgdlGiven   :1;
    unsigned  BSIM4v7ckappasGiven   :1;
    unsigned  BSIM4v7ckappadGiven   :1;
    unsigned  BSIM4v7cfGiven   :1;
    unsigned  BSIM4v7vfbcvGiven   :1;
    unsigned  BSIM4v7clcGiven   :1;
    unsigned  BSIM4v7cleGiven   :1;
    unsigned  BSIM4v7dwcGiven   :1;
    unsigned  BSIM4v7dlcGiven   :1;
    unsigned  BSIM4v7xwGiven    :1;
    unsigned  BSIM4v7xlGiven    :1;
    unsigned  BSIM4v7dlcigGiven   :1;
    unsigned  BSIM4v7dlcigdGiven   :1;
    unsigned  BSIM4v7dwjGiven   :1;
    unsigned  BSIM4v7noffGiven  :1;
    unsigned  BSIM4v7voffcvGiven :1;
    unsigned  BSIM4v7acdeGiven  :1;
    unsigned  BSIM4v7moinGiven  :1;
    unsigned  BSIM4v7tcjGiven   :1;
    unsigned  BSIM4v7tcjswGiven :1;
    unsigned  BSIM4v7tcjswgGiven :1;
    unsigned  BSIM4v7tpbGiven    :1;
    unsigned  BSIM4v7tpbswGiven  :1;
    unsigned  BSIM4v7tpbswgGiven :1;
    unsigned  BSIM4v7dmcgGiven :1;
    unsigned  BSIM4v7dmciGiven :1;
    unsigned  BSIM4v7dmdgGiven :1;
    unsigned  BSIM4v7dmcgtGiven :1;
    unsigned  BSIM4v7xgwGiven :1;
    unsigned  BSIM4v7xglGiven :1;
    unsigned  BSIM4v7rshgGiven :1;
    unsigned  BSIM4v7ngconGiven :1;


    /* Length dependence */
    unsigned  BSIM4v7lcdscGiven   :1;
    unsigned  BSIM4v7lcdscbGiven   :1;
    unsigned  BSIM4v7lcdscdGiven   :1;
    unsigned  BSIM4v7lcitGiven   :1;
    unsigned  BSIM4v7lnfactorGiven   :1;
    unsigned  BSIM4v7lxjGiven   :1;
    unsigned  BSIM4v7lvsatGiven   :1;
    unsigned  BSIM4v7latGiven   :1;
    unsigned  BSIM4v7la0Given   :1;
    unsigned  BSIM4v7lagsGiven   :1;
    unsigned  BSIM4v7la1Given   :1;
    unsigned  BSIM4v7la2Given   :1;
    unsigned  BSIM4v7lketaGiven   :1;
    unsigned  BSIM4v7lnsubGiven   :1;
    unsigned  BSIM4v7lndepGiven   :1;
    unsigned  BSIM4v7lnsdGiven    :1;
    unsigned  BSIM4v7lphinGiven   :1;
    unsigned  BSIM4v7lngateGiven   :1;
    unsigned  BSIM4v7lgamma1Given   :1;
    unsigned  BSIM4v7lgamma2Given   :1;
    unsigned  BSIM4v7lvbxGiven   :1;
    unsigned  BSIM4v7lvbmGiven   :1;
    unsigned  BSIM4v7lxtGiven   :1;
    unsigned  BSIM4v7lk1Given   :1;
    unsigned  BSIM4v7lkt1Given   :1;
    unsigned  BSIM4v7lkt1lGiven   :1;
    unsigned  BSIM4v7lkt2Given   :1;
    unsigned  BSIM4v7lk2Given   :1;
    unsigned  BSIM4v7lk3Given   :1;
    unsigned  BSIM4v7lk3bGiven   :1;
    unsigned  BSIM4v7lw0Given   :1;
    unsigned  BSIM4v7ldvtp0Given :1;
    unsigned  BSIM4v7ldvtp1Given :1;
    unsigned  BSIM4v7ldvtp2Given :1;        /* New DIBL/Rout */
    unsigned  BSIM4v7ldvtp3Given :1;
    unsigned  BSIM4v7ldvtp4Given :1;
    unsigned  BSIM4v7ldvtp5Given :1;
    unsigned  BSIM4v7llpe0Given   :1;
    unsigned  BSIM4v7llpebGiven   :1;
    unsigned  BSIM4v7ldvt0Given   :1;
    unsigned  BSIM4v7ldvt1Given   :1;
    unsigned  BSIM4v7ldvt2Given   :1;
    unsigned  BSIM4v7ldvt0wGiven   :1;
    unsigned  BSIM4v7ldvt1wGiven   :1;
    unsigned  BSIM4v7ldvt2wGiven   :1;
    unsigned  BSIM4v7ldroutGiven   :1;
    unsigned  BSIM4v7ldsubGiven   :1;
    unsigned  BSIM4v7lvth0Given   :1;
    unsigned  BSIM4v7luaGiven   :1;
    unsigned  BSIM4v7lua1Given   :1;
    unsigned  BSIM4v7lubGiven   :1;
    unsigned  BSIM4v7lub1Given   :1;
    unsigned  BSIM4v7lucGiven   :1;
    unsigned  BSIM4v7luc1Given   :1;
    unsigned  BSIM4v7ludGiven     :1;
    unsigned  BSIM4v7lud1Given     :1;
    unsigned  BSIM4v7lupGiven     :1;
    unsigned  BSIM4v7llpGiven     :1;
    unsigned  BSIM4v7lu0Given   :1;
    unsigned  BSIM4v7leuGiven   :1;
        unsigned  BSIM4v7lucsGiven   :1;
    unsigned  BSIM4v7luteGiven   :1;
        unsigned  BSIM4v7lucsteGiven  :1;
    unsigned  BSIM4v7lvoffGiven   :1;
    unsigned  BSIM4v7ltvoffGiven   :1;
    unsigned  BSIM4v7ltnfactorGiven  :1;         /* v4.7 Temp dep of leakage current */
    unsigned  BSIM4v7lteta0Given   :1;           /* v4.7 temp dep of leakage current */
    unsigned  BSIM4v7ltvoffcvGiven   :1;        /* v4.7 temp dep of leakage current */
    unsigned  BSIM4v7lminvGiven   :1;
    unsigned  BSIM4v7lminvcvGiven   :1;
    unsigned  BSIM4v7lrdswGiven   :1;
    unsigned  BSIM4v7lrswGiven   :1;
    unsigned  BSIM4v7lrdwGiven   :1;
    unsigned  BSIM4v7lprwgGiven   :1;
    unsigned  BSIM4v7lprwbGiven   :1;
    unsigned  BSIM4v7lprtGiven   :1;
    unsigned  BSIM4v7leta0Given   :1;
    unsigned  BSIM4v7letabGiven   :1;
    unsigned  BSIM4v7lpclmGiven   :1;
    unsigned  BSIM4v7lpdibl1Given   :1;
    unsigned  BSIM4v7lpdibl2Given   :1;
    unsigned  BSIM4v7lpdiblbGiven   :1;
    unsigned  BSIM4v7lfproutGiven   :1;
    unsigned  BSIM4v7lpditsGiven    :1;
    unsigned  BSIM4v7lpditsdGiven    :1;
    unsigned  BSIM4v7lpscbe1Given   :1;
    unsigned  BSIM4v7lpscbe2Given   :1;
    unsigned  BSIM4v7lpvagGiven   :1;
    unsigned  BSIM4v7ldeltaGiven  :1;
    unsigned  BSIM4v7lwrGiven   :1;
    unsigned  BSIM4v7ldwgGiven   :1;
    unsigned  BSIM4v7ldwbGiven   :1;
    unsigned  BSIM4v7lb0Given   :1;
    unsigned  BSIM4v7lb1Given   :1;
    unsigned  BSIM4v7lalpha0Given   :1;
    unsigned  BSIM4v7lalpha1Given   :1;
    unsigned  BSIM4v7lbeta0Given   :1;
    unsigned  BSIM4v7lvfbGiven   :1;
    unsigned  BSIM4v7lagidlGiven   :1;
    unsigned  BSIM4v7lbgidlGiven   :1;
    unsigned  BSIM4v7lcgidlGiven   :1;
    unsigned  BSIM4v7legidlGiven   :1;
    unsigned  BSIM4v7lfgidlGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7lkgidlGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7lrgidlGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7lagislGiven   :1;
    unsigned  BSIM4v7lbgislGiven   :1;
    unsigned  BSIM4v7lcgislGiven   :1;
    unsigned  BSIM4v7legislGiven   :1;
    unsigned  BSIM4v7lfgislGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7lkgislGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7lrgislGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7laigcGiven   :1;
    unsigned  BSIM4v7lbigcGiven   :1;
    unsigned  BSIM4v7lcigcGiven   :1;
    unsigned  BSIM4v7laigsdGiven   :1;
    unsigned  BSIM4v7lbigsdGiven   :1;
    unsigned  BSIM4v7lcigsdGiven   :1;
    unsigned  BSIM4v7laigsGiven   :1;
    unsigned  BSIM4v7lbigsGiven   :1;
    unsigned  BSIM4v7lcigsGiven   :1;
    unsigned  BSIM4v7laigdGiven   :1;
    unsigned  BSIM4v7lbigdGiven   :1;
    unsigned  BSIM4v7lcigdGiven   :1;
    unsigned  BSIM4v7laigbaccGiven   :1;
    unsigned  BSIM4v7lbigbaccGiven   :1;
    unsigned  BSIM4v7lcigbaccGiven   :1;
    unsigned  BSIM4v7laigbinvGiven   :1;
    unsigned  BSIM4v7lbigbinvGiven   :1;
    unsigned  BSIM4v7lcigbinvGiven   :1;
    unsigned  BSIM4v7lnigcGiven   :1;
    unsigned  BSIM4v7lnigbinvGiven   :1;
    unsigned  BSIM4v7lnigbaccGiven   :1;
    unsigned  BSIM4v7lntoxGiven   :1;
    unsigned  BSIM4v7leigbinvGiven   :1;
    unsigned  BSIM4v7lpigcdGiven   :1;
    unsigned  BSIM4v7lpoxedgeGiven   :1;
    unsigned  BSIM4v7lxrcrg1Given   :1;
    unsigned  BSIM4v7lxrcrg2Given   :1;
    unsigned  BSIM4v7llambdaGiven    :1;
    unsigned  BSIM4v7lvtlGiven    :1;
    unsigned  BSIM4v7lxnGiven    :1;
    unsigned  BSIM4v7lvfbsdoffGiven    :1;
    unsigned  BSIM4v7ltvfbsdoffGiven    :1;

    /* CV model */
    unsigned  BSIM4v7lcgslGiven   :1;
    unsigned  BSIM4v7lcgdlGiven   :1;
    unsigned  BSIM4v7lckappasGiven   :1;
    unsigned  BSIM4v7lckappadGiven   :1;
    unsigned  BSIM4v7lcfGiven   :1;
    unsigned  BSIM4v7lclcGiven   :1;
    unsigned  BSIM4v7lcleGiven   :1;
    unsigned  BSIM4v7lvfbcvGiven   :1;
    unsigned  BSIM4v7lnoffGiven   :1;
    unsigned  BSIM4v7lvoffcvGiven :1;
    unsigned  BSIM4v7lacdeGiven   :1;
    unsigned  BSIM4v7lmoinGiven   :1;

    /* Width dependence */
    unsigned  BSIM4v7wcdscGiven   :1;
    unsigned  BSIM4v7wcdscbGiven   :1;
    unsigned  BSIM4v7wcdscdGiven   :1;
    unsigned  BSIM4v7wcitGiven   :1;
    unsigned  BSIM4v7wnfactorGiven   :1;
    unsigned  BSIM4v7wxjGiven   :1;
    unsigned  BSIM4v7wvsatGiven   :1;
    unsigned  BSIM4v7watGiven   :1;
    unsigned  BSIM4v7wa0Given   :1;
    unsigned  BSIM4v7wagsGiven   :1;
    unsigned  BSIM4v7wa1Given   :1;
    unsigned  BSIM4v7wa2Given   :1;
    unsigned  BSIM4v7wketaGiven   :1;
    unsigned  BSIM4v7wnsubGiven   :1;
    unsigned  BSIM4v7wndepGiven   :1;
    unsigned  BSIM4v7wnsdGiven    :1;
    unsigned  BSIM4v7wphinGiven   :1;
    unsigned  BSIM4v7wngateGiven   :1;
    unsigned  BSIM4v7wgamma1Given   :1;
    unsigned  BSIM4v7wgamma2Given   :1;
    unsigned  BSIM4v7wvbxGiven   :1;
    unsigned  BSIM4v7wvbmGiven   :1;
    unsigned  BSIM4v7wxtGiven   :1;
    unsigned  BSIM4v7wk1Given   :1;
    unsigned  BSIM4v7wkt1Given   :1;
    unsigned  BSIM4v7wkt1lGiven   :1;
    unsigned  BSIM4v7wkt2Given   :1;
    unsigned  BSIM4v7wk2Given   :1;
    unsigned  BSIM4v7wk3Given   :1;
    unsigned  BSIM4v7wk3bGiven   :1;
    unsigned  BSIM4v7ww0Given   :1;
    unsigned  BSIM4v7wdvtp0Given :1;
    unsigned  BSIM4v7wdvtp1Given :1;
    unsigned  BSIM4v7wdvtp2Given :1;        /* New DIBL/Rout */
    unsigned  BSIM4v7wdvtp3Given :1;
    unsigned  BSIM4v7wdvtp4Given :1;
    unsigned  BSIM4v7wdvtp5Given :1;
    unsigned  BSIM4v7wlpe0Given   :1;
    unsigned  BSIM4v7wlpebGiven   :1;
    unsigned  BSIM4v7wdvt0Given   :1;
    unsigned  BSIM4v7wdvt1Given   :1;
    unsigned  BSIM4v7wdvt2Given   :1;
    unsigned  BSIM4v7wdvt0wGiven   :1;
    unsigned  BSIM4v7wdvt1wGiven   :1;
    unsigned  BSIM4v7wdvt2wGiven   :1;
    unsigned  BSIM4v7wdroutGiven   :1;
    unsigned  BSIM4v7wdsubGiven   :1;
    unsigned  BSIM4v7wvth0Given   :1;
    unsigned  BSIM4v7wuaGiven   :1;
    unsigned  BSIM4v7wua1Given   :1;
    unsigned  BSIM4v7wubGiven   :1;
    unsigned  BSIM4v7wub1Given   :1;
    unsigned  BSIM4v7wucGiven   :1;
    unsigned  BSIM4v7wuc1Given   :1;
    unsigned  BSIM4v7wudGiven     :1;
    unsigned  BSIM4v7wud1Given     :1;
    unsigned  BSIM4v7wupGiven     :1;
    unsigned  BSIM4v7wlpGiven     :1;
    unsigned  BSIM4v7wu0Given   :1;
    unsigned  BSIM4v7weuGiven   :1;
        unsigned  BSIM4v7wucsGiven  :1;
    unsigned  BSIM4v7wuteGiven   :1;
        unsigned  BSIM4v7wucsteGiven  :1;
    unsigned  BSIM4v7wvoffGiven   :1;
    unsigned  BSIM4v7wtvoffGiven   :1;
    unsigned  BSIM4v7wtnfactorGiven  :1;         /* v4.7 Temp dep of leakage current */
    unsigned  BSIM4v7wteta0Given   :1;           /* v4.7 temp dep of leakage current */
    unsigned  BSIM4v7wtvoffcvGiven   :1;        /* v4.7 temp dep of leakage current */
    unsigned  BSIM4v7wminvGiven   :1;
    unsigned  BSIM4v7wminvcvGiven   :1;
    unsigned  BSIM4v7wrdswGiven   :1;
    unsigned  BSIM4v7wrswGiven   :1;
    unsigned  BSIM4v7wrdwGiven   :1;
    unsigned  BSIM4v7wprwgGiven   :1;
    unsigned  BSIM4v7wprwbGiven   :1;
    unsigned  BSIM4v7wprtGiven   :1;
    unsigned  BSIM4v7weta0Given   :1;
    unsigned  BSIM4v7wetabGiven   :1;
    unsigned  BSIM4v7wpclmGiven   :1;
    unsigned  BSIM4v7wpdibl1Given   :1;
    unsigned  BSIM4v7wpdibl2Given   :1;
    unsigned  BSIM4v7wpdiblbGiven   :1;
    unsigned  BSIM4v7wfproutGiven   :1;
    unsigned  BSIM4v7wpditsGiven    :1;
    unsigned  BSIM4v7wpditsdGiven    :1;
    unsigned  BSIM4v7wpscbe1Given   :1;
    unsigned  BSIM4v7wpscbe2Given   :1;
    unsigned  BSIM4v7wpvagGiven   :1;
    unsigned  BSIM4v7wdeltaGiven  :1;
    unsigned  BSIM4v7wwrGiven   :1;
    unsigned  BSIM4v7wdwgGiven   :1;
    unsigned  BSIM4v7wdwbGiven   :1;
    unsigned  BSIM4v7wb0Given   :1;
    unsigned  BSIM4v7wb1Given   :1;
    unsigned  BSIM4v7walpha0Given   :1;
    unsigned  BSIM4v7walpha1Given   :1;
    unsigned  BSIM4v7wbeta0Given   :1;
    unsigned  BSIM4v7wvfbGiven   :1;
    unsigned  BSIM4v7wagidlGiven   :1;
    unsigned  BSIM4v7wbgidlGiven   :1;
    unsigned  BSIM4v7wcgidlGiven   :1;
    unsigned  BSIM4v7wegidlGiven   :1;
    unsigned  BSIM4v7wfgidlGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7wkgidlGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7wrgidlGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7wagislGiven   :1;
    unsigned  BSIM4v7wbgislGiven   :1;
    unsigned  BSIM4v7wcgislGiven   :1;
    unsigned  BSIM4v7wegislGiven   :1;
    unsigned  BSIM4v7wfgislGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7wkgislGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7wrgislGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7waigcGiven   :1;
    unsigned  BSIM4v7wbigcGiven   :1;
    unsigned  BSIM4v7wcigcGiven   :1;
    unsigned  BSIM4v7waigsdGiven   :1;
    unsigned  BSIM4v7wbigsdGiven   :1;
    unsigned  BSIM4v7wcigsdGiven   :1;
    unsigned  BSIM4v7waigsGiven   :1;
    unsigned  BSIM4v7wbigsGiven   :1;
    unsigned  BSIM4v7wcigsGiven   :1;
    unsigned  BSIM4v7waigdGiven   :1;
    unsigned  BSIM4v7wbigdGiven   :1;
    unsigned  BSIM4v7wcigdGiven   :1;
    unsigned  BSIM4v7waigbaccGiven   :1;
    unsigned  BSIM4v7wbigbaccGiven   :1;
    unsigned  BSIM4v7wcigbaccGiven   :1;
    unsigned  BSIM4v7waigbinvGiven   :1;
    unsigned  BSIM4v7wbigbinvGiven   :1;
    unsigned  BSIM4v7wcigbinvGiven   :1;
    unsigned  BSIM4v7wnigcGiven   :1;
    unsigned  BSIM4v7wnigbinvGiven   :1;
    unsigned  BSIM4v7wnigbaccGiven   :1;
    unsigned  BSIM4v7wntoxGiven   :1;
    unsigned  BSIM4v7weigbinvGiven   :1;
    unsigned  BSIM4v7wpigcdGiven   :1;
    unsigned  BSIM4v7wpoxedgeGiven   :1;
    unsigned  BSIM4v7wxrcrg1Given   :1;
    unsigned  BSIM4v7wxrcrg2Given   :1;
    unsigned  BSIM4v7wlambdaGiven    :1;
    unsigned  BSIM4v7wvtlGiven    :1;
    unsigned  BSIM4v7wxnGiven    :1;
    unsigned  BSIM4v7wvfbsdoffGiven    :1;
    unsigned  BSIM4v7wtvfbsdoffGiven    :1;

    /* CV model */
    unsigned  BSIM4v7wcgslGiven   :1;
    unsigned  BSIM4v7wcgdlGiven   :1;
    unsigned  BSIM4v7wckappasGiven   :1;
    unsigned  BSIM4v7wckappadGiven   :1;
    unsigned  BSIM4v7wcfGiven   :1;
    unsigned  BSIM4v7wclcGiven   :1;
    unsigned  BSIM4v7wcleGiven   :1;
    unsigned  BSIM4v7wvfbcvGiven   :1;
    unsigned  BSIM4v7wnoffGiven   :1;
    unsigned  BSIM4v7wvoffcvGiven :1;
    unsigned  BSIM4v7wacdeGiven   :1;
    unsigned  BSIM4v7wmoinGiven   :1;

    /* Cross-term dependence */
    unsigned  BSIM4v7pcdscGiven   :1;
    unsigned  BSIM4v7pcdscbGiven   :1;
    unsigned  BSIM4v7pcdscdGiven   :1;
    unsigned  BSIM4v7pcitGiven   :1;
    unsigned  BSIM4v7pnfactorGiven   :1;
    unsigned  BSIM4v7pxjGiven   :1;
    unsigned  BSIM4v7pvsatGiven   :1;
    unsigned  BSIM4v7patGiven   :1;
    unsigned  BSIM4v7pa0Given   :1;
    unsigned  BSIM4v7pagsGiven   :1;
    unsigned  BSIM4v7pa1Given   :1;
    unsigned  BSIM4v7pa2Given   :1;
    unsigned  BSIM4v7pketaGiven   :1;
    unsigned  BSIM4v7pnsubGiven   :1;
    unsigned  BSIM4v7pndepGiven   :1;
    unsigned  BSIM4v7pnsdGiven    :1;
    unsigned  BSIM4v7pphinGiven   :1;
    unsigned  BSIM4v7pngateGiven   :1;
    unsigned  BSIM4v7pgamma1Given   :1;
    unsigned  BSIM4v7pgamma2Given   :1;
    unsigned  BSIM4v7pvbxGiven   :1;
    unsigned  BSIM4v7pvbmGiven   :1;
    unsigned  BSIM4v7pxtGiven   :1;
    unsigned  BSIM4v7pk1Given   :1;
    unsigned  BSIM4v7pkt1Given   :1;
    unsigned  BSIM4v7pkt1lGiven   :1;
    unsigned  BSIM4v7pkt2Given   :1;
    unsigned  BSIM4v7pk2Given   :1;
    unsigned  BSIM4v7pk3Given   :1;
    unsigned  BSIM4v7pk3bGiven   :1;
    unsigned  BSIM4v7pw0Given   :1;
    unsigned  BSIM4v7pdvtp0Given :1;
    unsigned  BSIM4v7pdvtp1Given :1;
    unsigned  BSIM4v7pdvtp2Given :1;        /* New DIBL/Rout */
    unsigned  BSIM4v7pdvtp3Given :1;
    unsigned  BSIM4v7pdvtp4Given :1;
    unsigned  BSIM4v7pdvtp5Given :1;
    unsigned  BSIM4v7plpe0Given   :1;
    unsigned  BSIM4v7plpebGiven   :1;
    unsigned  BSIM4v7pdvt0Given   :1;
    unsigned  BSIM4v7pdvt1Given   :1;
    unsigned  BSIM4v7pdvt2Given   :1;
    unsigned  BSIM4v7pdvt0wGiven   :1;
    unsigned  BSIM4v7pdvt1wGiven   :1;
    unsigned  BSIM4v7pdvt2wGiven   :1;
    unsigned  BSIM4v7pdroutGiven   :1;
    unsigned  BSIM4v7pdsubGiven   :1;
    unsigned  BSIM4v7pvth0Given   :1;
    unsigned  BSIM4v7puaGiven   :1;
    unsigned  BSIM4v7pua1Given   :1;
    unsigned  BSIM4v7pubGiven   :1;
    unsigned  BSIM4v7pub1Given   :1;
    unsigned  BSIM4v7pucGiven   :1;
    unsigned  BSIM4v7puc1Given   :1;
    unsigned  BSIM4v7pudGiven     :1;
    unsigned  BSIM4v7pud1Given     :1;
    unsigned  BSIM4v7pupGiven     :1;
    unsigned  BSIM4v7plpGiven     :1;
    unsigned  BSIM4v7pu0Given   :1;
    unsigned  BSIM4v7peuGiven   :1;
        unsigned  BSIM4v7pucsGiven   :1;
    unsigned  BSIM4v7puteGiven   :1;
        unsigned  BSIM4v7pucsteGiven  :1;
    unsigned  BSIM4v7pvoffGiven   :1;
    unsigned  BSIM4v7ptvoffGiven   :1;
    unsigned  BSIM4v7ptnfactorGiven  :1;         /* v4.7 Temp dep of leakage current */
    unsigned  BSIM4v7pteta0Given   :1;           /* v4.7 temp dep of leakage current */
    unsigned  BSIM4v7ptvoffcvGiven   :1;        /* v4.7 temp dep of leakage current */
    unsigned  BSIM4v7pminvGiven   :1;
    unsigned  BSIM4v7pminvcvGiven   :1;
    unsigned  BSIM4v7prdswGiven   :1;
    unsigned  BSIM4v7prswGiven   :1;
    unsigned  BSIM4v7prdwGiven   :1;
    unsigned  BSIM4v7pprwgGiven   :1;
    unsigned  BSIM4v7pprwbGiven   :1;
    unsigned  BSIM4v7pprtGiven   :1;
    unsigned  BSIM4v7peta0Given   :1;
    unsigned  BSIM4v7petabGiven   :1;
    unsigned  BSIM4v7ppclmGiven   :1;
    unsigned  BSIM4v7ppdibl1Given   :1;
    unsigned  BSIM4v7ppdibl2Given   :1;
    unsigned  BSIM4v7ppdiblbGiven   :1;
    unsigned  BSIM4v7pfproutGiven   :1;
    unsigned  BSIM4v7ppditsGiven    :1;
    unsigned  BSIM4v7ppditsdGiven    :1;
    unsigned  BSIM4v7ppscbe1Given   :1;
    unsigned  BSIM4v7ppscbe2Given   :1;
    unsigned  BSIM4v7ppvagGiven   :1;
    unsigned  BSIM4v7pdeltaGiven  :1;
    unsigned  BSIM4v7pwrGiven   :1;
    unsigned  BSIM4v7pdwgGiven   :1;
    unsigned  BSIM4v7pdwbGiven   :1;
    unsigned  BSIM4v7pb0Given   :1;
    unsigned  BSIM4v7pb1Given   :1;
    unsigned  BSIM4v7palpha0Given   :1;
    unsigned  BSIM4v7palpha1Given   :1;
    unsigned  BSIM4v7pbeta0Given   :1;
    unsigned  BSIM4v7pvfbGiven   :1;
    unsigned  BSIM4v7pagidlGiven   :1;
    unsigned  BSIM4v7pbgidlGiven   :1;
    unsigned  BSIM4v7pcgidlGiven   :1;
    unsigned  BSIM4v7pegidlGiven   :1;
    unsigned  BSIM4v7pfgidlGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7pkgidlGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7prgidlGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7pagislGiven   :1;
    unsigned  BSIM4v7pbgislGiven   :1;
    unsigned  BSIM4v7pcgislGiven   :1;
    unsigned  BSIM4v7pegislGiven   :1;
    unsigned  BSIM4v7pfgislGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7pkgislGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7prgislGiven   :1;    /* v4.7 New GIDL/GISL */
    unsigned  BSIM4v7paigcGiven   :1;
    unsigned  BSIM4v7pbigcGiven   :1;
    unsigned  BSIM4v7pcigcGiven   :1;
    unsigned  BSIM4v7paigsdGiven   :1;
    unsigned  BSIM4v7pbigsdGiven   :1;
    unsigned  BSIM4v7pcigsdGiven   :1;
    unsigned  BSIM4v7paigsGiven   :1;
    unsigned  BSIM4v7pbigsGiven   :1;
    unsigned  BSIM4v7pcigsGiven   :1;
    unsigned  BSIM4v7paigdGiven   :1;
    unsigned  BSIM4v7pbigdGiven   :1;
    unsigned  BSIM4v7pcigdGiven   :1;
    unsigned  BSIM4v7paigbaccGiven   :1;
    unsigned  BSIM4v7pbigbaccGiven   :1;
    unsigned  BSIM4v7pcigbaccGiven   :1;
    unsigned  BSIM4v7paigbinvGiven   :1;
    unsigned  BSIM4v7pbigbinvGiven   :1;
    unsigned  BSIM4v7pcigbinvGiven   :1;
    unsigned  BSIM4v7pnigcGiven   :1;
    unsigned  BSIM4v7pnigbinvGiven   :1;
    unsigned  BSIM4v7pnigbaccGiven   :1;
    unsigned  BSIM4v7pntoxGiven   :1;
    unsigned  BSIM4v7peigbinvGiven   :1;
    unsigned  BSIM4v7ppigcdGiven   :1;
    unsigned  BSIM4v7ppoxedgeGiven   :1;
    unsigned  BSIM4v7pxrcrg1Given   :1;
    unsigned  BSIM4v7pxrcrg2Given   :1;
    unsigned  BSIM4v7plambdaGiven    :1;
    unsigned  BSIM4v7pvtlGiven    :1;
    unsigned  BSIM4v7pxnGiven    :1;
    unsigned  BSIM4v7pvfbsdoffGiven    :1;
    unsigned  BSIM4v7ptvfbsdoffGiven    :1;

    /* CV model */
    unsigned  BSIM4v7pcgslGiven   :1;
    unsigned  BSIM4v7pcgdlGiven   :1;
    unsigned  BSIM4v7pckappasGiven   :1;
    unsigned  BSIM4v7pckappadGiven   :1;
    unsigned  BSIM4v7pcfGiven   :1;
    unsigned  BSIM4v7pclcGiven   :1;
    unsigned  BSIM4v7pcleGiven   :1;
    unsigned  BSIM4v7pvfbcvGiven   :1;
    unsigned  BSIM4v7pnoffGiven   :1;
    unsigned  BSIM4v7pvoffcvGiven :1;
    unsigned  BSIM4v7pacdeGiven   :1;
    unsigned  BSIM4v7pmoinGiven   :1;

    unsigned  BSIM4v7useFringeGiven   :1;

    unsigned  BSIM4v7tnomGiven   :1;
    unsigned  BSIM4v7cgsoGiven   :1;
    unsigned  BSIM4v7cgdoGiven   :1;
    unsigned  BSIM4v7cgboGiven   :1;
    unsigned  BSIM4v7xpartGiven   :1;
    unsigned  BSIM4v7sheetResistanceGiven   :1;

    unsigned  BSIM4v7SjctSatCurDensityGiven   :1;
    unsigned  BSIM4v7SjctSidewallSatCurDensityGiven   :1;
    unsigned  BSIM4v7SjctGateSidewallSatCurDensityGiven   :1;
    unsigned  BSIM4v7SbulkJctPotentialGiven   :1;
    unsigned  BSIM4v7SbulkJctBotGradingCoeffGiven   :1;
    unsigned  BSIM4v7SsidewallJctPotentialGiven   :1;
    unsigned  BSIM4v7SGatesidewallJctPotentialGiven   :1;
    unsigned  BSIM4v7SbulkJctSideGradingCoeffGiven   :1;
    unsigned  BSIM4v7SunitAreaJctCapGiven   :1;
    unsigned  BSIM4v7SunitLengthSidewallJctCapGiven   :1;
    unsigned  BSIM4v7SbulkJctGateSideGradingCoeffGiven   :1;
    unsigned  BSIM4v7SunitLengthGateSidewallJctCapGiven   :1;
    unsigned  BSIM4v7SjctEmissionCoeffGiven :1;
    unsigned  BSIM4v7SjctTempExponentGiven        :1;

    unsigned  BSIM4v7DjctSatCurDensityGiven   :1;
    unsigned  BSIM4v7DjctSidewallSatCurDensityGiven   :1;
    unsigned  BSIM4v7DjctGateSidewallSatCurDensityGiven   :1;
    unsigned  BSIM4v7DbulkJctPotentialGiven   :1;
    unsigned  BSIM4v7DbulkJctBotGradingCoeffGiven   :1;
    unsigned  BSIM4v7DsidewallJctPotentialGiven   :1;
    unsigned  BSIM4v7DGatesidewallJctPotentialGiven   :1;
    unsigned  BSIM4v7DbulkJctSideGradingCoeffGiven   :1;
    unsigned  BSIM4v7DunitAreaJctCapGiven   :1;
    unsigned  BSIM4v7DunitLengthSidewallJctCapGiven   :1;
    unsigned  BSIM4v7DbulkJctGateSideGradingCoeffGiven   :1;
    unsigned  BSIM4v7DunitLengthGateSidewallJctCapGiven   :1;
    unsigned  BSIM4v7DjctEmissionCoeffGiven :1;
    unsigned  BSIM4v7DjctTempExponentGiven :1;

    unsigned  BSIM4v7oxideTrapDensityAGiven  :1;
    unsigned  BSIM4v7oxideTrapDensityBGiven  :1;
    unsigned  BSIM4v7oxideTrapDensityCGiven  :1;
    unsigned  BSIM4v7emGiven  :1;
    unsigned  BSIM4v7efGiven  :1;
    unsigned  BSIM4v7afGiven  :1;
    unsigned  BSIM4v7kfGiven  :1;

    unsigned  BSIM4v7vgsMaxGiven  :1;
    unsigned  BSIM4v7vgdMaxGiven  :1;
    unsigned  BSIM4v7vgbMaxGiven  :1;
    unsigned  BSIM4v7vdsMaxGiven  :1;
    unsigned  BSIM4v7vbsMaxGiven  :1;
    unsigned  BSIM4v7vbdMaxGiven  :1;

    unsigned  BSIM4v7LintGiven   :1;
    unsigned  BSIM4v7LlGiven   :1;
    unsigned  BSIM4v7LlcGiven   :1;
    unsigned  BSIM4v7LlnGiven   :1;
    unsigned  BSIM4v7LwGiven   :1;
    unsigned  BSIM4v7LwcGiven   :1;
    unsigned  BSIM4v7LwnGiven   :1;
    unsigned  BSIM4v7LwlGiven   :1;
    unsigned  BSIM4v7LwlcGiven   :1;
    unsigned  BSIM4v7LminGiven   :1;
    unsigned  BSIM4v7LmaxGiven   :1;

    unsigned  BSIM4v7WintGiven   :1;
    unsigned  BSIM4v7WlGiven   :1;
    unsigned  BSIM4v7WlcGiven   :1;
    unsigned  BSIM4v7WlnGiven   :1;
    unsigned  BSIM4v7WwGiven   :1;
    unsigned  BSIM4v7WwcGiven   :1;
    unsigned  BSIM4v7WwnGiven   :1;
    unsigned  BSIM4v7WwlGiven   :1;
    unsigned  BSIM4v7WwlcGiven   :1;
    unsigned  BSIM4v7WminGiven   :1;
    unsigned  BSIM4v7WmaxGiven   :1;

    /* added for stress effect */
    unsigned  BSIM4v7sarefGiven   :1;
    unsigned  BSIM4v7sbrefGiven   :1;
    unsigned  BSIM4v7wlodGiven  :1;
    unsigned  BSIM4v7ku0Given   :1;
    unsigned  BSIM4v7kvsatGiven  :1;
    unsigned  BSIM4v7kvth0Given  :1;
    unsigned  BSIM4v7tku0Given   :1;
    unsigned  BSIM4v7llodku0Given   :1;
    unsigned  BSIM4v7wlodku0Given   :1;
    unsigned  BSIM4v7llodvthGiven   :1;
    unsigned  BSIM4v7wlodvthGiven   :1;
    unsigned  BSIM4v7lku0Given   :1;
    unsigned  BSIM4v7wku0Given   :1;
    unsigned  BSIM4v7pku0Given   :1;
    unsigned  BSIM4v7lkvth0Given   :1;
    unsigned  BSIM4v7wkvth0Given   :1;
    unsigned  BSIM4v7pkvth0Given   :1;
    unsigned  BSIM4v7stk2Given   :1;
    unsigned  BSIM4v7lodk2Given  :1;
    unsigned  BSIM4v7steta0Given :1;
    unsigned  BSIM4v7lodeta0Given :1;

    unsigned  BSIM4v7webGiven   :1;
    unsigned  BSIM4v7wecGiven   :1;
    unsigned  BSIM4v7kvth0weGiven   :1;
    unsigned  BSIM4v7k2weGiven   :1;
    unsigned  BSIM4v7ku0weGiven   :1;
    unsigned  BSIM4v7screfGiven   :1;
    unsigned  BSIM4v7wpemodGiven   :1;
    unsigned  BSIM4v7lkvth0weGiven   :1;
    unsigned  BSIM4v7lk2weGiven   :1;
    unsigned  BSIM4v7lku0weGiven   :1;
    unsigned  BSIM4v7wkvth0weGiven   :1;
    unsigned  BSIM4v7wk2weGiven   :1;
    unsigned  BSIM4v7wku0weGiven   :1;
    unsigned  BSIM4v7pkvth0weGiven   :1;
    unsigned  BSIM4v7pk2weGiven   :1;
    unsigned  BSIM4v7pku0weGiven   :1;


} BSIM4v7model;


#ifndef NMOS
#define NMOS 1
#define PMOS -1
#endif /*NMOS*/


/* Instance parameters */
#define BSIM4v7_W                   1
#define BSIM4v7_L                   2
#define BSIM4v7_AS                  3
#define BSIM4v7_AD                  4
#define BSIM4v7_PS                  5
#define BSIM4v7_PD                  6
#define BSIM4v7_NRS                 7
#define BSIM4v7_NRD                 8
#define BSIM4v7_OFF                 9
#define BSIM4v7_IC                  10
#define BSIM4v7_IC_VDS              11
#define BSIM4v7_IC_VGS              12
#define BSIM4v7_IC_VBS              13
#define BSIM4v7_TRNQSMOD            14
#define BSIM4v7_RBODYMOD            15
#define BSIM4v7_RGATEMOD            16
#define BSIM4v7_GEOMOD              17
#define BSIM4v7_RGEOMOD             18
#define BSIM4v7_NF                  19
#define BSIM4v7_MIN                 20
#define BSIM4v7_ACNQSMOD            22
#define BSIM4v7_RBDB                23
#define BSIM4v7_RBSB                24
#define BSIM4v7_RBPB                25
#define BSIM4v7_RBPS                26
#define BSIM4v7_RBPD                27
#define BSIM4v7_SA                  28
#define BSIM4v7_SB                  29
#define BSIM4v7_SD                  30
#define BSIM4v7_DELVTO              31
#define BSIM4v7_XGW                 32
#define BSIM4v7_NGCON               33
#define BSIM4v7_SCA                 34
#define BSIM4v7_SCB                 35
#define BSIM4v7_SCC                 36
#define BSIM4v7_SC                  37
#define BSIM4v7_M                   38

/* Global parameters */
#define BSIM4v7_MOD_TEMPEOT         65
#define BSIM4v7_MOD_LEFFEOT         66
#define BSIM4v7_MOD_WEFFEOT         67
#define BSIM4v7_MOD_UCSTE           68
#define BSIM4v7_MOD_LUCSTE          69
#define BSIM4v7_MOD_WUCSTE          70
#define BSIM4v7_MOD_PUCSTE          71
#define BSIM4v7_MOD_UCS             72
#define BSIM4v7_MOD_LUCS            73
#define BSIM4v7_MOD_WUCS            74
#define BSIM4v7_MOD_PUCS            75
#define BSIM4v7_MOD_CVCHARGEMOD     76
#define BSIM4v7_MOD_ADOS            77
#define BSIM4v7_MOD_BDOS            78
#define BSIM4v7_MOD_TEMPMOD         79
#define BSIM4v7_MOD_MTRLMOD         80
#define BSIM4v7_MOD_IGCMOD          81
#define BSIM4v7_MOD_IGBMOD          82
#define BSIM4v7_MOD_ACNQSMOD        83
#define BSIM4v7_MOD_FNOIMOD         84
#define BSIM4v7_MOD_RDSMOD          85
#define BSIM4v7_MOD_DIOMOD          86
#define BSIM4v7_MOD_PERMOD          87
#define BSIM4v7_MOD_GEOMOD          88
#define BSIM4v7_MOD_RGEOMOD         89
#define BSIM4v7_MOD_RGATEMOD        90
#define BSIM4v7_MOD_RBODYMOD        91
#define BSIM4v7_MOD_CAPMOD          92
#define BSIM4v7_MOD_TRNQSMOD        93
#define BSIM4v7_MOD_MOBMOD          94
#define BSIM4v7_MOD_TNOIMOD         95
#define BSIM4v7_MOD_EOT             96
#define BSIM4v7_MOD_VDDEOT          97
#define BSIM4v7_MOD_TOXE            98
#define BSIM4v7_MOD_CDSC            99
#define BSIM4v7_MOD_CDSCB           100
#define BSIM4v7_MOD_CIT             101
#define BSIM4v7_MOD_NFACTOR         102
#define BSIM4v7_MOD_XJ              103
#define BSIM4v7_MOD_VSAT            104
#define BSIM4v7_MOD_AT              105
#define BSIM4v7_MOD_A0              106
#define BSIM4v7_MOD_A1              107
#define BSIM4v7_MOD_A2              108
#define BSIM4v7_MOD_KETA            109
#define BSIM4v7_MOD_NSUB            110
#define BSIM4v7_MOD_PHIG            111
#define BSIM4v7_MOD_EPSRGATE        112
#define BSIM4v7_MOD_EASUB           113
#define BSIM4v7_MOD_EPSRSUB         114
#define BSIM4v7_MOD_NI0SUB          115
#define BSIM4v7_MOD_BG0SUB          116
#define BSIM4v7_MOD_TBGASUB         117
#define BSIM4v7_MOD_TBGBSUB         118
#define BSIM4v7_MOD_NDEP            119
#define BSIM4v7_MOD_NGATE           120
#define BSIM4v7_MOD_GAMMA1          121
#define BSIM4v7_MOD_GAMMA2          122
#define BSIM4v7_MOD_VBX             123
#define BSIM4v7_MOD_BINUNIT         124
#define BSIM4v7_MOD_VBM             125
#define BSIM4v7_MOD_XT              126
#define BSIM4v7_MOD_K1              129
#define BSIM4v7_MOD_KT1             130
#define BSIM4v7_MOD_KT1L            131
#define BSIM4v7_MOD_K2              132
#define BSIM4v7_MOD_KT2             133
#define BSIM4v7_MOD_K3              134
#define BSIM4v7_MOD_K3B             135
#define BSIM4v7_MOD_W0              136
#define BSIM4v7_MOD_LPE0            137
#define BSIM4v7_MOD_DVT0            138
#define BSIM4v7_MOD_DVT1            139
#define BSIM4v7_MOD_DVT2            140
#define BSIM4v7_MOD_DVT0W           141
#define BSIM4v7_MOD_DVT1W           142
#define BSIM4v7_MOD_DVT2W           143
#define BSIM4v7_MOD_DROUT           144
#define BSIM4v7_MOD_DSUB            145
#define BSIM4v7_MOD_VTH0            146
#define BSIM4v7_MOD_UA              147
#define BSIM4v7_MOD_UA1             148
#define BSIM4v7_MOD_UB              149
#define BSIM4v7_MOD_UB1             150
#define BSIM4v7_MOD_UC              151
#define BSIM4v7_MOD_UC1             152
#define BSIM4v7_MOD_U0              153
#define BSIM4v7_MOD_UTE             154
#define BSIM4v7_MOD_VOFF            155
#define BSIM4v7_MOD_DELTA           156
#define BSIM4v7_MOD_RDSW            157
#define BSIM4v7_MOD_PRT             158
#define BSIM4v7_MOD_LDD             159
#define BSIM4v7_MOD_ETA             160
#define BSIM4v7_MOD_ETA0            161
#define BSIM4v7_MOD_ETAB            162
#define BSIM4v7_MOD_PCLM            163
#define BSIM4v7_MOD_PDIBL1          164
#define BSIM4v7_MOD_PDIBL2          165
#define BSIM4v7_MOD_PSCBE1          166
#define BSIM4v7_MOD_PSCBE2          167
#define BSIM4v7_MOD_PVAG            168
#define BSIM4v7_MOD_WR              169
#define BSIM4v7_MOD_DWG             170
#define BSIM4v7_MOD_DWB             171
#define BSIM4v7_MOD_B0              172
#define BSIM4v7_MOD_B1              173
#define BSIM4v7_MOD_ALPHA0          174
#define BSIM4v7_MOD_BETA0           175
#define BSIM4v7_MOD_PDIBLB          178
#define BSIM4v7_MOD_PRWG            179
#define BSIM4v7_MOD_PRWB            180
#define BSIM4v7_MOD_CDSCD           181
#define BSIM4v7_MOD_AGS             182
#define BSIM4v7_MOD_FRINGE          184
#define BSIM4v7_MOD_CGSL            186
#define BSIM4v7_MOD_CGDL            187
#define BSIM4v7_MOD_CKAPPAS         188
#define BSIM4v7_MOD_CF              189
#define BSIM4v7_MOD_CLC             190
#define BSIM4v7_MOD_CLE             191
#define BSIM4v7_MOD_PARAMCHK        192
#define BSIM4v7_MOD_VERSION         193
#define BSIM4v7_MOD_VFBCV           194
#define BSIM4v7_MOD_ACDE            195
#define BSIM4v7_MOD_MOIN            196
#define BSIM4v7_MOD_NOFF            197
#define BSIM4v7_MOD_IJTHDFWD        198
#define BSIM4v7_MOD_ALPHA1          199
#define BSIM4v7_MOD_VFB             200
#define BSIM4v7_MOD_TOXM            201
#define BSIM4v7_MOD_TCJ             202
#define BSIM4v7_MOD_TCJSW           203
#define BSIM4v7_MOD_TCJSWG          204
#define BSIM4v7_MOD_TPB             205
#define BSIM4v7_MOD_TPBSW           206
#define BSIM4v7_MOD_TPBSWG          207
#define BSIM4v7_MOD_VOFFCV          208
#define BSIM4v7_MOD_GBMIN           209
#define BSIM4v7_MOD_RBDB            210
#define BSIM4v7_MOD_RBSB            211
#define BSIM4v7_MOD_RBPB            212
#define BSIM4v7_MOD_RBPS            213
#define BSIM4v7_MOD_RBPD            214
#define BSIM4v7_MOD_DMCG            215
#define BSIM4v7_MOD_DMCI            216
#define BSIM4v7_MOD_DMDG            217
#define BSIM4v7_MOD_XGW             218
#define BSIM4v7_MOD_XGL             219
#define BSIM4v7_MOD_RSHG            220
#define BSIM4v7_MOD_NGCON           221
#define BSIM4v7_MOD_AGIDL           222
#define BSIM4v7_MOD_BGIDL           223
#define BSIM4v7_MOD_EGIDL           224
#define BSIM4v7_MOD_IJTHSFWD        225
#define BSIM4v7_MOD_XJBVD           226
#define BSIM4v7_MOD_XJBVS           227
#define BSIM4v7_MOD_BVD             228
#define BSIM4v7_MOD_BVS             229
#define BSIM4v7_MOD_TOXP            230
#define BSIM4v7_MOD_DTOX            231
#define BSIM4v7_MOD_XRCRG1          232
#define BSIM4v7_MOD_XRCRG2          233
#define BSIM4v7_MOD_EU              234
#define BSIM4v7_MOD_IJTHSREV        235
#define BSIM4v7_MOD_IJTHDREV        236
#define BSIM4v7_MOD_MINV            237
#define BSIM4v7_MOD_VOFFL           238
#define BSIM4v7_MOD_PDITS           239
#define BSIM4v7_MOD_PDITSD          240
#define BSIM4v7_MOD_PDITSL          241
#define BSIM4v7_MOD_TNOIA           242
#define BSIM4v7_MOD_TNOIB           243
#define BSIM4v7_MOD_NTNOI           244
#define BSIM4v7_MOD_FPROUT          245
#define BSIM4v7_MOD_LPEB            246
#define BSIM4v7_MOD_DVTP0           247
#define BSIM4v7_MOD_DVTP1           248
#define BSIM4v7_MOD_CGIDL           249
#define BSIM4v7_MOD_PHIN            250
#define BSIM4v7_MOD_RDSWMIN         251
#define BSIM4v7_MOD_RSW             252
#define BSIM4v7_MOD_RDW             253
#define BSIM4v7_MOD_RDWMIN          254
#define BSIM4v7_MOD_RSWMIN          255
#define BSIM4v7_MOD_NSD             256
#define BSIM4v7_MOD_CKAPPAD         257
#define BSIM4v7_MOD_DMCGT           258
#define BSIM4v7_MOD_AIGC            259
#define BSIM4v7_MOD_BIGC            260
#define BSIM4v7_MOD_CIGC            261
#define BSIM4v7_MOD_AIGBACC         262
#define BSIM4v7_MOD_BIGBACC         263
#define BSIM4v7_MOD_CIGBACC         264
#define BSIM4v7_MOD_AIGBINV         265
#define BSIM4v7_MOD_BIGBINV         266
#define BSIM4v7_MOD_CIGBINV         267
#define BSIM4v7_MOD_NIGC            268
#define BSIM4v7_MOD_NIGBACC         269
#define BSIM4v7_MOD_NIGBINV         270
#define BSIM4v7_MOD_NTOX            271
#define BSIM4v7_MOD_TOXREF          272
#define BSIM4v7_MOD_EIGBINV         273
#define BSIM4v7_MOD_PIGCD           274
#define BSIM4v7_MOD_POXEDGE         275
#define BSIM4v7_MOD_EPSROX          276
#define BSIM4v7_MOD_AIGSD           277
#define BSIM4v7_MOD_BIGSD           278
#define BSIM4v7_MOD_CIGSD           279
#define BSIM4v7_MOD_JSWGS           280
#define BSIM4v7_MOD_JSWGD           281
#define BSIM4v7_MOD_LAMBDA          282
#define BSIM4v7_MOD_VTL             283
#define BSIM4v7_MOD_LC              284
#define BSIM4v7_MOD_XN              285
#define BSIM4v7_MOD_RNOIA           286
#define BSIM4v7_MOD_RNOIB           287
#define BSIM4v7_MOD_VFBSDOFF        288
#define BSIM4v7_MOD_LINTNOI         289
#define BSIM4v7_MOD_UD              290
#define BSIM4v7_MOD_UD1             291
#define BSIM4v7_MOD_UP              292
#define BSIM4v7_MOD_LP              293
#define BSIM4v7_MOD_TVOFF           294
#define BSIM4v7_MOD_TVFBSDOFF       295
#define BSIM4v7_MOD_MINVCV          296
#define BSIM4v7_MOD_VOFFCVL         297
#define BSIM4v7_MOD_MTRLCOMPATMOD   380

/* Length dependence */
#define BSIM4v7_MOD_LCDSC            301
#define BSIM4v7_MOD_LCDSCB           302
#define BSIM4v7_MOD_LCIT             303
#define BSIM4v7_MOD_LNFACTOR         304
#define BSIM4v7_MOD_LXJ              305
#define BSIM4v7_MOD_LVSAT            306
#define BSIM4v7_MOD_LAT              307
#define BSIM4v7_MOD_LA0              308
#define BSIM4v7_MOD_LA1              309
#define BSIM4v7_MOD_LA2              310
#define BSIM4v7_MOD_LKETA            311
#define BSIM4v7_MOD_LNSUB            312
#define BSIM4v7_MOD_LNDEP            313
#define BSIM4v7_MOD_LNGATE           315
#define BSIM4v7_MOD_LGAMMA1          316
#define BSIM4v7_MOD_LGAMMA2          317
#define BSIM4v7_MOD_LVBX             318
#define BSIM4v7_MOD_LVBM             320
#define BSIM4v7_MOD_LXT              322
#define BSIM4v7_MOD_LK1              325
#define BSIM4v7_MOD_LKT1             326
#define BSIM4v7_MOD_LKT1L            327
#define BSIM4v7_MOD_LK2              328
#define BSIM4v7_MOD_LKT2             329
#define BSIM4v7_MOD_LK3              330
#define BSIM4v7_MOD_LK3B             331
#define BSIM4v7_MOD_LW0              332
#define BSIM4v7_MOD_LLPE0            333
#define BSIM4v7_MOD_LDVT0            334
#define BSIM4v7_MOD_LDVT1            335
#define BSIM4v7_MOD_LDVT2            336
#define BSIM4v7_MOD_LDVT0W           337
#define BSIM4v7_MOD_LDVT1W           338
#define BSIM4v7_MOD_LDVT2W           339
#define BSIM4v7_MOD_LDROUT           340
#define BSIM4v7_MOD_LDSUB            341
#define BSIM4v7_MOD_LVTH0            342
#define BSIM4v7_MOD_LUA              343
#define BSIM4v7_MOD_LUA1             344
#define BSIM4v7_MOD_LUB              345
#define BSIM4v7_MOD_LUB1             346
#define BSIM4v7_MOD_LUC              347
#define BSIM4v7_MOD_LUC1             348
#define BSIM4v7_MOD_LU0              349
#define BSIM4v7_MOD_LUTE             350
#define BSIM4v7_MOD_LVOFF            351
#define BSIM4v7_MOD_LDELTA           352
#define BSIM4v7_MOD_LRDSW            353
#define BSIM4v7_MOD_LPRT             354
#define BSIM4v7_MOD_LLDD             355
#define BSIM4v7_MOD_LETA             356
#define BSIM4v7_MOD_LETA0            357
#define BSIM4v7_MOD_LETAB            358
#define BSIM4v7_MOD_LPCLM            359
#define BSIM4v7_MOD_LPDIBL1          360
#define BSIM4v7_MOD_LPDIBL2          361
#define BSIM4v7_MOD_LPSCBE1          362
#define BSIM4v7_MOD_LPSCBE2          363
#define BSIM4v7_MOD_LPVAG            364
#define BSIM4v7_MOD_LWR              365
#define BSIM4v7_MOD_LDWG             366
#define BSIM4v7_MOD_LDWB             367
#define BSIM4v7_MOD_LB0              368
#define BSIM4v7_MOD_LB1              369
#define BSIM4v7_MOD_LALPHA0          370
#define BSIM4v7_MOD_LBETA0           371
#define BSIM4v7_MOD_LPDIBLB          374
#define BSIM4v7_MOD_LPRWG            375
#define BSIM4v7_MOD_LPRWB            376
#define BSIM4v7_MOD_LCDSCD           377
#define BSIM4v7_MOD_LAGS             378

#define BSIM4v7_MOD_LFRINGE          381
#define BSIM4v7_MOD_LCGSL            383
#define BSIM4v7_MOD_LCGDL            384
#define BSIM4v7_MOD_LCKAPPAS         385
#define BSIM4v7_MOD_LCF              386
#define BSIM4v7_MOD_LCLC             387
#define BSIM4v7_MOD_LCLE             388
#define BSIM4v7_MOD_LVFBCV           389
#define BSIM4v7_MOD_LACDE            390
#define BSIM4v7_MOD_LMOIN            391
#define BSIM4v7_MOD_LNOFF            392
#define BSIM4v7_MOD_LALPHA1          394
#define BSIM4v7_MOD_LVFB             395
#define BSIM4v7_MOD_LVOFFCV          396
#define BSIM4v7_MOD_LAGIDL           397
#define BSIM4v7_MOD_LBGIDL           398
#define BSIM4v7_MOD_LEGIDL           399
#define BSIM4v7_MOD_LXRCRG1          400
#define BSIM4v7_MOD_LXRCRG2          401
#define BSIM4v7_MOD_LEU              402
#define BSIM4v7_MOD_LMINV            403
#define BSIM4v7_MOD_LPDITS           404
#define BSIM4v7_MOD_LPDITSD          405
#define BSIM4v7_MOD_LFPROUT          406
#define BSIM4v7_MOD_LLPEB            407
#define BSIM4v7_MOD_LDVTP0           408
#define BSIM4v7_MOD_LDVTP1           409
#define BSIM4v7_MOD_LCGIDL           410
#define BSIM4v7_MOD_LPHIN            411
#define BSIM4v7_MOD_LRSW             412
#define BSIM4v7_MOD_LRDW             413
#define BSIM4v7_MOD_LNSD             414
#define BSIM4v7_MOD_LCKAPPAD         415
#define BSIM4v7_MOD_LAIGC            416
#define BSIM4v7_MOD_LBIGC            417
#define BSIM4v7_MOD_LCIGC            418
#define BSIM4v7_MOD_LAIGBACC         419
#define BSIM4v7_MOD_LBIGBACC         420
#define BSIM4v7_MOD_LCIGBACC         421
#define BSIM4v7_MOD_LAIGBINV         422
#define BSIM4v7_MOD_LBIGBINV         423
#define BSIM4v7_MOD_LCIGBINV         424
#define BSIM4v7_MOD_LNIGC            425
#define BSIM4v7_MOD_LNIGBACC         426
#define BSIM4v7_MOD_LNIGBINV         427
#define BSIM4v7_MOD_LNTOX            428
#define BSIM4v7_MOD_LEIGBINV         429
#define BSIM4v7_MOD_LPIGCD           430
#define BSIM4v7_MOD_LPOXEDGE         431
#define BSIM4v7_MOD_LAIGSD           432
#define BSIM4v7_MOD_LBIGSD           433
#define BSIM4v7_MOD_LCIGSD           434

#define BSIM4v7_MOD_LLAMBDA          435
#define BSIM4v7_MOD_LVTL             436
#define BSIM4v7_MOD_LXN              437
#define BSIM4v7_MOD_LVFBSDOFF        438
#define BSIM4v7_MOD_LUD              439
#define BSIM4v7_MOD_LUD1             440
#define BSIM4v7_MOD_LUP              441
#define BSIM4v7_MOD_LLP              442
#define BSIM4v7_MOD_LMINVCV          443

#define BSIM4v7_MOD_FGIDL            444                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_KGIDL            445                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_RGIDL            446                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_FGISL            447                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_KGISL            448                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_RGISL            449                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_LFGIDL           450                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_LKGIDL           451                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_LRGIDL           452                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_LFGISL           453                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_LKGISL           454                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_LRGISL           455                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_WFGIDL           456                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_WKGIDL           457                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_WRGIDL           458                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_WFGISL           459                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_WKGISL           460                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_WRGISL           461                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_PFGIDL           462                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_PKGIDL           463                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_PRGIDL           464                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_PFGISL           465                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_PKGISL           466                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_PRGISL           467                  /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_GIDLMOD          379                        /* v4.7 New GIDL/GISL*/
#define BSIM4v7_MOD_DVTP2           468                        /* v4.7 NEW DIBL/Rout*/
#define BSIM4v7_MOD_DVTP3           469                        /* v4.7 NEW DIBL/Rout*/
#define BSIM4v7_MOD_DVTP4           470                        /* v4.7 NEW DIBL/Rout*/
#define BSIM4v7_MOD_DVTP5           471                        /* v4.7 NEW DIBL/Rout*/
#define BSIM4v7_MOD_LDVTP2          472                        /* v4.7 NEW DIBL/Rout*/
#define BSIM4v7_MOD_LDVTP3          473                        /* v4.7 NEW DIBL/Rout*/
#define BSIM4v7_MOD_LDVTP4          474                        /* v4.7 NEW DIBL/Rout*/
#define BSIM4v7_MOD_LDVTP5          475                        /* v4.7 NEW DIBL/Rout*/
#define BSIM4v7_MOD_WDVTP2          476                        /* v4.7 NEW DIBL/Rout*/
#define BSIM4v7_MOD_WDVTP3          477                        /* v4.7 NEW DIBL/Rout*/
#define BSIM4v7_MOD_WDVTP4          478                        /* v4.7 NEW DIBL/Rout*/
#define BSIM4v7_MOD_WDVTP5          479                        /* v4.7 NEW DIBL/Rout*/
#define BSIM4v7_MOD_PDVTP2          480                        /* v4.7 NEW DIBL/Rout*/
#define BSIM4v7_MOD_PDVTP3          298                         /* v4.7 NEW DIBL/Rout*/
#define BSIM4v7_MOD_PDVTP4          299                         /* v4.7 NEW DIBL/Rout*/
#define BSIM4v7_MOD_PDVTP5          300                         /* v4.7 NEW DIBL/Rout*/

/* Width dependence */
#define BSIM4v7_MOD_WCDSC            481
#define BSIM4v7_MOD_WCDSCB           482
#define BSIM4v7_MOD_WCIT             483
#define BSIM4v7_MOD_WNFACTOR         484
#define BSIM4v7_MOD_WXJ              485
#define BSIM4v7_MOD_WVSAT            486
#define BSIM4v7_MOD_WAT              487
#define BSIM4v7_MOD_WA0              488
#define BSIM4v7_MOD_WA1              489
#define BSIM4v7_MOD_WA2              490
#define BSIM4v7_MOD_WKETA            491
#define BSIM4v7_MOD_WNSUB            492
#define BSIM4v7_MOD_WNDEP            493
#define BSIM4v7_MOD_WNGATE           495
#define BSIM4v7_MOD_WGAMMA1          496
#define BSIM4v7_MOD_WGAMMA2          497
#define BSIM4v7_MOD_WVBX             498
#define BSIM4v7_MOD_WVBM             500
#define BSIM4v7_MOD_WXT              502
#define BSIM4v7_MOD_WK1              505
#define BSIM4v7_MOD_WKT1             506
#define BSIM4v7_MOD_WKT1L            507
#define BSIM4v7_MOD_WK2              508
#define BSIM4v7_MOD_WKT2             509
#define BSIM4v7_MOD_WK3              510
#define BSIM4v7_MOD_WK3B             511
#define BSIM4v7_MOD_WW0              512
#define BSIM4v7_MOD_WLPE0            513
#define BSIM4v7_MOD_WDVT0            514
#define BSIM4v7_MOD_WDVT1            515
#define BSIM4v7_MOD_WDVT2            516
#define BSIM4v7_MOD_WDVT0W           517
#define BSIM4v7_MOD_WDVT1W           518
#define BSIM4v7_MOD_WDVT2W           519
#define BSIM4v7_MOD_WDROUT           520
#define BSIM4v7_MOD_WDSUB            521
#define BSIM4v7_MOD_WVTH0            522
#define BSIM4v7_MOD_WUA              523
#define BSIM4v7_MOD_WUA1             524
#define BSIM4v7_MOD_WUB              525
#define BSIM4v7_MOD_WUB1             526
#define BSIM4v7_MOD_WUC              527
#define BSIM4v7_MOD_WUC1             528
#define BSIM4v7_MOD_WU0              529
#define BSIM4v7_MOD_WUTE             530
#define BSIM4v7_MOD_WVOFF            531
#define BSIM4v7_MOD_WDELTA           532
#define BSIM4v7_MOD_WRDSW            533
#define BSIM4v7_MOD_WPRT             534
#define BSIM4v7_MOD_WLDD             535
#define BSIM4v7_MOD_WETA             536
#define BSIM4v7_MOD_WETA0            537
#define BSIM4v7_MOD_WETAB            538
#define BSIM4v7_MOD_WPCLM            539
#define BSIM4v7_MOD_WPDIBL1          540
#define BSIM4v7_MOD_WPDIBL2          541
#define BSIM4v7_MOD_WPSCBE1          542
#define BSIM4v7_MOD_WPSCBE2          543
#define BSIM4v7_MOD_WPVAG            544
#define BSIM4v7_MOD_WWR              545
#define BSIM4v7_MOD_WDWG             546
#define BSIM4v7_MOD_WDWB             547
#define BSIM4v7_MOD_WB0              548
#define BSIM4v7_MOD_WB1              549
#define BSIM4v7_MOD_WALPHA0          550
#define BSIM4v7_MOD_WBETA0           551
#define BSIM4v7_MOD_WPDIBLB          554
#define BSIM4v7_MOD_WPRWG            555
#define BSIM4v7_MOD_WPRWB            556
#define BSIM4v7_MOD_WCDSCD           557
#define BSIM4v7_MOD_WAGS             558

#define BSIM4v7_MOD_WFRINGE          561
#define BSIM4v7_MOD_WCGSL            563
#define BSIM4v7_MOD_WCGDL            564
#define BSIM4v7_MOD_WCKAPPAS         565
#define BSIM4v7_MOD_WCF              566
#define BSIM4v7_MOD_WCLC             567
#define BSIM4v7_MOD_WCLE             568
#define BSIM4v7_MOD_WVFBCV           569
#define BSIM4v7_MOD_WACDE            570
#define BSIM4v7_MOD_WMOIN            571
#define BSIM4v7_MOD_WNOFF            572
#define BSIM4v7_MOD_WALPHA1          574
#define BSIM4v7_MOD_WVFB             575
#define BSIM4v7_MOD_WVOFFCV          576
#define BSIM4v7_MOD_WAGIDL           577
#define BSIM4v7_MOD_WBGIDL           578
#define BSIM4v7_MOD_WEGIDL           579
#define BSIM4v7_MOD_WXRCRG1          580
#define BSIM4v7_MOD_WXRCRG2          581
#define BSIM4v7_MOD_WEU              582
#define BSIM4v7_MOD_WMINV            583
#define BSIM4v7_MOD_WPDITS           584
#define BSIM4v7_MOD_WPDITSD          585
#define BSIM4v7_MOD_WFPROUT          586
#define BSIM4v7_MOD_WLPEB            587
#define BSIM4v7_MOD_WDVTP0           588
#define BSIM4v7_MOD_WDVTP1           589
#define BSIM4v7_MOD_WCGIDL           590
#define BSIM4v7_MOD_WPHIN            591
#define BSIM4v7_MOD_WRSW             592
#define BSIM4v7_MOD_WRDW             593
#define BSIM4v7_MOD_WNSD             594
#define BSIM4v7_MOD_WCKAPPAD         595
#define BSIM4v7_MOD_WAIGC            596
#define BSIM4v7_MOD_WBIGC            597
#define BSIM4v7_MOD_WCIGC            598
#define BSIM4v7_MOD_WAIGBACC         599
#define BSIM4v7_MOD_WBIGBACC         600
#define BSIM4v7_MOD_WCIGBACC         601
#define BSIM4v7_MOD_WAIGBINV         602
#define BSIM4v7_MOD_WBIGBINV         603
#define BSIM4v7_MOD_WCIGBINV         604
#define BSIM4v7_MOD_WNIGC            605
#define BSIM4v7_MOD_WNIGBACC         606
#define BSIM4v7_MOD_WNIGBINV         607
#define BSIM4v7_MOD_WNTOX            608
#define BSIM4v7_MOD_WEIGBINV         609
#define BSIM4v7_MOD_WPIGCD           610
#define BSIM4v7_MOD_WPOXEDGE         611
#define BSIM4v7_MOD_WAIGSD           612
#define BSIM4v7_MOD_WBIGSD           613
#define BSIM4v7_MOD_WCIGSD           614
#define BSIM4v7_MOD_WLAMBDA          615
#define BSIM4v7_MOD_WVTL             616
#define BSIM4v7_MOD_WXN              617
#define BSIM4v7_MOD_WVFBSDOFF        618
#define BSIM4v7_MOD_WUD              619
#define BSIM4v7_MOD_WUD1             620
#define BSIM4v7_MOD_WUP              621
#define BSIM4v7_MOD_WLP              622
#define BSIM4v7_MOD_WMINVCV          623

/* Cross-term dependence */
#define BSIM4v7_MOD_PCDSC            661
#define BSIM4v7_MOD_PCDSCB           662
#define BSIM4v7_MOD_PCIT             663
#define BSIM4v7_MOD_PNFACTOR         664
#define BSIM4v7_MOD_PXJ              665
#define BSIM4v7_MOD_PVSAT            666
#define BSIM4v7_MOD_PAT              667
#define BSIM4v7_MOD_PA0              668
#define BSIM4v7_MOD_PA1              669
#define BSIM4v7_MOD_PA2              670
#define BSIM4v7_MOD_PKETA            671
#define BSIM4v7_MOD_PNSUB            672
#define BSIM4v7_MOD_PNDEP            673
#define BSIM4v7_MOD_PNGATE           675
#define BSIM4v7_MOD_PGAMMA1          676
#define BSIM4v7_MOD_PGAMMA2          677
#define BSIM4v7_MOD_PVBX             678

#define BSIM4v7_MOD_PVBM             680

#define BSIM4v7_MOD_PXT              682
#define BSIM4v7_MOD_PK1              685
#define BSIM4v7_MOD_PKT1             686
#define BSIM4v7_MOD_PKT1L            687
#define BSIM4v7_MOD_PK2              688
#define BSIM4v7_MOD_PKT2             689
#define BSIM4v7_MOD_PK3              690
#define BSIM4v7_MOD_PK3B             691
#define BSIM4v7_MOD_PW0              692
#define BSIM4v7_MOD_PLPE0            693

#define BSIM4v7_MOD_PDVT0            694
#define BSIM4v7_MOD_PDVT1            695
#define BSIM4v7_MOD_PDVT2            696

#define BSIM4v7_MOD_PDVT0W           697
#define BSIM4v7_MOD_PDVT1W           698
#define BSIM4v7_MOD_PDVT2W           699

#define BSIM4v7_MOD_PDROUT           700
#define BSIM4v7_MOD_PDSUB            701
#define BSIM4v7_MOD_PVTH0            702
#define BSIM4v7_MOD_PUA              703
#define BSIM4v7_MOD_PUA1             704
#define BSIM4v7_MOD_PUB              705
#define BSIM4v7_MOD_PUB1             706
#define BSIM4v7_MOD_PUC              707
#define BSIM4v7_MOD_PUC1             708
#define BSIM4v7_MOD_PU0              709
#define BSIM4v7_MOD_PUTE             710
#define BSIM4v7_MOD_PVOFF            711
#define BSIM4v7_MOD_PDELTA           712
#define BSIM4v7_MOD_PRDSW            713
#define BSIM4v7_MOD_PPRT             714
#define BSIM4v7_MOD_PLDD             715
#define BSIM4v7_MOD_PETA             716
#define BSIM4v7_MOD_PETA0            717
#define BSIM4v7_MOD_PETAB            718
#define BSIM4v7_MOD_PPCLM            719
#define BSIM4v7_MOD_PPDIBL1          720
#define BSIM4v7_MOD_PPDIBL2          721
#define BSIM4v7_MOD_PPSCBE1          722
#define BSIM4v7_MOD_PPSCBE2          723
#define BSIM4v7_MOD_PPVAG            724
#define BSIM4v7_MOD_PWR              725
#define BSIM4v7_MOD_PDWG             726
#define BSIM4v7_MOD_PDWB             727
#define BSIM4v7_MOD_PB0              728
#define BSIM4v7_MOD_PB1              729
#define BSIM4v7_MOD_PALPHA0          730
#define BSIM4v7_MOD_PBETA0           731
#define BSIM4v7_MOD_PPDIBLB          734

#define BSIM4v7_MOD_PPRWG            735
#define BSIM4v7_MOD_PPRWB            736

#define BSIM4v7_MOD_PCDSCD           737
#define BSIM4v7_MOD_PAGS             738

#define BSIM4v7_MOD_PFRINGE          741
#define BSIM4v7_MOD_PCGSL            743
#define BSIM4v7_MOD_PCGDL            744
#define BSIM4v7_MOD_PCKAPPAS         745
#define BSIM4v7_MOD_PCF              746
#define BSIM4v7_MOD_PCLC             747
#define BSIM4v7_MOD_PCLE             748
#define BSIM4v7_MOD_PVFBCV           749
#define BSIM4v7_MOD_PACDE            750
#define BSIM4v7_MOD_PMOIN            751
#define BSIM4v7_MOD_PNOFF            752
#define BSIM4v7_MOD_PALPHA1          754
#define BSIM4v7_MOD_PVFB             755
#define BSIM4v7_MOD_PVOFFCV          756
#define BSIM4v7_MOD_PAGIDL           757
#define BSIM4v7_MOD_PBGIDL           758
#define BSIM4v7_MOD_PEGIDL           759
#define BSIM4v7_MOD_PXRCRG1          760
#define BSIM4v7_MOD_PXRCRG2          761
#define BSIM4v7_MOD_PEU              762
#define BSIM4v7_MOD_PMINV            763
#define BSIM4v7_MOD_PPDITS           764
#define BSIM4v7_MOD_PPDITSD          765
#define BSIM4v7_MOD_PFPROUT          766
#define BSIM4v7_MOD_PLPEB            767
#define BSIM4v7_MOD_PDVTP0           768
#define BSIM4v7_MOD_PDVTP1           769
#define BSIM4v7_MOD_PCGIDL           770
#define BSIM4v7_MOD_PPHIN            771
#define BSIM4v7_MOD_PRSW             772
#define BSIM4v7_MOD_PRDW             773
#define BSIM4v7_MOD_PNSD             774
#define BSIM4v7_MOD_PCKAPPAD         775
#define BSIM4v7_MOD_PAIGC            776
#define BSIM4v7_MOD_PBIGC            777
#define BSIM4v7_MOD_PCIGC            778
#define BSIM4v7_MOD_PAIGBACC         779
#define BSIM4v7_MOD_PBIGBACC         780
#define BSIM4v7_MOD_PCIGBACC         781
#define BSIM4v7_MOD_PAIGBINV         782
#define BSIM4v7_MOD_PBIGBINV         783
#define BSIM4v7_MOD_PCIGBINV         784
#define BSIM4v7_MOD_PNIGC            785
#define BSIM4v7_MOD_PNIGBACC         786
#define BSIM4v7_MOD_PNIGBINV         787
#define BSIM4v7_MOD_PNTOX            788
#define BSIM4v7_MOD_PEIGBINV         789
#define BSIM4v7_MOD_PPIGCD           790
#define BSIM4v7_MOD_PPOXEDGE         791
#define BSIM4v7_MOD_PAIGSD           792
#define BSIM4v7_MOD_PBIGSD           793
#define BSIM4v7_MOD_PCIGSD           794

#define BSIM4v7_MOD_SAREF            795
#define BSIM4v7_MOD_SBREF            796
#define BSIM4v7_MOD_KU0              797
#define BSIM4v7_MOD_KVSAT            798
#define BSIM4v7_MOD_TKU0             799
#define BSIM4v7_MOD_LLODKU0          800
#define BSIM4v7_MOD_WLODKU0          801
#define BSIM4v7_MOD_LLODVTH          802
#define BSIM4v7_MOD_WLODVTH          803
#define BSIM4v7_MOD_LKU0             804
#define BSIM4v7_MOD_WKU0             805
#define BSIM4v7_MOD_PKU0             806
#define BSIM4v7_MOD_KVTH0            807
#define BSIM4v7_MOD_LKVTH0           808
#define BSIM4v7_MOD_WKVTH0           809
#define BSIM4v7_MOD_PKVTH0           810
#define BSIM4v7_MOD_WLOD                   811
#define BSIM4v7_MOD_STK2                   812
#define BSIM4v7_MOD_LODK2                   813
#define BSIM4v7_MOD_STETA0           814
#define BSIM4v7_MOD_LODETA0           815

#define BSIM4v7_MOD_WEB              816
#define BSIM4v7_MOD_WEC              817
#define BSIM4v7_MOD_KVTH0WE          818
#define BSIM4v7_MOD_K2WE             819
#define BSIM4v7_MOD_KU0WE            820
#define BSIM4v7_MOD_SCREF            821
#define BSIM4v7_MOD_WPEMOD           822
#define BSIM4v7_MOD_PMINVCV          823

#define BSIM4v7_MOD_PLAMBDA          825
#define BSIM4v7_MOD_PVTL             826
#define BSIM4v7_MOD_PXN              827
#define BSIM4v7_MOD_PVFBSDOFF        828

#define BSIM4v7_MOD_TNOM             831
#define BSIM4v7_MOD_CGSO             832
#define BSIM4v7_MOD_CGDO             833
#define BSIM4v7_MOD_CGBO             834
#define BSIM4v7_MOD_XPART            835
#define BSIM4v7_MOD_RSH              836
#define BSIM4v7_MOD_JSS              837
#define BSIM4v7_MOD_PBS              838
#define BSIM4v7_MOD_MJS              839
#define BSIM4v7_MOD_PBSWS            840
#define BSIM4v7_MOD_MJSWS            841
#define BSIM4v7_MOD_CJS              842
#define BSIM4v7_MOD_CJSWS            843
#define BSIM4v7_MOD_NMOS             844
#define BSIM4v7_MOD_PMOS             845
#define BSIM4v7_MOD_NOIA             846
#define BSIM4v7_MOD_NOIB             847
#define BSIM4v7_MOD_NOIC             848
#define BSIM4v7_MOD_LINT             849
#define BSIM4v7_MOD_LL               850
#define BSIM4v7_MOD_LLN              851
#define BSIM4v7_MOD_LW               852
#define BSIM4v7_MOD_LWN              853
#define BSIM4v7_MOD_LWL              854
#define BSIM4v7_MOD_LMIN             855
#define BSIM4v7_MOD_LMAX             856
#define BSIM4v7_MOD_WINT             857
#define BSIM4v7_MOD_WL               858
#define BSIM4v7_MOD_WLN              859
#define BSIM4v7_MOD_WW               860
#define BSIM4v7_MOD_WWN              861
#define BSIM4v7_MOD_WWL              862
#define BSIM4v7_MOD_WMIN             863
#define BSIM4v7_MOD_WMAX             864
#define BSIM4v7_MOD_DWC              865
#define BSIM4v7_MOD_DLC              866
#define BSIM4v7_MOD_XL               867
#define BSIM4v7_MOD_XW               868
#define BSIM4v7_MOD_EM               869
#define BSIM4v7_MOD_EF               870
#define BSIM4v7_MOD_AF               871
#define BSIM4v7_MOD_KF               872
#define BSIM4v7_MOD_NJS              873
#define BSIM4v7_MOD_XTIS             874
#define BSIM4v7_MOD_PBSWGS           875
#define BSIM4v7_MOD_MJSWGS           876
#define BSIM4v7_MOD_CJSWGS           877
#define BSIM4v7_MOD_JSWS             878
#define BSIM4v7_MOD_LLC              879
#define BSIM4v7_MOD_LWC              880
#define BSIM4v7_MOD_LWLC             881
#define BSIM4v7_MOD_WLC              882
#define BSIM4v7_MOD_WWC              883
#define BSIM4v7_MOD_WWLC             884
#define BSIM4v7_MOD_DWJ              885
#define BSIM4v7_MOD_JSD              886
#define BSIM4v7_MOD_PBD              887
#define BSIM4v7_MOD_MJD              888
#define BSIM4v7_MOD_PBSWD            889
#define BSIM4v7_MOD_MJSWD            890
#define BSIM4v7_MOD_CJD              891
#define BSIM4v7_MOD_CJSWD            892
#define BSIM4v7_MOD_NJD              893
#define BSIM4v7_MOD_XTID             894
#define BSIM4v7_MOD_PBSWGD           895
#define BSIM4v7_MOD_MJSWGD           896
#define BSIM4v7_MOD_CJSWGD           897
#define BSIM4v7_MOD_JSWD             898
#define BSIM4v7_MOD_DLCIG            899

/* trap-assisted tunneling */

#define BSIM4v7_MOD_JTSS             900
#define BSIM4v7_MOD_JTSD                   901
#define BSIM4v7_MOD_JTSSWS           902
#define BSIM4v7_MOD_JTSSWD           903
#define BSIM4v7_MOD_JTSSWGS           904
#define BSIM4v7_MOD_JTSSWGD           905
#define BSIM4v7_MOD_NJTS                    906
#define BSIM4v7_MOD_NJTSSW           907
#define BSIM4v7_MOD_NJTSSWG           908
#define BSIM4v7_MOD_XTSS                   909
#define BSIM4v7_MOD_XTSD                   910
#define BSIM4v7_MOD_XTSSWS           911
#define BSIM4v7_MOD_XTSSWD           912
#define BSIM4v7_MOD_XTSSWGS           913
#define BSIM4v7_MOD_XTSSWGD           914
#define BSIM4v7_MOD_TNJTS                   915
#define BSIM4v7_MOD_TNJTSSW           916
#define BSIM4v7_MOD_TNJTSSWG           917
#define BSIM4v7_MOD_VTSS             918
#define BSIM4v7_MOD_VTSD                   919
#define BSIM4v7_MOD_VTSSWS           920
#define BSIM4v7_MOD_VTSSWD           921
#define BSIM4v7_MOD_VTSSWGS           922
#define BSIM4v7_MOD_VTSSWGD           923
#define BSIM4v7_MOD_PUD              924
#define BSIM4v7_MOD_PUD1             925
#define BSIM4v7_MOD_PUP              926
#define BSIM4v7_MOD_PLP              927
#define BSIM4v7_MOD_JTWEFF           928

/* device questions */
#define BSIM4v7_DNODE                945
#define BSIM4v7_GNODEEXT             946
#define BSIM4v7_SNODE                947
#define BSIM4v7_BNODE                948
#define BSIM4v7_DNODEPRIME           949
#define BSIM4v7_GNODEPRIME           950
#define BSIM4v7_GNODEMIDE            951
#define BSIM4v7_GNODEMID             952
#define BSIM4v7_SNODEPRIME           953
#define BSIM4v7_BNODEPRIME           954
#define BSIM4v7_DBNODE               955
#define BSIM4v7_SBNODE               956
#define BSIM4v7_VBD                  957
#define BSIM4v7_VBS                  958
#define BSIM4v7_VGS                  959
#define BSIM4v7_VDS                  960
#define BSIM4v7_CD                   961
#define BSIM4v7_CBS                  962
#define BSIM4v7_CBD                  963
#define BSIM4v7_GM                   964
#define BSIM4v7_GDS                  965
#define BSIM4v7_GMBS                 966
#define BSIM4v7_GBD                  967
#define BSIM4v7_GBS                  968
#define BSIM4v7_QB                   969
#define BSIM4v7_CQB                  970
#define BSIM4v7_QG                   971
#define BSIM4v7_CQG                  972
#define BSIM4v7_QD                   973
#define BSIM4v7_CQD                  974
#define BSIM4v7_CGGB                 975
#define BSIM4v7_CGDB                 976
#define BSIM4v7_CGSB                 977
#define BSIM4v7_CBGB                 978
#define BSIM4v7_CAPBD                979
#define BSIM4v7_CQBD                 980
#define BSIM4v7_CAPBS                981
#define BSIM4v7_CQBS                 982
#define BSIM4v7_CDGB                 983
#define BSIM4v7_CDDB                 984
#define BSIM4v7_CDSB                 985
#define BSIM4v7_VON                  986
#define BSIM4v7_VDSAT                987
#define BSIM4v7_QBS                  988
#define BSIM4v7_QBD                  989
#define BSIM4v7_SOURCECONDUCT        990
#define BSIM4v7_DRAINCONDUCT         991
#define BSIM4v7_CBDB                 992
#define BSIM4v7_CBSB                 993
#define BSIM4v7_CSUB                   994
#define BSIM4v7_QINV                   995
#define BSIM4v7_IGIDL                   996
#define BSIM4v7_CSGB                 997
#define BSIM4v7_CSDB                 998
#define BSIM4v7_CSSB                 999
#define BSIM4v7_CGBB                 1000
#define BSIM4v7_CDBB                 1001
#define BSIM4v7_CSBB                 1002
#define BSIM4v7_CBBB                 1003
#define BSIM4v7_QS                   1004
#define BSIM4v7_IGISL                   1005
#define BSIM4v7_IGS                   1006
#define BSIM4v7_IGD                   1007
#define BSIM4v7_IGB                   1008
#define BSIM4v7_IGCS                   1009
#define BSIM4v7_IGCD                   1010
#define BSIM4v7_QDEF                   1011
#define BSIM4v7_DELVT0                   1012
#define BSIM4v7_GCRG                 1013
#define BSIM4v7_GTAU                 1014

#define BSIM4v7_MOD_LTVOFF           1051
#define BSIM4v7_MOD_LTVFBSDOFF       1052
#define BSIM4v7_MOD_WTVOFF           1053
#define BSIM4v7_MOD_WTVFBSDOFF       1054
#define BSIM4v7_MOD_PTVOFF           1055
#define BSIM4v7_MOD_PTVFBSDOFF       1056

#define BSIM4v7_MOD_LKVTH0WE          1061
#define BSIM4v7_MOD_LK2WE             1062
#define BSIM4v7_MOD_LKU0WE                1063
#define BSIM4v7_MOD_WKVTH0WE          1064
#define BSIM4v7_MOD_WK2WE             1065
#define BSIM4v7_MOD_WKU0WE                1066
#define BSIM4v7_MOD_PKVTH0WE          1067
#define BSIM4v7_MOD_PK2WE             1068
#define BSIM4v7_MOD_PKU0WE                1069

#define BSIM4v7_MOD_RBPS0               1101
#define BSIM4v7_MOD_RBPSL               1102
#define BSIM4v7_MOD_RBPSW               1103
#define BSIM4v7_MOD_RBPSNF              1104
#define BSIM4v7_MOD_RBPD0               1105
#define BSIM4v7_MOD_RBPDL               1106
#define BSIM4v7_MOD_RBPDW               1107
#define BSIM4v7_MOD_RBPDNF              1108

#define BSIM4v7_MOD_RBPBX0              1109
#define BSIM4v7_MOD_RBPBXL              1110
#define BSIM4v7_MOD_RBPBXW              1111
#define BSIM4v7_MOD_RBPBXNF             1112
#define BSIM4v7_MOD_RBPBY0              1113
#define BSIM4v7_MOD_RBPBYL              1114
#define BSIM4v7_MOD_RBPBYW              1115
#define BSIM4v7_MOD_RBPBYNF             1116

#define BSIM4v7_MOD_RBSBX0              1117
#define BSIM4v7_MOD_RBSBY0              1118
#define BSIM4v7_MOD_RBDBX0              1119
#define BSIM4v7_MOD_RBDBY0              1120

#define BSIM4v7_MOD_RBSDBXL             1121
#define BSIM4v7_MOD_RBSDBXW             1122
#define BSIM4v7_MOD_RBSDBXNF            1123
#define BSIM4v7_MOD_RBSDBYL             1124
#define BSIM4v7_MOD_RBSDBYW             1125
#define BSIM4v7_MOD_RBSDBYNF            1126

#define BSIM4v7_MOD_AGISL               1200
#define BSIM4v7_MOD_BGISL               1201
#define BSIM4v7_MOD_EGISL               1202
#define BSIM4v7_MOD_CGISL               1203
#define BSIM4v7_MOD_LAGISL              1204
#define BSIM4v7_MOD_LBGISL              1205
#define BSIM4v7_MOD_LEGISL              1206
#define BSIM4v7_MOD_LCGISL              1207
#define BSIM4v7_MOD_WAGISL              1208
#define BSIM4v7_MOD_WBGISL              1209
#define BSIM4v7_MOD_WEGISL              1210
#define BSIM4v7_MOD_WCGISL              1211
#define BSIM4v7_MOD_PAGISL              1212
#define BSIM4v7_MOD_PBGISL              1213
#define BSIM4v7_MOD_PEGISL              1214
#define BSIM4v7_MOD_PCGISL              1215

#define BSIM4v7_MOD_AIGS                1220
#define BSIM4v7_MOD_BIGS                1221
#define BSIM4v7_MOD_CIGS                1222
#define BSIM4v7_MOD_LAIGS               1223
#define BSIM4v7_MOD_LBIGS               1224
#define BSIM4v7_MOD_LCIGS               1225
#define BSIM4v7_MOD_WAIGS               1226
#define BSIM4v7_MOD_WBIGS               1227
#define BSIM4v7_MOD_WCIGS               1228
#define BSIM4v7_MOD_PAIGS               1229
#define BSIM4v7_MOD_PBIGS               1230
#define BSIM4v7_MOD_PCIGS               1231
#define BSIM4v7_MOD_AIGD                1232
#define BSIM4v7_MOD_BIGD                1233
#define BSIM4v7_MOD_CIGD                1234
#define BSIM4v7_MOD_LAIGD               1235
#define BSIM4v7_MOD_LBIGD               1236
#define BSIM4v7_MOD_LCIGD               1237
#define BSIM4v7_MOD_WAIGD               1238
#define BSIM4v7_MOD_WBIGD               1239
#define BSIM4v7_MOD_WCIGD               1240
#define BSIM4v7_MOD_PAIGD               1241
#define BSIM4v7_MOD_PBIGD               1242
#define BSIM4v7_MOD_PCIGD               1243
#define BSIM4v7_MOD_DLCIGD              1244

#define BSIM4v7_MOD_NJTSD               1250
#define BSIM4v7_MOD_NJTSSWD             1251
#define BSIM4v7_MOD_NJTSSWGD            1252
#define BSIM4v7_MOD_TNJTSD              1253
#define BSIM4v7_MOD_TNJTSSWD            1254
#define BSIM4v7_MOD_TNJTSSWGD           1255

/* v4.7 temp dep of leakage current  */

#define BSIM4v7_MOD_TNFACTOR          1256
#define BSIM4v7_MOD_TETA0                    1257
#define BSIM4v7_MOD_TVOFFCV           1258
#define BSIM4v7_MOD_LTNFACTOR         1260
#define BSIM4v7_MOD_LTETA0            1261
#define BSIM4v7_MOD_LTVOFFCV          1262
#define BSIM4v7_MOD_WTNFACTOR         1264
#define BSIM4v7_MOD_WTETA0            1265
#define BSIM4v7_MOD_WTVOFFCV          1266
#define BSIM4v7_MOD_PTNFACTOR         1268
#define BSIM4v7_MOD_PTETA0            1269
#define BSIM4v7_MOD_PTVOFFCV          1270

/* tnoiMod=2 (v4.7) */
#define BSIM4v7_MOD_TNOIC             1272
#define BSIM4v7_MOD_RNOIC             1273

#define BSIM4v7_MOD_VGS_MAX            1301
#define BSIM4v7_MOD_VGD_MAX            1302
#define BSIM4v7_MOD_VGB_MAX            1303
#define BSIM4v7_MOD_VDS_MAX            1304
#define BSIM4v7_MOD_VBS_MAX            1305
#define BSIM4v7_MOD_VBD_MAX            1306

#include "bsim4v7ext.h"

extern void BSIM4v7evaluate(double,double,double,BSIM4v7instance*,BSIM4v7model*,
        double*,double*,double*, double*, double*, double*, double*,
        double*, double*, double*, double*, double*, double*, double*,
        double*, double*, double*, double*, CKTcircuit*);
extern int BSIM4v7debug(BSIM4v7model*, BSIM4v7instance*, CKTcircuit*, int);
extern int BSIM4v7checkModel(BSIM4v7model*, BSIM4v7instance*, CKTcircuit*);
extern int BSIM4v7PAeffGeo(double, int, int, double, double, double, double, double *, double *, double *, double *);
extern int BSIM4v7RdseffGeo(double, int, int, int, double, double, double, double, double, int, double *);
extern int BSIM4v7RdsEndIso(double, double, double, double, double, double, int, int, double *);
extern int BSIM4v7RdsEndSha(double, double, double, double, double, double, int, int, double *);

#endif /*BSIM4v7*/
