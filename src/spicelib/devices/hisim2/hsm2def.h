/***********************************************************************

 HiSIM (Hiroshima University STARC IGFET Model)
 Copyright (C) 2012 Hiroshima University & STARC

 MODEL NAME : HiSIM
 ( VERSION : 2  SUBVERSION : 7  REVISION : 0 ) Beta
 
 FILE : hsm2def.h

 Date : 2012.10.25

 released by 
                Hiroshima University &
                Semiconductor Technology Academic Research Center (STARC)
***********************************************************************/

#ifndef HSM2
#define HSM2

#include "ngspice/ifsim.h"
#include "ngspice/gendefs.h"
#include "ngspice/cktdefs.h"
#include "ngspice/complex.h"
#include "ngspice/noisedef.h"

/* declarations for HiSIM2 MOSFETs */

/* unit-converted model parameters */
typedef struct sHSM2modelMKSParam {
  double HSM2_npext ;
  double HSM2_nsubcwpe ;
  double HSM2_nsubpwpe ;
  double HSM2_npextwpe ;
  double HSM2_ll ;
  double HSM2_wl ;
  double HSM2_svgsl ;
  double HSM2_svgsw ;
  double HSM2_svbsl ;
  double HSM2_slgl ;
  double HSM2_sub1l ;
  double HSM2_slg ;
  double HSM2_sub2l ;
  double HSM2_nsubcmax ;
  double HSM2_glksd3 ;
  double HSM2_gleak2 ;
  double HSM2_gleak4 ;
  double HSM2_gleak5 ;
  double HSM2_gleak7 ;
  double HSM2_cit ;
  double HSM2_ovslp ;
  double HSM2_dly3 ;
} HSM2modelMKSParam ;

/* binning parameters */
typedef struct sHSM2binningParam {
  double HSM2_vmax ;
  double HSM2_bgtmp1 ;
  double HSM2_bgtmp2 ;
  double HSM2_eg0 ;
  double HSM2_lover ;
  double HSM2_vfbover ;
  double HSM2_nover ;
  double HSM2_wl2 ;
  double HSM2_vfbc ;
  double HSM2_nsubc ;
  double HSM2_nsubp ;
  double HSM2_scp1 ;
  double HSM2_scp2 ;
  double HSM2_scp3 ;
  double HSM2_sc1 ;
  double HSM2_sc2 ;
  double HSM2_sc3 ;
  double HSM2_sc4 ;
  double HSM2_pgd1 ;
//double HSM2_pgd3 ;
  double HSM2_ndep ;
  double HSM2_ninv ;
  double HSM2_muecb0 ;
  double HSM2_muecb1 ;
  double HSM2_mueph1 ;
  double HSM2_vtmp ;
  double HSM2_wvth0 ;
  double HSM2_muesr1 ;
  double HSM2_muetmp ;
  double HSM2_sub1 ;
  double HSM2_sub2 ;
  double HSM2_svds ;
  double HSM2_svbs ;
  double HSM2_svgs ;
  double HSM2_nsti ;
  double HSM2_wsti ;
  double HSM2_scsti1 ;
  double HSM2_scsti2 ;
  double HSM2_vthsti ;
  double HSM2_muesti1 ;
  double HSM2_muesti2 ;
  double HSM2_muesti3 ;
  double HSM2_nsubpsti1 ;
  double HSM2_nsubpsti2 ;
  double HSM2_nsubpsti3 ;
  double HSM2_nsubcsti1;
  double HSM2_nsubcsti2;
  double HSM2_nsubcsti3;
  double HSM2_cgso ;
  double HSM2_cgdo ;
  double HSM2_js0 ;
  double HSM2_js0sw ;
  double HSM2_nj ;
  double HSM2_cisbk ;
  double HSM2_clm1 ;
  double HSM2_clm2 ;
  double HSM2_clm3 ;
  double HSM2_wfc ;
  double HSM2_gidl1 ;
  double HSM2_gidl2 ;
  double HSM2_gleak1 ;
  double HSM2_gleak2 ;
  double HSM2_gleak3 ;
  double HSM2_gleak6 ;
  double HSM2_glksd1 ;
  double HSM2_glksd2 ;
  double HSM2_glkb1 ;
  double HSM2_glkb2 ;
  double HSM2_nftrp ;
  double HSM2_nfalp ;
  double HSM2_vdiffj ;
  double HSM2_ibpc1 ;
  double HSM2_ibpc2 ;
} HSM2binningParam ;

/* unit-converted parameters for each instance */
typedef struct sHSM2hereMKSParam {
  double HSM2_nsubcdfm ;
} HSM2hereMKSParam ;

/* information needed for each instance */
typedef struct sHSM2instance {
  struct sHSM2model *HSM2modPtr;           /* pointer to model */
  struct sHSM2instance *HSM2nextInstance;  /* pointer to next instance of 
                                              current model*/
  IFuid HSM2name; /* pointer to character string naming this instance */
  int HSM2states; /* index into state table for this device */

  int HSM2dNode;      /* number of the drain node of the mosfet */
  int HSM2gNode;      /* number of the gate node of the mosfet */
  int HSM2sNode;      /* number of the source node of the mosfet */
  int HSM2bNode;      /* number of the bulk node of the mosfet */
  int HSM2dNodePrime; /* number od the inner drain node */
  int HSM2gNodePrime; /* number of the inner gate node */
  int HSM2sNodePrime; /* number od the inner source node */
  int HSM2bNodePrime;
  int HSM2dbNode;
  int HSM2sbNode;

  double HSM2_noiflick; /* for 1/f noise calc. */
  double HSM2_noithrml; /* for thermal noise calc. */
  double HSM2_noiigate; /* for induced gate noise */
  double HSM2_noicross; /* for induced gate noise */
  double HSM2_Qdrat;    /* for induced gate noise */

  /* instance */
  double HSM2_l;    /* the length of the channel region */
  double HSM2_w;    /* the width of the channel region */
  double HSM2_ad;   /* the area of the drain diffusion */
  double HSM2_as;   /* the area of the source diffusion */
  double HSM2_pd;   /* perimeter of drain junction [m] */
  double HSM2_ps;   /* perimeter of source junction [m] */
  double HSM2_nrd;  /* equivalent num of squares of drain [-] (unused) */
  double HSM2_nrs;  /* equivalent num of squares of source [-] (unused) */
  double HSM2_temp; /* lattice temperature [C] */
  double HSM2_dtemp;

  double HSM2_weff; /* the effective width of the channel region */
  double HSM2_weff_nf; /* Weff * NF */
  double HSM2_leff; /* the effective length of the channel region */

  int HSM2_corbnet  ;
  double HSM2_rbpb ;
  double HSM2_rbpd ;
  double HSM2_rbps ;
  double HSM2_rbdb ;
  double HSM2_rbsb ;

  int HSM2_corg ;
/*   double HSM2_rshg; */
  double HSM2_ngcon;
  double HSM2_xgw;
  double HSM2_xgl;
  double HSM2_nf;

  double HSM2_sa;
  double HSM2_sb;
  double HSM2_sd;
  double HSM2_nsubcdfm; /* DFM */
  double HSM2_mphdfm; /* DFM */
  double HSM2_m;
  
/* WPE */
  double HSM2_sca;  /* scc */
  double HSM2_scb;  /* scb */
  double HSM2_scc;  /* scc */

  int HSM2_called; /* flag to check the first call */
  /* previous values to evaluate initial guess */
  double HSM2_mode_prv;
  double HSM2_vbsc_prv;
  double HSM2_vdsc_prv;
  double HSM2_vgsc_prv;
  double HSM2_ps0_prv;
  double HSM2_ps0_dvbs_prv;
  double HSM2_ps0_dvds_prv;
  double HSM2_ps0_dvgs_prv;
  double HSM2_pds_prv;
  double HSM2_pds_dvbs_prv;
  double HSM2_pds_dvds_prv;
  double HSM2_pds_dvgs_prv;
  double HSM2_ids_prv;
  double HSM2_ids_dvbs_prv;
  double HSM2_ids_dvds_prv;
  double HSM2_ids_dvgs_prv;
  double HSM2_mode_prv2;
  double HSM2_vbsc_prv2;
  double HSM2_vdsc_prv2;
  double HSM2_vgsc_prv2;
  double HSM2_ps0_prv2;
  double HSM2_ps0_dvbs_prv2;
  double HSM2_ps0_dvds_prv2;
  double HSM2_ps0_dvgs_prv2;
  double HSM2_pds_prv2;
  double HSM2_pds_dvbs_prv2;
  double HSM2_pds_dvds_prv2;
  double HSM2_pds_dvgs_prv2;
  double HSM2_PS0Z_SCE_prv ;
  double HSM2_PS0Z_SCE_dvds_prv ;
  double HSM2_PS0Z_SCE_dvgs_prv ;
  double HSM2_PS0Z_SCE_dvbs_prv ;

  /* output */
  int    HSM2_capop;
  double HSM2_gd;
  double HSM2_gs;
  double HSM2_cgso;
  double HSM2_cgdo;
  double HSM2_cgbo;
  double HSM2_cdso;
  double HSM2_cddo;
  double HSM2_cdgo;
  double HSM2_csso;
  double HSM2_csdo;
  double HSM2_csgo;
  double HSM2_cqyd;
  double HSM2_cqyg;
  double HSM2_cqyb;
  double HSM2_von; /* vth */
  double HSM2_vdsat;
  double HSM2_ids; /* cdrain, HSM2_cd */
  double HSM2_gds;
  double HSM2_gm;
  double HSM2_gmbs;
  double HSM2_ibs; /* HSM2_cbs */
  double HSM2_ibd; /* HSM2_cbd */
  double HSM2_gbs;
  double HSM2_gbd;
  double HSM2_capbs;
  double HSM2_capbd;
  double HSM2_capgs;
  double HSM2_capgd;
  double HSM2_capgb;
  double HSM2_isub; /* HSM2_csub */
  double HSM2_gbgs;
  double HSM2_gbds;
  double HSM2_gbbs;
  double HSM2_qg;
  double HSM2_qd;
  double HSM2_qs;
  double HSM2_qb;  /* bulk charge qb = -(qg + qd + qs) */
  double HSM2_cggb;
  double HSM2_cgdb;
  double HSM2_cgsb;
  double HSM2_cbgb;
  double HSM2_cbdb;
  double HSM2_cbsb;
  double HSM2_cdgb;
  double HSM2_cddb;
  double HSM2_cdsb;

  double HSM2_mu; /* mobility */
  double HSM2_igidl; /* gate induced drain leakage */
  double HSM2_gigidlgs;
  double HSM2_gigidlds;
  double HSM2_gigidlbs;
  double HSM2_igisl; /* gate induced source leakage */
  double HSM2_gigislgd;
  double HSM2_gigislsd;
  double HSM2_gigislbd;
  double HSM2_igb; /* gate tunneling current (gate to bulk) */
  double HSM2_gigbg;
  double HSM2_gigbd;
  double HSM2_gigbb;
  double HSM2_gigbs;
  double HSM2_igs; /* gate tunneling current (gate to source) */
  double HSM2_gigsg;
  double HSM2_gigsd;
  double HSM2_gigsb;
  double HSM2_gigss;
  double HSM2_igd; /* gate tunneling current (gate to drain) */
  double HSM2_gigdg;
  double HSM2_gigdd;
  double HSM2_gigdb;
  double HSM2_gigds;

  /* NQS */
  double HSM2_tau ;
  double HSM2_tau_dVgs ;
  double HSM2_tau_dVds ;
  double HSM2_tau_dVbs ;
  double HSM2_Xd  ;
  double HSM2_Xd_dVgs  ;
  double HSM2_Xd_dVds  ;
  double HSM2_Xd_dVbs  ;
  double HSM2_Qi  ;
  double HSM2_Qi_dVgs  ;
  double HSM2_Qi_dVds  ;
  double HSM2_Qi_dVbs  ;
  double HSM2_taub  ;
  double HSM2_taub_dVgs  ;
  double HSM2_taub_dVds  ;
  double HSM2_taub_dVbs  ;
  double HSM2_Qb  ;
  double HSM2_Qb_dVgs  ;
  double HSM2_Qb_dVds  ;
  double HSM2_Qb_dVbs  ;
  double HSM2_alpha;

  /* internal variables */
  double HSM2_eg ;
  double HSM2_beta ;
  double HSM2_beta_inv ;
  double HSM2_beta2 ;
  double HSM2_betatnom ;
  double HSM2_nin ;
  double HSM2_egp12 ;
  double HSM2_egp32 ;
  double HSM2_lgate ;
  double HSM2_wg ;
  double HSM2_mueph ;
  double HSM2_mphn0 ;
  double HSM2_mphn1 ;
  double HSM2_muesr ;
  double HSM2_nsub ;
  double HSM2_qnsub ;
  double HSM2_qnsub_esi ;
  double HSM2_2qnsub_esi ;
  double HSM2_ptovr0 ;
  double HSM2_ptovr ;
  double HSM2_vmax0 ;
  double HSM2_vmax ;
  double HSM2_pb2 ;
  double HSM2_pb20 ;
  double HSM2_pb2c ;
  double HSM2_cnst0 ;
  double HSM2_cnst1 ;
  double HSM2_isbd ;
  double HSM2_isbd2 ;
  double HSM2_isbs ;
  double HSM2_isbs2 ;
  double HSM2_vbdt ;
  double HSM2_vbst ;
  double HSM2_exptemp ;
  double HSM2_wsti ;
  double HSM2_cnstpgd ;
  double HSM2_ninvp0 ;
  double HSM2_ninv0 ;
  double HSM2_grbpb ;
  double HSM2_grbpd ;
  double HSM2_grbps ;
  double HSM2_grbdb ;
  double HSM2_grbsb ;
  double HSM2_grg ;
  double HSM2_rs ;
  double HSM2_rd ;
  double HSM2_clmmod ;
  double HSM2_lgatesm ;
  double HSM2_dVthsm ;
  double HSM2_ddlt ;
  /* 2007.02.20--03.15 */
  double HSM2_xsub1 ;
  double HSM2_xsub2 ;
  double HSM2_xgate ;
  double HSM2_xvbs ;
  double HSM2_vg2const ;
  double HSM2_wdpl ;
  double HSM2_wdplp ;
  double HSM2_cfrng ;
  double HSM2_jd_nvtm_inv ;
  double HSM2_jd_expcd ;
  double HSM2_jd_expcs ;
  double HSM2_sqrt_eg ;

  double HSM2_egtnom ;
  double HSM2_cecox ;
  double HSM2_msc ;
  int HSM2_flg_pgd ;
  double HSM2_ndep_o_esi ;
  double HSM2_ninv_o_esi ;
  double HSM2_cqyb0 ;
  double HSM2_cnst0over ;
  double HSM2_costi00 ;
  double HSM2_nsti_p2 ;
  double HSM2_costi0 ;
  double HSM2_costi0_p2 ;
  double HSM2_costi1 ;
  double HSM2_pb2over ; /* for Qover model */
  double HSM2_ps0ldinib ;
  double HSM2_ptl0;
  double HSM2_pt40;
  double HSM2_gdl0;
  double HSM2_muecb0;
  double HSM2_muecb1;
  double HSM2_ktemp; /* lattice temperature [K] */
  double HSM2_mueph1 ;
  double HSM2_nsubp ;
  double HSM2_nsubc ;

  HSM2hereMKSParam hereMKS ; /* unit-converted parameters */
  HSM2binningParam pParam ; /* binning parameters */
  
  /* no use in SPICE3f5
      double HSM2drainSquares;       the length of the drain in squares
      double HSM2sourceSquares;      the length of the source in squares */
  double HSM2sourceConductance; /* cond. of source (or 0): set in setup */
  double HSM2drainConductance;  /* cond. of drain (or 0): set in setup */
  double HSM2internalGs; /* internal cond. of source for thermal noise calc. */
  double HSM2internalGd; /* internal cond. of drain for thermal noise calc. */

  double HSM2_icVBS; /* initial condition B-S voltage */
  double HSM2_icVDS; /* initial condition D-S voltage */
  double HSM2_icVGS; /* initial condition G-S voltage */
  int HSM2_off;      /* non-zero to indicate device is off for dc analysis */
  int HSM2_mode;     /* device mode : 1 = normal, -1 = inverse */

  unsigned HSM2_l_Given :1;
  unsigned HSM2_w_Given :1;
  unsigned HSM2_ad_Given :1;
  unsigned HSM2_as_Given    :1;
  /*  unsigned HSM2drainSquaresGiven  :1;
      unsigned HSM2sourceSquaresGiven :1;*/
  unsigned HSM2_pd_Given    :1;
  unsigned HSM2_ps_Given   :1;
  unsigned HSM2_nrd_Given  :1;
  unsigned HSM2_nrs_Given  :1;
  unsigned HSM2_temp_Given  :1;
  unsigned HSM2_dtemp_Given  :1;
  unsigned HSM2_icVBS_Given :1;
  unsigned HSM2_icVDS_Given :1;
  unsigned HSM2_icVGS_Given :1;
  unsigned HSM2_corbnet_Given  :1;
  unsigned HSM2_rbpb_Given :1;
  unsigned HSM2_rbpd_Given :1;
  unsigned HSM2_rbps_Given :1;
  unsigned HSM2_rbdb_Given :1;
  unsigned HSM2_rbsb_Given :1;
  unsigned HSM2_corg_Given  :1;
/*   unsigned HSM2_rshg_Given  :1; */
  unsigned HSM2_ngcon_Given  :1;
  unsigned HSM2_xgw_Given  :1;
  unsigned HSM2_xgl_Given  :1;
  unsigned HSM2_nf_Given  :1;
  unsigned HSM2_sa_Given  :1;
  unsigned HSM2_sb_Given  :1;
  unsigned HSM2_sd_Given  :1;
  unsigned HSM2_nsubcdfm_Given  :1; /* DFM */
  unsigned HSM2_mphdfm_Given  :1; /* DFM */
  unsigned HSM2_m_Given  :1;
  
 /* WPE */
  unsigned HSM2_sca_Given :1;	/* sca */
  unsigned HSM2_scb_Given :1;	/* scb */
  unsigned HSM2_scc_Given :1;	/* scc */
  /* pointer to sparse matrix */

  double *HSM2GgPtr;   /* pointer to sparse matrix element at (gate node,gate node) */
  double *HSM2GgpPtr;  /* pointer to sparse matrix element at (gate node,gate prime node) */
  double *HSM2GdpPtr;  /* pointer to sparse matrix element at (gate node,drain prime node) */
  double *HSM2GspPtr;  /* pointer to sparse matrix element at (gate node,source prime node) */
  double *HSM2GbpPtr;  /* pointer to sparse matrix element at (gate node,bulk prime node) */

  double *HSM2GPgPtr;  /* pointer to sparse matrix element at (gate prime node,gate node) */
  double *HSM2GPgpPtr;  /* pointer to sparse matrix element at (gate prime node,gate prime node) */
  double *HSM2GPdpPtr;  /* pointer to sparse matrix element at (gate prime node,drain prime node) */
  double *HSM2GPspPtr;  /* pointer to sparse matrix element at (gate prime node,source prime node) */
  double *HSM2GPbpPtr;  /* pointer to sparse matrix element at (gate prime node,bulk prime node) */

  double *HSM2DPdPtr;  /* pointer to sparse matrix element at (drain prime node,drain node) */
  double *HSM2DPdpPtr; /* pointer to sparse matrix element at (drain prime node,drain prime node) */
  double *HSM2DPgpPtr; /* pointer to sparse matrix element at (drain prime node,gate prime node) */
  double *HSM2DPspPtr; /* pointer to sparse matrix element at (drain prime node,source prime node) */
  double *HSM2DPbpPtr; /* pointer to sparse matrix element at (drain prime node,bulk prime node) */
  double *HSM2DPdbPtr; /* pointer to sparse matrix element at (drain prime node,drain body node) */

  double *HSM2DdPtr;   /* pointer to sparse matrix element at (Drain node,drain node) */
  double *HSM2DdpPtr;  /* pointer to sparse matrix element at (drain node,drain prime node) */

  double *HSM2SPsPtr;  /* pointer to sparse matrix element at (source prime node,source node) */
  double *HSM2SPspPtr; /* pointer to sparse matrix element at (source prime node,source prime node) */
  double *HSM2SPgpPtr; /* pointer to sparse matrix element at (source prime node,gate prime node) */
  double *HSM2SPdpPtr; /* pointer to sparse matrix element at (source prime node,drain prime node) */
  double *HSM2SPbpPtr; /* pointer to sparse matrix element at (source prime node,bulk prime node) */
  double *HSM2SPsbPtr; /* pointer to sparse matrix element at (source prime node,source body node) */

  double *HSM2SsPtr;   /* pointer to sparse matrix element at (source node,source node) */
  double *HSM2SspPtr;  /* pointer to sparse matrix element at (source node,source prime node) */

  double *HSM2BPgpPtr;  /* pointer to sparse matrix element at (bulk prime node,gate prime node) */
  double *HSM2BPbpPtr; /* pointer to sparse matrix element at (bulk prime node,bulk prime node) */
  double *HSM2BPdpPtr; /* pointer to sparse matrix element at (bulk prime node,drain prime node) */
  double *HSM2BPspPtr; /* pointer to sparse matrix element at (bulk prime node,source prime node) */
  double *HSM2BPbPtr;  /* pointer to sparse matrix element at (bulk prime node,bulk node) */
  double *HSM2BPdbPtr; /* pointer to sparse matrix element at (bulk prime node,source body node) */
  double *HSM2BPsbPtr; /* pointer to sparse matrix element at (bulk prime node,source body node) */

  double *HSM2DBdpPtr; /* pointer to sparse matrix element at (drain body node,drain prime node) */
  double *HSM2DBdbPtr; /* pointer to sparse matrix element at (drain body node,drain body node) */
  double *HSM2DBbpPtr; /* pointer to sparse matrix element at (drain body node,bulk prime node) */
  double *HSM2DBbPtr;  /* pointer to sparse matrix element at (drain body node,bulk node) */
  
  double *HSM2SBspPtr; /* pointer to sparse matrix element at (source body node,drain prime node) */
  double *HSM2SBbpPtr; /* pointer to sparse matrix element at (source body node,drain body node) */
  double *HSM2SBbPtr;  /* pointer to sparse matrix element at (source body node,bulk prime node) */
  double *HSM2SBsbPtr; /* pointer to sparse matrix element at (source body node,bulk node) */

  double *HSM2BsbPtr;  /* pointer to sparse matrix element at (bulk node,source body node) */
  double *HSM2BbpPtr;  /* pointer to sparse matrix element at (bulk node,bulk prime node) */
  double *HSM2BdbPtr;  /* pointer to sparse matrix element at (bulk node,drain body node) */
  double *HSM2BbPtr;   /* pointer to sparse matrix element at (bulk node,bulk node) */

#ifdef USE_OMP
    /* per instance storage of results, to update matrix and rhs at a later stage */
    double HSM2rhsdPrime;
    double HSM2rhsgPrime;
    double HSM2rhsbPrime;
    double HSM2rhssPrime;
    double HSM2rhsdb;
    double HSM2rhssb;

    double HSM2_1;
    double HSM2_2;
    double HSM2_3;
    double HSM2_4;
    double HSM2_5;
    double HSM2_6;
    double HSM2_7;
    double HSM2_8;
    double HSM2_9;
    double HSM2_10;
    double HSM2_11;
    double HSM2_12;
    double HSM2_13;
    double HSM2_14;
    double HSM2_15;
    double HSM2_16;
    double HSM2_17;
    double HSM2_18;
    double HSM2_19;
    double HSM2_20;
    double HSM2_21;
    double HSM2_22;
    double HSM2_23;
    double HSM2_24;
    double HSM2_25;
    double HSM2_26;
    double HSM2_27;
    double HSM2_28;
    double HSM2_29;
    double HSM2_30;
    double HSM2_31;
    double HSM2_32;
    double HSM2_33;
    double HSM2_34;
    double HSM2_35;
    double HSM2_36;
    double HSM2_37;
    double HSM2_38;
    double HSM2_39;
    double HSM2_40;
    double HSM2_41;
    double HSM2_42;
    double HSM2_43;
    double HSM2_44;
    double HSM2_45;
    double HSM2_46;
    double HSM2_47;
    double HSM2_48;
    double HSM2_49;
    double HSM2_50;
    double HSM2_51;
    double HSM2_52;
    double HSM2_53;
    double HSM2_54;
    double HSM2_55;
    double HSM2_56;
    double HSM2_57;
    double HSM2_58;
    double HSM2_59;
    double HSM2_60;
    double HSM2_61;
    double HSM2_62;
    double HSM2_63;
#endif

  /* common state values in hisim module */
#define HSM2vbd HSM2states+ 0
#define HSM2vbs HSM2states+ 1
#define HSM2vgs HSM2states+ 2
#define HSM2vds HSM2states+ 3
#define HSM2vdbs HSM2states+ 4
#define HSM2vdbd HSM2states+ 5
#define HSM2vsbs HSM2states+ 6
#define HSM2vges HSM2states+ 7

#define HSM2qb  HSM2states+ 8
#define HSM2cqb HSM2states+ 9
#define HSM2qg  HSM2states+ 10
#define HSM2cqg HSM2states+ 11
#define HSM2qd  HSM2states+ 12
#define HSM2cqd HSM2states+ 13

#define HSM2qbs HSM2states+ 14
#define HSM2cqbs HSM2states+ 15
#define HSM2qbd HSM2states+ 16
#define HSM2cqbd HSM2states+ 17

#define HSM2numStates 18

/* nqs charges */
#define HSM2qi_nqs HSM2states+ 19
#define HSM2qb_nqs HSM2states+ 20

#define HSM2numStatesNqs 21

  /* indices to the array of HiSIM2 NOISE SOURCES (the same as BSIM3) */
#define HSM2RDNOIZ       0
#define HSM2RSNOIZ       1
#define HSM2IDNOIZ       2
#define HSM2FLNOIZ       3
#define HSM2IGSNOIZ      4  /* shot noise */
#define HSM2IGDNOIZ      5  /* shot noise */
#define HSM2IGBNOIZ      6  /* shot noise */
#define HSM2IGNOIZ       7  /* induced gate noise */
#define HSM2TOTNOIZ      8

#define HSM2NSRCS        9  /* the number of HiSIM2 MOSFET noise sources */

#ifndef NONOISE
  double HSM2nVar[NSTATVARS][HSM2NSRCS];
#else /* NONOISE */
  double **HSM2nVar;
#endif /* NONOISE */

} HSM2instance ;


/* per model data */

typedef struct sHSM2model {       	/* model structure for a resistor */
  int HSM2modType;    		/* type index of this device type */
  struct sHSM2model *HSM2nextModel; /* pointer to next possible model 
					 in linked list */
  HSM2instance * HSM2instances;	/* pointer to list of instances 
				   that have this model */
  IFuid HSM2modName;       	/* pointer to the name of this model */

  /* --- end of generic struct GENmodel --- */

  int HSM2_type;      		/* device type: 1 = nmos,  -1 = pmos */
  int HSM2_level;               /* level */
  int HSM2_info;                /* information */
  int HSM2_noise;               /* noise model selecter see hsm2noi.c */
  int HSM2_version;             /* model version 200 */
  int HSM2_show;                /* show physical value 1, 2, ... , 11 */

  /* flags for initial guess */
  int HSM2_corsrd ;
  int HSM2_corg   ;
  int HSM2_coiprv ;
  int HSM2_copprv ;
  int HSM2_coadov ;
  int HSM2_coisub ;
  int HSM2_coiigs ;
  int HSM2_cogidl ;
  int HSM2_coovlp ;
  int HSM2_coflick ;
  int HSM2_coisti ;
  int HSM2_conqs  ; /* HiSIM2 */
  int HSM2_corbnet ;
  int HSM2_cothrml;
  int HSM2_coign; /* induced gate noise */
  int HSM2_codfm; /* DFM */
  int HSM2_corecip;
  int HSM2_coqy;
  int HSM2_coqovsm ;
  int HSM2_coerrrep;

  /* HiSIM original */
  double HSM2_vmax ;
  double HSM2_bgtmp1 ;
  double HSM2_bgtmp2 ;
  double HSM2_eg0 ;
  double HSM2_tox ;
  double HSM2_xld ;
  double HSM2_lover ;
  double HSM2_ddltmax ; /* Vdseff */
  double HSM2_ddltslp ; /* Vdseff */
  double HSM2_ddltict ; /* Vdseff */
  double HSM2_vfbover ;
  double HSM2_nover ;
  double HSM2_xwd ;
  double HSM2_xl ;
  double HSM2_xw ;
  double HSM2_saref ;
  double HSM2_sbref ;
  double HSM2_ll ;
  double HSM2_lld ;
  double HSM2_lln ;
  double HSM2_wl ;
  double HSM2_wl1 ;
  double HSM2_wl1p ;
  double HSM2_wl2 ;
  double HSM2_wl2p ;
  double HSM2_wld ;
  double HSM2_wln ;
  double HSM2_xqy ;
  double HSM2_xqy1 ;
  double HSM2_xqy2 ;
  double HSM2_qyrat ;
  double HSM2_rs;     /* source contact resistance */
  double HSM2_rd;     /* drain contact resistance */
  double HSM2_rsh;    /* source/drain diffusion sheet resistance */
  double HSM2_rshg;
/*   double HSM2_ngcon; */
/*   double HSM2_xgw; */
/*   double HSM2_xgl; */
/*   double HSM2_nf; */
  double HSM2_vfbc ;
  double HSM2_vbi ;
  double HSM2_vfbcl;
  double HSM2_vfbclp;
  double HSM2_nsubc ;
  double HSM2_parl2 ;
  double HSM2_lp ;
  double HSM2_nsubp ;
  double HSM2_nsubpl ;
  double HSM2_nsubpdlt;
  double HSM2_nsubpfac ;
  double HSM2_nsubpw ;
  double HSM2_nsubpwp ;
  double HSM2_scp1 ;
  double HSM2_scp2 ;
  double HSM2_scp3 ;
  double HSM2_sc1 ;
  double HSM2_sc2 ;
  double HSM2_sc3 ;
  double HSM2_sc4 ;
  double HSM2_pgd1 ;
  double HSM2_pgd2 ;
//double HSM2_pgd3 ;
  double HSM2_pgd4 ;
  double HSM2_ndep ;
  double HSM2_ndepl ;
  double HSM2_ndeplp ;
  double HSM2_ndepw ;
  double HSM2_ndepwp ;
  double HSM2_ninv ;
  double HSM2_ninvd ;
  double HSM2_muecb0 ;
  double HSM2_muecb1 ;
  double HSM2_mueph1 ;
  double HSM2_mueph0 ;
  double HSM2_muephw ;
  double HSM2_muepwp ;
  double HSM2_muepwd ;
  double HSM2_muephl ;
  double HSM2_mueplp ;
  double HSM2_muepld ;
  double HSM2_muephs ;
  double HSM2_muepsp ;
  double HSM2_vtmp ;
  double HSM2_wvth0 ;
  double HSM2_muesr1 ;
  double HSM2_muesr0 ;
  double HSM2_muesrw ;
  double HSM2_mueswp ;
  double HSM2_muesrl ;
  double HSM2_mueslp ;
  double HSM2_bb ;
  double HSM2_sub1 ;
  double HSM2_sub2 ;
  double HSM2_svgs ;
  double HSM2_svbs ;
  double HSM2_svbsl ;
  double HSM2_svds ;
  double HSM2_slg ;
  double HSM2_sub1l ;
  double HSM2_sub2l ;
  double HSM2_svgsl ;
  double HSM2_svgslp ;
  double HSM2_svgswp ;
  double HSM2_svgsw ;
  double HSM2_svbslp ;
  double HSM2_slgl ;
  double HSM2_slglp ;
  double HSM2_sub1lp ;
  double HSM2_nsti ;  
  double HSM2_wsti ;
  double HSM2_wstil ;
  double HSM2_wstilp ;
  double HSM2_wstiw ;
  double HSM2_wstiwp ;
  double HSM2_scsti1 ;
  double HSM2_scsti2 ;
  double HSM2_vthsti ;
  double HSM2_vdsti ;
  double HSM2_muesti1 ;
  double HSM2_muesti2 ;
  double HSM2_muesti3 ;
  double HSM2_nsubpsti1 ;
  double HSM2_nsubpsti2 ;
  double HSM2_nsubcsti1;
  double HSM2_nsubcsti2;
  double HSM2_nsubcsti3;
  double HSM2_nsubpsti3 ;
  double HSM2_lpext ;
  double HSM2_npext ;
  double HSM2_npextw ;
  double HSM2_npextwp ;
  double HSM2_scp22 ;
  double HSM2_scp21 ;
  double HSM2_bs1 ;
  double HSM2_bs2 ;
  double HSM2_cgso ;
  double HSM2_cgdo ;
  double HSM2_cgbo ;
  double HSM2_tpoly ;
  double HSM2_js0 ;
  double HSM2_js0sw ;
  double HSM2_nj ;
  double HSM2_njsw ;
  double HSM2_xti ;
  double HSM2_cj ;
  double HSM2_cjsw ;
  double HSM2_cjswg ;
  double HSM2_mj ;
  double HSM2_mjsw ;
  double HSM2_mjswg ;
  double HSM2_xti2 ;
  double HSM2_cisb ;
  double HSM2_cvb ;
  double HSM2_ctemp ;
  double HSM2_cisbk ;
  double HSM2_cvbk ;
  double HSM2_divx ;
  double HSM2_pb ;
  double HSM2_pbsw ;
  double HSM2_pbswg ;
  double HSM2_tcjbd ;
  double HSM2_tcjbs ;
  double HSM2_tcjbdsw ;
  double HSM2_tcjbssw ;
  double HSM2_tcjbdswg ;
  double HSM2_tcjbsswg ;

  double HSM2_clm1 ;
  double HSM2_clm2 ;
  double HSM2_clm3 ;
  double HSM2_clm5 ;
  double HSM2_clm6 ;
  double HSM2_muetmp ;
  double HSM2_vover ;
  double HSM2_voverp ;
  double HSM2_vovers ;
  double HSM2_voversp ;
  double HSM2_wfc ;
  double HSM2_nsubcw ;
  double HSM2_nsubcwp ;
  double HSM2_nsubcmax ;
  double HSM2_qme1 ;
  double HSM2_qme2 ;
  double HSM2_qme3 ;
  double HSM2_gidl1 ;
  double HSM2_gidl2 ;
  double HSM2_gidl3 ;
  double HSM2_gidl4 ;
  double HSM2_gidl6;
  double HSM2_gidl7;
  double HSM2_gidl5 ;
  double HSM2_gleak1 ;
  double HSM2_gleak2 ;
  double HSM2_gleak3 ;
  double HSM2_gleak4 ;
  double HSM2_gleak5 ;
  double HSM2_gleak6 ;
  double HSM2_gleak7 ;
  double HSM2_glksd1 ;
  double HSM2_glksd2 ;
  double HSM2_glksd3 ;
  double HSM2_glkb1 ;
  double HSM2_glkb2 ;
  double HSM2_glkb3 ;
  double HSM2_egig;
  double HSM2_igtemp2;
  double HSM2_igtemp3;
  double HSM2_vzadd0 ;
  double HSM2_pzadd0 ;
  double HSM2_nftrp ;
  double HSM2_nfalp ;
  double HSM2_falph ;
  double HSM2_cit ;
  double HSM2_kappa ;  
  double HSM2_vdiffj ; 
  double HSM2_dly1 ;
  double HSM2_dly2 ;
  double HSM2_dly3;
  double HSM2_tnom ;
  double HSM2_ovslp ;
  double HSM2_ovmag ;
  /* substrate resistances */
  double HSM2_gbmin;
  double HSM2_rbpb ;
  double HSM2_rbpd ;
  double HSM2_rbps ;
  double HSM2_rbdb ;
  double HSM2_rbsb ;
  /* IBPC */
  double HSM2_ibpc1 ;
  double HSM2_ibpc2 ;
  /* DFM */
  double HSM2_mphdfm ;

  double HSM2_ptl, HSM2_ptp, HSM2_pt2, HSM2_ptlp, HSM2_gdl, HSM2_gdlp  ;

  double HSM2_gdld ;
  double HSM2_pt4 ;
  double HSM2_pt4p ;
  double HSM2_muephl2 ;
  double HSM2_mueplp2 ;
  double HSM2_nsubcw2 ;
  double HSM2_nsubcwp2 ;
  double HSM2_muephw2 ;
  double HSM2_muepwp2 ;

  /* variables for WPE */
  double HSM2_web ;
  double HSM2_wec ;
  double HSM2_nsubcwpe ;	
  double HSM2_npextwpe ;	
  double HSM2_nsubpwpe ;	
  /* for Ps0_min */
  double HSM2_Vgsmin ; 
  double HSM2_sc3Vbs ; /* SC3 clamping  */
  double HSM2_byptol ; /* bypass control */
  double HSM2_muecb0lp;
  double HSM2_muecb1lp;


  /* binning parameters */
  double HSM2_lmin ;
  double HSM2_lmax ;
  double HSM2_wmin ;
  double HSM2_wmax ;
  double HSM2_lbinn ;
  double HSM2_wbinn ;

  /* Length dependence */
  double HSM2_lvmax ;
  double HSM2_lbgtmp1 ;
  double HSM2_lbgtmp2 ;
  double HSM2_leg0 ;
  double HSM2_llover ;
  double HSM2_lvfbover ;
  double HSM2_lnover ;
  double HSM2_lwl2 ;
  double HSM2_lvfbc ;
  double HSM2_lnsubc ;
  double HSM2_lnsubp ;
  double HSM2_lscp1 ;
  double HSM2_lscp2 ;
  double HSM2_lscp3 ;
  double HSM2_lsc1 ;
  double HSM2_lsc2 ;
  double HSM2_lsc3 ;
  double HSM2_lsc4 ;
  double HSM2_lpgd1 ;
//double HSM2_lpgd3 ;
  double HSM2_lndep ;
  double HSM2_lninv ;
  double HSM2_lmuecb0 ;
  double HSM2_lmuecb1 ;
  double HSM2_lmueph1 ;
  double HSM2_lvtmp ;
  double HSM2_lwvth0 ;
  double HSM2_lmuesr1 ;
  double HSM2_lmuetmp ;
  double HSM2_lsub1 ;
  double HSM2_lsub2 ;
  double HSM2_lsvds ;
  double HSM2_lsvbs ;
  double HSM2_lsvgs ;
  double HSM2_lnsti ;
  double HSM2_lwsti ;
  double HSM2_lscsti1 ;
  double HSM2_lscsti2 ;
  double HSM2_lvthsti ;
  double HSM2_lmuesti1 ;
  double HSM2_lmuesti2 ;
  double HSM2_lmuesti3 ;
  double HSM2_lnsubpsti1 ;
  double HSM2_lnsubpsti2 ;
  double HSM2_lnsubcsti1;
  double HSM2_lnsubcsti2;
  double HSM2_lnsubcsti3;
  double HSM2_lnsubpsti3 ;
  double HSM2_lcgso ;
  double HSM2_lcgdo ;
  double HSM2_ljs0 ;
  double HSM2_ljs0sw ;
  double HSM2_lnj ;
  double HSM2_lcisbk ;
  double HSM2_lclm1 ;
  double HSM2_lclm2 ;
  double HSM2_lclm3 ;
  double HSM2_lwfc ;
  double HSM2_lgidl1 ;
  double HSM2_lgidl2 ;
  double HSM2_lgleak1 ;
  double HSM2_lgleak2 ;
  double HSM2_lgleak3 ;
  double HSM2_lgleak6 ;
  double HSM2_lglksd1 ;
  double HSM2_lglksd2 ;
  double HSM2_lglkb1 ;
  double HSM2_lglkb2 ;
  double HSM2_lnftrp ;
  double HSM2_lnfalp ;
  double HSM2_lvdiffj ;
  double HSM2_libpc1 ;
  double HSM2_libpc2 ;

  /* Width dependence */
  double HSM2_wvmax ;
  double HSM2_wbgtmp1 ;
  double HSM2_wbgtmp2 ;
  double HSM2_weg0 ;
  double HSM2_wlover ;
  double HSM2_wvfbover ;
  double HSM2_wnover ;
  double HSM2_wwl2 ;
  double HSM2_wvfbc ;
  double HSM2_wnsubc ;
  double HSM2_wnsubp ;
  double HSM2_wscp1 ;
  double HSM2_wscp2 ;
  double HSM2_wscp3 ;
  double HSM2_wsc1 ;
  double HSM2_wsc2 ;
  double HSM2_wsc3 ;
  double HSM2_wsc4 ;
  double HSM2_wpgd1 ;
//double HSM2_wpgd3 ;
  double HSM2_wndep ;
  double HSM2_wninv ;
  double HSM2_wmuecb0 ;
  double HSM2_wmuecb1 ;
  double HSM2_wmueph1 ;
  double HSM2_wvtmp ;
  double HSM2_wwvth0 ;
  double HSM2_wmuesr1 ;
  double HSM2_wmuetmp ;
  double HSM2_wsub1 ;
  double HSM2_wsub2 ;
  double HSM2_wsvds ;
  double HSM2_wsvbs ;
  double HSM2_wsvgs ;
  double HSM2_wnsti ;
  double HSM2_wwsti ;
  double HSM2_wscsti1 ;
  double HSM2_wscsti2 ;
  double HSM2_wvthsti ;
  double HSM2_wmuesti1 ;
  double HSM2_wmuesti2 ;
  double HSM2_wmuesti3 ;
  double HSM2_wnsubpsti1 ;
  double HSM2_wnsubpsti2 ;
  double HSM2_wnsubcsti1;
  double HSM2_wnsubcsti2;
  double HSM2_wnsubcsti3;
  double HSM2_wnsubpsti3 ;
  double HSM2_wcgso ;
  double HSM2_wcgdo ;
  double HSM2_wjs0 ;
  double HSM2_wjs0sw ;
  double HSM2_wnj ;
  double HSM2_wcisbk ;
  double HSM2_wclm1 ;
  double HSM2_wclm2 ;
  double HSM2_wclm3 ;
  double HSM2_wwfc ;
  double HSM2_wgidl1 ;
  double HSM2_wgidl2 ;
  double HSM2_wgleak1 ;
  double HSM2_wgleak2 ;
  double HSM2_wgleak3 ;
  double HSM2_wgleak6 ;
  double HSM2_wglksd1 ;
  double HSM2_wglksd2 ;
  double HSM2_wglkb1 ;
  double HSM2_wglkb2 ;
  double HSM2_wnftrp ;
  double HSM2_wnfalp ;
  double HSM2_wvdiffj ;
  double HSM2_wibpc1 ;
  double HSM2_wibpc2 ;

  /* Cross-term dependence */
  double HSM2_pvmax ;
  double HSM2_pbgtmp1 ;
  double HSM2_pbgtmp2 ;
  double HSM2_peg0 ;
  double HSM2_plover ;
  double HSM2_pvfbover ;
  double HSM2_pnover ;
  double HSM2_pwl2 ;
  double HSM2_pvfbc ;
  double HSM2_pnsubc ;
  double HSM2_pnsubp ;
  double HSM2_pscp1 ;
  double HSM2_pscp2 ;
  double HSM2_pscp3 ;
  double HSM2_psc1 ;
  double HSM2_psc2 ;
  double HSM2_psc3 ;
  double HSM2_psc4 ;
  double HSM2_ppgd1 ;
//double HSM2_ppgd3 ;
  double HSM2_pndep ;
  double HSM2_pninv ;
  double HSM2_pmuecb0 ;
  double HSM2_pmuecb1 ;
  double HSM2_pmueph1 ;
  double HSM2_pvtmp ;
  double HSM2_pwvth0 ;
  double HSM2_pmuesr1 ;
  double HSM2_pmuetmp ;
  double HSM2_psub1 ;
  double HSM2_psub2 ;
  double HSM2_psvds ;
  double HSM2_psvbs ;
  double HSM2_psvgs ;
  double HSM2_pnsti ;
  double HSM2_pwsti ;
  double HSM2_pscsti1 ;
  double HSM2_pscsti2 ;
  double HSM2_pvthsti ;
  double HSM2_pmuesti1 ;
  double HSM2_pmuesti2 ;
  double HSM2_pmuesti3 ;
  double HSM2_pnsubpsti1 ;
  double HSM2_pnsubpsti2 ;
  double HSM2_pnsubcsti1;
  double HSM2_pnsubcsti2;
  double HSM2_pnsubcsti3;
  double HSM2_pnsubpsti3 ;
  double HSM2_pcgso ;
  double HSM2_pcgdo ;
  double HSM2_pjs0 ;
  double HSM2_pjs0sw ;
  double HSM2_pnj ;
  double HSM2_pcisbk ;
  double HSM2_pclm1 ;
  double HSM2_pclm2 ;
  double HSM2_pclm3 ;
  double HSM2_pwfc ;
  double HSM2_pgidl1 ;
  double HSM2_pgidl2 ;
  double HSM2_pgleak1 ;
  double HSM2_pgleak2 ;
  double HSM2_pgleak3 ;
  double HSM2_pgleak6 ;
  double HSM2_pglksd1 ;
  double HSM2_pglksd2 ;
  double HSM2_pglkb1 ;
  double HSM2_pglkb2 ;
  double HSM2_pnftrp ;
  double HSM2_pnfalp ;
  double HSM2_pvdiffj ;
  double HSM2_pibpc1 ;
  double HSM2_pibpc2 ;

  /* internal variables */
  double HSM2_vcrit ;
  int HSM2_flg_qme ;
  double HSM2_qme12 ;
  double HSM2_ktnom ;
  int HSM2_bypass_enable ;

  double HSM2vgsMax;
  double HSM2vgdMax;
  double HSM2vgbMax;
  double HSM2vdsMax;
  double HSM2vbsMax;
  double HSM2vbdMax;

  HSM2modelMKSParam modelMKS ; /* unit-converted parameters */

#ifdef USE_OMP
    int HSM2InstCount;
    struct sHSM2instance **HSM2InstanceArray;
#endif

  /* flag for model */
  unsigned HSM2_type_Given  :1;
  unsigned HSM2_level_Given  :1;
  unsigned HSM2_info_Given  :1;
  unsigned HSM2_noise_Given :1;
  unsigned HSM2_version_Given :1;
  unsigned HSM2_show_Given :1;
  unsigned HSM2_corsrd_Given  :1;
  unsigned HSM2_corg_Given    :1;
  unsigned HSM2_coiprv_Given  :1;
  unsigned HSM2_copprv_Given  :1;
  unsigned HSM2_coadov_Given  :1;
  unsigned HSM2_coisub_Given  :1;
  unsigned HSM2_coiigs_Given  :1;
  unsigned HSM2_cogidl_Given  :1;
  unsigned HSM2_coovlp_Given  :1;
  unsigned HSM2_coflick_Given  :1;
  unsigned HSM2_coisti_Given  :1;
  unsigned HSM2_conqs_Given  :1;
  unsigned HSM2_corbnet_Given  :1;
  unsigned HSM2_cothrml_Given  :1;
  unsigned HSM2_coign_Given  :1; /* induced gate noise */
  unsigned HSM2_codfm_Given  :1; /* DFM */
  unsigned HSM2_corecip_Given  :1;
  unsigned HSM2_coqy_Given  :1;
  unsigned HSM2_coqovsm_Given  :1;
  unsigned HSM2_coerrrep_Given :1;
  unsigned HSM2_kappa_Given :1;  
  unsigned HSM2_vdiffj_Given :1; 
  unsigned HSM2_vmax_Given  :1;
  unsigned HSM2_bgtmp1_Given  :1;
  unsigned HSM2_bgtmp2_Given  :1;
  unsigned HSM2_eg0_Given  :1;
  unsigned HSM2_tox_Given  :1;
  unsigned HSM2_xld_Given  :1;
  unsigned HSM2_lover_Given  :1;
  unsigned HSM2_ddltmax_Given  :1; /* Vdseff */
  unsigned HSM2_ddltslp_Given  :1; /* Vdseff */
  unsigned HSM2_ddltict_Given  :1; /* Vdseff */
  unsigned HSM2_vfbover_Given  :1;
  unsigned HSM2_nover_Given  :1;
  unsigned HSM2_xwd_Given  :1; 
  unsigned HSM2_xl_Given  :1;
  unsigned HSM2_xw_Given  :1;
  unsigned HSM2_saref_Given  :1;
  unsigned HSM2_sbref_Given  :1;
  unsigned HSM2_ll_Given  :1;
  unsigned HSM2_lld_Given  :1;
  unsigned HSM2_lln_Given  :1;
  unsigned HSM2_wl_Given  :1;
  unsigned HSM2_wl1_Given  :1;
  unsigned HSM2_wl1p_Given  :1;
  unsigned HSM2_wl2_Given  :1;
  unsigned HSM2_wl2p_Given  :1;
  unsigned HSM2_wld_Given  :1;
  unsigned HSM2_wln_Given  :1; 
  unsigned HSM2_xqy_Given  :1;   
  unsigned HSM2_xqy1_Given  :1;   
  unsigned HSM2_xqy2_Given  :1;   
  unsigned HSM2_qyrat_Given  :1;   
  unsigned HSM2_rs_Given  :1;
  unsigned HSM2_rd_Given  :1;
  unsigned HSM2_rsh_Given  :1;
  unsigned HSM2_rshg_Given  :1;
/*   unsigned HSM2_ngcon_Given  :1; */
/*   unsigned HSM2_xgw_Given  :1; */
/*   unsigned HSM2_xgl_Given  :1; */
/*   unsigned HSM2_nf_Given  :1; */
  unsigned HSM2_vfbc_Given  :1;
  unsigned HSM2_vbi_Given  :1;
  unsigned HSM2_vfbcl_Given :1;
  unsigned HSM2_vfbclp_Given :1;
  unsigned HSM2_nsubc_Given  :1;
  unsigned HSM2_parl2_Given  :1;
  unsigned HSM2_lp_Given  :1;
  unsigned HSM2_nsubp_Given  :1;
  unsigned HSM2_nsubpl_Given  :1;
  unsigned HSM2_nsubpdlt_Given :1;
  unsigned HSM2_nsubpfac_Given  :1;
  unsigned HSM2_nsubpw_Given  :1;
  unsigned HSM2_nsubpwp_Given  :1;
  unsigned HSM2_scp1_Given  :1;
  unsigned HSM2_scp2_Given  :1;
  unsigned HSM2_scp3_Given  :1;
  unsigned HSM2_sc1_Given  :1;
  unsigned HSM2_sc2_Given  :1;
  unsigned HSM2_sc3_Given  :1;
  unsigned HSM2_sc4_Given  :1;
  unsigned HSM2_pgd1_Given  :1;
  unsigned HSM2_pgd2_Given  :1;
//unsigned HSM2_pgd3_Given  :1;
  unsigned HSM2_pgd4_Given  :1;
  unsigned HSM2_ndep_Given  :1;
  unsigned HSM2_ndepl_Given  :1;
  unsigned HSM2_ndeplp_Given  :1;
  unsigned HSM2_ndepw_Given  :1;
  unsigned HSM2_ndepwp_Given  :1;
  unsigned HSM2_ninv_Given  :1;
  unsigned HSM2_ninvd_Given  :1;
  unsigned HSM2_muecb0_Given  :1;
  unsigned HSM2_muecb1_Given  :1;
  unsigned HSM2_mueph1_Given  :1;
  unsigned HSM2_mueph0_Given  :1;
  unsigned HSM2_muephw_Given  :1;
  unsigned HSM2_muepwp_Given  :1;
  unsigned HSM2_muepwd_Given  :1;
  unsigned HSM2_muephl_Given  :1;
  unsigned HSM2_mueplp_Given  :1;
  unsigned HSM2_muepld_Given  :1;
  unsigned HSM2_muephs_Given  :1;
  unsigned HSM2_muepsp_Given  :1;
  unsigned HSM2_vtmp_Given  :1;
  unsigned HSM2_wvth0_Given  :1;
  unsigned HSM2_muesr1_Given  :1;
  unsigned HSM2_muesr0_Given  :1;
  unsigned HSM2_muesrl_Given  :1;
  unsigned HSM2_mueslp_Given  :1;
  unsigned HSM2_muesrw_Given  :1;
  unsigned HSM2_mueswp_Given  :1;
  unsigned HSM2_bb_Given  :1;
  unsigned HSM2_sub1_Given  :1;
  unsigned HSM2_sub2_Given  :1;
  unsigned HSM2_svgs_Given  :1; 
  unsigned HSM2_svbs_Given  :1; 
  unsigned HSM2_svbsl_Given  :1; 
  unsigned HSM2_svds_Given  :1; 
  unsigned HSM2_slg_Given  :1; 
  unsigned HSM2_sub1l_Given  :1; 
  unsigned HSM2_sub2l_Given  :1;  
  unsigned HSM2_svgsl_Given  :1;
  unsigned HSM2_svgslp_Given  :1;
  unsigned HSM2_svgswp_Given  :1;
  unsigned HSM2_svgsw_Given  :1;
  unsigned HSM2_svbslp_Given  :1;
  unsigned HSM2_slgl_Given  :1;
  unsigned HSM2_slglp_Given  :1;
  unsigned HSM2_sub1lp_Given  :1;
  unsigned HSM2_nsti_Given  :1;  
  unsigned HSM2_wsti_Given  :1;
  unsigned HSM2_wstil_Given  :1;
  unsigned HSM2_wstilp_Given  :1;
  unsigned HSM2_wstiw_Given  :1;
  unsigned HSM2_wstiwp_Given  :1;
  unsigned HSM2_scsti1_Given  :1;
  unsigned HSM2_scsti2_Given  :1;
  unsigned HSM2_vthsti_Given  :1;
  unsigned HSM2_vdsti_Given  :1;
  unsigned HSM2_muesti1_Given  :1;
  unsigned HSM2_muesti2_Given  :1;
  unsigned HSM2_muesti3_Given  :1;
  unsigned HSM2_nsubpsti1_Given  :1;
  unsigned HSM2_nsubpsti2_Given  :1;
  unsigned HSM2_nsubcsti1_Given :1;
  unsigned HSM2_nsubcsti2_Given :1;
  unsigned HSM2_nsubcsti3_Given :1;
  unsigned HSM2_nsubpsti3_Given  :1;
  unsigned HSM2_lpext_Given  :1;
  unsigned HSM2_npext_Given  :1;
  unsigned HSM2_npextw_Given  :1;
  unsigned HSM2_npextwp_Given  :1;
  unsigned HSM2_scp22_Given  :1;
  unsigned HSM2_scp21_Given  :1;
  unsigned HSM2_bs1_Given  :1;
  unsigned HSM2_bs2_Given  :1;
  unsigned HSM2_cgso_Given  :1;
  unsigned HSM2_cgdo_Given  :1;
  unsigned HSM2_cgbo_Given  :1;
  unsigned HSM2_tpoly_Given  :1;
  unsigned HSM2_js0_Given  :1;
  unsigned HSM2_js0sw_Given  :1;
  unsigned HSM2_nj_Given  :1;
  unsigned HSM2_njsw_Given  :1;  
  unsigned HSM2_xti_Given  :1;
  unsigned HSM2_cj_Given  :1;
  unsigned HSM2_cjsw_Given  :1;
  unsigned HSM2_cjswg_Given  :1;
  unsigned HSM2_mj_Given  :1;
  unsigned HSM2_mjsw_Given  :1;
  unsigned HSM2_mjswg_Given  :1;
  unsigned HSM2_xti2_Given  :1;
  unsigned HSM2_cisb_Given  :1;
  unsigned HSM2_cvb_Given  :1;
  unsigned HSM2_ctemp_Given  :1;
  unsigned HSM2_cisbk_Given  :1;
  unsigned HSM2_cvbk_Given  :1;
  unsigned HSM2_divx_Given  :1;
  unsigned HSM2_pb_Given  :1;
  unsigned HSM2_pbsw_Given  :1;
  unsigned HSM2_pbswg_Given  :1;
  unsigned HSM2_tcjbd_Given :1;
  unsigned HSM2_tcjbs_Given :1;
  unsigned HSM2_tcjbdsw_Given :1;
  unsigned HSM2_tcjbssw_Given :1;
  unsigned HSM2_tcjbdswg_Given :1;
  unsigned HSM2_tcjbsswg_Given :1;

  unsigned HSM2_clm1_Given  :1;
  unsigned HSM2_clm2_Given  :1;
  unsigned HSM2_clm3_Given  :1;
  unsigned HSM2_clm5_Given  :1;
  unsigned HSM2_clm6_Given  :1;
  unsigned HSM2_muetmp_Given  :1;
  unsigned HSM2_vover_Given  :1;
  unsigned HSM2_voverp_Given  :1;
  unsigned HSM2_vovers_Given  :1;
  unsigned HSM2_voversp_Given  :1;
  unsigned HSM2_wfc_Given  :1;
  unsigned HSM2_nsubcw_Given  :1;
  unsigned HSM2_nsubcwp_Given  :1;
  unsigned HSM2_nsubcmax_Given  :1;
  unsigned HSM2_qme1_Given  :1;
  unsigned HSM2_qme2_Given  :1;
  unsigned HSM2_qme3_Given  :1;
  unsigned HSM2_gidl1_Given  :1;
  unsigned HSM2_gidl2_Given  :1;
  unsigned HSM2_gidl3_Given  :1;
  unsigned HSM2_gidl4_Given  :1;
  unsigned HSM2_gidl6_Given :1;
  unsigned HSM2_gidl7_Given :1;
  unsigned HSM2_gidl5_Given  :1;
  unsigned HSM2_gleak1_Given  :1;
  unsigned HSM2_gleak2_Given  :1;
  unsigned HSM2_gleak3_Given  :1;
  unsigned HSM2_gleak4_Given  :1;
  unsigned HSM2_gleak5_Given  :1;
  unsigned HSM2_gleak6_Given  :1;
  unsigned HSM2_gleak7_Given  :1;
  unsigned HSM2_glksd1_Given  :1;
  unsigned HSM2_glksd2_Given  :1;
  unsigned HSM2_glksd3_Given  :1;
  unsigned HSM2_glkb1_Given  :1;
  unsigned HSM2_glkb2_Given  :1;
  unsigned HSM2_glkb3_Given  :1;
  unsigned HSM2_egig_Given  :1;
  unsigned HSM2_igtemp2_Given  :1;
  unsigned HSM2_igtemp3_Given  :1;
  unsigned HSM2_vzadd0_Given  :1;
  unsigned HSM2_pzadd0_Given  :1;
  unsigned HSM2_nftrp_Given  :1;
  unsigned HSM2_nfalp_Given  :1;
  unsigned HSM2_cit_Given  :1;
  unsigned HSM2_falph_Given  :1;
  unsigned HSM2_dly1_Given :1;
  unsigned HSM2_dly2_Given :1;
  unsigned HSM2_dly3_Given :1;
  unsigned HSM2_tnom_Given :1;
  unsigned HSM2_ovslp_Given :1;
  unsigned HSM2_ovmag_Given :1;
  unsigned HSM2_gbmin_Given :1;
  unsigned HSM2_rbpb_Given :1;
  unsigned HSM2_rbpd_Given :1;
  unsigned HSM2_rbps_Given :1;
  unsigned HSM2_rbdb_Given :1;
  unsigned HSM2_rbsb_Given :1;
  unsigned HSM2_ibpc1_Given :1;
  unsigned HSM2_ibpc2_Given :1;
  unsigned HSM2_mphdfm_Given :1;

  unsigned HSM2_ptl_Given :1;
  unsigned HSM2_ptp_Given :1;
  unsigned HSM2_pt2_Given :1;
  unsigned HSM2_ptlp_Given :1;
  unsigned HSM2_gdl_Given :1;
  unsigned HSM2_gdlp_Given :1;

  unsigned HSM2_gdld_Given :1;
  unsigned HSM2_pt4_Given :1;
  unsigned HSM2_pt4p_Given :1;
  unsigned HSM2_muephl2_Given :1;
  unsigned HSM2_mueplp2_Given :1;
  unsigned HSM2_nsubcw2_Given :1;
  unsigned HSM2_nsubcwp2_Given :1;
  unsigned HSM2_muephw2_Given :1;
  unsigned HSM2_muepwp2_Given :1;

  /* val set flag for WPE */
  unsigned HSM2_web_Given :1;
  unsigned HSM2_wec_Given :1;
  unsigned HSM2_nsubcwpe_Given :1;
  unsigned HSM2_npextwpe_Given :1;
  unsigned HSM2_nsubpwpe_Given :1;
  unsigned HSM2_Vgsmin_Given :1;
  unsigned HSM2_sc3Vbs_Given :1;
  unsigned HSM2_byptol_Given :1;
  unsigned HSM2_muecb0lp_Given :1;
  unsigned HSM2_muecb1lp_Given :1;

  /* binning parameters */
  unsigned HSM2_lmin_Given :1;
  unsigned HSM2_lmax_Given :1;
  unsigned HSM2_wmin_Given :1;
  unsigned HSM2_wmax_Given :1;
  unsigned HSM2_lbinn_Given :1;
  unsigned HSM2_wbinn_Given :1;

  /* Length dependence */
  unsigned HSM2_lvmax_Given :1;
  unsigned HSM2_lbgtmp1_Given :1;
  unsigned HSM2_lbgtmp2_Given :1;
  unsigned HSM2_leg0_Given :1;
  unsigned HSM2_llover_Given :1;
  unsigned HSM2_lvfbover_Given :1;
  unsigned HSM2_lnover_Given :1;
  unsigned HSM2_lwl2_Given :1;
  unsigned HSM2_lvfbc_Given :1;
  unsigned HSM2_lnsubc_Given :1;
  unsigned HSM2_lnsubp_Given :1;
  unsigned HSM2_lscp1_Given :1;
  unsigned HSM2_lscp2_Given :1;
  unsigned HSM2_lscp3_Given :1;
  unsigned HSM2_lsc1_Given :1;
  unsigned HSM2_lsc2_Given :1;
  unsigned HSM2_lsc3_Given :1;
  unsigned HSM2_lsc4_Given :1;
  unsigned HSM2_lpgd1_Given :1;
//unsigned HSM2_lpgd3_Given :1;
  unsigned HSM2_lndep_Given :1;
  unsigned HSM2_lninv_Given :1;
  unsigned HSM2_lmuecb0_Given :1;
  unsigned HSM2_lmuecb1_Given :1;
  unsigned HSM2_lmueph1_Given :1;
  unsigned HSM2_lvtmp_Given :1;
  unsigned HSM2_lwvth0_Given :1;
  unsigned HSM2_lmuesr1_Given :1;
  unsigned HSM2_lmuetmp_Given :1;
  unsigned HSM2_lsub1_Given :1;
  unsigned HSM2_lsub2_Given :1;
  unsigned HSM2_lsvds_Given :1;
  unsigned HSM2_lsvbs_Given :1;
  unsigned HSM2_lsvgs_Given :1;
  unsigned HSM2_lnsti_Given :1;
  unsigned HSM2_lwsti_Given :1;
  unsigned HSM2_lscsti1_Given :1;
  unsigned HSM2_lscsti2_Given :1;
  unsigned HSM2_lvthsti_Given :1;
  unsigned HSM2_lmuesti1_Given :1;
  unsigned HSM2_lmuesti2_Given :1;
  unsigned HSM2_lmuesti3_Given :1;
  unsigned HSM2_lnsubpsti1_Given :1;
  unsigned HSM2_lnsubpsti2_Given :1;
  unsigned HSM2_lnsubcsti1_Given :1;
  unsigned HSM2_lnsubcsti2_Given :1;
  unsigned HSM2_lnsubcsti3_Given :1;
  unsigned HSM2_lnsubpsti3_Given :1;
  unsigned HSM2_lcgso_Given :1;
  unsigned HSM2_lcgdo_Given :1;
  unsigned HSM2_ljs0_Given :1;
  unsigned HSM2_ljs0sw_Given :1;
  unsigned HSM2_lnj_Given :1;
  unsigned HSM2_lcisbk_Given :1;
  unsigned HSM2_lclm1_Given :1;
  unsigned HSM2_lclm2_Given :1;
  unsigned HSM2_lclm3_Given :1;
  unsigned HSM2_lwfc_Given :1;
  unsigned HSM2_lgidl1_Given :1;
  unsigned HSM2_lgidl2_Given :1;
  unsigned HSM2_lgleak1_Given :1;
  unsigned HSM2_lgleak2_Given :1;
  unsigned HSM2_lgleak3_Given :1;
  unsigned HSM2_lgleak6_Given :1;
  unsigned HSM2_lglksd1_Given :1;
  unsigned HSM2_lglksd2_Given :1;
  unsigned HSM2_lglkb1_Given :1;
  unsigned HSM2_lglkb2_Given :1;
  unsigned HSM2_lnftrp_Given :1;
  unsigned HSM2_lnfalp_Given :1;
  unsigned HSM2_lvdiffj_Given :1;
  unsigned HSM2_libpc1_Given :1;
  unsigned HSM2_libpc2_Given :1;

  /* Width dependence */
  unsigned HSM2_wvmax_Given :1;
  unsigned HSM2_wbgtmp1_Given :1;
  unsigned HSM2_wbgtmp2_Given :1;
  unsigned HSM2_weg0_Given :1;
  unsigned HSM2_wlover_Given :1;
  unsigned HSM2_wvfbover_Given :1;
  unsigned HSM2_wnover_Given :1;
  unsigned HSM2_wwl2_Given :1;
  unsigned HSM2_wvfbc_Given :1;
  unsigned HSM2_wnsubc_Given :1;
  unsigned HSM2_wnsubp_Given :1;
  unsigned HSM2_wscp1_Given :1;
  unsigned HSM2_wscp2_Given :1;
  unsigned HSM2_wscp3_Given :1;
  unsigned HSM2_wsc1_Given :1;
  unsigned HSM2_wsc2_Given :1;
  unsigned HSM2_wsc3_Given :1;
  unsigned HSM2_wsc4_Given :1;
  unsigned HSM2_wpgd1_Given :1;
//unsigned HSM2_wpgd3_Given :1;
  unsigned HSM2_wndep_Given :1;
  unsigned HSM2_wninv_Given :1;
  unsigned HSM2_wmuecb0_Given :1;
  unsigned HSM2_wmuecb1_Given :1;
  unsigned HSM2_wmueph1_Given :1;
  unsigned HSM2_wvtmp_Given :1;
  unsigned HSM2_wwvth0_Given :1;
  unsigned HSM2_wmuesr1_Given :1;
  unsigned HSM2_wmuetmp_Given :1;
  unsigned HSM2_wsub1_Given :1;
  unsigned HSM2_wsub2_Given :1;
  unsigned HSM2_wsvds_Given :1;
  unsigned HSM2_wsvbs_Given :1;
  unsigned HSM2_wsvgs_Given :1;
  unsigned HSM2_wnsti_Given :1;
  unsigned HSM2_wwsti_Given :1;
  unsigned HSM2_wscsti1_Given :1;
  unsigned HSM2_wscsti2_Given :1;
  unsigned HSM2_wvthsti_Given :1;
  unsigned HSM2_wmuesti1_Given :1;
  unsigned HSM2_wmuesti2_Given :1;
  unsigned HSM2_wmuesti3_Given :1;
  unsigned HSM2_wnsubpsti1_Given :1;
  unsigned HSM2_wnsubpsti2_Given :1;
  unsigned HSM2_wnsubcsti1_Given :1;
  unsigned HSM2_wnsubcsti2_Given :1;
  unsigned HSM2_wnsubcsti3_Given :1;
  unsigned HSM2_wnsubpsti3_Given :1;
  unsigned HSM2_wcgso_Given :1;
  unsigned HSM2_wcgdo_Given :1;
  unsigned HSM2_wjs0_Given :1;
  unsigned HSM2_wjs0sw_Given :1;
  unsigned HSM2_wnj_Given :1;
  unsigned HSM2_wcisbk_Given :1;
  unsigned HSM2_wclm1_Given :1;
  unsigned HSM2_wclm2_Given :1;
  unsigned HSM2_wclm3_Given :1;
  unsigned HSM2_wwfc_Given :1;
  unsigned HSM2_wgidl1_Given :1;
  unsigned HSM2_wgidl2_Given :1;
  unsigned HSM2_wgleak1_Given :1;
  unsigned HSM2_wgleak2_Given :1;
  unsigned HSM2_wgleak3_Given :1;
  unsigned HSM2_wgleak6_Given :1;
  unsigned HSM2_wglksd1_Given :1;
  unsigned HSM2_wglksd2_Given :1;
  unsigned HSM2_wglkb1_Given :1;
  unsigned HSM2_wglkb2_Given :1;
  unsigned HSM2_wnftrp_Given :1;
  unsigned HSM2_wnfalp_Given :1;
  unsigned HSM2_wvdiffj_Given :1;
  unsigned HSM2_wibpc1_Given :1;
  unsigned HSM2_wibpc2_Given :1;

  /* Cross-term dependence */
  unsigned HSM2_pvmax_Given :1;
  unsigned HSM2_pbgtmp1_Given :1;
  unsigned HSM2_pbgtmp2_Given :1;
  unsigned HSM2_peg0_Given :1;
  unsigned HSM2_plover_Given :1;
  unsigned HSM2_pvfbover_Given :1;
  unsigned HSM2_pnover_Given :1;
  unsigned HSM2_pwl2_Given :1;
  unsigned HSM2_pvfbc_Given :1;
  unsigned HSM2_pnsubc_Given :1;
  unsigned HSM2_pnsubp_Given :1;
  unsigned HSM2_pscp1_Given :1;
  unsigned HSM2_pscp2_Given :1;
  unsigned HSM2_pscp3_Given :1;
  unsigned HSM2_psc1_Given :1;
  unsigned HSM2_psc2_Given :1;
  unsigned HSM2_psc3_Given :1;
  unsigned HSM2_psc4_Given :1;
  unsigned HSM2_ppgd1_Given :1;
//unsigned HSM2_ppgd3_Given :1;
  unsigned HSM2_pndep_Given :1;
  unsigned HSM2_pninv_Given :1;
  unsigned HSM2_pmuecb0_Given :1;
  unsigned HSM2_pmuecb1_Given :1;
  unsigned HSM2_pmueph1_Given :1;
  unsigned HSM2_pvtmp_Given :1;
  unsigned HSM2_pwvth0_Given :1;
  unsigned HSM2_pmuesr1_Given :1;
  unsigned HSM2_pmuetmp_Given :1;
  unsigned HSM2_psub1_Given :1;
  unsigned HSM2_psub2_Given :1;
  unsigned HSM2_psvds_Given :1;
  unsigned HSM2_psvbs_Given :1;
  unsigned HSM2_psvgs_Given :1;
  unsigned HSM2_pnsti_Given :1;
  unsigned HSM2_pwsti_Given :1;
  unsigned HSM2_pscsti1_Given :1;
  unsigned HSM2_pscsti2_Given :1;
  unsigned HSM2_pvthsti_Given :1;
  unsigned HSM2_pmuesti1_Given :1;
  unsigned HSM2_pmuesti2_Given :1;
  unsigned HSM2_pmuesti3_Given :1;
  unsigned HSM2_pnsubpsti1_Given :1;
  unsigned HSM2_pnsubpsti2_Given :1;
  unsigned HSM2_pnsubcsti1_Given :1;
  unsigned HSM2_pnsubcsti2_Given :1;
  unsigned HSM2_pnsubcsti3_Given :1;
  unsigned HSM2_pnsubpsti3_Given :1;
  unsigned HSM2_pcgso_Given :1;
  unsigned HSM2_pcgdo_Given :1;
  unsigned HSM2_pjs0_Given :1;
  unsigned HSM2_pjs0sw_Given :1;
  unsigned HSM2_pnj_Given :1;
  unsigned HSM2_pcisbk_Given :1;
  unsigned HSM2_pclm1_Given :1;
  unsigned HSM2_pclm2_Given :1;
  unsigned HSM2_pclm3_Given :1;
  unsigned HSM2_pwfc_Given :1;
  unsigned HSM2_pgidl1_Given :1;
  unsigned HSM2_pgidl2_Given :1;
  unsigned HSM2_pgleak1_Given :1;
  unsigned HSM2_pgleak2_Given :1;
  unsigned HSM2_pgleak3_Given :1;
  unsigned HSM2_pgleak6_Given :1;
  unsigned HSM2_pglksd1_Given :1;
  unsigned HSM2_pglksd2_Given :1;
  unsigned HSM2_pglkb1_Given :1;
  unsigned HSM2_pglkb2_Given :1;
  unsigned HSM2_pnftrp_Given :1;
  unsigned HSM2_pnfalp_Given :1;
  unsigned HSM2_pvdiffj_Given :1;
  unsigned HSM2_pibpc1_Given :1;
  unsigned HSM2_pibpc2_Given :1;

  unsigned  HSM2vgsMaxGiven  :1;
  unsigned  HSM2vgdMaxGiven  :1;
  unsigned  HSM2vgbMaxGiven  :1;
  unsigned  HSM2vdsMaxGiven  :1;
  unsigned  HSM2vbsMaxGiven  :1;
  unsigned  HSM2vbdMaxGiven  :1;

} HSM2model;

#ifndef NMOS
#define NMOS 1
#define PMOS -1
#endif /*NMOS*/

#define HSM2_BAD_PARAM -1

/* flags */
#define HSM2_MOD_NMOS     1
#define HSM2_MOD_PMOS     2
#define HSM2_MOD_LEVEL    3
#define HSM2_MOD_INFO     4
#define HSM2_MOD_NOISE    5
#define HSM2_MOD_VERSION  6
#define HSM2_MOD_SHOW     7
#define HSM2_MOD_CORSRD  11
#define HSM2_MOD_COIPRV  12
#define HSM2_MOD_COPPRV  13
#define HSM2_MOD_COADOV  17
#define HSM2_MOD_COISUB  21
#define HSM2_MOD_COIIGS    22
#define HSM2_MOD_COGIDL 23
#define HSM2_MOD_COOVLP  24
#define HSM2_MOD_COFLICK 25
#define HSM2_MOD_COISTI  26
#define HSM2_MOD_CONQS   29 /* HiSIM2 */
#define HSM2_MOD_COTHRML 30
#define HSM2_MOD_COIGN   31 /* induced gate noise */
#define HSM2_MOD_CORG    32
#define HSM2_MOD_CORBNET 33
#define HSM2_MOD_CODFM   36 /* DFM */
#define HSM2_MOD_CORECIP 37
#define HSM2_MOD_COQY    38
#define HSM2_MOD_COQOVSM 39
#define HSM2_MOD_COERRREP     153
/* device parameters */
#define HSM2_L           51
#define HSM2_W           52
#define HSM2_AD          53
#define HSM2_AS          54
#define HSM2_PD          55
#define HSM2_PS          56
#define HSM2_NRD         57
#define HSM2_NRS         58
#define HSM2_TEMP        59
#define HSM2_DTEMP       60
#define HSM2_OFF         61
#define HSM2_IC_VBS      62
#define HSM2_IC_VDS      63
#define HSM2_IC_VGS      64
#define HSM2_IC          65
#define HSM2_CORBNET     66
#define HSM2_RBPB        67
#define HSM2_RBPD        68
#define HSM2_RBPS        69
#define HSM2_RBDB        70
#define HSM2_RBSB        71
#define HSM2_CORG        72
/* #define HSM2_RSHG        73 */
#define HSM2_NGCON       74
#define HSM2_XGW         75 
#define HSM2_XGL         76
#define HSM2_NF          77
#define HSM2_SA          78
#define HSM2_SB          79
#define HSM2_SD          80
#define HSM2_NSUBCDFM    82
#define HSM2_MPHDFM      84
#define HSM2_M           83

/* val symbol for WPE */
#define HSM2_SCA	 85	/* sca */
#define HSM2_SCB	 86	/* scb */
#define HSM2_SCC	 87	/* scc */

/* model parameters */
#define HSM2_MOD_VMAX      100
#define HSM2_MOD_BGTMP1    101
#define HSM2_MOD_BGTMP2    102
#define HSM2_MOD_EG0       103
#define HSM2_MOD_TOX       104
#define HSM2_MOD_XLD       105
#define HSM2_MOD_LOVER     106
#define HSM2_MOD_DDLTMAX   421 /* Vdseff */
#define HSM2_MOD_DDLTSLP   422 /* Vdseff */
#define HSM2_MOD_DDLTICT   423 /* Vdseff */
#define HSM2_MOD_VFBOVER   428
#define HSM2_MOD_NOVER     430
#define HSM2_MOD_XWD       107
#define HSM2_MOD_XL        112
#define HSM2_MOD_XW        117
#define HSM2_MOD_SAREF     433
#define HSM2_MOD_SBREF     434
#define HSM2_MOD_LL        108
#define HSM2_MOD_LLD       109
#define HSM2_MOD_LLN       110
#define HSM2_MOD_WL        111
#define HSM2_MOD_WL1       113
#define HSM2_MOD_WL1P      114
#define HSM2_MOD_WL2       407
#define HSM2_MOD_WL2P      408
#define HSM2_MOD_WLD       115
#define HSM2_MOD_WLN       116

#define HSM2_MOD_XQY       178
#define HSM2_MOD_XQY1      118
#define HSM2_MOD_XQY2      120
#define HSM2_MOD_QYRAT     991

#define HSM2_MOD_RSH       119
#define HSM2_MOD_RSHG      384
/* #define HSM2_MOD_NGCON     385 */
/* #define HSM2_MOD_XGW       386 */
/* #define HSM2_MOD_XGL       387 */
/* #define HSM2_MOD_NF        388 */
#define HSM2_MOD_RS        398
#define HSM2_MOD_RD        399

#define HSM2_MOD_VFBC      121
#define HSM2_MOD_VBI       122
#define HSM2_MOD_NSUBC     123
#define HSM2_MOD_VFBCL     272
#define HSM2_MOD_VFBCLP    273
#define HSM2_MOD_TNOM      124
#define HSM2_MOD_PARL2     125
#define HSM2_MOD_SC1       126
#define HSM2_MOD_SC2       127
#define HSM2_MOD_SC3       128
#define HSM2_MOD_SC4       460
#define HSM2_MOD_NDEP      129
#define HSM2_MOD_NDEPL     419
#define HSM2_MOD_NDEPLP    420
#define HSM2_MOD_NDEPW     469
#define HSM2_MOD_NDEPWP    470
#define HSM2_MOD_NINV      130
#define HSM2_MOD_NINVD     300
#define HSM2_MOD_MUECB0    131
#define HSM2_MOD_MUECB1    132
#define HSM2_MOD_MUEPH1    133
#define HSM2_MOD_MUEPH0    134
#define HSM2_MOD_MUEPHW    135
#define HSM2_MOD_MUEPWP    136
#define HSM2_MOD_MUEPWD    333
#define HSM2_MOD_MUEPHL    137
#define HSM2_MOD_MUEPLP    138
#define HSM2_MOD_MUEPLD    150
#define HSM2_MOD_MUEPHS    139
#define HSM2_MOD_MUEPSP    140
#define HSM2_MOD_VTMP      141
#define HSM2_MOD_WVTH0 	   142
#define HSM2_MOD_MUESR1    143
#define HSM2_MOD_MUESR0    144
#define HSM2_MOD_MUESRL    145
#define HSM2_MOD_MUESLP    146
#define HSM2_MOD_MUESRW    147
#define HSM2_MOD_MUESWP    148
#define HSM2_MOD_BB        149

#define HSM2_MOD_SUB1      151
#define HSM2_MOD_SUB2      152
#define HSM2_MOD_CGSO      154
#define HSM2_MOD_CGDO      155
#define HSM2_MOD_CGBO      156
#define HSM2_MOD_JS0       157
#define HSM2_MOD_JS0SW     158
#define HSM2_MOD_NJ        159
#define HSM2_MOD_NJSW      160
#define HSM2_MOD_XTI       161
#define HSM2_MOD_CJ        162
#define HSM2_MOD_CJSW      163
#define HSM2_MOD_CJSWG     164
#define HSM2_MOD_MJ        165
#define HSM2_MOD_MJSW      166
#define HSM2_MOD_MJSWG     167
#define HSM2_MOD_XTI2      168
#define HSM2_MOD_CISB      169
#define HSM2_MOD_CVB       170
#define HSM2_MOD_CTEMP     171
#define HSM2_MOD_CISBK     172
#define HSM2_MOD_CVBK      173
#define HSM2_MOD_DIVX      174
#define HSM2_MOD_PB        175
#define HSM2_MOD_PBSW      176
#define HSM2_MOD_PBSWG     177
#define HSM2_MOD_TPOLY     179
#define HSM2_MOD_LP        180
#define HSM2_MOD_NSUBP     181
#define HSM2_MOD_NSUBPL    196
#define HSM2_MOD_NSUBPFAC  197
#define HSM2_MOD_NSUBPDLT     274
#define HSM2_MOD_NSUBPW    182
#define HSM2_MOD_NSUBPWP   183
#define HSM2_MOD_SCP1      184
#define HSM2_MOD_SCP2      185
#define HSM2_MOD_SCP3      186
#define HSM2_MOD_PGD1      187
#define HSM2_MOD_PGD2      188
//#define HSM2_MOD_PGD3      189
#define HSM2_MOD_PGD4      190
#define HSM2_MOD_CLM1      191
#define HSM2_MOD_CLM2      192
#define HSM2_MOD_CLM3      193
#define HSM2_MOD_CLM5      402
#define HSM2_MOD_CLM6      403
#define HSM2_MOD_MUETMP    195

#define HSM2_MOD_VOVER     199
#define HSM2_MOD_VOVERP    200
#define HSM2_MOD_WFC       201
#define HSM2_MOD_NSUBCW    249
#define HSM2_MOD_NSUBCWP   250
#define HSM2_MOD_NSUBCMAX  248
#define HSM2_MOD_QME1      202
#define HSM2_MOD_QME2      203
#define HSM2_MOD_QME3      204
#define HSM2_MOD_GIDL1     205
#define HSM2_MOD_GIDL2     206
#define HSM2_MOD_GIDL3     207
#define HSM2_MOD_GLEAK1    208
#define HSM2_MOD_GLEAK2    209
#define HSM2_MOD_GLEAK3    210
#define HSM2_MOD_GLEAK4    211
#define HSM2_MOD_GLEAK5    212
#define HSM2_MOD_GLEAK6    213
#define HSM2_MOD_GLEAK7    214
#define HSM2_MOD_GLKSD1    215
#define HSM2_MOD_GLKSD2    216
#define HSM2_MOD_GLKSD3    217
#define HSM2_MOD_GLKB1     218
#define HSM2_MOD_GLKB2     219
#define HSM2_MOD_GLKB3     429
#define HSM2_MOD_EGIG      220
#define HSM2_MOD_IGTEMP2   221
#define HSM2_MOD_IGTEMP3   222
#define HSM2_MOD_VZADD0    223
#define HSM2_MOD_PZADD0    224
#define HSM2_MOD_NSTI      225
#define HSM2_MOD_WSTI      226
#define HSM2_MOD_WSTIL     227
#define HSM2_MOD_WSTILP    231
#define HSM2_MOD_WSTIW     234
#define HSM2_MOD_WSTIWP    228
#define HSM2_MOD_SCSTI1    229
#define HSM2_MOD_SCSTI2    230
#define HSM2_MOD_VTHSTI    232
#define HSM2_MOD_VDSTI     233
#define HSM2_MOD_MUESTI1   235
#define HSM2_MOD_MUESTI2   236
#define HSM2_MOD_MUESTI3   237
#define HSM2_MOD_NSUBPSTI1 238
#define HSM2_MOD_NSUBPSTI2 239
#define HSM2_MOD_NSUBPSTI3 240
#define HSM2_MOD_NSUBCSTI1    198
#define HSM2_MOD_NSUBCSTI2    247
#define HSM2_MOD_NSUBCSTI3    252
#define HSM2_MOD_LPEXT     241
#define HSM2_MOD_NPEXT     242
#define HSM2_MOD_NPEXTW    471
#define HSM2_MOD_NPEXTWP   472
#define HSM2_MOD_SCP22     243
#define HSM2_MOD_SCP21     244
#define HSM2_MOD_BS1       245
#define HSM2_MOD_BS2       246
#define HSM2_MOD_KAPPA     251
#define HSM2_MOD_VDIFFJ    254
#define HSM2_MOD_DLY1      255
#define HSM2_MOD_DLY2      256
#define HSM2_MOD_DLY3      257
#define HSM2_MOD_NFTRP     258
#define HSM2_MOD_NFALP     259
#define HSM2_MOD_FALPH     263
#define HSM2_MOD_CIT       260
#define HSM2_MOD_OVSLP     261
#define HSM2_MOD_OVMAG     262
#define HSM2_MOD_GIDL4     281
#define HSM2_MOD_GIDL5     282
#define HSM2_MOD_GIDL6        189
#define HSM2_MOD_GIDL7        194
#define HSM2_MOD_SVGS      283
#define HSM2_MOD_SVBS      284
#define HSM2_MOD_SVBSL     285
#define HSM2_MOD_SVDS      286
#define HSM2_MOD_SLG       287
#define HSM2_MOD_SUB1L     290
#define HSM2_MOD_SUB2L     292
#define HSM2_MOD_VOVERS    303
#define HSM2_MOD_VOVERSP   304
#define HSM2_MOD_SVGSL     305
#define HSM2_MOD_SVGSLP    306
#define HSM2_MOD_SVGSWP    307
#define HSM2_MOD_SVGSW     308
#define HSM2_MOD_SVBSLP    309
#define HSM2_MOD_SLGL      310
#define HSM2_MOD_SLGLP     311
#define HSM2_MOD_SUB1LP    312
#define HSM2_MOD_IBPC1     404
#define HSM2_MOD_IBPC2     405
#define HSM2_MOD_MPHDFM    409

#define HSM2_MOD_PTL       450
#define HSM2_MOD_PTP       451
#define HSM2_MOD_PT2       452
#define HSM2_MOD_PTLP      455
#define HSM2_MOD_GDL       453
#define HSM2_MOD_GDLP      454

#define HSM2_MOD_GDLD      456
#define HSM2_MOD_PT4       457
#define HSM2_MOD_PT4P      465
#define HSM2_MOD_MUEPHL2   458
#define HSM2_MOD_MUEPLP2   459
#define HSM2_MOD_NSUBCW2   461
#define HSM2_MOD_NSUBCWP2  462
#define HSM2_MOD_MUEPHW2   463
#define HSM2_MOD_MUEPWP2   464

/* val symbol for WPE */
#define HSM2_MOD_WEB	    88
#define HSM2_MOD_WEC	    89
#define HSM2_MOD_NSUBCWPE  91
#define HSM2_MOD_NPEXTWPE  41
#define HSM2_MOD_NSUBPWPE  43

#define HSM2_MOD_VGSMIN    466
#define HSM2_MOD_SC3VBS    467
#define HSM2_MOD_BYPTOL    468
#define HSM2_MOD_MUECB0LP   473
#define HSM2_MOD_MUECB1LP   474

/* binning parameters */
#define HSM2_MOD_LMIN       1000
#define HSM2_MOD_LMAX       1001
#define HSM2_MOD_WMIN       1002
#define HSM2_MOD_WMAX       1003
#define HSM2_MOD_LBINN      1004
#define HSM2_MOD_WBINN      1005

/* Length dependence */
#define HSM2_MOD_LVMAX      1100
#define HSM2_MOD_LBGTMP1    1101
#define HSM2_MOD_LBGTMP2    1102
#define HSM2_MOD_LEG0       1103
#define HSM2_MOD_LLOVER     1106
#define HSM2_MOD_LVFBOVER   1428
#define HSM2_MOD_LNOVER     1430
#define HSM2_MOD_LWL2       1407
#define HSM2_MOD_LVFBC      1121
#define HSM2_MOD_LNSUBC     1123
#define HSM2_MOD_LNSUBP     1181
#define HSM2_MOD_LSCP1      1184
#define HSM2_MOD_LSCP2      1185
#define HSM2_MOD_LSCP3      1186
#define HSM2_MOD_LSC1       1126
#define HSM2_MOD_LSC2       1127
#define HSM2_MOD_LSC3       1128
#define HSM2_MOD_LSC4       1270
#define HSM2_MOD_LPGD1      1187
//#define HSM2_MOD_LPGD3      1189
#define HSM2_MOD_LNDEP      1129
#define HSM2_MOD_LNINV      1130
#define HSM2_MOD_LMUECB0    1131
#define HSM2_MOD_LMUECB1    1132
#define HSM2_MOD_LMUEPH1    1133
#define HSM2_MOD_LVTMP      1141
#define HSM2_MOD_LWVTH0     1142
#define HSM2_MOD_LMUESR1    1143
#define HSM2_MOD_LMUETMP    1195
#define HSM2_MOD_LSUB1      1151
#define HSM2_MOD_LSUB2      1152
#define HSM2_MOD_LSVDS      1286
#define HSM2_MOD_LSVBS      1284
#define HSM2_MOD_LSVGS      1283
#define HSM2_MOD_LNSTI      1225
#define HSM2_MOD_LWSTI      1226
#define HSM2_MOD_LSCSTI1    1229
#define HSM2_MOD_LSCSTI2    1230
#define HSM2_MOD_LVTHSTI    1232
#define HSM2_MOD_LMUESTI1   1235
#define HSM2_MOD_LMUESTI2   1236
#define HSM2_MOD_LMUESTI3   1237
#define HSM2_MOD_LNSUBPSTI1 1238
#define HSM2_MOD_LNSUBPSTI2 1239
#define HSM2_MOD_LNSUBPSTI3 1240
#define HSM2_MOD_LNSUBCSTI1   253
#define HSM2_MOD_LNSUBCSTI2   264
#define HSM2_MOD_LNSUBCSTI3   265
#define HSM2_MOD_LCGSO      1154
#define HSM2_MOD_LCGDO      1155
#define HSM2_MOD_LJS0       1157
#define HSM2_MOD_LJS0SW     1158
#define HSM2_MOD_LNJ        1159
#define HSM2_MOD_LCISBK     1172
#define HSM2_MOD_LCLM1      1191
#define HSM2_MOD_LCLM2      1192
#define HSM2_MOD_LCLM3      1193
#define HSM2_MOD_LWFC       1201
#define HSM2_MOD_LGIDL1     1205
#define HSM2_MOD_LGIDL2     1206
#define HSM2_MOD_LGLEAK1    1208
#define HSM2_MOD_LGLEAK2    1209
#define HSM2_MOD_LGLEAK3    1210
#define HSM2_MOD_LGLEAK6    1213
#define HSM2_MOD_LGLKSD1    1215
#define HSM2_MOD_LGLKSD2    1216
#define HSM2_MOD_LGLKB1     1218
#define HSM2_MOD_LGLKB2     1219
#define HSM2_MOD_LNFTRP     1258
#define HSM2_MOD_LNFALP     1259
#define HSM2_MOD_LVDIFFJ    1254
#define HSM2_MOD_LIBPC1     1404
#define HSM2_MOD_LIBPC2     1405

/* Width dependence */
#define HSM2_MOD_WVMAX      2100
#define HSM2_MOD_WBGTMP1    2101
#define HSM2_MOD_WBGTMP2    2102
#define HSM2_MOD_WEG0       2103
#define HSM2_MOD_WLOVER     2106
#define HSM2_MOD_WVFBOVER   2428
#define HSM2_MOD_WNOVER     2430
#define HSM2_MOD_WWL2       2407
#define HSM2_MOD_WVFBC      2121
#define HSM2_MOD_WNSUBC     2123
#define HSM2_MOD_WNSUBP     2181
#define HSM2_MOD_WSCP1      2184
#define HSM2_MOD_WSCP2      2185
#define HSM2_MOD_WSCP3      2186
#define HSM2_MOD_WSC1       2126
#define HSM2_MOD_WSC2       2127
#define HSM2_MOD_WSC3       2128
#define HSM2_MOD_WSC4       2270
#define HSM2_MOD_WPGD1      2187
//#define HSM2_MOD_WPGD3      2189
#define HSM2_MOD_WNDEP      2129
#define HSM2_MOD_WNINV      2130
#define HSM2_MOD_WMUECB0    2131
#define HSM2_MOD_WMUECB1    2132
#define HSM2_MOD_WMUEPH1    2133
#define HSM2_MOD_WVTMP      2141
#define HSM2_MOD_WWVTH0     2142
#define HSM2_MOD_WMUESR1    2143
#define HSM2_MOD_WMUETMP    2195
#define HSM2_MOD_WSUB1      2151
#define HSM2_MOD_WSUB2      2152
#define HSM2_MOD_WSVDS      2286
#define HSM2_MOD_WSVBS      2284
#define HSM2_MOD_WSVGS      2283
#define HSM2_MOD_WNSTI      2225
#define HSM2_MOD_WWSTI      2226
#define HSM2_MOD_WSCSTI1    2229
#define HSM2_MOD_WSCSTI2    2230
#define HSM2_MOD_WVTHSTI    2232
#define HSM2_MOD_WMUESTI1   2235
#define HSM2_MOD_WMUESTI2   2236
#define HSM2_MOD_WMUESTI3   2237
#define HSM2_MOD_WNSUBPSTI1 2238
#define HSM2_MOD_WNSUBPSTI2 2239
#define HSM2_MOD_WNSUBPSTI3 2240
#define HSM2_MOD_WNSUBCSTI1   266
#define HSM2_MOD_WNSUBCSTI2   267
#define HSM2_MOD_WNSUBCSTI3   268
#define HSM2_MOD_WCGSO      2154
#define HSM2_MOD_WCGDO      2155
#define HSM2_MOD_WJS0       2157
#define HSM2_MOD_WJS0SW     2158
#define HSM2_MOD_WNJ        2159
#define HSM2_MOD_WCISBK     2172
#define HSM2_MOD_WCLM1      2191
#define HSM2_MOD_WCLM2      2192
#define HSM2_MOD_WCLM3      2193
#define HSM2_MOD_WWFC       2201
#define HSM2_MOD_WGIDL1     2205
#define HSM2_MOD_WGIDL2     2206
#define HSM2_MOD_WGLEAK1    2208
#define HSM2_MOD_WGLEAK2    2209
#define HSM2_MOD_WGLEAK3    2210
#define HSM2_MOD_WGLEAK6    2213
#define HSM2_MOD_WGLKSD1    2215
#define HSM2_MOD_WGLKSD2    2216
#define HSM2_MOD_WGLKB1     2218
#define HSM2_MOD_WGLKB2     2219
#define HSM2_MOD_WNFTRP     2258
#define HSM2_MOD_WNFALP     2259
#define HSM2_MOD_WVDIFFJ    2254
#define HSM2_MOD_WIBPC1     2404
#define HSM2_MOD_WIBPC2     2405

/* Cross-term dependence */
#define HSM2_MOD_PVMAX      3100
#define HSM2_MOD_PBGTMP1    3101
#define HSM2_MOD_PBGTMP2    3102
#define HSM2_MOD_PEG0       3103
#define HSM2_MOD_PLOVER     3106
#define HSM2_MOD_PVFBOVER   3428
#define HSM2_MOD_PNOVER     3430
#define HSM2_MOD_PWL2       3407
#define HSM2_MOD_PVFBC      3121
#define HSM2_MOD_PNSUBC     3123
#define HSM2_MOD_PNSUBP     3181
#define HSM2_MOD_PSCP1      3184
#define HSM2_MOD_PSCP2      3185
#define HSM2_MOD_PSCP3      3186
#define HSM2_MOD_PSC1       3126
#define HSM2_MOD_PSC2       3127
#define HSM2_MOD_PSC3       3128
#define HSM2_MOD_PSC4       3270
#define HSM2_MOD_PPGD1      3187
//#define HSM2_MOD_PPGD3      3189
#define HSM2_MOD_PNDEP      3129
#define HSM2_MOD_PNINV      3130
#define HSM2_MOD_PMUECB0    3131
#define HSM2_MOD_PMUECB1    3132
#define HSM2_MOD_PMUEPH1    3133
#define HSM2_MOD_PVTMP      3141
#define HSM2_MOD_PWVTH0     3142
#define HSM2_MOD_PMUESR1    3143
#define HSM2_MOD_PMUETMP    3195
#define HSM2_MOD_PSUB1      3151
#define HSM2_MOD_PSUB2      3152
#define HSM2_MOD_PSVDS      3286
#define HSM2_MOD_PSVBS      3284
#define HSM2_MOD_PSVGS      3283
#define HSM2_MOD_PNSTI      3225
#define HSM2_MOD_PWSTI      3226
#define HSM2_MOD_PSCSTI1    3229
#define HSM2_MOD_PSCSTI2    3230
#define HSM2_MOD_PVTHSTI    3232
#define HSM2_MOD_PMUESTI1   3235
#define HSM2_MOD_PMUESTI2   3236
#define HSM2_MOD_PMUESTI3   3237
#define HSM2_MOD_PNSUBPSTI1 3238
#define HSM2_MOD_PNSUBPSTI2 3239
#define HSM2_MOD_PNSUBPSTI3 3240
#define HSM2_MOD_PNSUBCSTI1   269
#define HSM2_MOD_PNSUBCSTI2   270
#define HSM2_MOD_PNSUBCSTI3   271
#define HSM2_MOD_PCGSO      3154
#define HSM2_MOD_PCGDO      3155
#define HSM2_MOD_PJS0       3157
#define HSM2_MOD_PJS0SW     3158
#define HSM2_MOD_PNJ        3159
#define HSM2_MOD_PCISBK     3172
#define HSM2_MOD_PCLM1      3191
#define HSM2_MOD_PCLM2      3192
#define HSM2_MOD_PCLM3      3193
#define HSM2_MOD_PWFC       3201
#define HSM2_MOD_PGIDL1     3205
#define HSM2_MOD_PGIDL2     3206
#define HSM2_MOD_PGLEAK1    3208
#define HSM2_MOD_PGLEAK2    3209
#define HSM2_MOD_PGLEAK3    3210
#define HSM2_MOD_PGLEAK6    3213
#define HSM2_MOD_PGLKSD1    3215
#define HSM2_MOD_PGLKSD2    3216
#define HSM2_MOD_PGLKB1     3218
#define HSM2_MOD_PGLKB2     3219
#define HSM2_MOD_PNFTRP     3258
#define HSM2_MOD_PNFALP     3259
#define HSM2_MOD_PVDIFFJ    3254
#define HSM2_MOD_PIBPC1     3404
#define HSM2_MOD_PIBPC2     3405

/* device questions */
#define HSM2_DNODE          341
#define HSM2_GNODE          342
#define HSM2_SNODE          343
#define HSM2_BNODE          344
#define HSM2_DNODEPRIME     345
#define HSM2_SNODEPRIME     346
#define HSM2_BNODEPRIME     395
#define HSM2_DBNODE         396
#define HSM2_SBNODE         397
#define HSM2_VBD            347
#define HSM2_VBS            348
#define HSM2_VGS            349
#define HSM2_VDS            350
#define HSM2_CD             351
#define HSM2_CBS            352
#define HSM2_CBD            353
#define HSM2_GM             354
#define HSM2_GDS            355
#define HSM2_GMBS           356
#define HSM2_GBD            357
#define HSM2_GBS            358
#define HSM2_QB             359
#define HSM2_CQB            360
#define HSM2_QG             361
#define HSM2_CQG            362
#define HSM2_QD             363
#define HSM2_CQD            364
#define HSM2_CGG            365
#define HSM2_CGD            366
#define HSM2_CGS            367
#define HSM2_CBG            368
#define HSM2_CAPBD          369
#define HSM2_CQBD           370
#define HSM2_CAPBS          371
#define HSM2_CQBS           372
#define HSM2_CDG            373
#define HSM2_CDD            374
#define HSM2_CDS            375
#define HSM2_VON            376
#define HSM2_VDSAT          377
#define HSM2_QBS            378
#define HSM2_QBD            379
#define HSM2_SOURCECONDUCT  380
#define HSM2_DRAINCONDUCT   381
#define HSM2_CBDB           382
#define HSM2_CBSB           383
#define HSM2_MOD_RBPB       389
#define HSM2_MOD_RBPD       390
#define HSM2_MOD_RBPS       391
#define HSM2_MOD_RBDB       392
#define HSM2_MOD_RBSB       393
#define HSM2_MOD_GBMIN      394

#define HSM2_ISUB           410
#define HSM2_IGIDL          411
#define HSM2_IGISL          412
#define HSM2_IGD            413
#define HSM2_IGS            414
#define HSM2_IGB            415
#define HSM2_CGSO           416
#define HSM2_CGBO           417
#define HSM2_CGDO           418

#define HSM2_MOD_TCJBD         92
#define HSM2_MOD_TCJBS         93
#define HSM2_MOD_TCJBDSW       94 
#define HSM2_MOD_TCJBSSW       95  
#define HSM2_MOD_TCJBDSWG      96   
#define HSM2_MOD_TCJBSSWG      97   

#define HSM2_MOD_VGS_MAX    4001
#define HSM2_MOD_VGD_MAX    4002
#define HSM2_MOD_VGB_MAX    4003
#define HSM2_MOD_VDS_MAX    4004
#define HSM2_MOD_VBS_MAX    4005
#define HSM2_MOD_VBD_MAX    4006

#include "hsm2ext.h"

/*
extern void HSM2evaluate(double,double,double,HSM2instance*,HSM2model*,
        double*,double*,double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, CKTcircuit*);
*/

#endif /*HSM2*/

