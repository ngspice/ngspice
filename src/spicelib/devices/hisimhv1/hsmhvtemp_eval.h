/***********************************************************************

 HiSIM (Hiroshima University STARC IGFET Model)
 Copyright (C) 2012 Hiroshima University & STARC

 MODEL NAME : HiSIM_HV 
 ( VERSION : 1  SUBVERSION : 2  REVISION : 4 )
 Model Parameter VERSION : 1.23
 FILE : hsmhvtemp_eval.h

 DATE : 2013.04.30

 released by
                Hiroshima University &
                Semiconductor Technology Academic Research Center (STARC)
***********************************************************************/
#define C_rdtemp_min 5.0e-3
#define C_rdtemp_dlt 1.0e-2

    TTEMP = ckt->CKTtemp;
    if ( here->HSMHV_dtemp_Given ) { TTEMP = TTEMP + here->HSMHV_dtemp ; }
    TTEMP0 = TTEMP ; 
#ifdef HSMHVEVAL
    /* Self heating */
    TTEMP = TTEMP + deltemp ; 
#endif
/*Tdiff0 =   TTEMP0 - model->HSMHV_ktnom ;*/
/*Tdiff0_2 = TTEMP0 * TTEMP0 - model->HSMHV_ktnom * model->HSMHV_ktnom ;*/
Tdiff =    TTEMP  - model->HSMHV_ktnom ;
Tdiff_2 =  TTEMP  * TTEMP  - model->HSMHV_ktnom * model->HSMHV_ktnom ;

        /* Band gap */
        T1 = TTEMP - model->HSMHV_ktnom ;
        T2 = TTEMP * TTEMP - model->HSMHV_ktnom * model->HSMHV_ktnom ;
        here->HSMHV_eg = Eg = here->HSMHV_egtnom - pParam->HSMHV_bgtmp1 * T1
          - pParam->HSMHV_bgtmp2 * T2 ;
        here->HSMHV_sqrt_eg = sqrt( Eg ) ;
#ifdef HSMHVEVAL
        Eg_dT = -pParam->HSMHV_bgtmp1 - 2.0e0 * TTEMP * pParam->HSMHV_bgtmp2 ;
#endif

	T1 = 1.0 / TTEMP ;
	T2 = 1.0 / model->HSMHV_ktnom ;
	T3 = here->HSMHV_egtnom + model->HSMHV_egig
	  + model->HSMHV_igtemp2 * ( T1 - T2 )
	  + model->HSMHV_igtemp3 * ( T1 * T1 - T2 * T2 ) ;
 	here->HSMHV_egp12 = sqrt ( T3 ) ;
 	here->HSMHV_egp32 = T3 * here->HSMHV_egp12 ;

        
        /* Inverse of the thermal voltage */
        here->HSMHV_beta = beta = C_QE / (C_KB * TTEMP) ;
        here->HSMHV_beta_inv = 1.0 / beta ;
        here->HSMHV_beta2 = beta * beta ;
        here->HSMHV_betatnom = C_QE / (C_KB * model->HSMHV_ktnom) ;
#ifdef HSMHVEVAL
        beta_dT=-C_QE/(C_KB*TTEMP*TTEMP);
        beta_inv_dT = C_KB / C_QE ;
#endif

        /* Intrinsic carrier concentration */
        here->HSMHV_nin = Nin = C_Nin0 * Fn_Pow (TTEMP / model->HSMHV_ktnom, 1.5e0) 
          * exp (- Eg / 2.0e0 * beta + here->HSMHV_egtnom / 2.0e0 * here->HSMHV_betatnom) ;
#ifdef HSMHVEVAL
        Nin_dT = C_Nin0 * exp (- Eg / 2.0e0 * beta + here->HSMHV_egtnom / 2.0e0 * here->HSMHV_betatnom)
          * 1.5e0 * Fn_Pow ( TTEMP / model->HSMHV_ktnom , 0.5e0 ) / model->HSMHV_ktnom 
          + C_Nin0 * Fn_Pow (TTEMP / model->HSMHV_ktnom, 1.5e0) 
          * exp (- Eg / 2.0e0 * beta + here->HSMHV_egtnom / 2.0e0 * here->HSMHV_betatnom)
          * ( - Eg / 2.0e0 * beta_dT - beta / 2.0e0 * Eg_dT );
#endif

        /* Phonon Scattering (temperature-dependent part) */
        T1 =  Fn_Pow (TTEMP / model->HSMHV_ktnom, pParam->HSMHV_muetmp) ;
        here->HSMHV_mphn0 = T1 / here->HSMHV_mueph ;
        here->HSMHV_mphn1 = here->HSMHV_mphn0 * model->HSMHV_mueph0 ;
#ifdef HSMHVEVAL
        T1_dT = pParam->HSMHV_muetmp * Fn_Pow(TTEMP / model->HSMHV_ktnom, pParam->HSMHV_muetmp - 1.0 )
          / model->HSMHV_ktnom ;
        mphn0_dT = T1_dT / here->HSMHV_mueph ;
#endif

        /* Pocket Overlap (temperature-dependent part) */
        here->HSMHV_ptovr = here->HSMHV_ptovr0 / beta ;
#ifdef HSMHVEVAL
        ptovr_dT =  here->HSMHV_ptovr0 * beta_inv_dT ; 
#endif

        /* Velocity Temperature Dependence */
        T1 = TTEMP  / model->HSMHV_ktnom ;
        T3 = TTEMP0 - model->HSMHV_ktnom ;
        T4 = TTEMP0 * TTEMP0 - model->HSMHV_ktnom * model->HSMHV_ktnom ;
        T0 = 1.8 + 0.4 * T1 + 0.1 * T1 * T1 - pParam->HSMHV_vtmp * (1.0 - T1) ;
	 if ( model->HSMHV_cotemp != 2 ) { /* without deltemp (COTEMP=0,1,3) */
        here->HSMHV_vmax = here->HSMHV_vmax0 * pParam->HSMHV_vmax
          / T0
          * ( 1.0 + model->HSMHV_vmaxt1 * T3 + model->HSMHV_vmaxt2 * T4 ) ;
#ifdef HSMHVEVAL
        Vmax_dT=-here->HSMHV_vmax0 * pParam->HSMHV_vmax 
          / ( T0 * T0 ) * ( 1.0 + model->HSMHV_vmaxt1 * T3 + model->HSMHV_vmaxt2 * T4 )
          * 1/model->HSMHV_ktnom * (0.4 + 0.2 * T1 + pParam->HSMHV_vtmp) ;
#endif
	 } else { /* with deltemp (COTEMP=2) */
	   here->HSMHV_vmax = here->HSMHV_vmax0 * pParam->HSMHV_vmax
	     / T0 
	     * ( 1.0 + model->HSMHV_vmaxt1 * Tdiff + model->HSMHV_vmaxt2 * Tdiff_2 ) ;
#ifdef HSMHVEVAL
	   /* under development */
	   Vmax_dT = here->HSMHV_vmax0 * pParam->HSMHV_vmax 
	     / ( T0 * T0 )
	     * ( ( model->HSMHV_vmaxt1 + 2.0 * TTEMP * model->HSMHV_vmaxt2 ) * T0 
		- ( 1.0 + model->HSMHV_vmaxt1 * Tdiff + model->HSMHV_vmaxt2 * Tdiff_2 )
		* 1/model->HSMHV_ktnom * (0.4 + 0.2 * T1 + pParam->HSMHV_vtmp) ) ;
#endif
	 }
	 if ( model->HSMHV_cotemp != 2 ) { /* without deltemp (COTEMP=0,1,3) */
        here->HSMHV_ninvd = here->HSMHV_ninvd0 * ( 1.0 + model->HSMHV_ninvdt1 * T3 + model->HSMHV_ninvdt2 * T4 ) ;
#ifdef HSMHVEVAL
	ninvd_dT = 0.0 ;
#endif
	 } else { /* with deltemp (COTEMP=2) */
	   /* under development */
	   here->HSMHV_ninvd = here->HSMHV_ninvd0 * ( 1.0 + model->HSMHV_ninvdt1 * Tdiff + model->HSMHV_ninvdt2 * Tdiff_2 ) ;
#ifdef HSMHVEVAL
	   ninvd_dT = here->HSMHV_ninvd0 * ( model->HSMHV_ninvdt1 + 2.0 * TTEMP * model->HSMHV_ninvdt2 ) ;
#endif
	 }
      
	/* Temperature Dependence of RTH0 */
	pParam->HSMHV_rth = ( pParam->HSMHV_rth0 + model->HSMHV_rthtemp1 * T3 + model->HSMHV_rthtemp2 * T4  ) * here->HSMHV_rthtemp0 ;


	/* Temperature Dependence of POWRAT */
        T2 = pParam->HSMHV_powrat + model->HSMHV_prattemp1 * T3 + model->HSMHV_prattemp2 * T4  ;
	Fn_SL( T2 , T2 , 0 , 0.05 , T0 );
	Fn_SU( here->HSMHV_powratio , T2 , 1 , 0.05 , T0 );


        /* 2 phi_B (temperature-dependent) */
        /* @temp, with pocket */
        here->HSMHV_pb2 =  2.0e0 / beta * log (here->HSMHV_nsub / Nin) ;
#ifdef HSMHVEVAL
        Pb2_dT = - (here->HSMHV_pb2 * beta_dT  + 2.0e0 / Nin * Nin_dT ) / beta ;
#endif

        /* Depletion Width */
        T1 = 2.0e0 * C_ESI / C_QE ;
        here->HSMHV_wdpl = sqrt ( T1 / here->HSMHV_nsub ) ;
        here->HSMHV_wdplp = sqrt( T1 / ( here->HSMHV_nsubp ) ) ; 

        /* Coefficient of the F function for bulk charge */
        here->HSMHV_cnst0 = sqrt ( 2.0 * C_ESI * C_QE * here->HSMHV_nsub / beta ) ;

        /* cnst1: n_{p0} / p_{p0} */
        T1 = Nin / here->HSMHV_nsub ;
        here->HSMHV_cnst1 = T1 * T1 ;
#ifdef HSMHVEVAL
        cnst0_dT = 0.5e0 / here->HSMHV_cnst0 * 2.0 * C_ESI * C_QE * here->HSMHV_nsub * beta_inv_dT ;
        cnst1_dT = 2.0e0 * Nin * Nin_dT / here->HSMHV_nsub / here->HSMHV_nsub ;
#endif

	if ( pParam->HSMHV_nover != 0.0 ) {
	  here->HSMHV_cnst0over = here->HSMHV_cnst0 * sqrt( pParam->HSMHV_nover / here->HSMHV_nsub ) ;     
#ifdef HSMHVEVAL
         cnst0over_dT = cnst0_dT * sqrt( pParam->HSMHV_nover / here->HSMHV_nsub ) ; 
#endif
        /* ps0ldinib : Ps0_iniB for Ps0LD */
        T1 = here->HSMHV_cnst0over * model->HSMHV_tox / here->HSMHV_cecox  ;
        T2 = pParam->HSMHV_nover / Nin ;
        here->HSMHV_ps0ldinib = T2 * T2 / ( T1 * T1 );
#ifdef HSMHVEVAL
        T1_dT = cnst0over_dT  * model->HSMHV_tox / here->HSMHV_cecox  ;  
        T2_dT = - Nin_dT * T2 / Nin; 
        ps0ldinib_dT = 2.0 * here->HSMHV_ps0ldinib * ( T2_dT * T1 - T2 * T1_dT ) / ( T1 * T2 );
#endif
	}
	if ( pParam->HSMHV_novers != 0.0 ) {
	  here->HSMHV_cnst0overs = here->HSMHV_cnst0 * sqrt( pParam->HSMHV_novers / here->HSMHV_nsub ) ;     
#ifdef HSMHVEVAL
	  cnst0overs_dT = cnst0_dT * sqrt( pParam->HSMHV_novers / here->HSMHV_nsub ) ;
#endif
          /* ps0ldinib : Ps0_iniB for Ps0LD */
          T1 = here->HSMHV_cnst0overs * model->HSMHV_tox / here->HSMHV_cecox  ;
          T2 = pParam->HSMHV_novers / Nin ;
          here->HSMHV_ps0ldinibs = T2 * T2 / ( T1 * T1 );
#ifdef HSMHVEVAL
         T1_dT = cnst0overs_dT  * model->HSMHV_tox / here->HSMHV_cecox  ;  
         T2_dT = - Nin_dT * T2 / Nin; 
         ps0ldinibs_dT = 2.0 * here->HSMHV_ps0ldinibs * ( T2_dT * T1 - T2 * T1_dT ) / ( T1 * T2 );
#endif
 	}

	/* temperature-dependent resistance model */
        T3 = model->HSMHV_ktnom ;
        T1 = TTEMP0 - T3 ;
        T4 = TTEMP0 * TTEMP0 - T3 * T3 ;
	/* drain side */
	if ( pParam->HSMHV_rd > 0.0 ) {
         T2 = here->HSMHV_rdtemp0
	   * ( here->HSMHV_ldrift1 * pParam->HSMHV_rdslp1 * C_m2um   + pParam->HSMHV_rdict1 )
	   * ( here->HSMHV_ldrift2 * model->HSMHV_rdslp2 * C_m2um   + model->HSMHV_rdict2 ) ;

	 if ( model->HSMHV_cotemp == 1 ) { /* without deltemp (COTEMP=1) */
	 here->HSMHV_rd = ( pParam->HSMHV_rd + T1 * modelMKS->HSMHV_rdtemp1 + modelMKS->HSMHV_rdtemp2 * T4 ) * T2 ;
	 Fn_SL( here->HSMHV_rd, here->HSMHV_rd, C_rdtemp_min * pParam->HSMHV_rd, C_rdtemp_dlt * pParam->HSMHV_rd, T0 );
#ifdef HSMHVEVAL
	 Rd0_dT = 0.0 ;
#endif
	 } else { /* with deltemp (COTEMP=0,2,3) */
	   here->HSMHV_rd = ( pParam->HSMHV_rd + modelMKS->HSMHV_rdtemp1 * Tdiff + modelMKS->HSMHV_rdtemp2 * Tdiff_2 ) * T2 ;
	   Fn_SL( here->HSMHV_rd, here->HSMHV_rd, C_rdtemp_min * pParam->HSMHV_rd, C_rdtemp_dlt * pParam->HSMHV_rd, T0 );
#ifdef HSMHVEVAL
	   Rd0_dT = ( modelMKS->HSMHV_rdtemp1 + 2.0 * TTEMP * modelMKS->HSMHV_rdtemp2 ) * T2 * T0 ;
#endif
	 }

	} else {
	  here->HSMHV_rd = 0.0 ;
	}
	/* source side (asymmetric case) */
	if ( pParam->HSMHV_rs > 0.0 ) {
	    T2 = here->HSMHV_rdtemp0
	      * ( here->HSMHV_ldrift1s * pParam->HSMHV_rdslp1 * C_m2um   + pParam->HSMHV_rdict1 ) 
	      * ( here->HSMHV_ldrift2s * model->HSMHV_rdslp2 * C_m2um   + model->HSMHV_rdict2 ) ;
	    
	 if ( model->HSMHV_cotemp == 1 ) { /* without deltemp (COTEMP=1) */
	    here->HSMHV_rs = ( pParam->HSMHV_rs + T1 * modelMKS->HSMHV_rdtemp1 + modelMKS->HSMHV_rdtemp2 * T4 ) * T2 ;
	    Fn_SL( here->HSMHV_rs, here->HSMHV_rs, C_rdtemp_min * pParam->HSMHV_rs, C_rdtemp_dlt * pParam->HSMHV_rs, T0 );
#ifdef HSMHVEVAL
	    Rs0_dT = 0.0 ;
#endif
	 } else { /* with deltemp (COTEMP=0,2,3) */
	    here->HSMHV_rs = ( pParam->HSMHV_rs + modelMKS->HSMHV_rdtemp1 * Tdiff + modelMKS->HSMHV_rdtemp2 * Tdiff_2 ) * T2 ;
	    Fn_SL( here->HSMHV_rs, here->HSMHV_rs, C_rdtemp_min * pParam->HSMHV_rs, C_rdtemp_dlt * pParam->HSMHV_rs, T0 );
#ifdef HSMHVEVAL
	    Rs0_dT = ( modelMKS->HSMHV_rdtemp1 + 2.0 * TTEMP * modelMKS->HSMHV_rdtemp2 ) * T2 * T0 ;
#endif
	 }

	} else {
	  here->HSMHV_rs = 0.0 ;
	}
        if ( pParam->HSMHV_rdvd > 0.0 ) {
          T4 = here->HSMHV_rdvdtemp0 * ( here->HSMHV_ldrift1 * pParam->HSMHV_rdslp1 * C_m2um  + pParam->HSMHV_rdict1 )
	                              * ( here->HSMHV_ldrift2 * model->HSMHV_rdslp2 * C_m2um  + model->HSMHV_rdict2 ) ;

	  T1 = ( 1 -  pParam->HSMHV_rdov13 ) * here->HSMHV_loverld * C_m2um ; 
	  T0 = - model->HSMHV_rdov11 / ( model->HSMHV_rdov12 + small ) ;
	  T3 = ( T0 * here->HSMHV_loverld * C_m2um  + 1.0 + model->HSMHV_rdov11 ) ;
	  Fn_SL( T5 , T3 * T4 , T4 , 10.0e-3 , T6 ) ;
	  Fn_SU( T7  , T5    , T4 * ( model->HSMHV_rdov11 + 1.0) , 50.0e-6 , T6 ) ;
	  Fn_SL( T2  , T7 + T1 * T4  , 0, 50.0e-6 , T6 ) ;

          T3 = model->HSMHV_ktnom ;
          T1 = TTEMP0 - T3 ;
	 if ( model->HSMHV_cotemp == 0 || model->HSMHV_cotemp == 1 ) { /* without deltemp (COTEMP=0,1) */
	  here->HSMHV_rdvd = ( pParam->HSMHV_rdvd + T1 * modelMKS->HSMHV_rdvdtemp1 + modelMKS->HSMHV_rdvdtemp2 * ( TTEMP0 * TTEMP0 - T3 * T3 ) ) * T2 ;
	  Fn_SL( here->HSMHV_rdvd, here->HSMHV_rdvd, C_rdtemp_min * pParam->HSMHV_rdvd, C_rdtemp_dlt * pParam->HSMHV_rdvd, T0 );
#ifdef HSMHVEVAL
          Rdvd_dT = 0.0 ;
#endif
	 } else { /* with deltemp (COTEMP=2,3) */
	   here->HSMHV_rdvd = ( pParam->HSMHV_rdvd + modelMKS->HSMHV_rdvdtemp1 * Tdiff + modelMKS->HSMHV_rdvdtemp2 * Tdiff_2 ) * T2 ;
	   Fn_SL( here->HSMHV_rdvd, here->HSMHV_rdvd, C_rdtemp_min * pParam->HSMHV_rdvd, C_rdtemp_dlt * pParam->HSMHV_rdvd, T0 );
#ifdef HSMHVEVAL
	   Rdvd_dT = ( modelMKS->HSMHV_rdvdtemp1 + 2.0 * TTEMP * modelMKS->HSMHV_rdvdtemp2 ) * T2 * T0 ;
#endif
	 }

          T4 = here->HSMHV_rdvdtemp0 * ( here->HSMHV_ldrift1s * pParam->HSMHV_rdslp1 * C_m2um  + pParam->HSMHV_rdict1 )
	                              * ( here->HSMHV_ldrift2s * model->HSMHV_rdslp2 * C_m2um  + model->HSMHV_rdict2 ) ;

	  T1 = ( 1 -  pParam->HSMHV_rdov13 ) * here->HSMHV_lovers * C_m2um ; 
	  T0 = - model->HSMHV_rdov11 / ( model->HSMHV_rdov12 + small ) ;
	  T3 = ( T0 * here->HSMHV_lovers * C_m2um + 1.0 + model->HSMHV_rdov11 ) ;
	  Fn_SL( T5 , T3 * T4 , T4 , 10.0e-3 , T6 ) ;
	  Fn_SU( T7  , T5    , T4 * ( model->HSMHV_rdov11 + 1.0) , 50.0e-6 , T6 ) ;
	  Fn_SL( T2  , T7 + T1 * T4  , 0, 50.0e-6 , T6 ) ;

          T3 = model->HSMHV_ktnom ;
          T1 = TTEMP0 - T3 ;
	 if ( model->HSMHV_cotemp == 0 || model->HSMHV_cotemp == 1 ) { /* without deltemp (COTEMP=0,1) */
	  here->HSMHV_rsvd = ( pParam->HSMHV_rdvd + T1 * modelMKS->HSMHV_rdvdtemp1 + modelMKS->HSMHV_rdvdtemp2 * ( TTEMP0 * TTEMP0 - T3 * T3 ) ) * T2 ;
	  Fn_SL( here->HSMHV_rsvd, here->HSMHV_rsvd, C_rdtemp_min * pParam->HSMHV_rdvd, C_rdtemp_dlt * pParam->HSMHV_rdvd, T0 );
#ifdef HSMHVEVAL
          Rsvd_dT = 0.0 ;
#endif
	 } else { /* with deltemp (COTEMP=2,3) */
	   here->HSMHV_rsvd = ( pParam->HSMHV_rdvd + modelMKS->HSMHV_rdvdtemp1 * Tdiff + modelMKS->HSMHV_rdvdtemp2 * Tdiff_2 ) * T2 ;
	   Fn_SL( here->HSMHV_rsvd, here->HSMHV_rsvd, C_rdtemp_min * pParam->HSMHV_rdvd, C_rdtemp_dlt * pParam->HSMHV_rdvd, T0 );
#ifdef HSMHVEVAL
	   Rsvd_dT = ( modelMKS->HSMHV_rdvdtemp1 + 2.0 * TTEMP * modelMKS->HSMHV_rdvdtemp2 ) * T2 * T0 ;
#endif
	 }
	} else {
	  here->HSMHV_rdvd = 0.0 ;
	  here->HSMHV_rsvd = 0.0 ;
	}

        /* for substrate-source/drain junction diode. */
        js   = pParam->HSMHV_js0
          * exp ((here->HSMHV_egtnom * here->HSMHV_betatnom - Eg * beta
                  + model->HSMHV_xti * log (TTEMP / model->HSMHV_ktnom)) / pParam->HSMHV_nj) ;
        jssw = pParam->HSMHV_js0sw
          * exp ((here->HSMHV_egtnom * here->HSMHV_betatnom - Eg * beta 
                  + model->HSMHV_xti * log (TTEMP / model->HSMHV_ktnom)) / model->HSMHV_njsw) ;

        js2  = pParam->HSMHV_js0
          * exp ((here->HSMHV_egtnom * here->HSMHV_betatnom - Eg * beta
                  + model->HSMHV_xti2 * log (TTEMP / model->HSMHV_ktnom)) / pParam->HSMHV_nj) ;  
        jssw2 = pParam->HSMHV_js0sw
          * exp ((here->HSMHV_egtnom * here->HSMHV_betatnom - Eg * beta
                  + model->HSMHV_xti2 * log (TTEMP / model->HSMHV_ktnom)) / model->HSMHV_njsw) ; 
      
#ifdef HSMHVEVAL
	T0 = - Eg * beta_dT - Eg_dT * beta ; /* Self heating */
	T1 = T0 + model->HSMHV_xti  / TTEMP ; /* Self heating */
	T2 = T0 + model->HSMHV_xti2 / TTEMP ; /* Self heating */

	js_dT =    js    * T1  / pParam->HSMHV_nj; /* Self heating */
	jssw_dT =  jssw  * T1/ model->HSMHV_njsw ; /* Self heating */
	js2_dT =   js2   * T2  / pParam->HSMHV_nj; /* Self heating */
	jssw2_dT = jssw2 * T2 / model->HSMHV_njsw; /* Self heating */
#endif
      
        here->HSMHV_isbd = here->HSMHV_ad * js + here->HSMHV_pd * jssw ;
        here->HSMHV_isbd2 = here->HSMHV_ad * js2 + here->HSMHV_pd * jssw2 ;
        here->HSMHV_isbs = here->HSMHV_as * js + here->HSMHV_ps * jssw ;
        here->HSMHV_isbs2 = here->HSMHV_as * js2 + here->HSMHV_ps * jssw2 ;

#ifdef HSMHVEVAL
	isbd_dT =  here->HSMHV_ad * js_dT  + here->HSMHV_pd * jssw_dT  ; /* Self heating */
	isbd2_dT = here->HSMHV_ad * js2_dT + here->HSMHV_pd * jssw2_dT ; /* Self heating */
	isbs_dT =  here->HSMHV_as * js_dT  + here->HSMHV_ps * jssw_dT  ; /* Self heating */
	isbs2_dT = here->HSMHV_as * js2_dT + here->HSMHV_ps * jssw2_dT ; /* Self heating */
#endif

        T1 = TTEMP / model->HSMHV_ktnom ;
        T0 = T1 * T1 ;
        T2 = here->HSMHV_isbd + small ;
        T3 = here->HSMHV_isbs + small ;
#ifdef HSMHVEVAL
        T1_dT = 1.0 / model->HSMHV_ktnom ; /* Self heating */
        T0_dT = 2.0 * T1 * T1_dT ;       /* Self heating */
        T2_dT = isbd_dT ;                /* Self heating */
	T3_dT = isbs_dT ;                /* Self heating */
#endif

        here->HSMHV_vbdt = pParam->HSMHV_nj / beta 
          * log ( pParam->HSMHV_vdiffj * T0 / T2 + 1.0 ) ;
        here->HSMHV_vbst = pParam->HSMHV_nj / beta 
          * log ( pParam->HSMHV_vdiffj * T0 / T3 + 1.0 ) ;

        here->HSMHV_exptemp = exp (( T1 - 1.0 ) * model->HSMHV_ctemp ) ;

#ifdef HSMHVEVAL
	vbdt_dT = - beta_dT / beta * here->HSMHV_vbdt
	  + pParam->HSMHV_nj / beta * pParam->HSMHV_vdiffj / ( pParam->HSMHV_vdiffj * T0 / T2 + 1.0 ) 
	  * ( T0_dT / T2 - T0 / T2 / T2 * T2_dT ) ; /* Self heating */
	vbst_dT = - beta_dT / beta * here->HSMHV_vbst
	  + pParam->HSMHV_nj / beta * pParam->HSMHV_vdiffj / ( pParam->HSMHV_vdiffj * T0 / T3 + 1.0 ) 
	  * ( T0_dT / T3 - T0 / T3 / T3 * T3_dT ) ; /* Self heating */
#endif

        here->HSMHV_jd_nvtm_inv = 1.0 / ( pParam->HSMHV_nj / beta ) ;
        here->HSMHV_jd_expcd = exp (here->HSMHV_vbdt * here->HSMHV_jd_nvtm_inv ) ;
        here->HSMHV_jd_expcs = exp (here->HSMHV_vbst * here->HSMHV_jd_nvtm_inv ) ;

#ifdef HSMHVEVAL
       exptemp_dT = model->HSMHV_ctemp / model->HSMHV_ktnom * here->HSMHV_exptemp ;       /* Self heating */
	jd_nvtm_inv_dT = beta_dT / pParam->HSMHV_nj ;                                  /* Self heating */
	jd_expcd_dT = here->HSMHV_jd_expcd
	  * ( vbdt_dT * here->HSMHV_jd_nvtm_inv + here->HSMHV_vbdt * jd_nvtm_inv_dT ) ; /* Self heating */
	jd_expcs_dT = here->HSMHV_jd_expcs
	  * ( vbst_dT * here->HSMHV_jd_nvtm_inv + here->HSMHV_vbst * jd_nvtm_inv_dT ) ; /* Self heating */
#endif

       /* costi0 and costi1 for STI transistor model (temperature-dependent part) */
       here->HSMHV_costi0 = here->HSMHV_costi00 * sqrt(here->HSMHV_beta_inv) ;
       here->HSMHV_costi0_p2 = here->HSMHV_costi0 * here->HSMHV_costi0 ;
       here->HSMHV_costi1 = here->HSMHV_nin * here->HSMHV_nin * here->HSMHV_nsti_p2 ;

/* end of HSMHVtemp_eval.h */
