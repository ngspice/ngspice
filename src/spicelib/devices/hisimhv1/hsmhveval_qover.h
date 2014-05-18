/***********************************************************************

 HiSIM (Hiroshima University STARC IGFET Model)
 Copyright (C) 2012 Hiroshima University & STARC

 MODEL NAME : HiSIM_HV 
 ( VERSION : 1  SUBVERSION : 2  REVISION : 4 )
 Model Parameter VERSION : 1.23
 FILE : hsmhveval_qover.h

 DATE : 2013.04.30

 released by
                Hiroshima University &
                Semiconductor Technology Academic Research Center (STARC)
***********************************************************************/

/*  Begin HSMHVevalQover */

    /*---------------------------------------------------*
     * Clamp -Vxbgmt.
     *-----------------*/
    T0 = - Vxbgmt;
    if ( T0 > Vbs_bnd ) {
      T1 =    T0   - Vbs_bnd;
      T1_dT =      - Vbs_bnd_dT;
      T2 =    Vbs_max    - Vbs_bnd;
      T2_dT = Vbs_max_dT - Vbs_bnd_dT;

      Fn_SUPoly4m( TY, T1, T2, T11, T0 );
      TY_dT = T1_dT * T11 + T2_dT * T0;

      T10 = Vbs_bnd + TY ;
      T10_dT = Vbs_bnd_dT + TY_dT ;
    }  else {
      T10 = T0 ;
      T11 = 1.0 ;
      T10_dT = 0.0;
    }
    Vxbgmtcl = - T10 - small2 ;
    Vxbgmtcl_dVxbgmt = T11;
    Vxbgmtcl_dT = - T10_dT;

    fac1 = cnst0over_func * Cox0_inv ;
    fac1_dVbs = 0.0; fac1_dVds = 0.0; fac1_dVgs = 0.0;
    fac1_dT = cnst0over_func_dT * Cox0_inv ;

    fac1p2 = fac1 * fac1 ;
    fac1p2_dT = 2.0 * fac1 * fac1_dT ;
  
    VgpLD = - Vgbgmt + pParam->HSMHV_vfbover;  
    VgpLD_dVgb = - 1.0e0 ;

    T0 = Nover_func / here->HSMHV_nin ;
    Pb2over = 2.0 / beta * log( T0 ) ;
    T0_dT = - T0 / here->HSMHV_nin * Nin_dT ;
    Pb2over_dT = - Pb2over / beta * beta_dT + 2.0 / beta / T0 * T0_dT ;

    Vgb_fb_LD =  - Vxbgmtcl ;

    /*-----------------------------------*
     * QsuLD: total charge = Accumulation | Depletion+inversion
     *-----------------*/
    if (   VgpLD  < Vgb_fb_LD ){   
      /*---------------------------*
       * Accumulation
       *-----------------*/
      flg_ovzone = -1 ; 
      T1 = 1.0 / ( beta * cnst0over_func ) ;
      T1_dT = - T1 * T1 * ( beta_dT * cnst0over_func + beta * cnst0over_func_dT ) ;
      TY = T1 * Cox0 ;
      Ac41 = 2.0 + 3.0 * C_SQRT_2 * TY ;
      Ac4 = 8.0 * Ac41 * Ac41 * Ac41 ;
      TY_dT = T1_dT  * Cox0 ;
      Ac41_dT = 3.0 * C_SQRT_2 * TY_dT ;
      Ac4_dT = 8.0 * 3.0 * Ac41 * Ac41 * Ac41_dT ;
  
      Ps0_min = here->HSMHV_eg - Pb2over ;
      Ps0_min_dT = Eg_dT - Pb2over_dT ; 
  
      TX = beta * ( VgpLD + Vxbgmtcl ) ;
      TX_dVxb = beta * Vxbgmtcl_dVxbgmt ;
      TX_dVgb = beta * VgpLD_dVgb ;
      TX_dT = beta_dT * ( VgpLD + Vxbgmtcl ) + beta * Vxbgmtcl_dT;
  
      Ac31 = 7.0 * C_SQRT_2 - 9.0 * TY * ( TX - 2.0 ) ;
      Ac31_dVxb = - 9.0 * TY * TX_dVxb ;
      Ac31_dVgb = - 9.0 * TY * TX_dVgb ;
      Ac31_dT = - 9.0 * ( TY_dT * ( TX - 2.0 ) + TY * TX_dT ); 
  
      Ac3 = Ac31 * Ac31 ;
      T1 = 2.0 * Ac31 ;
      Ac3_dVxb = T1 * Ac31_dVxb ;
      Ac3_dVgb = T1 * Ac31_dVgb ;
      Ac3_dT = T1 * Ac31_dT ; 
  
      Ac2 = sqrt( Ac4 + Ac3 ) ;
      T1 = 0.5 / Ac2 ;
      Ac2_dVxb = T1 *  Ac3_dVxb ;
      Ac2_dVgb = T1 *  Ac3_dVgb ;
      Ac2_dT = T1 *  ( Ac4_dT + Ac3_dT ); 
    
      Ac1 = -7.0 * C_SQRT_2 + Ac2 + 9.0 * TY * ( TX - 2.0 ) ;
      Ac1_dVxb = Ac2_dVxb + 9.0 * TY * TX_dVxb ;
      Ac1_dVgb = Ac2_dVgb + 9.0 * TY * TX_dVgb ;
      Ac1_dT = Ac2_dT + 9.0 * ( TY_dT * ( TX - 2.0 ) + TY * TX_dT ) ; 
  
      Acd = pow( Ac1 , C_1o3 ) ;
      T1 = C_1o3 / ( Acd * Acd ) ;
      Acd_dVxb = Ac1_dVxb * T1 ;
      Acd_dVgb = Ac1_dVgb * T1 ;
      Acd_dT = Ac1_dT * T1 ; 
  
      Acn = -4.0 * C_SQRT_2 - 12.0 * TY + 2.0 * Acd + C_SQRT_2 * Acd * Acd ;
      T1 = 2.0 + 2.0 * C_SQRT_2 * Acd ;
      Acn_dVxb = T1 * Acd_dVxb ;
      Acn_dVgb = T1 * Acd_dVgb ;
      Acn_dT = - 12.0 * TY_dT + T1 * Acd_dT ; 
   
      Chi = Acn / Acd ;
      T1 = 1.0 / ( Acd * Acd ) ;
      Chi_dVxb = ( Acn_dVxb * Acd - Acn * Acd_dVxb ) * T1 ;
      Chi_dVgb = ( Acn_dVgb * Acd - Acn * Acd_dVgb ) * T1 ;
      Chi_dT = ( Acn_dT * Acd - Acn * Acd_dT ) * T1 ; 
  
      Psa = Chi * beta_inv - Vxbgmtcl ;
      Psa_dVxb = Chi_dVxb * beta_inv - Vxbgmtcl_dVxbgmt ;
      Psa_dVgb = Chi_dVgb * beta_inv ;
      Psa_dT = Chi_dT * beta_inv + Chi * beta_inv_dT - Vxbgmtcl_dT ;
  
      T1 = Psa + Vxbgmtcl ;
      T1_dT = Psa_dT + Vxbgmtcl_dT ;
      T2 = T1 / Ps0_min ;
      T2_dT = ( T1_dT * Ps0_min - T1 * Ps0_min_dT ) / ( Ps0_min * Ps0_min ) ;
      T3 = sqrt( 1.0 + ( T2 * T2 ) ) ;
  
      T9 = T2 / T3 / Ps0_min ;
      T3_dVd = T9 * ( Psa_dVxb + Vxbgmtcl_dVxbgmt ) ;
      T3_dVg = T9 * Psa_dVgb ;
      T3_dT =  T2_dT * T2 / T3;

      Ps0LD = T1 / T3 - Vxbgmtcl ;
      T9 = 1.0 / ( T3 * T3 ) ;
      Ps0LD_dVxb = T9 * ( ( Psa_dVxb + Vxbgmtcl_dVxbgmt ) * T3 - T1 * T3_dVd ) - Vxbgmtcl_dVxbgmt ;
      Ps0LD_dVgb = T9 * ( Psa_dVgb * T3 - T1 * T3_dVg );
      Ps0LD_dT = T9 * ( T1_dT * T3 - T1 * T3_dT ) - Vxbgmtcl_dT;
     
      T2 = ( VgpLD - Ps0LD ) ;
      QsuLD = Cox0 * T2 ;
      QsuLD_dVxb = - Cox0 * Ps0LD_dVxb ;
      QsuLD_dVgb = Cox0 * ( VgpLD_dVgb - Ps0LD_dVgb ) ;
      QsuLD_dT = Cox0 * ( - Ps0LD_dT ) ; 
  
      QbuLD = QsuLD ;
      QbuLD_dVxb = QsuLD_dVxb ;
      QbuLD_dVgb = QsuLD_dVgb ;
      QbuLD_dT = QsuLD_dT ; 
  
    } else {
  
      /*---------------------------*
       * Depletion and inversion
       *-----------------*/

      /* initial value for a few fixpoint iterations
         to get Ps0_iniA from simplified Poisson equation: */
       flg_ovzone = 2 ;
       Chi = znbd3 ;
       Chi_dVxb = 0.0 ; Chi_dVgb = 0.0 ; Chi_dT = 0.0 ;

       Ps0_iniA= Chi/beta - Vxbgmtcl ;
       Ps0_iniA_dVxb = Chi_dVxb/beta - Vxbgmtcl_dVxbgmt ;
       Ps0_iniA_dVgb = Chi_dVgb/beta ;
       Ps0_iniA_dT   = Chi_dT/beta - Chi*beta_dT/(beta*beta) - Vxbgmtcl_dT;
      
      /* 1 .. 2 relaxation steps should be sufficient */
      for ( lp_ld = 1; lp_ld <= 2; lp_ld ++ ) {
        TY = exp(-Chi);
        TY_dVxb = -Chi_dVxb * TY;
        TY_dVgb = -Chi_dVgb * TY;
        TY_dT   = - Chi_dT  * TY;
        TX = 1.0e0 + 4.0e0 
           * ( beta * ( VgpLD + Vxbgmtcl ) - 1.0e0 + TY ) / ( fac1p2 * beta2 ) ;
        TX_dVxb = 4.0e0 * ( beta * ( Vxbgmtcl_dVxbgmt ) + TY_dVxb ) / ( fac1p2 * beta2 );
        TX_dVgb = 4.0e0 * ( beta * ( VgpLD_dVgb       ) + TY_dVgb ) / ( fac1p2 * beta2 );
        T1 = ( beta * ( VgpLD + Vxbgmtcl ) - 1.0e0 + TY );
        T1_dT = beta_dT * ( VgpLD + Vxbgmtcl ) + beta * Vxbgmtcl_dT + TY_dT;
        T3 = fac1p2 * beta2 ;
        T3_dT = fac1p2_dT * beta2 + fac1p2 * ( 2 * beta * beta_dT ) ;
        TX_dT = 4 * ( T1_dT * T3 - T1 * T3_dT ) / ( T3 * T3 );
        if ( TX < epsm10) {
          TX = epsm10;
          TX_dVxb = TX_dVgb = TX_dT = 0.0;
        }

        Ps0_iniA = VgpLD + fac1p2 * beta / 2.0e0 * ( 1.0e0 - sqrt( TX ) ) ;
        Ps0_iniA_dVxb =            - fac1p2 * beta / 2.0e0 * TX_dVxb * 0.5 / sqrt( TX );
        Ps0_iniA_dVgb = VgpLD_dVgb - fac1p2 * beta / 2.0e0 * TX_dVgb * 0.5 / sqrt( TX );
        T1 = fac1p2 * beta ;
        T1_dT = fac1p2_dT * beta + fac1p2 * beta_dT ;
        T2 = 1.0 - sqrt( TX );
        T2_dT = - 1.0e0 / ( 2.0e0 * sqrt( TX ) ) * TX_dT ;
        Ps0_iniA_dT = ( T1_dT * T2 + T1 * T2_dT ) / 2.0e0 ;
  
        Chi = beta * ( Ps0_iniA + Vxbgmtcl ) ;
        Chi_dVxb = beta * ( Ps0_iniA_dVxb + Vxbgmtcl_dVxbgmt ) ;
        Chi_dVgb = beta * ( Ps0_iniA_dVgb ) ;
        Chi_dT = beta_dT * ( Ps0_iniA + Vxbgmtcl ) + beta * ( Ps0_iniA_dT + Vxbgmtcl_dT );
      } /* End of iteration */

      if ( Chi < znbd3 ) { 

        flg_ovzone = 1 ; 

        /*-----------------------------------*
         * zone-D1
         * - Ps0_iniA is the analytical solution of QovLD=Qb0 with
         *   Qb0 being approximated by 3-degree polynomial.
         *
         *   new: Inclusion of exp(-Chi) term at right border
         *-----------------*/
        Ta =  1.0/(9.0*sqrt(2.0)) - (5.0+7.0*exp(-3.0)) / (54.0*sqrt(2.0+exp(-3.0)));
        Tb = (1.0+exp(-3.0)) / (2.0*sqrt(2.0+exp(-3.0))) - sqrt(2.0) / 3.0;
        Tc =  1.0/sqrt(2.0) + 1.0/(beta*fac1);
        Tc_dT = - (beta_dT*fac1 + beta*fac1_dT)/(beta2*fac1p2);
        Td = - (VgpLD + Vxbgmtcl) / fac1;
        Td_dVxb = - Vxbgmtcl_dVxbgmt / fac1;
        Td_dVgb = - VgpLD_dVgb / fac1;
        Td_dT   = - (Vxbgmtcl_dT*fac1 - (VgpLD+Vxbgmtcl)*fac1_dT)/fac1p2;
        Tq = Tb*Tb*Tb / (27.0*Ta*Ta*Ta) - Tb*Tc/(6.0*Ta*Ta) + Td/(2.0*Ta);
        Tq_dVxb = Td_dVxb/(2.0*Ta);
        Tq_dVgb = Td_dVgb / (2.0*Ta);
        Tq_dT   = - Tb/(6.0*Ta*Ta)*Tc_dT + Td_dT/(2.0*Ta);
        Tp = (3.0*Ta*Tc-Tb*Tb)/(9.0*Ta*Ta);
        Tp_dT = Tc_dT/(3.0*Ta);
        T5      = sqrt(Tq*Tq + Tp*Tp*Tp);
        T5_dVxb = 2.0*Tq*Tq_dVxb / (2.0*T5);
        T5_dVgb = 2.0*Tq*Tq_dVgb / (2.0*T5);
        T5_dT   = (2.0*Tq*Tq_dT + 3.0*Tp*Tp*Tp_dT) / (2.0*T5);
        Tu = pow(-Tq + T5,C_1o3);
        Tu_dVxb = Tu / (3.0 * (-Tq + T5)) * (-Tq_dVxb + T5_dVxb);
        Tu_dVgb = Tu / (3.0 * (-Tq + T5)) * (-Tq_dVgb + T5_dVgb);
        Tu_dT   = Tu / (3.0 * (-Tq + T5)) * (-Tq_dT   + T5_dT);
        Tv = -pow(Tq + T5,C_1o3);
        Tv_dVxb = Tv / (3.0 * (-Tq - T5)) * (-Tq_dVxb - T5_dVxb);
        Tv_dVgb = Tv / (3.0 * (-Tq - T5)) * (-Tq_dVgb - T5_dVgb);
        Tv_dT   = Tv / (3.0 * (-Tq - T5)) * (-Tq_dT   - T5_dT );
        TX      = Tu + Tv - Tb/(3.0*Ta);
        TX_dVxb = Tu_dVxb + Tv_dVxb;
        TX_dVgb = Tu_dVgb + Tv_dVgb;
        TX_dT   = Tu_dT   + Tv_dT  ;
        
        Ps0_iniA = TX * beta_inv - Vxbgmtcl ;
        Ps0_iniA_dVxb = TX_dVxb * beta_inv - Vxbgmtcl_dVxbgmt;
        Ps0_iniA_dVgb = TX_dVgb * beta_inv;
	Ps0_iniA_dT = TX_dT * beta_inv + TX * beta_inv_dT - Vxbgmtcl_dT;

        Chi = beta * ( Ps0_iniA + Vxbgmtcl ) ;
        Chi_dVxb = beta * ( Ps0_iniA_dVxb + Vxbgmtcl_dVxbgmt ) ;
        Chi_dVgb = beta * ( Ps0_iniA_dVgb ) ;
        Chi_dT = beta_dT * ( Ps0_iniA + Vxbgmtcl ) + beta * ( Ps0_iniA_dT + Vxbgmtcl_dT );
      }

      if ( model->HSMHV_coqovsm > 0 ) {
	  /*-----------------------------------*
	   * - Ps0_iniB : upper bound.
	   *-----------------*/
        flg_ovzone += 2;

        VgpLD_shift = VgpLD + Vxbgmtcl + 0.1;
        VgpLD_shift_dVgb = VgpLD_dVgb;
        VgpLD_shift_dVxb = Vxbgmtcl_dVxbgmt;
        VgpLD_shift_dT   = Vxbgmtcl_dT;
        exp_bVbs = exp( beta * - Vxbgmtcl ) + small ;
        exp_bVbs_dVxb = - exp_bVbs * beta * Vxbgmtcl_dVxbgmt;
        exp_bVbs_dT   = - exp_bVbs * (beta_dT*Vxbgmtcl + beta*Vxbgmtcl_dT);
        T0 = here->HSMHV_nin / Nover_func;
        T0_dT = Nin_dT / Nover_func;
        cnst1over = T0 * T0;
        cnst1over_dT = 2.0 * T0 * T0_dT;
        gamma = cnst1over * exp_bVbs ;
        gamma_dVxb = cnst1over * exp_bVbs_dVxb;
        gamma_dT   = cnst1over_dT * exp_bVbs + cnst1over * exp_bVbs_dT;
     
        T0    = beta2 * fac1p2;
        T0_dT = 2.0 * beta * fac1 * (beta_dT*fac1+beta*fac1_dT);

        psi = beta*VgpLD_shift;
        psi_dVgb = beta*VgpLD_shift_dVgb;
        psi_dVxb = beta*VgpLD_shift_dVxb;
        psi_dT   = beta_dT*VgpLD_shift + beta*VgpLD_shift_dT;
        Chi_1      = log(gamma*T0 + psi*psi) - log(cnst1over*T0) + beta*Vxbgmtcl;
        Chi_1_dVgb = 2.0*psi*psi_dVgb/ (gamma*T0 + psi*psi);
        Chi_1_dVxb = (gamma_dVxb*T0+2.0*psi*psi_dVxb)/(gamma*T0+psi*psi)
                            + beta*Vxbgmtcl_dVxbgmt;    
        Chi_1_dT   = (gamma_dT*T0+gamma*T0_dT+2.0*psi*psi_dT)/(gamma*T0+psi*psi)
                            - (cnst1over_dT*T0 + cnst1over*T0_dT)/(cnst1over*T0)
                            + beta_dT*Vxbgmtcl + beta*Vxbgmtcl_dT;

        Fn_SU2( Chi_1, Chi_1, psi, 1.0, T1, T2 );
        Chi_1_dVgb = Chi_1_dVgb*T1 + psi_dVgb*T2;
        Chi_1_dVxb = Chi_1_dVxb*T1 + psi_dVxb*T2;
        Chi_1_dT   = Chi_1_dT  *T1 + psi_dT  *T2;

     /* 1 fixpoint step for getting more accurate Chi_B */
        psi      -= Chi_1 ;
        psi_dVgb -= Chi_1_dVgb ;
        psi_dVxb -= Chi_1_dVxb ;
        psi_dT   -= Chi_1_dT ;
     
        psi      += beta*0.1 ;
        psi_dT   += beta_dT*0.1 ;

/*        psi_B = psi;*/
/*        arg_B = psi*psi/(gamma*T0);*/
        Chi_B = log(gamma*T0 + psi*psi) - log(cnst1over*T0) + beta*Vxbgmtcl;
        Chi_B_dVgb = 2.0*psi*psi_dVgb/ (gamma*T0 + psi*psi);
        Chi_B_dVxb = (gamma_dVxb*T0+2.0*psi*psi_dVxb)/(gamma*T0+psi*psi)
                            + beta*Vxbgmtcl_dVxbgmt;    
        Chi_B_dT   = (gamma_dT*T0+gamma*T0_dT+2.0*psi*psi_dT)/(gamma*T0+psi*psi)
                            - (cnst1over_dT*T0 + cnst1over*T0_dT)/(cnst1over*T0)
                            + beta_dT*Vxbgmtcl + beta*Vxbgmtcl_dT;
        Ps0_iniB      = Chi_B/beta - Vxbgmtcl ;
/*        Ps0_iniB_dVgb = Chi_B_dVgb/beta;
        Ps0_iniB_dVxb = Chi_B_dVxb/beta- Vxbgmtcl_dVxbgmt;
        Ps0_iniB_dT   = Chi_B_dT/beta - Chi_B/(beta*beta)*beta_dT - Vxbgmtcl_dT;
*/
        
        /* construction of Ps0LD by taking Ps0_iniB as an upper limit of Ps0_iniA
         *
         * Limiting is done for Chi rather than for Ps0LD, to avoid shifting
         * for Fn_SU2 */

        Chi_A = Chi;
        Chi_A_dVxb = Chi_dVxb;
        Chi_A_dVgb = Chi_dVgb;
        Chi_A_dT   = Chi_dT;

        Fn_SU2( Chi, Chi_A, Chi_B, c_ps0ini_2*75.00, T1, T2 ); /* org: 50 */
        Chi_dVgb = Chi_A_dVgb * T1 + Chi_B_dVgb * T2;
        Chi_dVxb = Chi_A_dVxb * T1 + Chi_B_dVxb * T2;
        Chi_dT   = Chi_A_dT   * T1 + Chi_B_dT   * T2;

      }

        /* updating Ps0LD */
        Ps0LD= Chi/beta - Vxbgmtcl ;
        Ps0LD_dVgb = Chi_dVgb/beta;
        Ps0LD_dVxb = Chi_dVxb/beta- Vxbgmtcl_dVxbgmt;
        Ps0LD_dT   = Chi_dT/beta - Chi/(beta*beta)*beta_dT - Vxbgmtcl_dT;

      T1      = Chi - 1.0 + exp(-Chi);
      T1_dVxb = (1.0 - exp(-Chi)) * Chi_dVxb ;
      T1_dVgb = (1.0 - exp(-Chi)) * Chi_dVgb ;
      T1_dT   = (1.0 - exp(-Chi)) * Chi_dT   ;
      if (T1 < epsm10) {
         T1 = epsm10 ;
         T1_dVxb = 0.0 ;
         T1_dVgb = 0.0 ;
         T1_dT   = 0.0 ;
      }
      T2 = sqrt(T1);
      QbuLD = cnst0over_func * T2 ;
      T3 = cnst0over_func * 0.5 / T2 ;
      QbuLD_dVxb = T3 * T1_dVxb ;
      QbuLD_dVgb = T3 * T1_dVgb ;
      QbuLD_dT = T3 * T1_dT + cnst0over_func_dT * T2 ; 
     
      /*-----------------------------------------------------------*
       * QsuLD : Qovs or Qovd in unit area.
       * note: QsuLD = Qdep+Qinv. 
       *-----------------*/
      QsuLD = Cox0 * ( VgpLD - Ps0LD ) ;
      QsuLD_dVxb = Cox0 * ( - Ps0LD_dVxb ) ;
      QsuLD_dVgb = Cox0 * ( VgpLD_dVgb - Ps0LD_dVgb ) ;
      QsuLD_dT = Cox0 * ( - Ps0LD_dT ) ;

      if ( model->HSMHV_coqovsm == 1 ) { /* take initial values from analytical model */ 

   
        /*---------------------------------------------------*
         * Calculation of Ps0LD. (beginning of Newton loop) 
         * - Fs0 : Fs0 = 0 is the equation to be solved. 
         * - dPs0 : correction value. 
         *-----------------*/

        /* initial value too close to flat band should not be used */
//      if (Ps0LD+Vxbgmtcl < 1.0e-2) Ps0LD = 1.0e-2 - Vxbgmtcl;
        exp_bVbs = exp( beta * - Vxbgmtcl ) ;
        T0 = here->HSMHV_nin / Nover_func;
        cnst1over = T0 * T0;
        cnst1over_dT = 2.0 * T0 * ( Nin_dT / Nover_func );
        cfs1 = cnst1over * exp_bVbs ;
    
        flg_conv = 0 ;
        for ( lp_s0 = 1 ; lp_s0 <= 2*lp_s0_max + 1 ; lp_s0 ++ ) { 

            Chi = beta * ( Ps0LD + Vxbgmtcl ) ;
   
            if ( Chi < znbd5 ) { 
              /*-------------------------------------------*
               * zone-D1/D2. (Ps0LD)
               *-----------------*/
              fi = Chi * Chi * Chi 
                * ( cn_im53 + Chi * ( cn_im54 + Chi * cn_im55 ) ) ;
              fi_dChi = Chi * Chi 
                * ( 3 * cn_im53 + Chi * ( 4 * cn_im54 + Chi * 5 * cn_im55 ) ) ;
      
              fs01 = cfs1 * fi * fi ;
              fs01_dPs0 = cfs1 * beta * 2 * fi * fi_dChi ;

              fb = Chi * ( cn_nc51 
                 + Chi * ( cn_nc52 
                 + Chi * ( cn_nc53 
                 + Chi * ( cn_nc54 + Chi * cn_nc55 ) ) ) ) ;
              fb_dChi = cn_nc51 
                 + Chi * ( 2 * cn_nc52 
                 + Chi * ( 3 * cn_nc53
                 + Chi * ( 4 * cn_nc54 + Chi * 5 * cn_nc55 ) ) ) ;

              fs02 = sqrt( fb * fb + fs01 + small ) ;
              fs02_dPs0 = ( beta * fb_dChi * 2 * fb + fs01_dPs0 ) / ( fs02 + fs02 ) ;

            } else {
             /*-------------------------------------------*
              * zone-D3. (Ps0LD)
              *-----------------*/
             if ( Chi < large_arg ) { /* avoid exp_Chi to become extremely large */
	        exp_Chi = exp( Chi ) ;
	        fs01 = cfs1 * ( exp_Chi - 1.0e0 ) ;
	        fs01_dPs0 = cfs1 * beta * ( exp_Chi ) ;
             } else {
                exp_bPs0 = exp( beta*Ps0LD ) ;
                fs01     = cnst1over * ( exp_bPs0 - exp_bVbs ) ;
                fs01_dPs0 = cnst1over * beta * exp_bPs0 ;
             }
             fs02 = sqrt( Chi - 1.0 + fs01 ) ;
             fs02_dPs0 = ( beta + fs01_dPs0 ) / fs02 * 0.5 ;
   
            } /* end of if ( Chi  ... ) block */
            /*-----------------------------------------------------------*
             * Fs0
             *-----------------*/
            Fs0 = VgpLD - Ps0LD - fac1 * fs02 ;
            Fs0_dPs0 = - 1.0e0 - fac1 * fs02_dPs0 ;

            if ( flg_conv == 1 ) break ;

            dPs0 = - Fs0 / Fs0_dPs0 ;

            /*-------------------------------------------*
             * Update Ps0LD .
             *-----------------*/
            dPlim = 0.5*dP_max*(1.0 + Fn_Max(1.e0,fabs(Ps0LD))) ;
            if ( fabs( dPs0 ) > dPlim ) dPs0 = dPlim * Fn_Sgn( dPs0 ) ;

            Ps0LD = Ps0LD + dPs0 ;

            TX = -Vxbgmtcl + ps_conv / 2 ;
            if ( Ps0LD < TX ) Ps0LD = TX ;
      
            /*-------------------------------------------*
             * Check convergence. 
             *-----------------*/
            if ( fabs( dPs0 ) <= ps_conv && fabs( Fs0 ) <= gs_conv ) {
              flg_conv = 1 ;
            }
      
        } /* end of Ps0LD Newton loop */

        /*-------------------------------------------*
         * Procedure for diverged case.
         *-----------------*/
        if ( flg_conv == 0 ) { 
          fprintf( stderr , 
                   "*** warning(HiSIM_HV): Went Over Iteration Maximum (Ps0LD)\n" ) ;
          fprintf( stderr , " -Vxbgmtcl = %e   Vgbgmt = %e\n" , -Vxbgmtcl , Vgbgmt ) ;
        } 

        /*---------------------------------------------------*
         * Evaluate derivatives of Ps0LD. 
         *-----------------*/
        Chi_dT = beta_dT * ( Ps0LD + Vxbgmtcl ) + beta * Vxbgmtcl_dT;
        exp_bVbs_dT = - ( beta_dT * Vxbgmtcl + beta * Vxbgmtcl_dT ) * exp_bVbs;
        cfs1_dT = exp_bVbs * cnst1over_dT + exp_bVbs_dT * cnst1over;

        if ( Chi < znbd5 ) { 
          fs01_dVbs = cfs1 * beta * fi * ( - fi + 2 * fi_dChi ) ; /* fs01_dVxbgmtcl */
          fs01_dT = cfs1 * 2 * fi * fi_dChi * Chi_dT + fi * fi * cfs1_dT ;
          T2 = 1.0e0 / ( fs02 + fs02 ) ;
          fs02_dVbs = ( + beta * fb_dChi * 2 * fb + fs01_dVbs ) * T2 ; /* fs02_dVxbgmtcl */
          fs02_dT = ( 2 * fb * fb_dChi * Chi_dT + fs01_dT ) * T2 ;
        } else {
          if ( Chi < large_arg ) {
            fs01_dVbs = + cfs1 * beta ; /* fs01_dVxbgmtcl */
            exp_Chi_dT  = exp_Chi * Chi_dT ;
            fs01_dT     = ( exp_Chi - 1.0e0 ) * cfs1_dT + cfs1 * exp_Chi_dT ;
          } else {
            fs01_dVbs   = + cfs1 * beta ;
            exp_bPs0_dT = exp_bPs0 * Ps0LD * beta_dT ;
            fs01_dT     = cnst1over_dT*(exp_bPs0-exp_bVbs) + cnst1over*(exp_bPs0_dT-exp_bVbs_dT) ;
          }
          T2 = 0.5e0 / fs02 ;
          fs02_dVbs = ( + beta + fs01_dVbs ) * T2 ; /* fs02_dVxbgmtcl */
          fs02_dT = T2 * ( Chi_dT + fs01_dT ) ;
        }

        T1 = 1.0 / Fs0_dPs0 ;
        Ps0LD_dVxb = - ( - fac1 * fs02_dVbs ) * T1 ;
        Ps0LD_dVds = 0.0 ;
        Ps0LD_dVgb = - ( VgpLD_dVgb - fac1_dVgs * fs02 ) * T1 ;
        Ps0LD_dT = - ( - ( fac1 * fs02_dT + fac1_dT * fs02 ) ) * T1;

        Chi_dT = beta_dT * ( Ps0LD + Vxbgmtcl ) + beta * ( Ps0LD_dT + Vxbgmtcl_dT );

        if ( Chi < znbd5 ) { 
          /*-------------------------------------------*
           * zone-D1/D2. (Ps0LD)
           *-----------------*/
          if ( Chi < znbd3 ) { flg_ovzone = 1; }
                        else { flg_ovzone = 2; }

          Xi0 = fb * fb + epsm10 ;
          T1 = 2 * fb * fb_dChi * beta ;
          Xi0_dVbs = T1 * ( Ps0LD_dVxb + 1.0 ) ; /* Xi0_dVxbgmtcl */
          Xi0_dVds = T1 * Ps0LD_dVds ;
          Xi0_dVgs = T1 * Ps0LD_dVgb ;
          Xi0_dT = 2 * fb * fb_dChi * Chi_dT ;

          Xi0p12 = fb + epsm10 ;
          T1 = fb_dChi * beta ;
          Xi0p12_dVbs = T1 * ( Ps0LD_dVxb + 1.0 ) ; /* Xi0p12_dVxbgmtcl */
          Xi0p12_dVds = T1 * Ps0LD_dVds ;
          Xi0p12_dVgs = T1 * Ps0LD_dVgb ;
          Xi0p12_dT = fb_dChi * Chi_dT ;

          Xi0p32 = fb * fb * fb + epsm10 ;
          T1 = 3 * fb * fb * fb_dChi * beta ;
          Xi0p32_dVbs = T1 * ( Ps0LD_dVxb + 1.0 ) ; /* Xi0p32_dVxbgmtcl */
          Xi0p32_dVds = T1 * Ps0LD_dVds ;
          Xi0p32_dVgs = T1 * Ps0LD_dVgb ;
          Xi0p32_dT = 3 * fb * fb * fb_dChi * Chi_dT ;
 
        } else { 
          /*-------------------------------------------*
           * zone-D3. (Ps0LD)
           *-----------------*/
          flg_ovzone = 3 ;

          Xi0 = Chi - 1.0e0 ;
          Xi0_dVbs = beta * ( Ps0LD_dVxb + 1.0e0 ) ; /* Xi0_dVxbgmtcl */
          Xi0_dVds = beta * Ps0LD_dVds ;
          Xi0_dVgs = beta * Ps0LD_dVgb ;
          Xi0_dT = Chi_dT ;
 
          Xi0p12 = sqrt( Xi0 ) ;
          T1 = 0.5e0 / Xi0p12 ;
          Xi0p12_dVbs = T1 * Xi0_dVbs ;
          Xi0p12_dVds = T1 * Xi0_dVds ;
          Xi0p12_dVgs = T1 * Xi0_dVgs ;
          Xi0p12_dT = T1 * Xi0_dT ;

          Xi0p32 = Xi0 * Xi0p12 ;
          T1 = 1.5e0 * Xi0p12 ;
          Xi0p32_dVbs = T1 * Xi0_dVbs ;
          Xi0p32_dVds = T1 * Xi0_dVds ;
          Xi0p32_dVgs = T1 * Xi0_dVgs ;
          Xi0p32_dT = T1 * Xi0_dT ;

        } /* end of if ( Chi  ... ) block */
    
        /*-----------------------------------------------------------*
         * - Recalculate the derivatives of fs01 and fs02.
         *-----------------*/
        fs01_dVbs = Ps0LD_dVxb * fs01_dPs0 + fs01_dVbs ;
        fs01_dVds = Ps0LD_dVds * fs01_dPs0 ;
        fs01_dVgs = Ps0LD_dVgb * fs01_dPs0 ;
        fs01_dT   = Ps0LD_dT * fs01_dPs0 + fs01_dT;
        fs02_dVbs = Ps0LD_dVxb * fs02_dPs0 + fs02_dVbs ;
/*        fs02_dVxb = Ps0LD_dVds * fs02_dPs0 ;*/
        fs02_dVgb = Ps0LD_dVgb * fs02_dPs0 ;
        fs02_dT   = Ps0LD_dT * fs02_dPs0 + fs02_dT;

        /*-----------------------------------------------------------*
         * QbuLD and QiuLD
         *-----------------*/
        QbuLD = cnst0over_func * Xi0p12 ;
        QbuLD_dVxb = cnst0over_func * Xi0p12_dVbs ;
        QbuLD_dVgb = cnst0over_func * Xi0p12_dVgs ;
        QbuLD_dT =   cnst0over_func * Xi0p12_dT + cnst0over_func_dT * Xi0p12;

        T1 = 1.0 / ( fs02 + Xi0p12 ) ;
        QiuLD = cnst0over_func * fs01 * T1 ;
        T2 = 1.0 / ( fs01 + epsm10 ) ;
        QiuLD_dVbs = QiuLD * ( fs01_dVbs * T2 - ( fs02_dVbs + Xi0p12_dVbs ) * T1 ) ;
        QiuLD_dVgs = QiuLD * ( fs01_dVgs * T2 - ( fs02_dVgb + Xi0p12_dVgs ) * T1 ) ;
        T1_dT = - T1 * T1 * ( fs02_dT + Xi0p12_dT );
        QiuLD_dT = cnst0over_func * ( fs01 * T1_dT + T1 * fs01_dT ) + fs01 * T1 * cnst0over_func_dT;

        /*-----------------------------------------------------------*
         * Extrapolation: X_dVxbgmt = X_dVxbgmtcl * Vxbgmtcl_dVxbgmt
         *-----------------*/
        QbuLD_dVxb *= Vxbgmtcl_dVxbgmt ;
        QiuLD_dVbs *= Vxbgmtcl_dVxbgmt ;

        /*-----------------------------------------------------------*
         * Total overlap charge
         *-----------------*/
        QsuLD = QbuLD + QiuLD;
        QsuLD_dVxb = QbuLD_dVxb + QiuLD_dVbs;
        QsuLD_dVgb = QbuLD_dVgb + QiuLD_dVgs;
        QsuLD_dT =   QbuLD_dT   + QiuLD_dT;

      } /* end of COQOVSM branches */

    } /* end of Vgbgmt region blocks */
  
    /* convert to source ref. */
/*    Ps0LD_dVbs = Ps0LD_dVxb * Vxbgmt_dVbs + Ps0LD_dVgb * Vgbgmt_dVbs ;*/
    Ps0LD_dVds = Ps0LD_dVxb * Vxbgmt_dVds + Ps0LD_dVgb * Vgbgmt_dVds ;
/*    Ps0LD_dVgs = Ps0LD_dVxb * Vxbgmt_dVgs + Ps0LD_dVgb * Vgbgmt_dVgs ;*/

    QsuLD_dVbs = QsuLD_dVxb * Vxbgmt_dVbs + QsuLD_dVgb * Vgbgmt_dVbs ;
    QsuLD_dVds = QsuLD_dVxb * Vxbgmt_dVds + QsuLD_dVgb * Vgbgmt_dVds ;
    QsuLD_dVgs = QsuLD_dVxb * Vxbgmt_dVgs + QsuLD_dVgb * Vgbgmt_dVgs ;

    QbuLD_dVbs = QbuLD_dVxb * Vxbgmt_dVbs + QbuLD_dVgb * Vgbgmt_dVbs ;
    QbuLD_dVds = QbuLD_dVxb * Vxbgmt_dVds + QbuLD_dVgb * Vgbgmt_dVds ;
    QbuLD_dVgs = QbuLD_dVxb * Vxbgmt_dVgs + QbuLD_dVgb * Vgbgmt_dVgs ;
  
    /* inversion charge = total - depletion */
    QiuLD = QsuLD - QbuLD  ;
    QiuLD_dVbs = QsuLD_dVbs - QbuLD_dVbs ;
    QiuLD_dVds = QsuLD_dVds - QbuLD_dVds ;
    QiuLD_dVgs = QsuLD_dVgs - QbuLD_dVgs ;
    QiuLD_dT = QsuLD_dT - QbuLD_dT ; 

/*  End HSMHVevalQover */
