/**** BSIM4v7.0 Released by Darsen Lu 04/08/2011 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4check.c of BSIM4v7.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
 * Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
 * Project Director: Prof. Chenming Hu.
 * Modified by Xuemei Xi, 04/06/2001.
 * Modified by Xuemei Xi, 10/05/2001.
 * Modified by Xuemei Xi, 11/15/2002.
 * Modified by Xuemei Xi, 05/09/2003.
 * Modified by Xuemei Xi, 03/04/2004.
 * Modified by Xuemei Xi, 07/29/2005.
 * Modified by Mohan Dunga, 12/13/2006
 * Modified by Mohan Dunga, Wenwei Yang, 05/18/2007.
 * Modified by  Wenwei Yang, 07/31/2008 .
 * Modified by Tanvir Morshed, Darsen Lu 03/27/2011
 **********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "bsim4v7def.h."
#include "ngspice/trandefs.h"
#include "ngspice/const.h"
#include "ngspice/sperror.h"
#include "ngspice/devdefs.h"
#include "ngspice/suffix.h"

int
BSIM4v7checkModel(
BSIM4v7model *model,
BSIM4v7instance *here,
CKTcircuit *ckt)
{
struct bsim4SizeDependParam *pParam;
int Fatal_Flag = 0;
FILE *fplog;
    
    if ((fplog = fopen("bsim4.out", "w")) != NULL)
    {   pParam = here->pParam;
        fprintf(fplog, "BSIM4v7: Berkeley Short Channel IGFET Model-4\n");
        fprintf(fplog, "Developed by Xuemei (Jane) Xi, Mohan Dunga, Prof. Ali Niknejad and Prof. Chenming Hu in 2003.\n");
        fprintf(fplog, "\n");
        fprintf(fplog, "++++++++++ BSIM4v7 PARAMETER CHECKING BELOW ++++++++++\n");

        if ((strcmp(model->BSIM4v7version, "4.7.0")) && (strcmp(model->BSIM4v7version, "4.70")) && (strcmp(model->BSIM4v7version, "4.7")))
        {  fprintf(fplog, "Warning: This model is BSIM4v7.0; you specified a wrong version number.\n");
           printf("Warning: This model is BSIM4v7.0; you specified a wrong version number.\n");
        }
        fprintf(fplog, "Model = %s\n", model->BSIM4v7modName);


        if ((here->BSIM4v7rgateMod == 2) || (here->BSIM4v7rgateMod == 3))
        {   if ((here->BSIM4v7trnqsMod == 1) || (here->BSIM4v7acnqsMod == 1))
            {   fprintf(fplog, "Warning: You've selected both Rg and charge deficit NQS; select one only.\n");
                printf("Warning: You've selected both Rg and charge deficit NQS; select one only.\n");
            }
        }

        if (model->BSIM4v7toxe <= 0.0)
        {   fprintf(fplog, "Fatal: Toxe = %g is not positive.\n",
                    model->BSIM4v7toxe);
            printf("Fatal: Toxe = %g is not positive.\n", model->BSIM4v7toxe);
            Fatal_Flag = 1;
        }
        if (model->BSIM4v7toxp <= 0.0)
        {   fprintf(fplog, "Fatal: Toxp = %g is not positive.\n",
                    model->BSIM4v7toxp);
            printf("Fatal: Toxp = %g is not positive.\n", model->BSIM4v7toxp);
            Fatal_Flag = 1;
        }
        if (model->BSIM4v7eot <= 0.0)
        {   fprintf(fplog, "Fatal: EOT = %g is not positive.\n",
                    model->BSIM4v7eot);
            printf("Fatal: EOT = %g is not positive.\n", model->BSIM4v7eot);
            Fatal_Flag = 1;
        }
        if (model->BSIM4v7epsrgate < 0.0)
        {   fprintf(fplog, "Fatal: Epsrgate = %g is not positive.\n",
                    model->BSIM4v7epsrgate);
            printf("Fatal: Epsrgate = %g is not positive.\n", model->BSIM4v7epsrgate);
            Fatal_Flag = 1;
        }
        if (model->BSIM4v7epsrsub < 0.0)
        {   fprintf(fplog, "Fatal: Epsrsub = %g is not positive.\n",
                    model->BSIM4v7epsrsub);
            printf("Fatal: Epsrsub = %g is not positive.\n", model->BSIM4v7epsrsub);
            Fatal_Flag = 1;
        }
        if (model->BSIM4v7easub < 0.0)
        {   fprintf(fplog, "Fatal: Easub = %g is not positive.\n",
                    model->BSIM4v7easub);
            printf("Fatal: Easub = %g is not positive.\n", model->BSIM4v7easub);
            Fatal_Flag = 1;
        }
        if (model->BSIM4v7ni0sub <= 0.0)
        {   fprintf(fplog, "Fatal: Ni0sub = %g is not positive.\n",
                    model->BSIM4v7ni0sub);
            printf("Fatal: Easub = %g is not positive.\n", model->BSIM4v7ni0sub);
            Fatal_Flag = 1;
        }

        if (model->BSIM4v7toxm <= 0.0)
        {   fprintf(fplog, "Fatal: Toxm = %g is not positive.\n",
                    model->BSIM4v7toxm);
            printf("Fatal: Toxm = %g is not positive.\n", model->BSIM4v7toxm);
            Fatal_Flag = 1;
        }

        if (model->BSIM4v7toxref <= 0.0)
        {   fprintf(fplog, "Fatal: Toxref = %g is not positive.\n",
                    model->BSIM4v7toxref);
            printf("Fatal: Toxref = %g is not positive.\n", model->BSIM4v7toxref);
            Fatal_Flag = 1;
        }

        if (pParam->BSIM4v7lpe0 < -pParam->BSIM4v7leff)
        {   fprintf(fplog, "Fatal: Lpe0 = %g is less than -Leff.\n",
                    pParam->BSIM4v7lpe0);
            printf("Fatal: Lpe0 = %g is less than -Leff.\n",
                        pParam->BSIM4v7lpe0);
            Fatal_Flag = 1;
        }
        if (model->BSIM4v7lintnoi > pParam->BSIM4v7leff/2)
        {   fprintf(fplog, "Fatal: Lintnoi = %g is too large - Leff for noise is negative.\n",
                    model->BSIM4v7lintnoi);
            printf("Fatal: Lintnoi = %g is too large - Leff for noise is negative.\n",
                    model->BSIM4v7lintnoi);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4v7lpeb < -pParam->BSIM4v7leff)
        {   fprintf(fplog, "Fatal: Lpeb = %g is less than -Leff.\n",
                    pParam->BSIM4v7lpeb);
            printf("Fatal: Lpeb = %g is less than -Leff.\n",
                        pParam->BSIM4v7lpeb);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4v7ndep <= 0.0)
        {   fprintf(fplog, "Fatal: Ndep = %g is not positive.\n",
                    pParam->BSIM4v7ndep);
            printf("Fatal: Ndep = %g is not positive.\n",
                   pParam->BSIM4v7ndep);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4v7phi <= 0.0)
        {   fprintf(fplog, "Fatal: Phi = %g is not positive. Please check Phin and Ndep\n",
                    pParam->BSIM4v7phi);
            fprintf(fplog, "           Phin = %g  Ndep = %g \n", 
                        pParam->BSIM4v7phin, pParam->BSIM4v7ndep);
            printf("Fatal: Phi = %g is not positive. Please check Phin and Ndep\n",
                    pParam->BSIM4v7phi);
            printf("           Phin = %g  Ndep = %g \n", 
                        pParam->BSIM4v7phin, pParam->BSIM4v7ndep);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4v7nsub <= 0.0)
        {   fprintf(fplog, "Fatal: Nsub = %g is not positive.\n",
                    pParam->BSIM4v7nsub);
            printf("Fatal: Nsub = %g is not positive.\n",
                   pParam->BSIM4v7nsub);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4v7ngate < 0.0)
        {   fprintf(fplog, "Fatal: Ngate = %g is not positive.\n",
                    pParam->BSIM4v7ngate);
            printf("Fatal: Ngate = %g Ngate is not positive.\n",
                   pParam->BSIM4v7ngate);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4v7ngate > 1.e25)
        {   fprintf(fplog, "Fatal: Ngate = %g is too high.\n",
                    pParam->BSIM4v7ngate);
            printf("Fatal: Ngate = %g Ngate is too high\n",
                   pParam->BSIM4v7ngate);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4v7xj <= 0.0)
        {   fprintf(fplog, "Fatal: Xj = %g is not positive.\n",
                    pParam->BSIM4v7xj);
            printf("Fatal: Xj = %g is not positive.\n", pParam->BSIM4v7xj);
            Fatal_Flag = 1;
        }

        if (pParam->BSIM4v7dvt1 < 0.0)
        {   fprintf(fplog, "Fatal: Dvt1 = %g is negative.\n",
                    pParam->BSIM4v7dvt1);   
            printf("Fatal: Dvt1 = %g is negative.\n", pParam->BSIM4v7dvt1);   
            Fatal_Flag = 1;
        }
            
        if (pParam->BSIM4v7dvt1w < 0.0)
        {   fprintf(fplog, "Fatal: Dvt1w = %g is negative.\n",
                    pParam->BSIM4v7dvt1w);
            printf("Fatal: Dvt1w = %g is negative.\n", pParam->BSIM4v7dvt1w);
            Fatal_Flag = 1;
        }
            
        if (pParam->BSIM4v7w0 == -pParam->BSIM4v7weff)
        {   fprintf(fplog, "Fatal: (W0 + Weff) = 0 causing divided-by-zero.\n");
            printf("Fatal: (W0 + Weff) = 0 causing divided-by-zero.\n");
            Fatal_Flag = 1;
        }   

        if (pParam->BSIM4v7dsub < 0.0)
        {   fprintf(fplog, "Fatal: Dsub = %g is negative.\n", pParam->BSIM4v7dsub);
            printf("Fatal: Dsub = %g is negative.\n", pParam->BSIM4v7dsub);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4v7b1 == -pParam->BSIM4v7weff)
        {   fprintf(fplog, "Fatal: (B1 + Weff) = 0 causing divided-by-zero.\n");
            printf("Fatal: (B1 + Weff) = 0 causing divided-by-zero.\n");
            Fatal_Flag = 1;
        }  
        if (here->BSIM4v7u0temp <= 0.0)
        {   fprintf(fplog, "Fatal: u0 at current temperature = %g is not positive.\n", here->BSIM4v7u0temp);
            printf("Fatal: u0 at current temperature = %g is not positive.\n",
                   here->BSIM4v7u0temp);
            Fatal_Flag = 1;
        }
    
        if (pParam->BSIM4v7delta < 0.0)
        {   fprintf(fplog, "Fatal: Delta = %g is less than zero.\n",
                    pParam->BSIM4v7delta);
            printf("Fatal: Delta = %g is less than zero.\n", pParam->BSIM4v7delta);
            Fatal_Flag = 1;
        }      

        if (here->BSIM4v7vsattemp <= 0.0)
        {   fprintf(fplog, "Fatal: Vsat at current temperature = %g is not positive.\n", here->BSIM4v7vsattemp);
            printf("Fatal: Vsat at current temperature = %g is not positive.\n",
                   here->BSIM4v7vsattemp);
            Fatal_Flag = 1;
        }

        if (pParam->BSIM4v7pclm <= 0.0)
        {   fprintf(fplog, "Fatal: Pclm = %g is not positive.\n", pParam->BSIM4v7pclm);
            printf("Fatal: Pclm = %g is not positive.\n", pParam->BSIM4v7pclm);
            Fatal_Flag = 1;
        }

        if (pParam->BSIM4v7drout < 0.0)
        {   fprintf(fplog, "Fatal: Drout = %g is negative.\n", pParam->BSIM4v7drout);
            printf("Fatal: Drout = %g is negative.\n", pParam->BSIM4v7drout);
            Fatal_Flag = 1;
        }

        if (here->BSIM4v7m < 1.0)
        {   fprintf(fplog, "Fatal: Number of multiplier = %g is smaller than one.\n", here->BSIM4v7m);
            printf("Fatal: Number of multiplier = %g is smaller than one.\n", here->BSIM4v7m);
            Fatal_Flag = 1;
        }

        if (here->BSIM4v7nf < 1.0)
        {   fprintf(fplog, "Fatal: Number of finger = %g is smaller than one.\n", here->BSIM4v7nf);
            printf("Fatal: Number of finger = %g is smaller than one.\n", here->BSIM4v7nf);
            Fatal_Flag = 1;
        }

        if((here->BSIM4v7sa > 0.0) && (here->BSIM4v7sb > 0.0) && 
               ((here->BSIM4v7nf == 1.0) || ((here->BSIM4v7nf > 1.0) && (here->BSIM4v7sd > 0.0))) )
        {   if (model->BSIM4v7saref <= 0.0)
            {   fprintf(fplog, "Fatal: SAref = %g is not positive.\n",model->BSIM4v7saref);
                     printf("Fatal: SAref = %g is not positive.\n",model->BSIM4v7saref);
                     Fatal_Flag = 1;
            }
            if (model->BSIM4v7sbref <= 0.0)
            {   fprintf(fplog, "Fatal: SBref = %g is not positive.\n",model->BSIM4v7sbref);
                    printf("Fatal: SBref = %g is not positive.\n",model->BSIM4v7sbref);
                    Fatal_Flag = 1;
            }
         }

        if ((here->BSIM4v7l + model->BSIM4v7xl) <= model->BSIM4v7xgl)
        {   fprintf(fplog, "Fatal: The parameter xgl must be smaller than Ldrawn+XL.\n");
            printf("Fatal: The parameter xgl must be smaller than Ldrawn+XL.\n");
            Fatal_Flag = 1;
        }
        if (here->BSIM4v7ngcon < 1.0)
        {   fprintf(fplog, "Fatal: The parameter ngcon cannot be smaller than one.\n");
            printf("Fatal: The parameter ngcon cannot be smaller than one.\n");
            Fatal_Flag = 1;
        }
        if ((here->BSIM4v7ngcon != 1.0) && (here->BSIM4v7ngcon != 2.0))
        {   here->BSIM4v7ngcon = 1.0;
            fprintf(fplog, "Warning: Ngcon must be equal to one or two; reset to 1.0.\n");
            printf("Warning: Ngcon must be equal to one or two; reset to 1.0.\n");
        }

        if (model->BSIM4v7gbmin < 1.0e-20)
        {   fprintf(fplog, "Warning: Gbmin = %g is too small.\n",
                    model->BSIM4v7gbmin);
            printf("Warning: Gbmin = %g is too small.\n", model->BSIM4v7gbmin);
        }

        /* Check saturation parameters */
        if (pParam->BSIM4v7fprout < 0.0)
        {   fprintf(fplog, "Fatal: fprout = %g is negative.\n",
                    pParam->BSIM4v7fprout);
            printf("Fatal: fprout = %g is negative.\n", pParam->BSIM4v7fprout);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4v7pdits < 0.0)
        {   fprintf(fplog, "Fatal: pdits = %g is negative.\n",
                    pParam->BSIM4v7pdits);
            printf("Fatal: pdits = %g is negative.\n", pParam->BSIM4v7pdits);
            Fatal_Flag = 1;
        }
        if (model->BSIM4v7pditsl < 0.0)
        {   fprintf(fplog, "Fatal: pditsl = %g is negative.\n",
                    model->BSIM4v7pditsl);
            printf("Fatal: pditsl = %g is negative.\n", model->BSIM4v7pditsl);
            Fatal_Flag = 1;
        }

        /* Check gate current parameters */
      if (model->BSIM4v7igbMod) {
        if (pParam->BSIM4v7nigbinv <= 0.0)
        {   fprintf(fplog, "Fatal: nigbinv = %g is non-positive.\n",
                    pParam->BSIM4v7nigbinv);
            printf("Fatal: nigbinv = %g is non-positive.\n", pParam->BSIM4v7nigbinv);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4v7nigbacc <= 0.0)
        {   fprintf(fplog, "Fatal: nigbacc = %g is non-positive.\n",
                    pParam->BSIM4v7nigbacc);
            printf("Fatal: nigbacc = %g is non-positive.\n", pParam->BSIM4v7nigbacc);
            Fatal_Flag = 1;
        }
      }
      if (model->BSIM4v7igcMod) {
        if (pParam->BSIM4v7nigc <= 0.0)
        {   fprintf(fplog, "Fatal: nigc = %g is non-positive.\n",
                    pParam->BSIM4v7nigc);
            printf("Fatal: nigc = %g is non-positive.\n", pParam->BSIM4v7nigc);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4v7poxedge <= 0.0)
        {   fprintf(fplog, "Fatal: poxedge = %g is non-positive.\n",
                    pParam->BSIM4v7poxedge);
            printf("Fatal: poxedge = %g is non-positive.\n", pParam->BSIM4v7poxedge);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4v7pigcd <= 0.0)
        {   fprintf(fplog, "Fatal: pigcd = %g is non-positive.\n",
                    pParam->BSIM4v7pigcd);
            printf("Fatal: pigcd = %g is non-positive.\n", pParam->BSIM4v7pigcd);
            Fatal_Flag = 1;
        }
      }

        /* Check capacitance parameters */
        if (pParam->BSIM4v7clc < 0.0)
        {   fprintf(fplog, "Fatal: Clc = %g is negative.\n", pParam->BSIM4v7clc);
            printf("Fatal: Clc = %g is negative.\n", pParam->BSIM4v7clc);
            Fatal_Flag = 1;
        }      

        /* Check overlap capacitance parameters */
        if (pParam->BSIM4v7ckappas < 0.02)
        {   fprintf(fplog, "Warning: ckappas = %g is too small. Set to 0.02\n",
                    pParam->BSIM4v7ckappas);
            printf("Warning: ckappas = %g is too small.\n", pParam->BSIM4v7ckappas);
            pParam->BSIM4v7ckappas = 0.02;
       }
        if (pParam->BSIM4v7ckappad < 0.02)
        {   fprintf(fplog, "Warning: ckappad = %g is too small. Set to 0.02\n",
                    pParam->BSIM4v7ckappad);
            printf("Warning: ckappad = %g is too small.\n", pParam->BSIM4v7ckappad);
            pParam->BSIM4v7ckappad = 0.02;
        }

        if (model->BSIM4v7vtss < 0.0)
        {   fprintf(fplog, "Fatal: Vtss = %g is negative.\n",
                        model->BSIM4v7vtss);
            printf("Fatal: Vtss = %g is negative.\n",
                        model->BSIM4v7vtss);
            Fatal_Flag = 1;
        }
        if (model->BSIM4v7vtsd < 0.0)
        {   fprintf(fplog, "Fatal: Vtsd = %g is negative.\n",
                        model->BSIM4v7vtsd);
            printf("Fatal: Vtsd = %g is negative.\n",
                    model->BSIM4v7vtsd);
            Fatal_Flag = 1;
        }
        if (model->BSIM4v7vtssws < 0.0)
        {   fprintf(fplog, "Fatal: Vtssws = %g is negative.\n",
                    model->BSIM4v7vtssws);
            printf("Fatal: Vtssws = %g is negative.\n",
                    model->BSIM4v7vtssws);
            Fatal_Flag = 1;
        }
        if (model->BSIM4v7vtsswd < 0.0)
        {   fprintf(fplog, "Fatal: Vtsswd = %g is negative.\n",
                    model->BSIM4v7vtsswd);
            printf("Fatal: Vtsswd = %g is negative.\n",
                    model->BSIM4v7vtsswd);
            Fatal_Flag = 1;
        }
        if (model->BSIM4v7vtsswgs < 0.0)
        {   fprintf(fplog, "Fatal: Vtsswgs = %g is negative.\n",
                    model->BSIM4v7vtsswgs);
            printf("Fatal: Vtsswgs = %g is negative.\n",
                    model->BSIM4v7vtsswgs);
            Fatal_Flag = 1;
        }
        if (model->BSIM4v7vtsswgd < 0.0)
        {   fprintf(fplog, "Fatal: Vtsswgd = %g is negative.\n",
                    model->BSIM4v7vtsswgd);
            printf("Fatal: Vtsswgd = %g is negative.\n",
                        model->BSIM4v7vtsswgd);
            Fatal_Flag = 1;
        }


      if (model->BSIM4v7paramChk ==1)
      {
/* Check L and W parameters */ 
        if (pParam->BSIM4v7leff <= 1.0e-9)
        {   fprintf(fplog, "Warning: Leff = %g <= 1.0e-9. Recommended Leff >= 1e-8 \n", 
                    pParam->BSIM4v7leff);
            printf("Warning: Leff = %g <= 1.0e-9. Recommended Leff >= 1e-8 \n",
                    pParam->BSIM4v7leff);
        }    
        
        if (pParam->BSIM4v7leffCV <= 1.0e-9)
        {   fprintf(fplog, "Warning: Leff for CV = %g <= 1.0e-9. Recommended LeffCV >=1e-8 \n",
                    pParam->BSIM4v7leffCV);
            printf("Warning: Leff for CV = %g <= 1.0e-9. Recommended LeffCV >=1e-8 \n",
                    pParam->BSIM4v7leffCV);
        }  
        
        if (pParam->BSIM4v7weff <= 1.0e-9)
        {   fprintf(fplog, "Warning: Weff = %g <= 1.0e-9. Recommended Weff >=1e-7 \n",
                    pParam->BSIM4v7weff);
            printf("Warning: Weff = %g <= 1.0e-9. Recommended Weff >=1e-7 \n",
                   pParam->BSIM4v7weff);
        }             
        
        if (pParam->BSIM4v7weffCV <= 1.0e-9)
        {   fprintf(fplog, "Warning: Weff for CV = %g <= 1.0e-9. Recommended WeffCV >= 1e-7 \n",
                    pParam->BSIM4v7weffCV);
            printf("Warning: Weff for CV = %g <= 1.0e-9. Recommended WeffCV >= 1e-7 \n",
                    pParam->BSIM4v7weffCV);
        }        
        
        /* Check threshold voltage parameters */
        if (model->BSIM4v7toxe < 1.0e-10)
        {   fprintf(fplog, "Warning: Toxe = %g is less than 1A. Recommended Toxe >= 5A\n",
                    model->BSIM4v7toxe);
            printf("Warning: Toxe = %g is less than 1A. Recommended Toxe >= 5A\n", model->BSIM4v7toxe);
        }
        if (model->BSIM4v7toxp < 1.0e-10)
        {   fprintf(fplog, "Warning: Toxp = %g is less than 1A. Recommended Toxp >= 5A\n",
                    model->BSIM4v7toxp);
            printf("Warning: Toxp = %g is less than 1A. Recommended Toxp >= 5A\n", model->BSIM4v7toxp);
        }
        if (model->BSIM4v7toxm < 1.0e-10)
        {   fprintf(fplog, "Warning: Toxm = %g is less than 1A. Recommended Toxm >= 5A\n",
                    model->BSIM4v7toxm);
            printf("Warning: Toxm = %g is less than 1A. Recommended Toxm >= 5A\n", model->BSIM4v7toxm);
        }

        if (pParam->BSIM4v7ndep <= 1.0e12)
        {   fprintf(fplog, "Warning: Ndep = %g may be too small.\n",
                    pParam->BSIM4v7ndep);
            printf("Warning: Ndep = %g may be too small.\n",
                   pParam->BSIM4v7ndep);
        }
        else if (pParam->BSIM4v7ndep >= 1.0e21)
        {   fprintf(fplog, "Warning: Ndep = %g may be too large.\n",
                    pParam->BSIM4v7ndep);
            printf("Warning: Ndep = %g may be too large.\n",
                   pParam->BSIM4v7ndep);
        }

         if (pParam->BSIM4v7nsub <= 1.0e14)
        {   fprintf(fplog, "Warning: Nsub = %g may be too small.\n",
                    pParam->BSIM4v7nsub);
            printf("Warning: Nsub = %g may be too small.\n",
                   pParam->BSIM4v7nsub);
        }
        else if (pParam->BSIM4v7nsub >= 1.0e21)
        {   fprintf(fplog, "Warning: Nsub = %g may be too large.\n",
                    pParam->BSIM4v7nsub);
            printf("Warning: Nsub = %g may be too large.\n",
                   pParam->BSIM4v7nsub);
        }

        if ((pParam->BSIM4v7ngate > 0.0) &&
            (pParam->BSIM4v7ngate <= 1.e18))
        {   fprintf(fplog, "Warning: Ngate = %g is less than 1.E18cm^-3.\n",
                    pParam->BSIM4v7ngate);
            printf("Warning: Ngate = %g is less than 1.E18cm^-3.\n",
                   pParam->BSIM4v7ngate);
        }
       
        if (pParam->BSIM4v7dvt0 < 0.0)
        {   fprintf(fplog, "Warning: Dvt0 = %g is negative.\n",
                    pParam->BSIM4v7dvt0);   
            printf("Warning: Dvt0 = %g is negative.\n", pParam->BSIM4v7dvt0);   
        }
            
        if (fabs(1.0e-8 / (pParam->BSIM4v7w0 + pParam->BSIM4v7weff)) > 10.0)
        {   fprintf(fplog, "Warning: (W0 + Weff) may be too small.\n");
            printf("Warning: (W0 + Weff) may be too small.\n");
        }

/* Check subthreshold parameters */
        if (pParam->BSIM4v7nfactor < 0.0)
        {   fprintf(fplog, "Warning: Nfactor = %g is negative.\n",
                    pParam->BSIM4v7nfactor);
            printf("Warning: Nfactor = %g is negative.\n", pParam->BSIM4v7nfactor);
        }
        if (pParam->BSIM4v7cdsc < 0.0)
        {   fprintf(fplog, "Warning: Cdsc = %g is negative.\n",
                    pParam->BSIM4v7cdsc);
            printf("Warning: Cdsc = %g is negative.\n", pParam->BSIM4v7cdsc);
        }
        if (pParam->BSIM4v7cdscd < 0.0)
        {   fprintf(fplog, "Warning: Cdscd = %g is negative.\n",
                    pParam->BSIM4v7cdscd);
            printf("Warning: Cdscd = %g is negative.\n", pParam->BSIM4v7cdscd);
        }
/* Check DIBL parameters */
        if (here->BSIM4v7eta0 < 0.0)
        {   fprintf(fplog, "Warning: Eta0 = %g is negative.\n",
                    here->BSIM4v7eta0); 
            printf("Warning: Eta0 = %g is negative.\n", here->BSIM4v7eta0); 
        }
              
/* Check Abulk parameters */            
         if (fabs(1.0e-8 / (pParam->BSIM4v7b1 + pParam->BSIM4v7weff)) > 10.0)
               {   fprintf(fplog, "Warning: (B1 + Weff) may be too small.\n");
                   printf("Warning: (B1 + Weff) may be too small.\n");
        }    
    

/* Check Saturation parameters */
             if (pParam->BSIM4v7a2 < 0.01)
        {   fprintf(fplog, "Warning: A2 = %g is too small. Set to 0.01.\n", pParam->BSIM4v7a2);
            printf("Warning: A2 = %g is too small. Set to 0.01.\n",
                   pParam->BSIM4v7a2);
            pParam->BSIM4v7a2 = 0.01;
        }
        else if (pParam->BSIM4v7a2 > 1.0)
        {   fprintf(fplog, "Warning: A2 = %g is larger than 1. A2 is set to 1 and A1 is set to 0.\n",
                    pParam->BSIM4v7a2);
            printf("Warning: A2 = %g is larger than 1. A2 is set to 1 and A1 is set to 0.\n",
                   pParam->BSIM4v7a2);
            pParam->BSIM4v7a2 = 1.0;
            pParam->BSIM4v7a1 = 0.0;
        }

        if (pParam->BSIM4v7prwg < 0.0)
        {   fprintf(fplog, "Warning: Prwg = %g is negative. Set to zero.\n",
                    pParam->BSIM4v7prwg);
            printf("Warning: Prwg = %g is negative. Set to zero.\n",
                   pParam->BSIM4v7prwg);
            pParam->BSIM4v7prwg = 0.0;
        }

        if (pParam->BSIM4v7rdsw < 0.0)
        {   fprintf(fplog, "Warning: Rdsw = %g is negative. Set to zero.\n",
                    pParam->BSIM4v7rdsw);
            printf("Warning: Rdsw = %g is negative. Set to zero.\n",
                   pParam->BSIM4v7rdsw);
            pParam->BSIM4v7rdsw = 0.0;
            pParam->BSIM4v7rds0 = 0.0;
        }

        if (pParam->BSIM4v7rds0 < 0.0)
        {   fprintf(fplog, "Warning: Rds at current temperature = %g is negative. Set to zero.\n",
                    pParam->BSIM4v7rds0);
            printf("Warning: Rds at current temperature = %g is negative. Set to zero.\n",
                   pParam->BSIM4v7rds0);
            pParam->BSIM4v7rds0 = 0.0;
        }

        if (pParam->BSIM4v7rdswmin < 0.0)
        {   fprintf(fplog, "Warning: Rdswmin at current temperature = %g is negative. Set to zero.\n",
                    pParam->BSIM4v7rdswmin);
            printf("Warning: Rdswmin at current temperature = %g is negative. Set to zero.\n",
                   pParam->BSIM4v7rdswmin);
            pParam->BSIM4v7rdswmin = 0.0;
        }

        if (pParam->BSIM4v7pscbe2 <= 0.0)
        {   fprintf(fplog, "Warning: Pscbe2 = %g is not positive.\n",
                    pParam->BSIM4v7pscbe2);
            printf("Warning: Pscbe2 = %g is not positive.\n", pParam->BSIM4v7pscbe2);
        }

        if (pParam->BSIM4v7vsattemp < 1.0e3)
        {   fprintf(fplog, "Warning: Vsat at current temperature = %g may be too small.\n", pParam->BSIM4v7vsattemp);
           printf("Warning: Vsat at current temperature = %g may be too small.\n", pParam->BSIM4v7vsattemp);
        }

      if((model->BSIM4v7lambdaGiven) && (pParam->BSIM4v7lambda > 0.0) ) 
      {
        if (pParam->BSIM4v7lambda > 1.0e-9)
        {   fprintf(fplog, "Warning: Lambda = %g may be too large.\n", pParam->BSIM4v7lambda);
           printf("Warning: Lambda = %g may be too large.\n", pParam->BSIM4v7lambda);
        }
      }

      if((model->BSIM4v7vtlGiven) && (pParam->BSIM4v7vtl > 0.0) )
      {  
        if (pParam->BSIM4v7vtl < 6.0e4)
        {   fprintf(fplog, "Warning: Thermal velocity vtl = %g may be too small.\n", pParam->BSIM4v7vtl);
           printf("Warning: Thermal velocity vtl = %g may be too small.\n", pParam->BSIM4v7vtl);
        }

        if (pParam->BSIM4v7xn < 3.0)
        {   fprintf(fplog, "Warning: back scattering coeff xn = %g is too small.\n", pParam->BSIM4v7xn);
            printf("Warning: back scattering coeff xn = %g is too small. Reset to 3.0 \n", pParam->BSIM4v7xn);
            pParam->BSIM4v7xn = 3.0;
        }

        if (model->BSIM4v7lc < 0.0)
        {   fprintf(fplog, "Warning: back scattering coeff lc = %g is too small.\n", model->BSIM4v7lc);
            printf("Warning: back scattering coeff lc = %g is too small. Reset to 0.0\n", model->BSIM4v7lc);
            pParam->BSIM4v7lc = 0.0;
        }
      }

        if (pParam->BSIM4v7pdibl1 < 0.0)
        {   fprintf(fplog, "Warning: Pdibl1 = %g is negative.\n",
                    pParam->BSIM4v7pdibl1);
            printf("Warning: Pdibl1 = %g is negative.\n", pParam->BSIM4v7pdibl1);
        }
        if (pParam->BSIM4v7pdibl2 < 0.0)
        {   fprintf(fplog, "Warning: Pdibl2 = %g is negative.\n",
                    pParam->BSIM4v7pdibl2);
            printf("Warning: Pdibl2 = %g is negative.\n", pParam->BSIM4v7pdibl2);
        }

/* Check stress effect parameters */        
        if((here->BSIM4v7sa > 0.0) && (here->BSIM4v7sb > 0.0) && 
               ((here->BSIM4v7nf == 1.0) || ((here->BSIM4v7nf > 1.0) && (here->BSIM4v7sd > 0.0))) )
        {   if (model->BSIM4v7lodk2 <= 0.0)
            {   fprintf(fplog, "Warning: LODK2 = %g is not positive.\n",model->BSIM4v7lodk2);
                printf("Warning: LODK2 = %g is not positive.\n",model->BSIM4v7lodk2);
            }
            if (model->BSIM4v7lodeta0 <= 0.0)
            {   fprintf(fplog, "Warning: LODETA0 = %g is not positive.\n",model->BSIM4v7lodeta0);
                    printf("Warning: LODETA0 = %g is not positive.\n",model->BSIM4v7lodeta0);
            }
         }

/* Check gate resistance parameters */        
        if (here->BSIM4v7rgateMod == 1)
        {   if (model->BSIM4v7rshg <= 0.0)
                printf("Warning: rshg should be positive for rgateMod = 1.\n");
        }
        else if (here->BSIM4v7rgateMod == 2)
        {   if (model->BSIM4v7rshg <= 0.0)
                printf("Warning: rshg <= 0.0 for rgateMod = 2.\n");
            else if (pParam->BSIM4v7xrcrg1 <= 0.0)
                     printf("Warning: xrcrg1 <= 0.0 for rgateMod = 2.\n");
         }
        if (here->BSIM4v7rgateMod == 3)
        {   if (model->BSIM4v7rshg <= 0.0)
                printf("Warning: rshg should be positive for rgateMod = 3.\n");
            else if (pParam->BSIM4v7xrcrg1 <= 0.0)
                     printf("Warning: xrcrg1 should be positive for rgateMod = 3.\n");
         }
         
/* Check capacitance parameters */
        if (pParam->BSIM4v7noff < 0.1)
        {   fprintf(fplog, "Warning: Noff = %g is too small.\n",
                    pParam->BSIM4v7noff);
            printf("Warning: Noff = %g is too small.\n", pParam->BSIM4v7noff);
        }

        if (pParam->BSIM4v7voffcv < -0.5)
        {   fprintf(fplog, "Warning: Voffcv = %g is too small.\n",
                    pParam->BSIM4v7voffcv);
            printf("Warning: Voffcv = %g is too small.\n", pParam->BSIM4v7voffcv);
        }
        if (pParam->BSIM4v7moin < 5.0)
        {   fprintf(fplog, "Warning: Moin = %g is too small.\n",
                    pParam->BSIM4v7moin);
            printf("Warning: Moin = %g is too small.\n", pParam->BSIM4v7moin);
        }
        if (pParam->BSIM4v7moin > 25.0)
        {   fprintf(fplog, "Warning: Moin = %g is too large.\n",
                    pParam->BSIM4v7moin);
            printf("Warning: Moin = %g is too large.\n", pParam->BSIM4v7moin);
        }
        if(model->BSIM4v7capMod ==2) {
                if (pParam->BSIM4v7acde < 0.1)
                {   fprintf(fplog, "Warning:  Acde = %g is too small.\n",
                            pParam->BSIM4v7acde);
                    printf("Warning: Acde = %g is too small.\n", pParam->BSIM4v7acde);
                }
                if (pParam->BSIM4v7acde > 1.6)
                {   fprintf(fplog, "Warning:  Acde = %g is too large.\n",
                            pParam->BSIM4v7acde);
                    printf("Warning: Acde = %g is too large.\n", pParam->BSIM4v7acde);
                }
        }

/* Check overlap capacitance parameters */
        if (model->BSIM4v7cgdo < 0.0)
        {   fprintf(fplog, "Warning: cgdo = %g is negative. Set to zero.\n", model->BSIM4v7cgdo);
            printf("Warning: cgdo = %g is negative. Set to zero.\n", model->BSIM4v7cgdo);
            model->BSIM4v7cgdo = 0.0;
        }      
        if (model->BSIM4v7cgso < 0.0)
        {   fprintf(fplog, "Warning: cgso = %g is negative. Set to zero.\n", model->BSIM4v7cgso);
            printf("Warning: cgso = %g is negative. Set to zero.\n", model->BSIM4v7cgso);
            model->BSIM4v7cgso = 0.0;
        }      
        if (model->BSIM4v7cgbo < 0.0)
        {   fprintf(fplog, "Warning: cgbo = %g is negative. Set to zero.\n", model->BSIM4v7cgbo);
            printf("Warning: cgbo = %g is negative. Set to zero.\n", model->BSIM4v7cgbo);
            model->BSIM4v7cgbo = 0.0;
        }      

        /* v4.7 */
        if (model->BSIM4v7tnoiMod == 1 || model->BSIM4v7tnoiMod == 2) { 
            if (model->BSIM4v7tnoia < 0.0) {   
                fprintf(fplog, "Warning: tnoia = %g is negative. Set to zero.\n", model->BSIM4v7tnoia);
                printf("Warning: tnoia = %g is negative. Set to zero.\n", model->BSIM4v7tnoia);
                model->BSIM4v7tnoia = 0.0;
            }
            if (model->BSIM4v7tnoib < 0.0) {   
                fprintf(fplog, "Warning: tnoib = %g is negative. Set to zero.\n", model->BSIM4v7tnoib);
                printf("Warning: tnoib = %g is negative. Set to zero.\n", model->BSIM4v7tnoib);
                model->BSIM4v7tnoib = 0.0;
            }
            if (model->BSIM4v7rnoia < 0.0) { 
                fprintf(fplog, "Warning: rnoia = %g is negative. Set to zero.\n", model->BSIM4v7rnoia);
                printf("Warning: rnoia = %g is negative. Set to zero.\n", model->BSIM4v7rnoia);
                model->BSIM4v7rnoia = 0.0;
            }
            if (model->BSIM4v7rnoib < 0.0) {
                fprintf(fplog, "Warning: rnoib = %g is negative. Set to zero.\n", model->BSIM4v7rnoib);
                printf("Warning: rnoib = %g is negative. Set to zero.\n", model->BSIM4v7rnoib);
                model->BSIM4v7rnoib = 0.0;
            }
        }

        /* v4.7 */
        if (model->BSIM4v7tnoiMod == 2) { 
            if (model->BSIM4v7tnoic < 0.0) {
                fprintf(fplog, "Warning: tnoic = %g is negative. Set to zero.\n", model->BSIM4v7tnoic);
                printf("Warning: tnoic = %g is negative. Set to zero.\n", model->BSIM4v7tnoic);
                model->BSIM4v7tnoic = 0.0;
            }
            if (model->BSIM4v7rnoic < 0.0) {
                fprintf(fplog, "Warning: rnoic = %g is negative. Set to zero.\n", model->BSIM4v7rnoic);
                printf("Warning: rnoic = %g is negative. Set to zero.\n", model->BSIM4v7rnoic);
                model->BSIM4v7rnoic = 0.0;
            }
        }

        /* Limits of Njs and Njd modified in BSIM4v7.7 */
        if (model->BSIM4v7SjctEmissionCoeff < 0.1) {
            fprintf(fplog, "Warning: Njs = %g is less than 0.1. Setting Njs to 0.1.\n", model->BSIM4v7SjctEmissionCoeff);
            printf("Warning: Njs = %g is less than 0.1. Setting Njs to 0.1.\n", model->BSIM4v7SjctEmissionCoeff);
            model->BSIM4v7SjctEmissionCoeff = 0.1;
        } 
        else if (model->BSIM4v7SjctEmissionCoeff < 0.7) { 
            fprintf(fplog, "Warning: Njs = %g is less than 0.7.\n", model->BSIM4v7SjctEmissionCoeff);
            printf("Warning: Njs = %g is less than 0.7.\n", model->BSIM4v7SjctEmissionCoeff);
        }
        if (model->BSIM4v7DjctEmissionCoeff < 0.1) { 
            fprintf(fplog, "Warning: Njd = %g is less than 0.1. Setting Njd to 0.1.\n", model->BSIM4v7DjctEmissionCoeff);
            printf("Warning: Njd = %g is less than 0.1. Setting Njd to 0.1.\n", model->BSIM4v7DjctEmissionCoeff);
            model->BSIM4v7DjctEmissionCoeff = 0.1;
        }
        else if (model->BSIM4v7DjctEmissionCoeff < 0.7) {
            fprintf(fplog, "Warning: Njd = %g is less than 0.7.\n", model->BSIM4v7DjctEmissionCoeff);
            printf("Warning: Njd = %g is less than 0.7.\n", model->BSIM4v7DjctEmissionCoeff);
        }

        if (model->BSIM4v7njtsstemp < 0.0)
        {   fprintf(fplog, "Warning: Njts = %g is negative at temperature = %g.\n",
                    model->BSIM4v7njtsstemp, ckt->CKTtemp);
            printf("Warning: Njts = %g is negative at temperature = %g.\n",
                   model->BSIM4v7njtsstemp, ckt->CKTtemp);
        }
        if (model->BSIM4v7njtsswstemp < 0.0)
        {   fprintf(fplog, "Warning: Njtssw = %g is negative at temperature = %g.\n",
                    model->BSIM4v7njtsswstemp, ckt->CKTtemp);
            printf("Warning: Njtssw = %g is negative at temperature = %g.\n",
                    model->BSIM4v7njtsswstemp, ckt->CKTtemp);
        }
        if (model->BSIM4v7njtsswgstemp < 0.0)
        {   fprintf(fplog, "Warning: Njtsswg = %g is negative at temperature = %g.\n",
                    model->BSIM4v7njtsswgstemp, ckt->CKTtemp);
            printf("Warning: Njtsswg = %g is negative at temperature = %g.\n",
                    model->BSIM4v7njtsswgstemp, ckt->CKTtemp);
        }

        if (model->BSIM4v7njtsdGiven && model->BSIM4v7njtsdtemp < 0.0)
        {   fprintf(fplog, "Warning: Njtsd = %g is negative at temperature = %g.\n",
                    model->BSIM4v7njtsdtemp, ckt->CKTtemp);
            printf("Warning: Njtsd = %g is negative at temperature = %g.\n",
                    model->BSIM4v7njtsdtemp, ckt->CKTtemp);
        }
        if (model->BSIM4v7njtsswdGiven && model->BSIM4v7njtsswdtemp < 0.0)
        {   fprintf(fplog, "Warning: Njtsswd = %g is negative at temperature = %g.\n",
                    model->BSIM4v7njtsswdtemp, ckt->CKTtemp);
            printf("Warning: Njtsswd = %g is negative at temperature = %g.\n",
                    model->BSIM4v7njtsswdtemp, ckt->CKTtemp);
        }
        if (model->BSIM4v7njtsswgdGiven && model->BSIM4v7njtsswgdtemp < 0.0)
        {   fprintf(fplog, "Warning: Njtsswgd = %g is negative at temperature = %g.\n",
                    model->BSIM4v7njtsswgdtemp, ckt->CKTtemp);
            printf("Warning: Njtsswgd = %g is negative at temperature = %g.\n",
                    model->BSIM4v7njtsswgdtemp, ckt->CKTtemp);
        }

       if (model->BSIM4v7ntnoi < 0.0)
        {   fprintf(fplog, "Warning: ntnoi = %g is negative. Set to zero.\n", model->BSIM4v7ntnoi);
            printf("Warning: ntnoi = %g is negative. Set to zero.\n", model->BSIM4v7ntnoi);
            model->BSIM4v7ntnoi = 0.0;
        }

        /* diode model */
       if (model->BSIM4v7SbulkJctBotGradingCoeff >= 0.99)
        {   fprintf(fplog, "Warning: MJS = %g is too big. Set to 0.99.\n", model->BSIM4v7SbulkJctBotGradingCoeff);
            printf("Warning: MJS = %g is too big. Set to 0.99.\n", model->BSIM4v7SbulkJctBotGradingCoeff);
            model->BSIM4v7SbulkJctBotGradingCoeff = 0.99;
        }
       if (model->BSIM4v7SbulkJctSideGradingCoeff >= 0.99)
        {   fprintf(fplog, "Warning: MJSWS = %g is too big. Set to 0.99.\n", model->BSIM4v7SbulkJctSideGradingCoeff);
            printf("Warning: MJSWS = %g is too big. Set to 0.99.\n", model->BSIM4v7SbulkJctSideGradingCoeff);
            model->BSIM4v7SbulkJctSideGradingCoeff = 0.99;
        }
       if (model->BSIM4v7SbulkJctGateSideGradingCoeff >= 0.99)
        {   fprintf(fplog, "Warning: MJSWGS = %g is too big. Set to 0.99.\n", model->BSIM4v7SbulkJctGateSideGradingCoeff);
            printf("Warning: MJSWGS = %g is too big. Set to 0.99.\n", model->BSIM4v7SbulkJctGateSideGradingCoeff);
            model->BSIM4v7SbulkJctGateSideGradingCoeff = 0.99;
        }

       if (model->BSIM4v7DbulkJctBotGradingCoeff >= 0.99)
        {   fprintf(fplog, "Warning: MJD = %g is too big. Set to 0.99.\n", model->BSIM4v7DbulkJctBotGradingCoeff);
            printf("Warning: MJD = %g is too big. Set to 0.99.\n", model->BSIM4v7DbulkJctBotGradingCoeff);
            model->BSIM4v7DbulkJctBotGradingCoeff = 0.99;
        }
       if (model->BSIM4v7DbulkJctSideGradingCoeff >= 0.99)
        {   fprintf(fplog, "Warning: MJSWD = %g is too big. Set to 0.99.\n", model->BSIM4v7DbulkJctSideGradingCoeff);
            printf("Warning: MJSWD = %g is too big. Set to 0.99.\n", model->BSIM4v7DbulkJctSideGradingCoeff);
            model->BSIM4v7DbulkJctSideGradingCoeff = 0.99;
        }
       if (model->BSIM4v7DbulkJctGateSideGradingCoeff >= 0.99)
        {   fprintf(fplog, "Warning: MJSWGD = %g is too big. Set to 0.99.\n", model->BSIM4v7DbulkJctGateSideGradingCoeff);
            printf("Warning: MJSWGD = %g is too big. Set to 0.99.\n", model->BSIM4v7DbulkJctGateSideGradingCoeff);
            model->BSIM4v7DbulkJctGateSideGradingCoeff = 0.99;
        }
        if (model->BSIM4v7wpemod == 1)
        {
                if (model->BSIM4v7scref <= 0.0)
                {   fprintf(fplog, "Warning: SCREF = %g is not positive. Set to 1e-6.\n", model->BSIM4v7scref);
                    printf("Warning: SCREF = %g is not positive. Set to 1e-6.\n", model->BSIM4v7scref);
                    model->BSIM4v7scref = 1e-6;
                }
                /*Move these checks to temp.c for sceff calculation*/
                        /*
                        if (here->BSIM4v7sca < 0.0)
                {   fprintf(fplog, "Warning: SCA = %g is negative. Set to 0.0.\n", here->BSIM4v7sca);
                    printf("Warning: SCA = %g is negative. Set to 0.0.\n", here->BSIM4v7sca);
                    here->BSIM4v7sca = 0.0;
                }
                if (here->BSIM4v7scb < 0.0)
                {   fprintf(fplog, "Warning: SCB = %g is negative. Set to 0.0.\n", here->BSIM4v7scb);
                    printf("Warning: SCB = %g is negative. Set to 0.0.\n", here->BSIM4v7scb);
                    here->BSIM4v7scb = 0.0;
                }
                if (here->BSIM4v7scc < 0.0)
                {   fprintf(fplog, "Warning: SCC = %g is negative. Set to 0.0.\n", here->BSIM4v7scc);
                    printf("Warning: SCC = %g is negative. Set to 0.0.\n", here->BSIM4v7scc);
                    here->BSIM4v7scc = 0.0;
                }
                if (here->BSIM4v7sc < 0.0)
                {   fprintf(fplog, "Warning: SC = %g is negative. Set to 0.0.\n", here->BSIM4v7sc);
                    printf("Warning: SC = %g is negative. Set to 0.0.\n", here->BSIM4v7sc);
                    here->BSIM4v7sc = 0.0;
                }
                                */
                                
        }
     }/* loop for the parameter check for warning messages */      
        fclose(fplog);
    }
    else
    {   fprintf(stderr, "Warning: Can't open log file. Parameter checking skipped.\n");
    }

    return(Fatal_Flag);
}

