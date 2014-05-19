/**** BSIM4.8.0 Released by Navid Paydavosi 11/01/2013 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4check.c of BSIM4.8.0.
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
 * Modified by Pankaj Kumar Thakur, 07/23/2012
 * Modified by Navid Paydavosi, 08/21/2013
 **********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "bsim4def.h"
#include "ngspice/trandefs.h"
#include "ngspice/const.h"
#include "ngspice/sperror.h"
#include "ngspice/devdefs.h"
#include "ngspice/suffix.h"

int
BSIM4checkModel(
BSIM4model *model,
BSIM4instance *here,
CKTcircuit *ckt)
{
struct bsim4SizeDependParam *pParam;
int Fatal_Flag = 0;
FILE *fplog;
    
    if ((fplog = fopen("bsim4.out", "w")) != NULL)
    {   pParam = here->pParam;
        fprintf(fplog, "BSIM4: Berkeley Short Channel IGFET Model-4\n");
        fprintf(fplog, "Developed by Xuemei (Jane) Xi, Mohan Dunga, Prof. Ali Niknejad and Prof. Chenming Hu in 2003.\n");
        fprintf(fplog, "\n");
        fprintf(fplog, "++++++++++ BSIM4 PARAMETER CHECKING BELOW ++++++++++\n");

        if ((strcmp(model->BSIM4version, "4.8.0")) && (strcmp(model->BSIM4version, "4.80")) && (strcmp(model->BSIM4version, "4.8")))
        {  fprintf(fplog, "Warning: This model is BSIM4.8.0; you specified a wrong version number.\n");
           printf("Warning: This model is BSIM4.8.0; you specified a wrong version number.\n");
        }
        fprintf(fplog, "Model = %s\n", model->BSIM4modName);


        if ((here->BSIM4rgateMod == 2) || (here->BSIM4rgateMod == 3))
        {   if ((here->BSIM4trnqsMod == 1) || (here->BSIM4acnqsMod == 1))
            {   fprintf(fplog, "Warning: You've selected both Rg and charge deficit NQS; select one only.\n");
                printf("Warning: You've selected both Rg and charge deficit NQS; select one only.\n");
            }
        }

        if (model->BSIM4toxe <= 0.0)
        {   fprintf(fplog, "Fatal: Toxe = %g is not positive.\n",
                    model->BSIM4toxe);
            printf("Fatal: Toxe = %g is not positive.\n", model->BSIM4toxe);
            Fatal_Flag = 1;
        }
        if (here->BSIM4toxp <= 0.0)
        {   fprintf(fplog, "Fatal: Toxp = %g is not positive.\n",
                    here->BSIM4toxp);
            printf("Fatal: Toxp = %g is not positive.\n", here->BSIM4toxp);
            Fatal_Flag = 1;
        }
        if (model->BSIM4eot <= 0.0)
        {   fprintf(fplog, "Fatal: EOT = %g is not positive.\n",
                    model->BSIM4eot);
            printf("Fatal: EOT = %g is not positive.\n", model->BSIM4eot);
            Fatal_Flag = 1;
        }
        if (model->BSIM4epsrgate < 0.0)
        {   fprintf(fplog, "Fatal: Epsrgate = %g is not positive.\n",
                    model->BSIM4epsrgate);
            printf("Fatal: Epsrgate = %g is not positive.\n", model->BSIM4epsrgate);
            Fatal_Flag = 1;
        }
        if (model->BSIM4epsrsub < 0.0)
        {   fprintf(fplog, "Fatal: Epsrsub = %g is not positive.\n",
                    model->BSIM4epsrsub);
            printf("Fatal: Epsrsub = %g is not positive.\n", model->BSIM4epsrsub);
            Fatal_Flag = 1;
        }
        if (model->BSIM4easub < 0.0)
        {   fprintf(fplog, "Fatal: Easub = %g is not positive.\n",
                    model->BSIM4easub);
            printf("Fatal: Easub = %g is not positive.\n", model->BSIM4easub);
            Fatal_Flag = 1;
        }
        if (model->BSIM4ni0sub <= 0.0)
        {   fprintf(fplog, "Fatal: Ni0sub = %g is not positive.\n",
                    model->BSIM4ni0sub);
            printf("Fatal: Easub = %g is not positive.\n", model->BSIM4ni0sub);
            Fatal_Flag = 1;
        }

        if (model->BSIM4toxm <= 0.0)
        {   fprintf(fplog, "Fatal: Toxm = %g is not positive.\n",
                    model->BSIM4toxm);
            printf("Fatal: Toxm = %g is not positive.\n", model->BSIM4toxm);
            Fatal_Flag = 1;
        }

        if (model->BSIM4toxref <= 0.0)
        {   fprintf(fplog, "Fatal: Toxref = %g is not positive.\n",
                    model->BSIM4toxref);
            printf("Fatal: Toxref = %g is not positive.\n", model->BSIM4toxref);
            Fatal_Flag = 1;
        }

        if (pParam->BSIM4lpe0 < -pParam->BSIM4leff)
        {   fprintf(fplog, "Fatal: Lpe0 = %g is less than -Leff.\n",
                    pParam->BSIM4lpe0);
            printf("Fatal: Lpe0 = %g is less than -Leff.\n",
                        pParam->BSIM4lpe0);
            Fatal_Flag = 1;
        }
        if (model->BSIM4lintnoi > pParam->BSIM4leff/2)
        {   fprintf(fplog, "Fatal: Lintnoi = %g is too large - Leff for noise is negative.\n",
                    model->BSIM4lintnoi);
            printf("Fatal: Lintnoi = %g is too large - Leff for noise is negative.\n",
                    model->BSIM4lintnoi);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4lpeb < -pParam->BSIM4leff)
        {   fprintf(fplog, "Fatal: Lpeb = %g is less than -Leff.\n",
                    pParam->BSIM4lpeb);
            printf("Fatal: Lpeb = %g is less than -Leff.\n",
                        pParam->BSIM4lpeb);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4ndep <= 0.0)
        {   fprintf(fplog, "Fatal: Ndep = %g is not positive.\n",
                    pParam->BSIM4ndep);
            printf("Fatal: Ndep = %g is not positive.\n",
                   pParam->BSIM4ndep);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4phi <= 0.0)
        {   fprintf(fplog, "Fatal: Phi = %g is not positive. Please check Phin and Ndep\n",
                    pParam->BSIM4phi);
            fprintf(fplog, "           Phin = %g  Ndep = %g \n", 
                        pParam->BSIM4phin, pParam->BSIM4ndep);
            printf("Fatal: Phi = %g is not positive. Please check Phin and Ndep\n",
                    pParam->BSIM4phi);
            printf("           Phin = %g  Ndep = %g \n", 
                        pParam->BSIM4phin, pParam->BSIM4ndep);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4nsub <= 0.0)
        {   fprintf(fplog, "Fatal: Nsub = %g is not positive.\n",
                    pParam->BSIM4nsub);
            printf("Fatal: Nsub = %g is not positive.\n",
                   pParam->BSIM4nsub);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4ngate < 0.0)
        {   fprintf(fplog, "Fatal: Ngate = %g is not positive.\n",
                    pParam->BSIM4ngate);
            printf("Fatal: Ngate = %g Ngate is not positive.\n",
                   pParam->BSIM4ngate);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4ngate > 1.e25)
        {   fprintf(fplog, "Fatal: Ngate = %g is too high.\n",
                    pParam->BSIM4ngate);
            printf("Fatal: Ngate = %g Ngate is too high\n",
                   pParam->BSIM4ngate);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4xj <= 0.0)
        {   fprintf(fplog, "Fatal: Xj = %g is not positive.\n",
                    pParam->BSIM4xj);
            printf("Fatal: Xj = %g is not positive.\n", pParam->BSIM4xj);
            Fatal_Flag = 1;
        }

        if (pParam->BSIM4dvt1 < 0.0)
        {   fprintf(fplog, "Fatal: Dvt1 = %g is negative.\n",
                    pParam->BSIM4dvt1);   
            printf("Fatal: Dvt1 = %g is negative.\n", pParam->BSIM4dvt1);   
            Fatal_Flag = 1;
        }
            
        if (pParam->BSIM4dvt1w < 0.0)
        {   fprintf(fplog, "Fatal: Dvt1w = %g is negative.\n",
                    pParam->BSIM4dvt1w);
            printf("Fatal: Dvt1w = %g is negative.\n", pParam->BSIM4dvt1w);
            Fatal_Flag = 1;
        }
            
        if (pParam->BSIM4w0 == -pParam->BSIM4weff)
        {   fprintf(fplog, "Fatal: (W0 + Weff) = 0 causing divided-by-zero.\n");
            printf("Fatal: (W0 + Weff) = 0 causing divided-by-zero.\n");
            Fatal_Flag = 1;
        }   

        if (pParam->BSIM4dsub < 0.0)
        {   fprintf(fplog, "Fatal: Dsub = %g is negative.\n", pParam->BSIM4dsub);
            printf("Fatal: Dsub = %g is negative.\n", pParam->BSIM4dsub);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4b1 == -pParam->BSIM4weff)
        {   fprintf(fplog, "Fatal: (B1 + Weff) = 0 causing divided-by-zero.\n");
            printf("Fatal: (B1 + Weff) = 0 causing divided-by-zero.\n");
            Fatal_Flag = 1;
        }  
        if (here->BSIM4u0temp <= 0.0)
        {   fprintf(fplog, "Fatal: u0 at current temperature = %g is not positive.\n", here->BSIM4u0temp);
            printf("Fatal: u0 at current temperature = %g is not positive.\n",
                   here->BSIM4u0temp);
            Fatal_Flag = 1;
        }
    
        if (pParam->BSIM4delta < 0.0)
        {   fprintf(fplog, "Fatal: Delta = %g is less than zero.\n",
                    pParam->BSIM4delta);
            printf("Fatal: Delta = %g is less than zero.\n", pParam->BSIM4delta);
            Fatal_Flag = 1;
        }      

        if (here->BSIM4vsattemp <= 0.0)
        {   fprintf(fplog, "Fatal: Vsat at current temperature = %g is not positive.\n", here->BSIM4vsattemp);
            printf("Fatal: Vsat at current temperature = %g is not positive.\n",
                   here->BSIM4vsattemp);
            Fatal_Flag = 1;
        }

        if (pParam->BSIM4pclm <= 0.0)
        {   fprintf(fplog, "Fatal: Pclm = %g is not positive.\n", pParam->BSIM4pclm);
            printf("Fatal: Pclm = %g is not positive.\n", pParam->BSIM4pclm);
            Fatal_Flag = 1;
        }

        if (pParam->BSIM4drout < 0.0)
        {   fprintf(fplog, "Fatal: Drout = %g is negative.\n", pParam->BSIM4drout);
            printf("Fatal: Drout = %g is negative.\n", pParam->BSIM4drout);
            Fatal_Flag = 1;
        }

        if (here->BSIM4m < 1.0)
        {   fprintf(fplog, "Fatal: Number of multiplier = %g is smaller than one.\n", here->BSIM4m);
            printf("Fatal: Number of multiplier = %g is smaller than one.\n", here->BSIM4m);
            Fatal_Flag = 1;
        }

        if (here->BSIM4nf < 1.0)
        {   fprintf(fplog, "Fatal: Number of finger = %g is smaller than one.\n", here->BSIM4nf);
            printf("Fatal: Number of finger = %g is smaller than one.\n", here->BSIM4nf);
            Fatal_Flag = 1;
        }

        if((here->BSIM4sa > 0.0) && (here->BSIM4sb > 0.0) && 
               ((here->BSIM4nf == 1.0) || ((here->BSIM4nf > 1.0) && (here->BSIM4sd > 0.0))) )
        {   if (model->BSIM4saref <= 0.0)
            {   fprintf(fplog, "Fatal: SAref = %g is not positive.\n",model->BSIM4saref);
                     printf("Fatal: SAref = %g is not positive.\n",model->BSIM4saref);
                     Fatal_Flag = 1;
            }
            if (model->BSIM4sbref <= 0.0)
            {   fprintf(fplog, "Fatal: SBref = %g is not positive.\n",model->BSIM4sbref);
                    printf("Fatal: SBref = %g is not positive.\n",model->BSIM4sbref);
                    Fatal_Flag = 1;
            }
         }

        if ((here->BSIM4l + model->BSIM4xl) <= model->BSIM4xgl)
        {   fprintf(fplog, "Fatal: The parameter xgl must be smaller than Ldrawn+XL.\n");
            printf("Fatal: The parameter xgl must be smaller than Ldrawn+XL.\n");
            Fatal_Flag = 1;
        }
        if (here->BSIM4ngcon < 1.0)
        {   fprintf(fplog, "Fatal: The parameter ngcon cannot be smaller than one.\n");
            printf("Fatal: The parameter ngcon cannot be smaller than one.\n");
            Fatal_Flag = 1;
        }
        if ((here->BSIM4ngcon != 1.0) && (here->BSIM4ngcon != 2.0))
        {   here->BSIM4ngcon = 1.0;
            fprintf(fplog, "Warning: Ngcon must be equal to one or two; reset to 1.0.\n");
            printf("Warning: Ngcon must be equal to one or two; reset to 1.0.\n");
        }

        if (model->BSIM4gbmin < 1.0e-20)
        {   fprintf(fplog, "Warning: Gbmin = %g is too small.\n",
                    model->BSIM4gbmin);
            printf("Warning: Gbmin = %g is too small.\n", model->BSIM4gbmin);
        }

        /* Check saturation parameters */
        if (pParam->BSIM4fprout < 0.0)
        {   fprintf(fplog, "Fatal: fprout = %g is negative.\n",
                    pParam->BSIM4fprout);
            printf("Fatal: fprout = %g is negative.\n", pParam->BSIM4fprout);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4pdits < 0.0)
        {   fprintf(fplog, "Fatal: pdits = %g is negative.\n",
                    pParam->BSIM4pdits);
            printf("Fatal: pdits = %g is negative.\n", pParam->BSIM4pdits);
            Fatal_Flag = 1;
        }
        if (model->BSIM4pditsl < 0.0)
        {   fprintf(fplog, "Fatal: pditsl = %g is negative.\n",
                    model->BSIM4pditsl);
            printf("Fatal: pditsl = %g is negative.\n", model->BSIM4pditsl);
            Fatal_Flag = 1;
        }

        /* Check gate current parameters */
      if (model->BSIM4igbMod) {
        if (pParam->BSIM4nigbinv <= 0.0)
        {   fprintf(fplog, "Fatal: nigbinv = %g is non-positive.\n",
                    pParam->BSIM4nigbinv);
            printf("Fatal: nigbinv = %g is non-positive.\n", pParam->BSIM4nigbinv);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4nigbacc <= 0.0)
        {   fprintf(fplog, "Fatal: nigbacc = %g is non-positive.\n",
                    pParam->BSIM4nigbacc);
            printf("Fatal: nigbacc = %g is non-positive.\n", pParam->BSIM4nigbacc);
            Fatal_Flag = 1;
        }
      }
      if (model->BSIM4igcMod) {
        if (pParam->BSIM4nigc <= 0.0)
        {   fprintf(fplog, "Fatal: nigc = %g is non-positive.\n",
                    pParam->BSIM4nigc);
            printf("Fatal: nigc = %g is non-positive.\n", pParam->BSIM4nigc);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4poxedge <= 0.0)
        {   fprintf(fplog, "Fatal: poxedge = %g is non-positive.\n",
                    pParam->BSIM4poxedge);
            printf("Fatal: poxedge = %g is non-positive.\n", pParam->BSIM4poxedge);
            Fatal_Flag = 1;
        }
        if (pParam->BSIM4pigcd <= 0.0)
        {   fprintf(fplog, "Fatal: pigcd = %g is non-positive.\n",
                    pParam->BSIM4pigcd);
            printf("Fatal: pigcd = %g is non-positive.\n", pParam->BSIM4pigcd);
            Fatal_Flag = 1;
        }
      }

        /* Check capacitance parameters */
        if (pParam->BSIM4clc < 0.0)
        {   fprintf(fplog, "Fatal: Clc = %g is negative.\n", pParam->BSIM4clc);
            printf("Fatal: Clc = %g is negative.\n", pParam->BSIM4clc);
            Fatal_Flag = 1;
        }      

        /* Check overlap capacitance parameters */
        if (pParam->BSIM4ckappas < 0.02)
        {   fprintf(fplog, "Warning: ckappas = %g is too small. Set to 0.02\n",
                    pParam->BSIM4ckappas);
            printf("Warning: ckappas = %g is too small.\n", pParam->BSIM4ckappas);
            pParam->BSIM4ckappas = 0.02;
       }
        if (pParam->BSIM4ckappad < 0.02)
        {   fprintf(fplog, "Warning: ckappad = %g is too small. Set to 0.02\n",
                    pParam->BSIM4ckappad);
            printf("Warning: ckappad = %g is too small.\n", pParam->BSIM4ckappad);
            pParam->BSIM4ckappad = 0.02;
        }

        if (model->BSIM4vtss < 0.0)
        {   fprintf(fplog, "Fatal: Vtss = %g is negative.\n",
                        model->BSIM4vtss);
            printf("Fatal: Vtss = %g is negative.\n",
                        model->BSIM4vtss);
            Fatal_Flag = 1;
        }
        if (model->BSIM4vtsd < 0.0)
        {   fprintf(fplog, "Fatal: Vtsd = %g is negative.\n",
                        model->BSIM4vtsd);
            printf("Fatal: Vtsd = %g is negative.\n",
                    model->BSIM4vtsd);
            Fatal_Flag = 1;
        }
        if (model->BSIM4vtssws < 0.0)
        {   fprintf(fplog, "Fatal: Vtssws = %g is negative.\n",
                    model->BSIM4vtssws);
            printf("Fatal: Vtssws = %g is negative.\n",
                    model->BSIM4vtssws);
            Fatal_Flag = 1;
        }
        if (model->BSIM4vtsswd < 0.0)
        {   fprintf(fplog, "Fatal: Vtsswd = %g is negative.\n",
                    model->BSIM4vtsswd);
            printf("Fatal: Vtsswd = %g is negative.\n",
                    model->BSIM4vtsswd);
            Fatal_Flag = 1;
        }
        if (model->BSIM4vtsswgs < 0.0)
        {   fprintf(fplog, "Fatal: Vtsswgs = %g is negative.\n",
                    model->BSIM4vtsswgs);
            printf("Fatal: Vtsswgs = %g is negative.\n",
                    model->BSIM4vtsswgs);
            Fatal_Flag = 1;
        }
        if (model->BSIM4vtsswgd < 0.0)
        {   fprintf(fplog, "Fatal: Vtsswgd = %g is negative.\n",
                    model->BSIM4vtsswgd);
            printf("Fatal: Vtsswgd = %g is negative.\n",
                        model->BSIM4vtsswgd);
            Fatal_Flag = 1;
        }


      if (model->BSIM4paramChk ==1)
      {
/* Check L and W parameters */ 
        if (pParam->BSIM4leff <= 1.0e-9)
        {   fprintf(fplog, "Warning: Leff = %g <= 1.0e-9. Recommended Leff >= 1e-8 \n", 
                    pParam->BSIM4leff);
            printf("Warning: Leff = %g <= 1.0e-9. Recommended Leff >= 1e-8 \n",
                    pParam->BSIM4leff);
        }    
        
        if (pParam->BSIM4leffCV <= 1.0e-9)
        {   fprintf(fplog, "Warning: Leff for CV = %g <= 1.0e-9. Recommended LeffCV >=1e-8 \n",
                    pParam->BSIM4leffCV);
            printf("Warning: Leff for CV = %g <= 1.0e-9. Recommended LeffCV >=1e-8 \n",
                    pParam->BSIM4leffCV);
        }  
        
        if (pParam->BSIM4weff <= 1.0e-9)
        {   fprintf(fplog, "Warning: Weff = %g <= 1.0e-9. Recommended Weff >=1e-7 \n",
                    pParam->BSIM4weff);
            printf("Warning: Weff = %g <= 1.0e-9. Recommended Weff >=1e-7 \n",
                   pParam->BSIM4weff);
        }             
        
        if (pParam->BSIM4weffCV <= 1.0e-9)
        {   fprintf(fplog, "Warning: Weff for CV = %g <= 1.0e-9. Recommended WeffCV >= 1e-7 \n",
                    pParam->BSIM4weffCV);
            printf("Warning: Weff for CV = %g <= 1.0e-9. Recommended WeffCV >= 1e-7 \n",
                    pParam->BSIM4weffCV);
        }        
        
        /* Check threshold voltage parameters */
        if (model->BSIM4toxe < 1.0e-10)
        {   fprintf(fplog, "Warning: Toxe = %g is less than 1A. Recommended Toxe >= 5A\n",
                    model->BSIM4toxe);
            printf("Warning: Toxe = %g is less than 1A. Recommended Toxe >= 5A\n", model->BSIM4toxe);
        }
        if (here->BSIM4toxp < 1.0e-10)
        {   fprintf(fplog, "Warning: Toxp = %g is less than 1A. Recommended Toxp >= 5A\n",
                    here->BSIM4toxp);
            printf("Warning: Toxp = %g is less than 1A. Recommended Toxp >= 5A\n", here->BSIM4toxp);
        }
        if (model->BSIM4toxm < 1.0e-10)
        {   fprintf(fplog, "Warning: Toxm = %g is less than 1A. Recommended Toxm >= 5A\n",
                    model->BSIM4toxm);
            printf("Warning: Toxm = %g is less than 1A. Recommended Toxm >= 5A\n", model->BSIM4toxm);
        }

        if (pParam->BSIM4ndep <= 1.0e12)
        {   fprintf(fplog, "Warning: Ndep = %g may be too small.\n",
                    pParam->BSIM4ndep);
            printf("Warning: Ndep = %g may be too small.\n",
                   pParam->BSIM4ndep);
        }
        else if (pParam->BSIM4ndep >= 1.0e21)
        {   fprintf(fplog, "Warning: Ndep = %g may be too large.\n",
                    pParam->BSIM4ndep);
            printf("Warning: Ndep = %g may be too large.\n",
                   pParam->BSIM4ndep);
        }

         if (pParam->BSIM4nsub <= 1.0e14)
        {   fprintf(fplog, "Warning: Nsub = %g may be too small.\n",
                    pParam->BSIM4nsub);
            printf("Warning: Nsub = %g may be too small.\n",
                   pParam->BSIM4nsub);
        }
        else if (pParam->BSIM4nsub >= 1.0e21)
        {   fprintf(fplog, "Warning: Nsub = %g may be too large.\n",
                    pParam->BSIM4nsub);
            printf("Warning: Nsub = %g may be too large.\n",
                   pParam->BSIM4nsub);
        }

        if ((pParam->BSIM4ngate > 0.0) &&
            (pParam->BSIM4ngate <= 1.e18))
        {   fprintf(fplog, "Warning: Ngate = %g is less than 1.E18cm^-3.\n",
                    pParam->BSIM4ngate);
            printf("Warning: Ngate = %g is less than 1.E18cm^-3.\n",
                   pParam->BSIM4ngate);
        }
       
        if (pParam->BSIM4dvt0 < 0.0)
        {   fprintf(fplog, "Warning: Dvt0 = %g is negative.\n",
                    pParam->BSIM4dvt0);   
            printf("Warning: Dvt0 = %g is negative.\n", pParam->BSIM4dvt0);   
        }
            
        if (fabs(1.0e-8 / (pParam->BSIM4w0 + pParam->BSIM4weff)) > 10.0)
        {   fprintf(fplog, "Warning: (W0 + Weff) may be too small.\n");
            printf("Warning: (W0 + Weff) may be too small.\n");
        }

/* Check subthreshold parameters */
        if (pParam->BSIM4nfactor < 0.0)
        {   fprintf(fplog, "Warning: Nfactor = %g is negative.\n",
                    pParam->BSIM4nfactor);
            printf("Warning: Nfactor = %g is negative.\n", pParam->BSIM4nfactor);
        }
        if (pParam->BSIM4cdsc < 0.0)
        {   fprintf(fplog, "Warning: Cdsc = %g is negative.\n",
                    pParam->BSIM4cdsc);
            printf("Warning: Cdsc = %g is negative.\n", pParam->BSIM4cdsc);
        }
        if (pParam->BSIM4cdscd < 0.0)
        {   fprintf(fplog, "Warning: Cdscd = %g is negative.\n",
                    pParam->BSIM4cdscd);
            printf("Warning: Cdscd = %g is negative.\n", pParam->BSIM4cdscd);
        }
/* Check DIBL parameters */
        if (here->BSIM4eta0 < 0.0)
        {   fprintf(fplog, "Warning: Eta0 = %g is negative.\n",
                    here->BSIM4eta0); 
            printf("Warning: Eta0 = %g is negative.\n", here->BSIM4eta0); 
        }
              
/* Check Abulk parameters */            
         if (fabs(1.0e-8 / (pParam->BSIM4b1 + pParam->BSIM4weff)) > 10.0)
               {   fprintf(fplog, "Warning: (B1 + Weff) may be too small.\n");
                   printf("Warning: (B1 + Weff) may be too small.\n");
        }    
    

/* Check Saturation parameters */
             if (pParam->BSIM4a2 < 0.01)
        {   fprintf(fplog, "Warning: A2 = %g is too small. Set to 0.01.\n", pParam->BSIM4a2);
            printf("Warning: A2 = %g is too small. Set to 0.01.\n",
                   pParam->BSIM4a2);
            pParam->BSIM4a2 = 0.01;
        }
        else if (pParam->BSIM4a2 > 1.0)
        {   fprintf(fplog, "Warning: A2 = %g is larger than 1. A2 is set to 1 and A1 is set to 0.\n",
                    pParam->BSIM4a2);
            printf("Warning: A2 = %g is larger than 1. A2 is set to 1 and A1 is set to 0.\n",
                   pParam->BSIM4a2);
            pParam->BSIM4a2 = 1.0;
            pParam->BSIM4a1 = 0.0;
        }

        if (pParam->BSIM4prwg < 0.0)
        {   fprintf(fplog, "Warning: Prwg = %g is negative. Set to zero.\n",
                    pParam->BSIM4prwg);
            printf("Warning: Prwg = %g is negative. Set to zero.\n",
                   pParam->BSIM4prwg);
            pParam->BSIM4prwg = 0.0;
        }

        if (pParam->BSIM4rdsw < 0.0)
        {   fprintf(fplog, "Warning: Rdsw = %g is negative. Set to zero.\n",
                    pParam->BSIM4rdsw);
            printf("Warning: Rdsw = %g is negative. Set to zero.\n",
                   pParam->BSIM4rdsw);
            pParam->BSIM4rdsw = 0.0;
            pParam->BSIM4rds0 = 0.0;
        }

        if (pParam->BSIM4rds0 < 0.0)
        {   fprintf(fplog, "Warning: Rds at current temperature = %g is negative. Set to zero.\n",
                    pParam->BSIM4rds0);
            printf("Warning: Rds at current temperature = %g is negative. Set to zero.\n",
                   pParam->BSIM4rds0);
            pParam->BSIM4rds0 = 0.0;
        }

        if (pParam->BSIM4rdswmin < 0.0)
        {   fprintf(fplog, "Warning: Rdswmin at current temperature = %g is negative. Set to zero.\n",
                    pParam->BSIM4rdswmin);
            printf("Warning: Rdswmin at current temperature = %g is negative. Set to zero.\n",
                   pParam->BSIM4rdswmin);
            pParam->BSIM4rdswmin = 0.0;
        }

        if (pParam->BSIM4pscbe2 <= 0.0)
        {   fprintf(fplog, "Warning: Pscbe2 = %g is not positive.\n",
                    pParam->BSIM4pscbe2);
            printf("Warning: Pscbe2 = %g is not positive.\n", pParam->BSIM4pscbe2);
        }

        if (pParam->BSIM4vsattemp < 1.0e3)
        {   fprintf(fplog, "Warning: Vsat at current temperature = %g may be too small.\n", pParam->BSIM4vsattemp);
           printf("Warning: Vsat at current temperature = %g may be too small.\n", pParam->BSIM4vsattemp);
        }

      if((model->BSIM4lambdaGiven) && (pParam->BSIM4lambda > 0.0) ) 
      {
        if (pParam->BSIM4lambda > 1.0e-9)
        {   fprintf(fplog, "Warning: Lambda = %g may be too large.\n", pParam->BSIM4lambda);
           printf("Warning: Lambda = %g may be too large.\n", pParam->BSIM4lambda);
        }
      }

      if((model->BSIM4vtlGiven) && (pParam->BSIM4vtl > 0.0) )
      {  
        if (pParam->BSIM4vtl < 6.0e4)
        {   fprintf(fplog, "Warning: Thermal velocity vtl = %g may be too small.\n", pParam->BSIM4vtl);
           printf("Warning: Thermal velocity vtl = %g may be too small.\n", pParam->BSIM4vtl);
        }

        if (pParam->BSIM4xn < 3.0)
        {   fprintf(fplog, "Warning: back scattering coeff xn = %g is too small.\n", pParam->BSIM4xn);
            printf("Warning: back scattering coeff xn = %g is too small. Reset to 3.0 \n", pParam->BSIM4xn);
            pParam->BSIM4xn = 3.0;
        }

        if (model->BSIM4lc < 0.0)
        {   fprintf(fplog, "Warning: back scattering coeff lc = %g is too small.\n", model->BSIM4lc);
            printf("Warning: back scattering coeff lc = %g is too small. Reset to 0.0\n", model->BSIM4lc);
            pParam->BSIM4lc = 0.0;
        }
      }

        if (pParam->BSIM4pdibl1 < 0.0)
        {   fprintf(fplog, "Warning: Pdibl1 = %g is negative.\n",
                    pParam->BSIM4pdibl1);
            printf("Warning: Pdibl1 = %g is negative.\n", pParam->BSIM4pdibl1);
        }
        if (pParam->BSIM4pdibl2 < 0.0)
        {   fprintf(fplog, "Warning: Pdibl2 = %g is negative.\n",
                    pParam->BSIM4pdibl2);
            printf("Warning: Pdibl2 = %g is negative.\n", pParam->BSIM4pdibl2);
        }

/* Check stress effect parameters */        
        if((here->BSIM4sa > 0.0) && (here->BSIM4sb > 0.0) && 
               ((here->BSIM4nf == 1.0) || ((here->BSIM4nf > 1.0) && (here->BSIM4sd > 0.0))) )
        {   if (model->BSIM4lodk2 <= 0.0)
            {   fprintf(fplog, "Warning: LODK2 = %g is not positive.\n",model->BSIM4lodk2);
                printf("Warning: LODK2 = %g is not positive.\n",model->BSIM4lodk2);
            }
            if (model->BSIM4lodeta0 <= 0.0)
            {   fprintf(fplog, "Warning: LODETA0 = %g is not positive.\n",model->BSIM4lodeta0);
                    printf("Warning: LODETA0 = %g is not positive.\n",model->BSIM4lodeta0);
            }
         }

/* Check gate resistance parameters */        
        if (here->BSIM4rgateMod == 1)
        {   if (model->BSIM4rshg <= 0.0)
                printf("Warning: rshg should be positive for rgateMod = 1.\n");
        }
        else if (here->BSIM4rgateMod == 2)
        {   if (model->BSIM4rshg <= 0.0)
                printf("Warning: rshg <= 0.0 for rgateMod = 2.\n");
            else if (pParam->BSIM4xrcrg1 <= 0.0)
                     printf("Warning: xrcrg1 <= 0.0 for rgateMod = 2.\n");
         }
        if (here->BSIM4rgateMod == 3)
        {   if (model->BSIM4rshg <= 0.0)
                printf("Warning: rshg should be positive for rgateMod = 3.\n");
            else if (pParam->BSIM4xrcrg1 <= 0.0)
                     printf("Warning: xrcrg1 should be positive for rgateMod = 3.\n");
         }
         
         /* Check body resistance parameters */ 
        if (model->BSIM4rbps0 <= 0.0)
        {   fprintf(fplog, "Fatal: RBPS0 = %g is not positive.\n",
                model->BSIM4rbps0 );
            printf("Fatal: RBPS0 = %g is not positive.\n", model->BSIM4rbps0);
            Fatal_Flag = 1;
        }
            if (model->BSIM4rbpd0 <= 0.0)
        {   fprintf(fplog, "Fatal: RBPD0 = %g is not positive.\n",
                model->BSIM4rbpd0 );
            printf("Fatal: RBPD0 = %g is not positive.\n", model->BSIM4rbpd0);
            Fatal_Flag = 1;
        }
            if (model->BSIM4rbpbx0 <= 0.0)
        {   fprintf(fplog, "Fatal: RBPBX0 = %g is not positive.\n",
                model->BSIM4rbpbx0);
            printf("Fatal: RBPBX0 = %g is not positive.\n", model->BSIM4rbpbx0);
            Fatal_Flag = 1;
        }
            if (model->BSIM4rbpby0 <= 0.0)
        {   fprintf(fplog, "Fatal: RBPBY0 = %g is not positive.\n",
                model->BSIM4rbpby0);
            printf("Fatal: RBPBY0 = %g is not positive.\n", model->BSIM4rbpby0);
            Fatal_Flag = 1;
        }
            if (model->BSIM4rbdbx0 <= 0.0)
        {   fprintf(fplog, "Fatal: RBDBX0 = %g is not positive.\n",
                model->BSIM4rbdbx0);
            printf("Fatal: RBDBX0 = %g is not positive.\n", model->BSIM4rbdbx0);
            Fatal_Flag = 1;
        }
            if (model->BSIM4rbdby0 <= 0.0)
        {   fprintf(fplog, "Fatal: RBDBY0 = %g is not positive.\n",
                model->BSIM4rbdby0);
            printf("Fatal: RBDBY0 = %g is not positive.\n", model->BSIM4rbdby0);
            Fatal_Flag = 1;
        }
            if (model->BSIM4rbsbx0 <= 0.0)
        {   fprintf(fplog, "Fatal: RBSBX0 = %g is not positive.\n",
                model->BSIM4rbsbx0);
            printf("Fatal: RBSBX0 = %g is not positive.\n", model->BSIM4rbsbx0);
            Fatal_Flag = 1;
        }
            if (model->BSIM4rbsby0 <= 0.0)
        {   fprintf(fplog, "Fatal: RBSBY0 = %g is not positive.\n",
                model->BSIM4rbsby0);
            printf("Fatal: RBSBY0 = %g is not positive.\n", model->BSIM4rbsby0);
            Fatal_Flag = 1;
        }

/* Check capacitance parameters */
        if (pParam->BSIM4noff < 0.1)
        {   fprintf(fplog, "Warning: Noff = %g is too small.\n",
                    pParam->BSIM4noff);
            printf("Warning: Noff = %g is too small.\n", pParam->BSIM4noff);
        }

        if (pParam->BSIM4voffcv < -0.5)
        {   fprintf(fplog, "Warning: Voffcv = %g is too small.\n",
                    pParam->BSIM4voffcv);
            printf("Warning: Voffcv = %g is too small.\n", pParam->BSIM4voffcv);
        }
        if (pParam->BSIM4moin < 5.0)
        {   fprintf(fplog, "Warning: Moin = %g is too small.\n",
                    pParam->BSIM4moin);
            printf("Warning: Moin = %g is too small.\n", pParam->BSIM4moin);
        }
        if (pParam->BSIM4moin > 25.0)
        {   fprintf(fplog, "Warning: Moin = %g is too large.\n",
                    pParam->BSIM4moin);
            printf("Warning: Moin = %g is too large.\n", pParam->BSIM4moin);
        }
        if(model->BSIM4capMod ==2) {
                if (pParam->BSIM4acde < 0.1)
                {   fprintf(fplog, "Warning:  Acde = %g is too small.\n",
                            pParam->BSIM4acde);
                    printf("Warning: Acde = %g is too small.\n", pParam->BSIM4acde);
                }
                if (pParam->BSIM4acde > 1.6)
                {   fprintf(fplog, "Warning:  Acde = %g is too large.\n",
                            pParam->BSIM4acde);
                    printf("Warning: Acde = %g is too large.\n", pParam->BSIM4acde);
                }
        }

/* Check overlap capacitance parameters */
        if (model->BSIM4cgdo < 0.0)
        {   fprintf(fplog, "Warning: cgdo = %g is negative. Set to zero.\n", model->BSIM4cgdo);
            printf("Warning: cgdo = %g is negative. Set to zero.\n", model->BSIM4cgdo);
            model->BSIM4cgdo = 0.0;
        }      
        if (model->BSIM4cgso < 0.0)
        {   fprintf(fplog, "Warning: cgso = %g is negative. Set to zero.\n", model->BSIM4cgso);
            printf("Warning: cgso = %g is negative. Set to zero.\n", model->BSIM4cgso);
            model->BSIM4cgso = 0.0;
        }      
        if (model->BSIM4cgbo < 0.0)
        {   fprintf(fplog, "Warning: cgbo = %g is negative. Set to zero.\n", model->BSIM4cgbo);
            printf("Warning: cgbo = %g is negative. Set to zero.\n", model->BSIM4cgbo);
            model->BSIM4cgbo = 0.0;
        }      

        /* v4.7 */
        if (model->BSIM4tnoiMod == 1 || model->BSIM4tnoiMod == 2) { 
            if (model->BSIM4tnoia < 0.0) {   
                fprintf(fplog, "Warning: tnoia = %g is negative. Set to zero.\n", model->BSIM4tnoia);
                printf("Warning: tnoia = %g is negative. Set to zero.\n", model->BSIM4tnoia);
                model->BSIM4tnoia = 0.0;
            }
            if (model->BSIM4tnoib < 0.0) {   
                fprintf(fplog, "Warning: tnoib = %g is negative. Set to zero.\n", model->BSIM4tnoib);
                printf("Warning: tnoib = %g is negative. Set to zero.\n", model->BSIM4tnoib);
                model->BSIM4tnoib = 0.0;
            }
            if (model->BSIM4rnoia < 0.0) { 
                fprintf(fplog, "Warning: rnoia = %g is negative. Set to zero.\n", model->BSIM4rnoia);
                printf("Warning: rnoia = %g is negative. Set to zero.\n", model->BSIM4rnoia);
                model->BSIM4rnoia = 0.0;
            }
            if (model->BSIM4rnoib < 0.0) {
                fprintf(fplog, "Warning: rnoib = %g is negative. Set to zero.\n", model->BSIM4rnoib);
                printf("Warning: rnoib = %g is negative. Set to zero.\n", model->BSIM4rnoib);
                model->BSIM4rnoib = 0.0;
            }
        }

        /* v4.7 */
        if (model->BSIM4tnoiMod == 2) { 
            if (model->BSIM4tnoic < 0.0) {
                fprintf(fplog, "Warning: tnoic = %g is negative. Set to zero.\n", model->BSIM4tnoic);
                printf("Warning: tnoic = %g is negative. Set to zero.\n", model->BSIM4tnoic);
                model->BSIM4tnoic = 0.0;
            }
            if (model->BSIM4rnoic < 0.0) {
                fprintf(fplog, "Warning: rnoic = %g is negative. Set to zero.\n", model->BSIM4rnoic);
                printf("Warning: rnoic = %g is negative. Set to zero.\n", model->BSIM4rnoic);
                model->BSIM4rnoic = 0.0;
            }
        }

        /* Limits of Njs and Njd modified in BSIM4.7 */
        if (model->BSIM4SjctEmissionCoeff < 0.1) {
            fprintf(fplog, "Warning: Njs = %g is less than 0.1. Setting Njs to 0.1.\n", model->BSIM4SjctEmissionCoeff);
            printf("Warning: Njs = %g is less than 0.1. Setting Njs to 0.1.\n", model->BSIM4SjctEmissionCoeff);
            model->BSIM4SjctEmissionCoeff = 0.1;
        } 
        else if (model->BSIM4SjctEmissionCoeff < 0.7) { 
            fprintf(fplog, "Warning: Njs = %g is less than 0.7.\n", model->BSIM4SjctEmissionCoeff);
            printf("Warning: Njs = %g is less than 0.7.\n", model->BSIM4SjctEmissionCoeff);
        }
        if (model->BSIM4DjctEmissionCoeff < 0.1) { 
            fprintf(fplog, "Warning: Njd = %g is less than 0.1. Setting Njd to 0.1.\n", model->BSIM4DjctEmissionCoeff);
            printf("Warning: Njd = %g is less than 0.1. Setting Njd to 0.1.\n", model->BSIM4DjctEmissionCoeff);
            model->BSIM4DjctEmissionCoeff = 0.1;
        }
        else if (model->BSIM4DjctEmissionCoeff < 0.7) {
            fprintf(fplog, "Warning: Njd = %g is less than 0.7.\n", model->BSIM4DjctEmissionCoeff);
            printf("Warning: Njd = %g is less than 0.7.\n", model->BSIM4DjctEmissionCoeff);
        }

        if (model->BSIM4njtsstemp < 0.0)
        {   fprintf(fplog, "Warning: Njts = %g is negative at temperature = %g.\n",
                    model->BSIM4njtsstemp, ckt->CKTtemp);
            printf("Warning: Njts = %g is negative at temperature = %g.\n",
                   model->BSIM4njtsstemp, ckt->CKTtemp);
        }
        if (model->BSIM4njtsswstemp < 0.0)
        {   fprintf(fplog, "Warning: Njtssw = %g is negative at temperature = %g.\n",
                    model->BSIM4njtsswstemp, ckt->CKTtemp);
            printf("Warning: Njtssw = %g is negative at temperature = %g.\n",
                    model->BSIM4njtsswstemp, ckt->CKTtemp);
        }
        if (model->BSIM4njtsswgstemp < 0.0)
        {   fprintf(fplog, "Warning: Njtsswg = %g is negative at temperature = %g.\n",
                    model->BSIM4njtsswgstemp, ckt->CKTtemp);
            printf("Warning: Njtsswg = %g is negative at temperature = %g.\n",
                    model->BSIM4njtsswgstemp, ckt->CKTtemp);
        }

        if (model->BSIM4njtsdGiven && model->BSIM4njtsdtemp < 0.0)
        {   fprintf(fplog, "Warning: Njtsd = %g is negative at temperature = %g.\n",
                    model->BSIM4njtsdtemp, ckt->CKTtemp);
            printf("Warning: Njtsd = %g is negative at temperature = %g.\n",
                    model->BSIM4njtsdtemp, ckt->CKTtemp);
        }
        if (model->BSIM4njtsswdGiven && model->BSIM4njtsswdtemp < 0.0)
        {   fprintf(fplog, "Warning: Njtsswd = %g is negative at temperature = %g.\n",
                    model->BSIM4njtsswdtemp, ckt->CKTtemp);
            printf("Warning: Njtsswd = %g is negative at temperature = %g.\n",
                    model->BSIM4njtsswdtemp, ckt->CKTtemp);
        }
        if (model->BSIM4njtsswgdGiven && model->BSIM4njtsswgdtemp < 0.0)
        {   fprintf(fplog, "Warning: Njtsswgd = %g is negative at temperature = %g.\n",
                    model->BSIM4njtsswgdtemp, ckt->CKTtemp);
            printf("Warning: Njtsswgd = %g is negative at temperature = %g.\n",
                    model->BSIM4njtsswgdtemp, ckt->CKTtemp);
        }

       if (model->BSIM4ntnoi < 0.0)
        {   fprintf(fplog, "Warning: ntnoi = %g is negative. Set to zero.\n", model->BSIM4ntnoi);
            printf("Warning: ntnoi = %g is negative. Set to zero.\n", model->BSIM4ntnoi);
            model->BSIM4ntnoi = 0.0;
        }

        /* diode model */
       if (model->BSIM4SbulkJctBotGradingCoeff >= 0.99)
        {   fprintf(fplog, "Warning: MJS = %g is too big. Set to 0.99.\n", model->BSIM4SbulkJctBotGradingCoeff);
            printf("Warning: MJS = %g is too big. Set to 0.99.\n", model->BSIM4SbulkJctBotGradingCoeff);
            model->BSIM4SbulkJctBotGradingCoeff = 0.99;
        }
       if (model->BSIM4SbulkJctSideGradingCoeff >= 0.99)
        {   fprintf(fplog, "Warning: MJSWS = %g is too big. Set to 0.99.\n", model->BSIM4SbulkJctSideGradingCoeff);
            printf("Warning: MJSWS = %g is too big. Set to 0.99.\n", model->BSIM4SbulkJctSideGradingCoeff);
            model->BSIM4SbulkJctSideGradingCoeff = 0.99;
        }
       if (model->BSIM4SbulkJctGateSideGradingCoeff >= 0.99)
        {   fprintf(fplog, "Warning: MJSWGS = %g is too big. Set to 0.99.\n", model->BSIM4SbulkJctGateSideGradingCoeff);
            printf("Warning: MJSWGS = %g is too big. Set to 0.99.\n", model->BSIM4SbulkJctGateSideGradingCoeff);
            model->BSIM4SbulkJctGateSideGradingCoeff = 0.99;
        }

       if (model->BSIM4DbulkJctBotGradingCoeff >= 0.99)
        {   fprintf(fplog, "Warning: MJD = %g is too big. Set to 0.99.\n", model->BSIM4DbulkJctBotGradingCoeff);
            printf("Warning: MJD = %g is too big. Set to 0.99.\n", model->BSIM4DbulkJctBotGradingCoeff);
            model->BSIM4DbulkJctBotGradingCoeff = 0.99;
        }
       if (model->BSIM4DbulkJctSideGradingCoeff >= 0.99)
        {   fprintf(fplog, "Warning: MJSWD = %g is too big. Set to 0.99.\n", model->BSIM4DbulkJctSideGradingCoeff);
            printf("Warning: MJSWD = %g is too big. Set to 0.99.\n", model->BSIM4DbulkJctSideGradingCoeff);
            model->BSIM4DbulkJctSideGradingCoeff = 0.99;
        }
       if (model->BSIM4DbulkJctGateSideGradingCoeff >= 0.99)
        {   fprintf(fplog, "Warning: MJSWGD = %g is too big. Set to 0.99.\n", model->BSIM4DbulkJctGateSideGradingCoeff);
            printf("Warning: MJSWGD = %g is too big. Set to 0.99.\n", model->BSIM4DbulkJctGateSideGradingCoeff);
            model->BSIM4DbulkJctGateSideGradingCoeff = 0.99;
        }
        if (model->BSIM4wpemod == 1)
        {
                if (model->BSIM4scref <= 0.0)
                {   fprintf(fplog, "Warning: SCREF = %g is not positive. Set to 1e-6.\n", model->BSIM4scref);
                    printf("Warning: SCREF = %g is not positive. Set to 1e-6.\n", model->BSIM4scref);
                    model->BSIM4scref = 1e-6;
                }
                /*Move these checks to temp.c for sceff calculation*/
                        /*
                        if (here->BSIM4sca < 0.0)
                {   fprintf(fplog, "Warning: SCA = %g is negative. Set to 0.0.\n", here->BSIM4sca);
                    printf("Warning: SCA = %g is negative. Set to 0.0.\n", here->BSIM4sca);
                    here->BSIM4sca = 0.0;
                }
                if (here->BSIM4scb < 0.0)
                {   fprintf(fplog, "Warning: SCB = %g is negative. Set to 0.0.\n", here->BSIM4scb);
                    printf("Warning: SCB = %g is negative. Set to 0.0.\n", here->BSIM4scb);
                    here->BSIM4scb = 0.0;
                }
                if (here->BSIM4scc < 0.0)
                {   fprintf(fplog, "Warning: SCC = %g is negative. Set to 0.0.\n", here->BSIM4scc);
                    printf("Warning: SCC = %g is negative. Set to 0.0.\n", here->BSIM4scc);
                    here->BSIM4scc = 0.0;
                }
                if (here->BSIM4sc < 0.0)
                {   fprintf(fplog, "Warning: SC = %g is negative. Set to 0.0.\n", here->BSIM4sc);
                    printf("Warning: SC = %g is negative. Set to 0.0.\n", here->BSIM4sc);
                    here->BSIM4sc = 0.0;
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

