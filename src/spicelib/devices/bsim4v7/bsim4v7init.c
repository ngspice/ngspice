#include "ngspice/config.h"

#include "ngspice/devdefs.h"

#include "bsim4v7itf.h."
#include "bsim4v7ext.h."
#include "bsim4v7init.h."


SPICEdev BSIM4v7info = {
    {
        "BSIM4v7",
        "Berkeley Short Channel IGFET Model-4",

        &BSIM4v7nSize,
        &BSIM4v7nSize,
        BSIM4v7names,

        &BSIM4v7pTSize,
        BSIM4v7pTable,

        &BSIM4v7mPTSize,
        BSIM4v7mPTable,

#ifdef XSPICE
/*----  Fixed by SDB 5.2.2003 to enable XSPICE/tclspice integration  -----*/
        NULL,  /* This is a SPICE device, it has no MIF info data */

        0,     /* This is a SPICE device, it has no MIF info data */
        NULL,  /* This is a SPICE device, it has no MIF info data */

        0,     /* This is a SPICE device, it has no MIF info data */
        NULL,  /* This is a SPICE device, it has no MIF info data */

        0,     /* This is a SPICE device, it has no MIF info data */
        NULL,  /* This is a SPICE device, it has no MIF info data */
/*---------------------------  End of SDB fix   -------------------------*/
#endif

        DEV_DEFAULT
    },

    BSIM4v7param,    /* DEVparam       */
    BSIM4v7mParam,   /* DEVmodParam    */
    BSIM4v7load,     /* DEVload        */
    BSIM4v7setup,    /* DEVsetup       */
    BSIM4v7unsetup,  /* DEVunsetup     */
    BSIM4v7setup,    /* DEVpzSetup     */
    BSIM4v7temp,     /* DEVtemperature */
    BSIM4v7trunc,    /* DEVtrunc       */
    NULL,          /* DEVfindBranch  */
    BSIM4v7acLoad,   /* DEVacLoad      */
    NULL,          /* DEVaccept      */
    BSIM4v7destroy,  /* DEVdestroy     */
    BSIM4v7mDelete,  /* DEVmodDelete   */
    BSIM4v7delete,   /* DEVdelete      */
    BSIM4v7getic,    /* DEVsetic       */
    BSIM4v7ask,      /* DEVask         */
    BSIM4v7mAsk,     /* DEVmodAsk      */
    BSIM4v7pzLoad,   /* DEVpzLoad      */
    BSIM4v7convTest, /* DEVconvTest    */
    NULL,          /* DEVsenSetup    */
    NULL,          /* DEVsenLoad     */
    NULL,          /* DEVsenUpdate   */
    NULL,          /* DEVsenAcLoad   */
    NULL,          /* DEVsenPrint    */
    NULL,          /* DEVsenTrunc    */
    NULL,          /* DEVdisto       */
    BSIM4v7noise,    /* DEVnoise       */
    BSIM4v7soaCheck, /* DEVsoaCheck    */
#ifdef CIDER
    NULL,          /* DEVdump        */
    NULL,          /* DEVacct        */
#endif
    &BSIM4v7iSize,   /* DEVinstSize    */
    &BSIM4v7mSize    /* DEVmodSize     */
};


SPICEdev *
get_bsim4v7_info(void)
{
    return &BSIM4v7info;
}
