/* NG-SPICE -- An electrical circuit simulator
 *
 * Copyright (c) 1990 University of California
 * Copyright (c) 2000 Arno W. Peters
 *
 * Permission to use, copy, modify, and distribute this software and
 * its documentation without fee, and without a written agreement is
 * hereby granted, provided that the above copyright notice, this
 * paragraph and the following three paragraphs appear in all copies.
 *
 * This software program and documentation are copyrighted by their
 * authors. The software program and documentation are supplied "as
 * is", without any accompanying services from the authors. The
 * authors do not warrant that the operation of the program will be
 * uninterrupted or error-free. The end-user understands that the
 * program was developed for research purposes and is advised not to
 * rely exclusively on the program for any reason.
 * 
 * IN NO EVENT SHALL THE AUTHORS BE LIABLE TO ANY PARTY FOR DIRECT,
 * INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
 * LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
 * DOCUMENTATION, EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE. THE AUTHORS SPECIFICALLY DISCLAIMS ANY
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
 * SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE AUTHORS
 * HAVE NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
 * ENHANCEMENTS, OR MODIFICATIONS. */

#include "ngspice/ngspice.h"
#include "assert.h"

#include "ngspice/devdefs.h"
#include "ngspice/ifsim.h"

#include "dev.h"
#include "ngspice/memory.h" /* to alloc, realloc devices*/


#ifdef XSPICE
/*saj headers for xspice*/
#include <string.h> /* for strcpy, strcat*/
#if (!defined HAS_WINGUI) && (!defined __MINGW32__) && (!defined _MSC_VER)
#include <dlfcn.h> /* to load libraries*/
typedef void *  funptr_t;
#else /* ifdef HAS_WINGUI */
#undef BOOLEAN
#include <windows.h>
typedef FARPROC funptr_t;
void *dlopen (const char *, int);
funptr_t dlsym (void *, const char *);
int dlclose (void *);
char *dlerror (void);
#define RTLD_LAZY	1	/* lazy function call binding */
#define RTLD_NOW	2	/* immediate function call binding */
#define RTLD_GLOBAL	4	/* symbols in this dlopen'ed obj are visible to other dlopen'ed objs */
static char errstr[128];
#endif /* ifndef HAS_WINGUI */
#include "ngspice/dllitf.h" /* the coreInfo Structure*/
#include "ngspice/evtudn.h" /*Use defined nodes */

Evt_Udn_Info_t  **g_evt_udn_info = NULL;
int g_evt_num_udn_types = 0;

/*The digital node type */
extern Evt_Udn_Info_t idn_digital_info;
int add_device(int n, SPICEdev **devs, int flag);
int add_udn(int,Evt_Udn_Info_t **);
/*saj*/
#endif


#include "asrc/asrcitf.h"
#include "bjt/bjtitf.h"
#include "bsim1/bsim1itf.h"
#include "bsim2/bsim2itf.h"
#include "bsim3/bsim3itf.h"
#include "bsim3v0/bsim3v0itf.h"
#include "bsim3v1/bsim3v1itf.h"
#include "bsim3v32/bsim3v32itf.h"
#include "bsim4/bsim4itf.h"
#include "bsim4v5/bsim4v5itf.h"
#include "bsim4v6/bsim4v6itf.h"
#include "bsim4v7/bsim4v7itf.h"
#include "bsim3soi_pd/b3soipditf.h"
#include "bsim3soi_fd/b3soifditf.h"
#include "bsim3soi_dd/b3soidditf.h"
#include "bsimsoi/b4soiitf.h"
#include "cap/capitf.h"
#include "cccs/cccsitf.h"
#include "ccvs/ccvsitf.h"
#include "csw/cswitf.h"
#include "dio/dioitf.h"
#include "hfet1/hfetitf.h"
#include "hfet2/hfet2itf.h"
#include "hisim2/hsm2itf.h"
#include "hisimhv1/hsmhvitf.h"
#include "ind/inditf.h"
#include "isrc/isrcitf.h"
#include "jfet/jfetitf.h"
#include "jfet2/jfet2itf.h"
#include "ltra/ltraitf.h"
#include "mes/mesitf.h"
#include "mesa/mesaitf.h"
#include "mos1/mos1itf.h"
#include "mos2/mos2itf.h"
#include "mos3/mos3itf.h"
#include "mos6/mos6itf.h"
#include "mos9/mos9itf.h"
#include "cpl/cplitf.h"
#include "res/resitf.h"
#include "soi3/soi3itf.h"
#include "sw/switf.h"
#include "tra/traitf.h"
#include "txl/txlitf.h"
#include "urc/urcitf.h"
#include "vbic/vbicitf.h"
#include "vccs/vccsitf.h"
#include "vcvs/vcvsitf.h"
#include "vsrc/vsrcitf.h"
#ifdef ADMS
#include "adms/hicum0/hicum0itf.h"
#include "adms/hicum2/hicum2itf.h"
#include "adms/mextram/bjt504titf.h"
#include "adms/ekv/ekvitf.h"
#include "adms/psp102/psp102itf.h"
#endif
#ifdef CIDER
/* Numerical devices (Cider integration) */
#include "nbjt/nbjtitf.h"
#include "nbjt2/nbjt2itf.h"
#include "numd/numditf.h"
#include "numd2/numd2itf.h"
#include "numos/numositf.h"
#endif
#ifdef NDEV
#include "ndev/ndevitf.h"
#endif

/*saj in xspice the DEVices size can be varied so DEVNUM is an int*/
#if defined XSPICE
   static int DEVNUM = 63;
#else
   #define DEVNUM 63
#endif

/*Make this dynamic for later attempt to make all devices dynamic*/
SPICEdev **DEVices=NULL;

/*Flag to indicate that device type it is,
 *0 = normal spice device
 *1 = xspice device
 */
#ifdef XSPICE
int *DEVicesfl=NULL;
int DEVflag(int type){
  if(type < DEVNUM && type >= 0)
    return DEVicesfl[type];
  else
    return -1;
}
#endif



void
spice_destroy_devices(void)
{
#ifdef XSPICE
    tfree(g_evt_udn_info);
    tfree(DEVicesfl);
#endif
    tfree(DEVices);
}


void
spice_init_devices(void)
{
#ifdef XSPICE
    /*Initilise the structs and add digital node type */
    g_evt_udn_info = TMALLOC(Evt_Udn_Info_t  *, 1);
    g_evt_num_udn_types = 1;
    g_evt_udn_info[0] =  &idn_digital_info;
    
    DEVicesfl = TMALLOC(int, DEVNUM);
    /* tmalloc should automatically zero the array! */
#endif

    DEVices = TMALLOC(SPICEdev *, DEVNUM);
    /* URC device MUST precede both resistors and capacitors */
    DEVices[ 0] = get_urc_info();
    DEVices[ 1] = get_asrc_info();
    DEVices[ 2] = get_bjt_info();
    DEVices[ 3] = get_bsim1_info();
    DEVices[ 4] = get_bsim2_info();
    DEVices[ 5] = get_bsim3_info();
    DEVices[ 6] = get_bsim3v0_info();
    DEVices[ 7] = get_bsim3v1_info();
    DEVices[ 8] = get_bsim3v32_info();
    DEVices[ 9] = get_b4soi_info();
    DEVices[10] = get_bsim4_info();
    DEVices[11] = get_bsim4v5_info();
    DEVices[12] = get_bsim4v6_info();
    DEVices[13] = get_bsim4v7_info();
    DEVices[14] = get_b3soipd_info();
    DEVices[15] = get_b3soifd_info();
    DEVices[16] = get_b3soidd_info();
    DEVices[17] = get_cap_info();
    DEVices[18] = get_cccs_info();
    DEVices[19] = get_ccvs_info();
    DEVices[20] = get_cpl_info();
    DEVices[21] = get_csw_info();
    DEVices[22] = get_dio_info();
    DEVices[23] = get_hfeta_info();
    DEVices[24] = get_hfet2_info();
    DEVices[25] = get_hsm2_info();  
    DEVices[26] = get_hsmhv_info();  
    DEVices[27] = get_ind_info();
    DEVices[28] = get_mut_info();
    DEVices[29] = get_isrc_info();
    DEVices[30] = get_jfet_info();
    DEVices[31] = get_jfet2_info();
    DEVices[32] = get_ltra_info();
    DEVices[33] = get_mes_info();
    DEVices[34] = get_mesa_info();
    DEVices[35] = get_mos1_info();
    DEVices[36] = get_mos2_info();
    DEVices[37] = get_mos3_info();
    DEVices[38] = get_mos6_info();
    DEVices[39] = get_mos9_info();
    DEVices[40] = get_res_info();
    DEVices[41] = get_soi3_info();
    DEVices[42] = get_sw_info();
    DEVices[43] = get_tra_info();
    DEVices[44] = get_txl_info();
    DEVices[45] = get_vbic_info();
    DEVices[46] = get_vccs_info();
    DEVices[47] = get_vcvs_info();
    DEVices[48] = get_vsrc_info();
#ifdef CIDER
    DEVices[49] = get_nbjt_info();
    DEVices[50] = get_nbjt2_info();
    DEVices[51] = get_numd_info();
    DEVices[52] = get_numd2_info();
    DEVices[53] = get_numos_info();   
#else
    DEVices[49] = NULL;
    DEVices[50] = NULL;
    DEVices[51] = NULL;
    DEVices[52] = NULL;
    DEVices[53] = NULL;  
#endif
          
#ifdef ADMS
    DEVices[54] = (SPICEdev*)get_hicum0_info();
    DEVices[55] = (SPICEdev*)get_hicum2_info();
    DEVices[56] = (SPICEdev*)get_bjt504t_info();
    DEVices[57] = (SPICEdev*)get_ekv_info();
    DEVices[58] = (SPICEdev*)get_psp102_info();
#else
    DEVices[54] = NULL;
    DEVices[55] = NULL;
    DEVices[56] = NULL;
    DEVices[57] = NULL;
    DEVices[58] = NULL;
#endif
   
#ifdef NDEV    /* NDEV */
   DEVices[59] = get_ndev_info();
#else
   DEVices[59] = NULL;
#endif
   DEVices[60] = NULL;
   DEVices[61] = NULL;


   return;
}

int num_devices(void)
{
    return DEVNUM;
}

IFdevice ** devices_ptr(void)
{
    return (IFdevice **) DEVices;
}

SPICEdev ** devices(void)
{
    return DEVices;
}


#ifdef DEVLIB
/*not yet usable*/

#ifdef ADMS
#define DEVICES_USED {"asrc", "bjt", "vbic", "bsim1", "bsim2", "bsim3", "bsim3v32", "bsim3v2", "bsim3v1", "bsim4", "bsim4v5", "bsim4v6", "bsim4v7", \
                      "bsim4soi", "bsim3soipd", "bsim3soifd", "bsim3soidd", "hisim2", "hisimhv1", \
                      "cap", "cccs", "ccvs", "csw", "dio", "hfet", "hfet2", "ind", "isrc", "jfet", "ltra", "mes", "mesa" ,"mos1", "mos2", "mos3", \
                      "mos6", "mos9", "res", "soi3", "sw", "tra", "urc", "vccs", "vcvs", "vsrc", "hicum0", "hicum2", "bjt504t", "ekv", "psp102"}
#else
#define DEVICES_USED {"asrc", "bjt", "vbic", "bsim1", "bsim2", "bsim3", "bsim3v32", "bsim3v2", "bsim3v1", "bsim4", "bsim4v5", "bsim4v6", "bsim4v7", \
                      "bsim4soi", "bsim3soipd", "bsim3soifd", "bsim3soidd", "hisim2", "hisimhv1", \
                      "cap", "cccs", "ccvs", "csw", "dio", "hfet", "hfet2", "ind", "isrc", "jfet", "ltra", "mes", "mesa" ,"mos1", "mos2", "mos3", \
                      "mos6", "mos9", "res", "soi3", "sw", "tra", "urc", "vccs", "vcvs", "vsrc"}
#endif
int load_dev(char *name) {
  char *msg;
  char libname[50];
  void *lib;
  funptr_t fetch;
  SPICEdev *device;

  strcpy(libname, "lib");
  strcat(libname,name);
  strcat(libname,".so");

  lib = dlopen(libname,RTLD_NOW);
  if(!lib){
    msg = dlerror();
    printf("%s\n", msg);
    return 1;
  }

  strcpy(libname, "get_");
  strcat(libname,name);
  strcat(libname,"_info");
  fetch = dlsym(lib,libname);

  if(!fetch){
    msg = dlerror();
    printf("%s\n", msg);
    return 1;
  }
  device = ((SPICEdev * (*)(void)) fetch) ();
  add_device(1,&device,0);
  return 0;
}

void load_alldevs(void){
  char *devs[] = DEVICES_USED;
  int num = NUMELEMS(devs);
  int i;
  for(i=0; i< num;i++)
    load_dev(devs[i]);
  return;
}
#endif

/*--------------------   XSPICE additions below  ----------------------*/
#ifdef XSPICE
#include "ngspice/mif.h"
#include "ngspice/cm.h"
#include "ngspice/cpextern.h"
#include "ngspice/fteext.h"  /* for ft_sim */
#include "ngspice/cktdefs.h" /* for DEVmaxnum */

static void relink(void) {

/*
 * This replacement done by SDB on 6.11.2003
 *
 * ft_sim->numDevices = num_devices();
 * DEVmaxnum = num_devices();
 */
  ft_sim->numDevices = DEVNUM;
  DEVmaxnum = DEVNUM;

  ft_sim->devices = devices_ptr();
  return;
}

int add_device(int n, SPICEdev **devs, int flag){
  int i;
  DEVices = TREALLOC(SPICEdev *, DEVices, DEVNUM + n);
  DEVicesfl = TREALLOC(int, DEVicesfl, DEVNUM + n);
  for(i = 0; i < n;i++){
#ifdef TRACE
      printf("Added device: %s\n",devs[i]->DEVpublic.name);
#endif
    DEVices[DEVNUM+i] = devs[i];

    /* added by SDB on 6.20.2003 */
    DEVices[DEVNUM+i]->DEVinstSize = &MIFiSize;

    DEVicesfl[DEVNUM+i] = flag;
  }
  DEVNUM += n;
  relink();
  return 0;
}

int add_udn(int n,Evt_Udn_Info_t **udns){
  int i;
  g_evt_udn_info = TREALLOC(Evt_Udn_Info_t  *, g_evt_udn_info, g_evt_num_udn_types + n);
  for(i = 0; i < n;i++){
#ifdef TRACE
      printf("Added udn: %s\n",udns[i]->name);
#endif
    g_evt_udn_info[g_evt_num_udn_types+i] = udns[i];
  }
  g_evt_num_udn_types += n;
  return 0;
}

extern struct coreInfo_t  coreInfo;

int load_opus(char *name){
  void *lib;
  const char *msg;
  int *num=NULL;
  struct coreInfo_t **core;
  SPICEdev **devs;
  Evt_Udn_Info_t  **udns;
  funptr_t fetch;

  lib = dlopen(name,RTLD_NOW);
  if(!lib){
    msg = dlerror();
    printf("%s\n", msg);
    return 1;
  }
  
  fetch = dlsym(lib,"CMdevNum");
  if(fetch){
    num = ((int * (*)(void)) fetch) ();
#ifdef TRACE
    printf("Got %u devices.\n",*num);
#endif
  }else{
    msg = dlerror();
    printf("%s\n", msg);
    return 1;
  }

  fetch = dlsym(lib,"CMdevs");
  if(fetch){
    devs = ((SPICEdev ** (*)(void)) fetch) ();
  }else{
    msg = dlerror();
    printf("%s\n", msg);
    return 1;
  }

  fetch = dlsym(lib,"CMgetCoreItfPtr");
  if(fetch){
    core = ((struct coreInfo_t ** (*)(void)) fetch) ();
    *core = &coreInfo;
  }else{
    msg = dlerror();
    printf("%s\n", msg);
    return 1;
  }
  add_device(*num,devs,1);

  fetch = dlsym(lib,"CMudnNum");
  if(fetch){
    num = ((int * (*)(void)) fetch) ();
#ifdef TRACE
    printf("Got %u udns.\n",*num);
#endif
  }else{
    msg = dlerror();
    printf("%s\n", msg);
    return 1;
  }

  fetch = dlsym(lib,"CMudns");
  if(fetch){
    udns = ((Evt_Udn_Info_t  ** (*)(void)) fetch) ();
  }else{
    msg = dlerror();
    printf("%s\n", msg);
    return 1;
  }

  add_udn(*num,udns);

  return 0;
}

#if defined(__MINGW32__) || defined(HAS_WINGUI) || defined(_MSC_VER)

void *dlopen(const char *name,int type)
{
	NG_IGNORE(type);
	return LoadLibrary(name);
}

funptr_t dlsym(void *hDll, const char *funcname)
{
	return GetProcAddress(hDll, funcname);
}

char *dlerror(void)
{
	LPVOID lpMsgBuf;

	FormatMessage(
		FORMAT_MESSAGE_ALLOCATE_BUFFER |
		FORMAT_MESSAGE_FROM_SYSTEM |
		FORMAT_MESSAGE_IGNORE_INSERTS,
		NULL,
		GetLastError(),
		MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		(LPTSTR) &lpMsgBuf,
		0,
		NULL
	);
	strcpy(errstr,lpMsgBuf);
	LocalFree(lpMsgBuf);
	return errstr;
}
#endif

#endif
/*--------------------   end of XSPICE additions  ----------------------*/
