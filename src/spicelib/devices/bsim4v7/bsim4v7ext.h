/**********
Copyright 2004 Regents of the University of California.  All rights reserved.
Author: 2000 Weidong Liu
Author: 2001- Xuemei Xi
File: bsim4v7ext.h
**********/

extern int BSIM4v7acLoad(GENmodel *,CKTcircuit*);
extern int BSIM4v7ask(CKTcircuit *,GENinstance*,int,IFvalue*,IFvalue*);
extern int BSIM4v7convTest(GENmodel *,CKTcircuit*);
extern int BSIM4v7delete(GENmodel*,IFuid,GENinstance**);
extern void BSIM4v7destroy(GENmodel**);
extern int BSIM4v7getic(GENmodel*,CKTcircuit*);
extern int BSIM4v7load(GENmodel*,CKTcircuit*);
extern int BSIM4v7mAsk(CKTcircuit*,GENmodel *,int, IFvalue*);
extern int BSIM4v7mDelete(GENmodel**,IFuid,GENmodel*);
extern int BSIM4v7mParam(int,IFvalue*,GENmodel*);
extern void BSIM4v7mosCap(CKTcircuit*, double, double, double, double,
        double, double, double, double, double, double, double,
        double, double, double, double, double, double, double*,
        double*, double*, double*, double*, double*, double*, double*,
        double*, double*, double*, double*, double*, double*, double*, 
        double*);
extern int BSIM4v7param(int,IFvalue*,GENinstance*,IFvalue*);
extern int BSIM4v7pzLoad(GENmodel*,CKTcircuit*,SPcomplex*);
extern int BSIM4v7setup(SMPmatrix*,GENmodel*,CKTcircuit*,int*);
extern int BSIM4v7temp(GENmodel*,CKTcircuit*);
extern int BSIM4v7trunc(GENmodel*,CKTcircuit*,double*);
extern int BSIM4v7noise(int,int,GENmodel*,CKTcircuit*,Ndata*,double*);
extern int BSIM4v7unsetup(GENmodel*,CKTcircuit*);
extern int BSIM4v7soaCheck(CKTcircuit *, GENmodel *);
