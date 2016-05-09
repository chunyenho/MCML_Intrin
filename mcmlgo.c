/***********************************************************
 *  Copyright Univ. of Texas M.D. Anderson Cancer Center
 *  1992.
 *
 *	Launch, move, and record photon weight.
 ****/

#include "mcml.h"
#include <mkl.h>
#define STANDARDTEST 0
/* testing program using fixed rnd seed. */

#define PARTIALREFLECTION 0
/* 1=split photon, 0=statistical reflection. */

#define COSZERO (1.0-1.0E-12)
/* cosine of about 1e-6 rad. */

#define COS90D  1.0E-6
/* cosine of about 1.57 - 1e-6 rad. */

/* This algorithm is mentioned in the ISO C standard, here extended
 *    for 32 bits.  */
int rand_r (unsigned int *seed)
{
    unsigned int next = *seed;
    int result;

    next *= 1103515245;
    next += 12345;
    result = (unsigned int) (next >>16) & 0x7ff;

    next *= 1103515245;
    next += 12345;
    result <<= 10;
    result ^= (unsigned int) (next >>16) & 0x3ff;

    next *= 1103515245;
    next += 12345;
    result <<= 10;
    result ^= (unsigned int) (next >>16) & 0x3ff;

    *seed = next;

    return result;
}

/* Vectorized random function for MIC. */
double mkl_rand (VSLStreamStatePtr stream, double * result, short * count)
{
	if ((*count) == 0)
	{
		vdRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, 1024, result, 0.0, 1.0 );	
		return result[0];
	}
	return result[(*count)];
}

/***********************************************************
 *	Compute the specular reflection.
 *
 *	If the first layer is a turbid medium, use the Fresnel
 *	reflection from the boundary of the first layer as the
 *	specular reflectance.
 *
 *	If the first layer is glass, multiple reflections in
 *	the first layer is considered to get the specular
 *	reflectance.
 *
 *	The subroutine assumes the Layerspecs array is correctly
 *	initialized.
 ****/
double Rspecular(LayerStruct * Layerspecs_Ptr)
{
    double r1, r2;
    /* direct reflections from the 1st and 2nd layers. */
    double temp;

    temp =(Layerspecs_Ptr[0].n - Layerspecs_Ptr[1].n)
          /(Layerspecs_Ptr[0].n + Layerspecs_Ptr[1].n);
    r1 = temp*temp;

    if((Layerspecs_Ptr[1].mua == 0.0)
            && (Layerspecs_Ptr[1].mus == 0.0))  { /* glass layer. */
        temp = (Layerspecs_Ptr[1].n - Layerspecs_Ptr[2].n)
               /(Layerspecs_Ptr[1].n + Layerspecs_Ptr[2].n);
        r2 = temp*temp;
        r1 = r1 + (1-r1)*(1-r1)*r2/(1-r1*r2);
    }

    return (r1);
}

/***********************************************************
 *	Initialize a photon packet.
 ****/
void Launch8Photon(double Rspecular,
                  LayerStruct  * Layerspecs_Ptr,
                  Photon8Struct * Photon_Ptr)
{
    int i;
#pragma simd 
    for(i = 0; i<8 ; i++)
    {
        Photon_Ptr->w[i]	 	= 1.0 - Rspecular;
        Photon_Ptr->dead[i] 	= 0;
        Photon_Ptr->layer[i] = 1;
        Photon_Ptr->s[i]	= 0;
        Photon_Ptr->sleft[i]= 0;

        Photon_Ptr->x[i] 	= 0.0;
        Photon_Ptr->y[i]	 	= 0.0;
        Photon_Ptr->z[i]	 	= 0.0;
        Photon_Ptr->ux[i]	= 0.0;
        Photon_Ptr->uy[i]	= 0.0;
        Photon_Ptr->uz[i]	= 1.0;

        if((Layerspecs_Ptr[1].mua == 0.0)
                && (Layerspecs_Ptr[1].mus == 0.0))  { /* glass layer. */
            Photon_Ptr->layer[i] 	= 2;
            Photon_Ptr->z[i]	= Layerspecs_Ptr[2].z0;
        }
    }
}

/***********************************************************
 *	Choose (sample) a new theta angle for photon propagation
 *	according to the anisotropy.
 *
 *	If anisotropy g is 0, then
 *		cos(theta) = 2*rand-1.
 *	otherwise
 *		sample according to the Henyey-Greenstein function.
 *
 *	Returns the cosine of the polar deflection angle theta.
 ****/
double SpinTheta(double g, VSLStreamStatePtr stream ,double  *  result,
                 short   *  count)
{
    double cost;

    if(g == 0.0)
    {
        cost = 2*mkl_rand(stream, result, count) -1;
        (*count) ++;
	(*count) %= 1024;
    }
    else {
        double temp = (1-g*g)/(1-g+2*g*mkl_rand(stream, result, count));
	(*count) ++;
        (*count) %= 1024;
        cost = (1+g*g - temp*temp)/(2*g);
        if(cost < -1) cost = -1;
        else if(cost > 1) cost = 1;
    }
    return(cost);
}


/***********************************************************
 *	Choose a new direction for photon propagation by
 *	sampling the polar deflection angle theta and the
 *	azimuthal angle psi.
 *
 *	Note:
 *  	theta: 0 - pi so sin(theta) is always positive
 *  	feel free to use sqrt() for cos(theta).
 *
 *  	psi:   0 - 2pi
 *  	for 0-pi  sin(psi) is +
 *  	for pi-2pi sin(psi) is -
 ****/
void Spin(double g,
          Photon8Struct * Photon_Ptr,
	  VSLStreamStatePtr  stream, double * result, short * count, int vec_num)
{
    double cost, sint;	/* cosine and sine of the */
    /* polar deflection angle theta. */
    double cosp, sinp;	/* cosine and sine of the */
    /* azimuthal angle psi. */
    double ux = Photon_Ptr->ux[vec_num];
    double uy = Photon_Ptr->uy[vec_num];
    double uz = Photon_Ptr->uz[vec_num];
    double psi;

    cost = SpinTheta(g, stream, result, count);
    sint = sqrt(1.0 - cost*cost);
    /* sqrt() is faster than sin(). */

    psi = 2.0*PI*mkl_rand(stream, result, count); /* spin psi 0-2pi. */
    (*count) ++;
    (*count) %= 1024;
    cosp = cos(psi);
    if(psi<PI)
        sinp = sqrt(1.0 - cosp*cosp);
    /* sqrt() is faster than sin(). */
    else
        sinp = - sqrt(1.0 - cosp*cosp);

    if(fabs(uz) > COSZERO)  { 	/* normal incident. */
        Photon_Ptr->ux[vec_num] = sint*cosp;
        Photon_Ptr->uy[vec_num] = sint*sinp;
        Photon_Ptr->uz[vec_num] = cost*SIGN(uz);
        /* SIGN() is faster than division. */
    } else  {		/* regular incident. */
        double temp = sqrt(1.0 - uz*uz);
        Photon_Ptr->ux[vec_num] = sint*(ux*uz*cosp - uy*sinp)
                         /temp + ux*cost;
        Photon_Ptr->uy[vec_num] = sint*(uy*uz*cosp + ux*sinp)
                         /temp + uy*cost;
        Photon_Ptr->uz[vec_num] = -sint*cosp*temp + uz*cost;
    }
}

/***********************************************************
 *	Move the photon s away in the current layer of medium.
 ****/
void Hop(Photon8Struct *	Photon_Ptr, int vec_num)
{
    double s = Photon_Ptr->s[vec_num];

    Photon_Ptr->x[vec_num] += s*Photon_Ptr->ux[vec_num];
    Photon_Ptr->y[vec_num] += s*Photon_Ptr->uy[vec_num];
    Photon_Ptr->z[vec_num] += s*Photon_Ptr->uz[vec_num];
}

/***********************************************************
 *	If uz != 0, return the photon step size in glass,
 *	Otherwise, return 0.
 *
 *	The step size is the distance between the current
 *	position and the boundary in the photon direction.
 *
 *	Make sure uz !=0 before calling this function.
 ****/
void StepSizeInGlass(Photon8Struct *  Photon_Ptr,
                     InputStruct  *  In_Ptr,
                     int vec_num)
{
    double dl_b;	/* step size to boundary. */
    short  layer = Photon_Ptr->layer[vec_num];
    double uz = Photon_Ptr->uz[vec_num];

    /* Stepsize to the boundary. */
    if(uz>0.0)
        dl_b = (In_Ptr->layerspecs[layer].z1 - Photon_Ptr->z[vec_num])
               /uz;
    else if(uz<0.0)
        dl_b = (In_Ptr->layerspecs[layer].z0 - Photon_Ptr->z[vec_num])
               /uz;
    else
        dl_b = 0.0;

    Photon_Ptr->s[vec_num] = dl_b;
}

/***********************************************************
 *	Pick a step size for a photon packet when it is in
 *	tissue.
 *	If the member sleft is zero, make a new step size
 *	with: -log(rnd)/(mua+mus).
 *	Otherwise, pick up the leftover in sleft.
 *
 *	Layer is the index to layer.
 *	In_Ptr is the input parameters.
 ****/
void StepSizeInTissue(Photon8Struct * Photon_Ptr,
                      InputStruct  * In_Ptr,
			VSLStreamStatePtr  stream, double * result, short * count, int vec_num)
{
    short  layer = Photon_Ptr->layer[vec_num];
    double mua = In_Ptr->layerspecs[layer].mua;
    double mus = In_Ptr->layerspecs[layer].mus;

    if(Photon_Ptr->sleft[vec_num] == 0.0) {  /* make a new step. */
        double rnd;

        do{
		rnd = mkl_rand(stream, result, count);
		(*count) ++;
        	(*count) %= 1024;
	}
        while( rnd <= 0.0 );    /* avoid zero. */
        Photon_Ptr->s[vec_num] = -log(rnd)/(mua+mus);
    } else {	/* take the leftover. */
        Photon_Ptr->s[vec_num] = Photon_Ptr->sleft[vec_num]/(mua+mus);
        Photon_Ptr->sleft[vec_num] = 0.0;
    }
}

/***********************************************************
 *	Check if the step will hit the boundary.
 *	Return 1 if hit boundary.
 *	Return 0 otherwise.
 *
 * 	If the projected step hits the boundary, the members
 *	s and sleft of Photon_Ptr are updated.
 ****/
Boolean HitBoundary(Photon8Struct *  Photon_Ptr,
                    InputStruct  *  In_Ptr, int vec_num)
{
    double dl_b;  /* length to boundary. */
    short  layer = Photon_Ptr->layer[vec_num];
    double uz = Photon_Ptr->uz[vec_num];
    Boolean hit;

    /* Distance to the boundary. */
    if(uz>0.0)
        dl_b = (In_Ptr->layerspecs[layer].z1
                - Photon_Ptr->z[vec_num])/uz;	/* dl_b>0. */
    else if(uz<0.0)
        dl_b = (In_Ptr->layerspecs[layer].z0
                - Photon_Ptr->z[vec_num])/uz;	/* dl_b>0. */

    if(uz != 0.0 && Photon_Ptr->s[vec_num] > dl_b) {
        /* not horizontal & crossing. */
        double mut = In_Ptr->layerspecs[layer].mua
                     + In_Ptr->layerspecs[layer].mus;

        Photon_Ptr->sleft[vec_num] = (Photon_Ptr->s[vec_num] - dl_b)*mut;
        Photon_Ptr->s[vec_num]    = dl_b;
        hit = 1;
    } else
        hit = 0;

    return(hit);
}

/***********************************************************
 *	Drop photon weight inside the tissue (not glass).
 *
 *  The photon is assumed not dead.
 *
 *	The weight drop is dw = w*mua/(mua+mus).
 *
 *	The dropped weight is assigned to the absorption array
 *	elements.
 ****/
void Drop(InputStruct  *	In_Ptr,
          Photon8Struct *	Photon_Ptr,
          OutStruct *		Out_Ptr, int vec_num)
{
    double dwa;		/* absorbed weight.*/
    double x = Photon_Ptr->x[vec_num];
    double y = Photon_Ptr->y[vec_num];
    short  iz, ir;	/* index to z & r. */
    short  layer = Photon_Ptr->layer[vec_num];
    double mua, mus;

    /* compute array indices. */
    iz = (short)(Photon_Ptr->z[vec_num]/In_Ptr->dz);
    if(iz>In_Ptr->nz-1) iz=In_Ptr->nz-1;

    ir = (short)(sqrt(x*x+y*y)/In_Ptr->dr);
    if(ir>In_Ptr->nr-1) ir=In_Ptr->nr-1;

    /* update photon weight. */
    mua = In_Ptr->layerspecs[layer].mua;
    mus = In_Ptr->layerspecs[layer].mus;
    dwa = Photon_Ptr->w[vec_num] * mua/(mua+mus);
    Photon_Ptr->w[vec_num] -= dwa;

    /* assign dwa to the absorption array element. */
    Out_Ptr->A_rz[ir][iz]	+= dwa;
}

/***********************************************************
 *	The photon weight is small, and the photon packet tries
 *	to survive a roulette.
 ****/
void Roulette(Photon8Struct * Photon_Ptr, VSLStreamStatePtr  stream, double * result, short * count, int vec_num) 
{
    if(Photon_Ptr->w[vec_num] == 0.0)
        Photon_Ptr->dead[vec_num] = 1;
    else {
    	if(mkl_rand(stream, result, count) < CHANCE) /* survived the roulette.*/
        	Photon_Ptr->w[vec_num] /= CHANCE;
    	else
        	Photon_Ptr->dead[vec_num] = 1;
	    (*count) ++;
        (*count) %= 1024;
    }
}

/***********************************************************
 *	Compute the Fresnel reflectance.
 *
 *	Make sure that the cosine of the incident angle a1
 *	is positive, and the case when the angle is greater
 *	than the critical angle is ruled out.
 *
 * 	Avoid trigonometric function operations as much as
 *	possible, because they are computation-intensive.
 ****/
double RFresnel(double n1,	/* incident refractive index.*/
                double n2,	/* transmit refractive index.*/
                double ca1,	/* cosine of the incident */
                /* angle. 0<a1<90 degrees. */
                double * ca2_Ptr)  /* pointer to the */
/* cosine of the transmission */
/* angle. a2>0. */
{
    double r;

    if(n1==n2) {			  	/** matched boundary. **/
        *ca2_Ptr = ca1;
        r = 0.0;
    } else if(ca1>COSZERO) {	/** normal incident. **/
        *ca2_Ptr = ca1;
        r = (n2-n1)/(n2+n1);
        r *= r;
    } else if(ca1<COS90D)  {	/** very slant. **/
        *ca2_Ptr = 0.0;
        r = 1.0;
    } else  {			  		/** general. **/
        double sa1, sa2;
        /* sine of the incident and transmission angles. */
        double ca2;

        sa1 = sqrt(1-ca1*ca1);
        sa2 = n1*sa1/n2;
        if(sa2>=1.0) {
            /* double check for total internal reflection. */
            *ca2_Ptr = 0.0;
            r = 1.0;
        } else  {
            double cap, cam;	/* cosines of the sum ap or */
            /* difference am of the two */
            /* angles. ap = a1+a2 */
            /* am = a1 - a2. */
            double sap, sam;	/* sines. */

            *ca2_Ptr = ca2 = sqrt(1-sa2*sa2);

            cap = ca1*ca2 - sa1*sa2; /* c+ = cc - ss. */
            cam = ca1*ca2 + sa1*sa2; /* c- = cc + ss. */
            sap = sa1*ca2 + ca1*sa2; /* s+ = sc + cs. */
            sam = sa1*ca2 - ca1*sa2; /* s- = sc - cs. */
            r = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam);
            /* rearranged for speed. */
        }
    }
    return(r);
}

/***********************************************************
 *	Record the photon weight exiting the first layer(uz<0),
 *	no matter whether the layer is glass or not, to the
 *	reflection array.
 *
 *	Update the photon weight as well.
 ****/
void RecordR(double			Refl,	/* reflectance. */
             InputStruct  *	In_Ptr,
             Photon8Struct *	Photon_Ptr,
             OutStruct *	Out_Ptr, int vec_num)
{
    double x = Photon_Ptr->x[vec_num];
    double y = Photon_Ptr->y[vec_num];
    short  ir, ia;	/* index to r & angle. */

    ir = (short)(sqrt(x*x+y*y)/In_Ptr->dr);
    if(ir>In_Ptr->nr-1) ir=In_Ptr->nr-1;

    ia = (short)(acos(-Photon_Ptr->uz[vec_num])/In_Ptr->da);
    if(ia>In_Ptr->na-1) ia=In_Ptr->na-1;

    /* assign photon to the reflection array element. */
    Out_Ptr->Rd_ra[ir][ia] += Photon_Ptr->w[vec_num]*(1.0-Refl);

    Photon_Ptr->w[vec_num] *= Refl;
}

/***********************************************************
 *	Record the photon weight exiting the last layer(uz>0),
 *	no matter whether the layer is glass or not, to the
 *	transmittance array.
 *
 *	Update the photon weight as well.
 ****/
void RecordT(double 		Refl,
             InputStruct  *	In_Ptr,
             Photon8Struct *	Photon_Ptr,
             OutStruct *	Out_Ptr, int vec_num)
{
    double x = Photon_Ptr->x[vec_num];
    double y = Photon_Ptr->y[vec_num];
    short  ir, ia;	/* index to r & angle. */

    ir = (short)(sqrt(x*x+y*y)/In_Ptr->dr);
    if(ir>In_Ptr->nr-1) ir=In_Ptr->nr-1;

    ia = (short)(acos(Photon_Ptr->uz[vec_num])/In_Ptr->da);
    if(ia>In_Ptr->na-1) ia=In_Ptr->na-1;

    /* assign photon to the transmittance array element. */
    Out_Ptr->Tt_ra[ir][ia] += Photon_Ptr->w[vec_num]*(1.0-Refl);

    Photon_Ptr->w[vec_num] *= Refl;
}

/***********************************************************
 *	Decide whether the photon will be transmitted or
 *	reflected on the upper boundary (uz<0) of the current
 *	layer.
 *
 *	If "layer" is the first layer, the photon packet will
 *	be partially transmitted and partially reflected if
 *	PARTIALREFLECTION is set to 1,
 *	or the photon packet will be either transmitted or
 *	reflected determined statistically if PARTIALREFLECTION
 *	is set to 0.
 *
 *	Record the transmitted photon weight as reflection.
 *
 *	If the "layer" is not the first layer and the photon
 *	packet is transmitted, move the photon to "layer-1".
 *
 *	Update the photon parmameters.
 ****/
void CrossUpOrNot(InputStruct  *	In_Ptr,
                  Photon8Struct *	Photon_Ptr,
                  OutStruct *		Out_Ptr,
		VSLStreamStatePtr  stream, double * result, short * count, int vec_num)
{
    double uz = Photon_Ptr->uz[vec_num]; /* z directional cosine. */
    double uz1;	/* cosines of transmission alpha. always */
    /* positive. */
    double r=0.0;	/* reflectance */
    short  layer = Photon_Ptr->layer[vec_num];
    double ni = In_Ptr->layerspecs[layer].n;
    double nt = In_Ptr->layerspecs[layer-1].n;

    /* Get r. */
    if( - uz <= In_Ptr->layerspecs[layer].cos_crit0)
        r=1.0;		      /* total internal reflection. */
    else r = RFresnel(ni, nt, -uz, &uz1);

#if PARTIALREFLECTION
    if(layer == 1 && r<1.0) {	/* partially transmitted. */
        Photon_Ptr->uz[vec_num] = -uz1;	/* transmitted photon. */
        RecordR(r, In_Ptr, Photon_Ptr, Out_Ptr, vec_num);
        Photon_Ptr->uz[vec_num] = -uz;	/* reflected photon. */
    } 
    else 
    {
	if(mkl_rand(stream, result, count) > r) { /* transmitted to layer-1. */
        	Photon_Ptr->layer[vec_num]--;
	        Photon_Ptr->ux[vec_num] *= ni/nt;
	        Photon_Ptr->uy[vec_num] *= ni/nt;
	        Photon_Ptr->uz[vec_num] = -uz1;
        } 
        else			      		/* reflected. */
        	Photon_Ptr->uz[vec_num] = -uz;
	(*count) ++;
        (*count) %= 1024;
    }
#else
    if(mkl_rand(stream, result, count) > r) {		/* transmitted to layer-1. */
        if(layer==1)  {
            Photon_Ptr->uz[vec_num] = -uz1;
            RecordR(0.0, In_Ptr, Photon_Ptr, Out_Ptr, vec_num);
            Photon_Ptr->dead[vec_num] = 1;
        } else {
            Photon_Ptr->layer[vec_num]--;
            Photon_Ptr->ux[vec_num] *= ni/nt;
            Photon_Ptr->uy[vec_num] *= ni/nt;
            Photon_Ptr->uz[vec_num] = -uz1;
        }
    } else 						/* reflected. */
        Photon_Ptr->uz[vec_num] = -uz;
    (*count) ++;
    (*count) %= 1024;
#endif
}

/***********************************************************
 *	Decide whether the photon will be transmitted  or be
 *	reflected on the bottom boundary (uz>0) of the current
 *	layer.
 *
 *	If the photon is transmitted, move the photon to
 *	"layer+1". If "layer" is the last layer, record the
 *	transmitted weight as transmittance. See comments for
 *	CrossUpOrNot.
 *
 *	Update the photon parmameters.
 ****/
void CrossDnOrNot(InputStruct  *	In_Ptr,
                  Photon8Struct *	Photon_Ptr,
                  OutStruct *		Out_Ptr,
		 VSLStreamStatePtr  stream, double * result, short * count, int vec_num)
{
    double uz = Photon_Ptr->uz[vec_num]; /* z directional cosine. */
    double uz1;	/* cosines of transmission alpha. */
    double r=0.0;	/* reflectance */
    short  layer = Photon_Ptr->layer[vec_num];
    double ni = In_Ptr->layerspecs[layer].n;
    double nt = In_Ptr->layerspecs[layer+1].n;

    /* Get r. */
    if( uz <= In_Ptr->layerspecs[layer].cos_crit1)
        r=1.0;		/* total internal reflection. */
    else r = RFresnel(ni, nt, uz, &uz1);

#if PARTIALREFLECTION
    if(layer == In_Ptr->num_layers && r<1.0) {
        Photon_Ptr->uz[vec_num] = uz1;
        RecordT(r, In_Ptr, Photon_Ptr, Out_Ptr, vec_num);
        Photon_Ptr->uz[vec_num] = -uz;
    } else {
	 if(mkl_rand(stream, result, count) > r) { /* transmitted to layer+1. */
        	Photon_Ptr->layer[vec_num]++;
	        Photon_Ptr->ux[vec_num] *= ni/nt;
	        Photon_Ptr->uy[vec_num] *= ni/nt;
	        Photon_Ptr->uz[vec_num] = uz1;
    	} 
	else 						/* reflected. */
        	Photon_Ptr->uz[vec_num] = -uz;
	(*count) ++;
        (*count) %= 1024;
   }
#else
    if(mkl_rand(stream, result, count) > r) {		/* transmitted to layer+1. */
        if(layer == In_Ptr->num_layers) {
            Photon_Ptr->uz[vec_num] = uz1;
            RecordT(0.0, In_Ptr, Photon_Ptr, Out_Ptr, vec_num);
            Photon_Ptr->dead[vec_num] = 1;
        } else {
            Photon_Ptr->layer[vec_num]++;
            Photon_Ptr->ux[vec_num] *= ni/nt;
            Photon_Ptr->uy[vec_num] *= ni/nt;
            Photon_Ptr->uz[vec_num] = uz1;
        }
    } else 						/* reflected. */
        Photon_Ptr->uz[vec_num] = -uz;
    (*count) ++;
    (*count) %= 1024;
#endif
}

/***********************************************************
 ****/
void CrossOrNot(InputStruct  *	In_Ptr,
                Photon8Struct *	Photon_Ptr,
                OutStruct    *	Out_Ptr,
		VSLStreamStatePtr  stream, double * result, short * count, int vec_num)
{
    if(Photon_Ptr->uz[vec_num] < 0.0)
        CrossUpOrNot(In_Ptr, Photon_Ptr, Out_Ptr, stream, result, count, vec_num);
    else
        CrossDnOrNot(In_Ptr, Photon_Ptr, Out_Ptr, stream, result, count, vec_num);
}

/***********************************************************
 *	Move the photon packet in glass layer.
 *	Horizontal photons are killed because they will
 *	never interact with tissue again.
 ****/
void HopInGlass(InputStruct  * In_Ptr,
                Photon8Struct * Photon_Ptr,
                OutStruct    * Out_Ptr,
		VSLStreamStatePtr * stream, double * result, short * count, int vec_num)
{
    double dl;     /* step size. 1/cm */

    if(Photon_Ptr->uz[vec_num] == 0.0) {
        /* horizontal photon in glass is killed. */
        Photon_Ptr->dead[vec_num] = 1;
    } else {
        StepSizeInGlass(Photon_Ptr, In_Ptr, vec_num);
        Hop(Photon_Ptr, vec_num);
        CrossOrNot(In_Ptr, Photon_Ptr, Out_Ptr, stream, result, count, vec_num);
    }
}

/***********************************************************
 *	Set a step size, move the photon, drop some weight,
 *	choose a new photon direction for propagation.
 *
 *	When a step size is long enough for the photon to
 *	hit an interface, this step is divided into two steps.
 *	First, move the photon to the boundary free of
 *	absorption or scattering, then decide whether the
 *	photon is reflected or transmitted.
 *	Then move the photon in the current or transmission
 *	medium with the unfinished stepsize to interaction
 *	site.  If the unfinished stepsize is still too long,
 *	repeat the above process.
 ****/
void HopDropSpinInTissue(InputStruct  *  In_Ptr,
                         Photon8Struct *  Photon_Ptr,
                         OutStruct    *  Out_Ptr, VSLStreamStatePtr stream, double * result, short * count,                          int vec_num)
{
    StepSizeInTissue(Photon_Ptr, In_Ptr, stream, result, count, vec_num);

    if(HitBoundary(Photon_Ptr, In_Ptr, vec_num)) {
        Hop(Photon_Ptr, vec_num);	/* move to boundary plane. */
        CrossOrNot(In_Ptr, Photon_Ptr, Out_Ptr, stream, result, count, vec_num);
    } else {
        Hop(Photon_Ptr,vec_num);
        Drop(In_Ptr, Photon_Ptr, Out_Ptr, vec_num);
        Spin(In_Ptr->layerspecs[Photon_Ptr->layer[vec_num]].g,
             Photon_Ptr, stream, result, count, vec_num);
    }
}

/***********************************************************
 ****/
void HopDropSpin(InputStruct  *  In_Ptr,
                 Photon8Struct *  Photon_Ptr,
                 OutStruct    *  Out_Ptr,
		VSLStreamStatePtr  stream, double * result, short * count, int vec_num)
{
    short layer = Photon_Ptr->layer[vec_num];

    if((In_Ptr->layerspecs[layer].mua == 0.0)
            && (In_Ptr->layerspecs[layer].mus == 0.0))
        /* glass layer. */
        HopInGlass(In_Ptr, Photon_Ptr, Out_Ptr, stream, result, count, vec_num);
    else
        HopDropSpinInTissue(In_Ptr, Photon_Ptr, Out_Ptr, stream, result, count, vec_num);

    if( Photon_Ptr->w[vec_num] < In_Ptr->Wth && !Photon_Ptr->dead[vec_num])
        Roulette(Photon_Ptr, stream, result, count, vec_num);
}
