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

void LaunchPhoton(double Rspecular,
        LayerStruct  * Layerspecs_Ptr,
        Photon8Struct * Photon_Ptr, int vec_num)
{
    Photon_Ptr->w[vec_num]       = 1.0 - Rspecular;
    Photon_Ptr->dead[vec_num]    = 0;
    Photon_Ptr->layer[vec_num] = 1;
    Photon_Ptr->s[vec_num]   = 0;
    Photon_Ptr->sleft[vec_num]= 0;

    Photon_Ptr->x[vec_num]   = 0.0;
    Photon_Ptr->y[vec_num]       = 0.0;
    Photon_Ptr->z[vec_num]       = 0.0;
    Photon_Ptr->ux[vec_num]  = 0.0;
    Photon_Ptr->uy[vec_num]  = 0.0;
    Photon_Ptr->uz[vec_num]  = 1.0;

    if((Layerspecs_Ptr[1].mua == 0.0)
            && (Layerspecs_Ptr[1].mus == 0.0))  { /* glass layer. */
        Photon_Ptr->layer[vec_num]   = 2;
        Photon_Ptr->z[vec_num]   = Layerspecs_Ptr[2].z0;
    }   
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
 ****/
void HopDropSpin(InputStruct  *  In_Ptr,
        Photon8Struct *  Photon_Ptr,
        OutStruct *  Out_Ptr,
        VSLStreamStatePtr  stream)
{
    int i;
    double s[8];
    short  layer[8];// = Photon_Ptr->layer;
    double mua[8];// = In_Ptr->layerspecs[layer].mua;
    double mus[8];// = In_Ptr->layerspecs[layer].mus;

    double dl_b[8];  /* length to boundary. */
    double uz[8];// = Photon_Ptr->uz;
    double mut[8];
    Boolean hit[8];
    double dwa[8];		/* absorbed weight.*/
    double x[8];// = Photon_Ptr->x;
    double y[8];// = Photon_Ptr->y;
    short  iz[8];	/* index to z & r. */

    double cost[8], sint[8];	/* cosine and sine of the */
    /* polar deflection angle theta. */
    double cosp[8], sinp[8];	/* cosine and sine of the */
    /* azimuthal angle psi. */
    double ux[8] ;
    double uy[8] ;//= Photon_Ptr->uy;
    double psi[8];
    double g[8];
    double temp[8];
    double cap[8], cam[8];	/* cosines of the sum ap or */
    /* difference am of the two */
    /* angles. ap = a1+a2 */
    /* am = a1 - a2. */
    double sap[8], sam[8];	/* sines. */
    double sa1[8], sa2[8];
    /* sine of the incident and transmission angles. */
    double ca2[8];
    double n1_temp[8] , n2_temp[8] , ca1_temp[8],* ca2_Ptr_temp[8] ;
    double uz1[8];	/* cosines of transmission alpha. always */
    /* positive. */
    double r[8];//=0.0;	/* reflectance */
    double ni[8];// = In_Ptr->layerspecs[layer].n;
    double nt[8];
    short  ir[8], ia[8];	/* index to r & angle. */
    double    result_t[8*7];
    vdRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, 8*9 , &result_t, 0.0, 1.0 );
    /////////////////////////////////////////////////////////////////////////////////////////////
#pragma simd
    for(i = 0 ; i<8 ; i++)
    {    
        layer[i] = Photon_Ptr->layer[i];
        mua[i] = In_Ptr->layerspecs[layer[i]].mua;
        mus[i] = In_Ptr->layerspecs[layer[i]].mus;
        if(Photon_Ptr->sleft[i] == 0.0) {  /* make a new step. */

            Photon_Ptr->s[i] = -log(result_t[i])/(mua[i]+mus[i]);
        } else {	/* take the leftover. */
            Photon_Ptr->s[i] = Photon_Ptr->sleft[i]/(mua[i]+mus[i]);
            Photon_Ptr->sleft[i] = 0.0;
        }
        /////////////////////////////////////////////////////////////////////////////////////////////
        layer[i] = Photon_Ptr->layer[i];
        uz[i] = Photon_Ptr->uz[i];
        /* Distance to the boundary. */
        if(uz[i]>0.0)
            dl_b[i] = (In_Ptr->layerspecs[layer[i]].z1
                    - Photon_Ptr->z[i])/uz[i];	/* dl_b>0. */
        else if(uz[i]<0.0)
            dl_b[i] = (In_Ptr->layerspecs[layer[i]].z0
                    - Photon_Ptr->z[i])/uz[i];	/* dl_b>0. */

        if(uz[i] != 0.0 && Photon_Ptr->s[i] > dl_b[i]) {
            /* not horizontal & crossing. */
            mut[i] = In_Ptr->layerspecs[layer[i]].mua
                + In_Ptr->layerspecs[layer[i]].mus;

            Photon_Ptr->sleft[i] = (Photon_Ptr->s[i] - dl_b[i])*mut[i];
            Photon_Ptr->s[i]    = dl_b[i];
            hit[i] = 1;
        } else
            hit[i] = 0;


        if(hit[i]) {

            s[i] = Photon_Ptr->s[i];
            Photon_Ptr->x[i] += s[i]*Photon_Ptr->ux[i];
            Photon_Ptr->y[i] += s[i]*Photon_Ptr->uy[i];
            Photon_Ptr->z[i] += s[i]*Photon_Ptr->uz[i];
            //        CrossOrNot(In_Ptr, Photon_Ptr, tmpOut_Ptr, rand_seed);
            if(Photon_Ptr->uz[i] < 0.0)
            {
                uz[i] = Photon_Ptr->uz[i]; /* z directional cosine. */
                r[i]=0.0;	/* reflectance */
                layer[i] = Photon_Ptr->layer[i];
                ni[i] = In_Ptr->layerspecs[layer[i]].n;
                nt[i] = In_Ptr->layerspecs[layer[i]-1].n;

                /* Get r. */
                if( - uz[i] <= In_Ptr->layerspecs[layer[i]].cos_crit0)
                    r[i]=1.0;		      /* total internal reflection. */
                else{
                    //     r = RFresnel(ni, nt, -uz, &uz1);
                    n1_temp[i] = ni[i];
                    n2_temp[i] = nt[i];
                    ca1_temp[i] = -uz[i];
                    ca2_Ptr_temp[i] = &uz1[i];
                    //    double r;

                    if(n1_temp[i]==n2_temp[i]) {			  	/** matched boundary. **/
                        *ca2_Ptr_temp[i] = ca1_temp[i];
                        r[i] = 0.0;
                    } else if(ca1_temp[i]>COSZERO) {	/** normal incident. **/
                        *(ca2_Ptr_temp[i]) = ca1_temp[i];
                        r[i] = (n2_temp[i]-n1_temp[i])/(n2_temp[i]+n1_temp[i]);
                        r[i] *= r[i];
                    } else if(ca1_temp[i]<COS90D)  {	/** very slant. **/
                        *(ca2_Ptr_temp[i]) = 0.0;
                        r[i] = 1.0;
                    } else  {			  		/** general. **/

                        sa1[i] = sqrt(1-ca1_temp[i]*ca1_temp[i]);
                        sa2[i] = n1_temp[i]*sa1[i]/n2_temp[i];
                        if(sa2[i]>=1.0) {
                            /* double check for total internal reflection. */
                            *ca2_Ptr_temp[i] = 0.0;
                            r[i] = 1.0;
                        } else  {

                            *(ca2_Ptr_temp[i]) = ca2[i] = sqrt(1-sa2[i]*sa2[i]);

                            cap[i] = ca1_temp[i]*ca2[i] - sa1[i]*sa2[i]; /* c+ = cc - ss. */
                            cam[i] = ca1_temp[i]*ca2[i] + sa1[i]*sa2[i]; /* c- = cc + ss. */
                            sap[i] = sa1[i]*ca2[i] + ca1_temp[i]*sa2[i]; /* s+ = sc + cs. */
                            sam[i] = sa1[i]*ca2[i] - ca1_temp[i]*sa2[i]; /* s- = sc - cs. */
                            r[i] = 0.5*sam[i]*sam[i]*(cam[i]*cam[i]+cap[i]*cap[i])/(sap[i]*sap[i]*cam[i]*cam[i]);
                            /* rearranged for speed. */
                        }
                    }

                }

                if(result_t[i+8] > r[i]) {		/* transmitted to layer-1. */
                    if(layer[i]==1)  {
                        Photon_Ptr->uz[i] = -uz1[i];
                        //            RecordR(0.0, In_Ptr, Photon_Ptr, tmpOut_Ptr);
                        x[i] = Photon_Ptr->x[i];
                        y[i]  = Photon_Ptr->y[i];

                        ir[i] = (short)(sqrt(x[i]*x[i]+y[i]*y[i])/In_Ptr->dr);
                        if(ir[i]>In_Ptr->nr-1) ir[i]=In_Ptr->nr-1;

                        ia[i] = (short)(acos(-Photon_Ptr->uz[i])/In_Ptr->da);
                        if(ia[i]>In_Ptr->na-1) ia[i]=In_Ptr->na-1;

                        /* assign photon to the reflection array element. */
                        //    tmpOut_Ptr->Rd_ra[ir][ia] += Photon_Ptr->w*(1.0-Refl);
                        Out_Ptr->Rd_ra[ir[i]][ia[i]] += Photon_Ptr->w[i]*(1.0);

                        Photon_Ptr->w[i] *= 0.0;

                        Photon_Ptr->dead[i] = 1;
                    } else {
                        Photon_Ptr->layer[i]--;
                        Photon_Ptr->ux[i] *= ni[i]/nt[i];
                        Photon_Ptr->uy[i] *= ni[i]/nt[i];
                        Photon_Ptr->uz[i] = -uz1[i];
                    }
                } else 						/* reflected. */
                    Photon_Ptr->uz[i] = -uz[i];

            }
            else
            {
                //        CrossDnOrNot(In_Ptr, Photon_Ptr, tmpOut_Ptr, rand_seed);
                uz[i] = Photon_Ptr->uz[i]; /* z directional cosine. */
                uz1[i];	/* cosines of transmission alpha. */
                r[i]=0.0;	/* reflectance */
                layer[i] = Photon_Ptr->layer[i];
                ni[i] = In_Ptr->layerspecs[layer[i]].n;
                nt[i] = In_Ptr->layerspecs[layer[i]+1].n;

                /* Get r. */
                if( uz[i] <= In_Ptr->layerspecs[layer[i]].cos_crit1)
                    r[i]=1.0;		/* total internal reflection. */
                else{
                    //         r = RFresnel(ni, nt, uz, &uz1);
                    //       double r;
                    n1_temp[i] = ni[i];
                    n2_temp[i] = nt[i];
                    ca1_temp[i] = uz[i];
                    ca2_Ptr_temp[i] = &uz1[i];

                    if(n1_temp[i]==n2_temp[i]) {			  	/** matched boundary. **/
                        *ca2_Ptr_temp[i] = ca1_temp[i];
                        r[i] = 0.0;
                    } else if(ca1_temp[i]>COSZERO) {	/** normal incident. **/
                        *ca2_Ptr_temp[i] = ca1_temp[i];
                        r[i] = (n2_temp[i]-n1_temp[i])/(n2_temp[i]+n1_temp[i]);
                        r[i] *= r[i];
                    } else if(ca1_temp[i]<COS90D)  {	/** very slant. **/
                        *(ca2_Ptr_temp[i]) = 0.0;
                        r[i] = 1.0;
                    } else  {			  		/** general. **/

                        sa1[i] = sqrt(1-ca1_temp[i]*ca1_temp[i]);
                        sa2[i] = n1_temp[i]*sa1[i]/n2_temp[i];
                        if(sa2[i]>=1.0) {
                            /* double check for total internal reflection. */
                            *(ca2_Ptr_temp[i]) = 0.0;
                            r[i] = 1.0;
                        } else  {

                            *(ca2_Ptr_temp[i]) = ca2[i] = sqrt(1-sa2[i]*sa2[i]);

                            cap[i] = ca1_temp[i]*ca2[i] - sa1[i]*sa2[i]; /* c+ = cc - ss. */
                            cam[i] = ca1_temp[i]*ca2[i] + sa1[i]*sa2[i]; /* c- = cc + ss. */
                            sap[i] = sa1[i]*ca2[i] + ca1_temp[i]*sa2[i]; /* s+ = sc + cs. */
                            sam[i] = sa1[i]*ca2[i] - ca1_temp[i]*sa2[i]; /* s- = sc - cs. */
                            r[i] = 0.5*sam[i]*sam[i]*(cam[i]*cam[i]+cap[i]*cap[i])/(sap[i]*sap[i]*cam[i]*cam[i]);
                            /* rearranged for speed. */
                        }
                    }

                }

                if(result_t[i+16] > r[i]) {		/* transmitted to layer+1. */
                    if(layer[i] == In_Ptr->num_layers) {
                        Photon_Ptr->uz[i] = uz1[i];
                        //            RecordT(0.0, In_Ptr, Photon_Ptr, tmpOut_Ptr);
                        x[i] = Photon_Ptr->x[i];
                        y[i] = Photon_Ptr->y[i];

                        ir[i] = (short)(sqrt(x[i]*x[i]+y[i]*y[i])/In_Ptr->dr);
                        if(ir[i]>In_Ptr->nr-1) ir[i]=In_Ptr->nr-1;

                        ia[i] = (short)(acos(Photon_Ptr->uz[i])/In_Ptr->da);
                        if(ia[i]>In_Ptr->na-1) ia[i]=In_Ptr->na-1;

                        /* assign photon to the transmittance array element. */
                        //    Out_Ptr->Tt_ra[ir][ia] += Photon_Ptr->w*(1.0-Refl);
                        Out_Ptr->Tt_ra[ir[i]][ia[i]] += Photon_Ptr->w[i]*(1.0);     

                        Photon_Ptr->w[i] *= 0.0;

                        Photon_Ptr->dead[i] = 1;
                    } else {
                        Photon_Ptr->layer[i]++;
                        Photon_Ptr->ux[i] *= ni[i]/nt[i];
                        Photon_Ptr->uy[i] *= ni[i]/nt[i];
                        Photon_Ptr->uz[i] = uz1[i];
                    }
                } else 						/* reflected. */
                    Photon_Ptr->uz[i] = -uz[i];

            }
        } else {
            s[i] = Photon_Ptr->s[i];
            Photon_Ptr->x[i] += s[i]*Photon_Ptr->ux[i];
            Photon_Ptr->y[i] += s[i]*Photon_Ptr->uy[i];
            Photon_Ptr->z[i] += s[i]*Photon_Ptr->uz[i];
            //        Drop(In_Ptr, Photon_Ptr, tmpOut_Ptr);
            /* compute array indices. */
            x[i] = Photon_Ptr->x[i];
            y[i] = Photon_Ptr->y[i];
            layer[i] = Photon_Ptr->layer[i];
            iz[i] = (short)(Photon_Ptr->z[i]/In_Ptr->dz);
            if(iz[i]>In_Ptr->nz-1) iz[i]=In_Ptr->nz-1;

            ir[i] = (short)(sqrt(x[i]*x[i]+y[i]*y[i])/In_Ptr->dr);
            if(ir[i]>In_Ptr->nr-1) ir[i]=In_Ptr->nr-1;

            /* update photon weight. */
            mua[i] = In_Ptr->layerspecs[layer[i]].mua;
            mus[i] = In_Ptr->layerspecs[layer[i]].mus;
            dwa[i] = Photon_Ptr->w[i] * mua[i]/(mua[i]+mus[i]);
            Photon_Ptr->w[i] -= dwa[i];

            /* assign dwa to the absorption array element. */
            //    Out_Ptr->A_rz[ir][iz]	+= dwa;

            Out_Ptr->A_rz[ir[i]][iz[i]] += dwa[i]; 
            //        Spin(In_Ptr->layerspecs[Photon_Ptr->layer].g,
            //             Photon_Ptr, rand_seed);
            g[i] = In_Ptr->layerspecs[Photon_Ptr->layer[i]].g;
            ux[i] = Photon_Ptr->ux[i];
            uy[i] = Photon_Ptr->uy[i];
            uz[i] = Photon_Ptr->uz[i];
            //    cost = SpinTheta(g, rand_seed);
            if(g[i] == 0.0)
                cost[i] = 2*result_t[i+24] -1;
            else {
                temp[i] = (1-g[i]*g[i])/(1-g[i]+2*g[i]*result_t[i+32]);
                cost[i] = (1+g[i]*g[i] - temp[i]*temp[i])/(2*g[i]);
                if(cost[i] < -1) cost[i] = -1;
                else if(cost[i] > 1) cost[i] = 1;
            }

            sint[i] = sqrt(1.0 - cost[i]*cost[i]);
            /* sqrt() is faster than sin(). */

            psi[i] = 2.0*PI*result_t[i+40]; /* spin psi 0-2pi. */
            cosp[i] = cos(psi[i]);
            if(psi[i]<PI)
                sinp[i] = sqrt(1.0 - cosp[i]*cosp[i]);
            /* sqrt() is faster than sin(). */
            else
                sinp[i] = - sqrt(1.0 - cosp[i]*cosp[i]);

            if(fabs(uz[i]) > COSZERO)  { 	/* normal incident. */
                Photon_Ptr->ux[i] = sint[i]*cosp[i];
                Photon_Ptr->uy[i] = sint[i]*sinp[i];
                Photon_Ptr->uz[i] = cost[i]*SIGN(uz[i]);
                /* SIGN() is faster than division. */
            } else  {		/* regular incident. */
                temp[i] = sqrt(1.0 - uz[i]*uz[i]);
                Photon_Ptr->ux[i] = sint[i]*(ux[i]*uz[i]*cosp[i] - uy[i]*sinp[i])
                    /temp[i] + ux[i]*cost[i];
                Photon_Ptr->uy[i] = sint[i]*(uy[i]*uz[i]*cosp[i] + ux[i]*sinp[i])
                    /temp[i] + uy[i]*cost[i];
                Photon_Ptr->uz[i] = -sint[i]*cosp[i]*temp[i] + uz[i]*cost[i];
            }

        }


        if( Photon_Ptr->w[i] < In_Ptr->Wth && !Photon_Ptr->dead[i])
        {

            if(Photon_Ptr->w[i] == 0.0)
                Photon_Ptr->dead[i] = 1;
            else if(result_t[i+48] < CHANCE) /* survived the roulette.*/
                Photon_Ptr->w[i] /= CHANCE;
            else
                Photon_Ptr->dead[i] = 1;

        }
    }
}
