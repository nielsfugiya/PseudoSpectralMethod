//
//  main.cpp
//  pstd2d
//
//  Created by Jianing Zhang on 1/1/14.
//  Copyright (c) 2014 Jianing Zhang. All rights reserved.
//

#include <cmath>
#include "cassert"
#include <iostream>
#include "fstream"
#include "fftw3.h"
#include "complex"
#include <cstdlib>
#include "ctime"
using namespace std;

double grid_xsize, grid_ysize, dl, dt;
int grid_nx,grid_ny, grid_nz;
double sizeparameter;
double resolution;
double wavelength;
double epsR;
double f;
double wavenumber;
double cdtds =0.9*2.0/M_PI/sqrt(2.0);


int qtime, maxTime = 600;
int npml;

double kx_hat, ky_hat;


double **q_e_xz, **q_e_yz, **q_h_xy, **q_h_yx;
double *be_x, *ce_x, *bh_x, *ch_x, *kedx, *khdx, *be_y, *ce_y, *bh_y, *ch_y, *kedy, *khdy;

//---------------------------------------------------------
// FIELDS & MATERIALS
//---------------------------------------------------------

double **ez, **hx, **hy, **epsilon;
double **ezx, **ezy, **hxy, **hyx;

double *xpos, *ypos;


double phi_inc = 0.0;

fftw_complex *xdiff, *ydiff;



void ReadParameters();
void InitializeGrid();
void InitializePML();
void InitializeScatter();
void PerformPSTD();
void UpdateE();
void UpdateH();
void FreeMemory();
double gauss(double t);
void CurlE();
void CurlH();
void PseudoSpecDiff(fftw_complex *data, int N);
double J_inc(int i, int j, double n);
void CalculateCoeff(double *kedl, double *be_l, double *ce_l, double *khdl, double *bh_l, double *ch_l, int grid_nl,double grid_dl);
double* ChebyshevGaussLobattoNodes(int N);
void FastChebyshevTransform(double *f, int N, int s);
void ChebyshevDerivativeCoefficients(double* f_hat, int N);
void FastChebyshevDerivative(double* f,int N);
void FourierTransformFields();



int main(int argc, const char * argv[])
{
    
    ReadParameters();
    InitializeGrid();
    InitializePML();
    InitializeScatter();
    PerformPSTD();
    FreeMemory();
    return 0;
}


void ReadParameters()
{
    
    ifstream fin("parameters");
    assert(fin.is_open());
    
    
    // COMPUTATIONAL GRID
    
    fin >> sizeparameter;
    fin >> wavelength;
    fin >> epsR;
    fin >> resolution;
    
    fin >> phi_inc;
    fin >> npml;
    
    
    fin.close();
    
    return;
}

void InitializeGrid()
{
    
    //=========================================================
    // SET-UP GRID
    //=========================================================
    
    dl = wavelength/resolution;
    grid_xsize = sizeparameter*wavelength/M_PI;
    grid_ysize = grid_xsize;
    
    // GET NUMBER OF GRID POINTS
    grid_nx = static_cast<int>(grid_xsize/dl) + 3;
    grid_ny = static_cast<int>(grid_ysize/dl) + 3;
    grid_nx = grid_nx + 2*npml + 8;
    grid_ny = grid_nx;
    
    printf("Size Parameter:\t%5.2f\n",sizeparameter);
    printf("Wave Resolution:\t%5.2f\n",resolution);
    printf("Grid Number of CPML:\t%d\n",npml);
    printf("Grid Number in X direction:\t%d\n",grid_nx);
    printf("Grid Number in Y direction:\t%d\n",grid_ny);
    
    // ADJUST GRID SIZES
    grid_xsize = (static_cast<double>(grid_nx))*dl;
    grid_ysize = (static_cast<double>(grid_ny))*dl;
    
    dt = cdtds*dl;
    
    //=========================================================
    // FIELDS
    //=========================================================
    //
    // ALOCATE MEMORY & INITIALIZE
    //---------------------------------------------------------
    

    ez = new double*[grid_nx];
    hx = new double*[grid_nx];
    hy = new double*[grid_nx];
    
    
    epsilon= new double*[grid_nx];
    
    
    for (int i = 0;i < grid_nx; i++)
    {
        

        ez[i] = new double[grid_ny];
        hx[i] = new double[grid_ny];
        hy[i] = new double[grid_ny];
        
        epsilon[i] = new double[grid_ny];
        
        //cout<< "Hello Hello\n";
        for (int j = 0;j < grid_ny; j++)
        {

            ez[i][j] = 0.0;
            
            
            hx[i][j] = 0.0;
            hy[i][j] = 0.0;
            
            epsilon[i][j] = 1.0;


        }
        
    }

    ezx = new double*[grid_nx];
    hyx = new double*[grid_nx];
    ezy = new double*[grid_nx];
    hxy = new double*[grid_nx];

    
    
    
    for(int i = 0; i < grid_nx; ++i)
    {
        ezx[i] = new double[grid_ny];
        ezy[i] = new double[grid_ny];
        hyx[i] = new double[grid_ny];
        hxy[i] = new double[grid_ny];
        
        for(int j = 0; j < grid_ny; ++j)
        {

            ezx[i][j] = 0.0;
            hyx[i][j] = 0.0;
            ezy[i][j] = 0.0;
            hxy[i][j] = 0.0;
            
        }
    }
    
    f =1.0/wavelength;
    wavenumber = 2*M_PI*f;
    
}


void InitializeScatter()
{
    
    double radius;
    double x, y;
    
    radius = sizeparameter*wavelength/2.0/M_PI;
    
    phi_inc = phi_inc*M_PI/180.0;
    
    kx_hat = cos(phi_inc);
    ky_hat = sin(phi_inc);
/*
    xpos = new double[grid_nx];
    ypos = new double[grid_ny];
    
    xpos = ChebyshevGaussLobattoNodes(grid_nx);
    ypos = ChebyshevGaussLobattoNodes(grid_ny);*/

    for (int i = 0; i < grid_nx; i++)
    {
        x = ((double)(i - grid_nx/2))*dl;
        for (int j =0; j < grid_ny; j++)
        {
            y = ((double)(j - grid_ny/2))*dl;
            
            if(x*x + y*y < radius*radius)
                
            {
                
                epsilon[i][j] = epsR;
                
            }
            else
            {
                epsilon[i][j] = 1.0;
            }
        }
    }
    
/*    for (int i = 0; i < grid_nx; i++)
    {
        x = xpos[i];
        for (int j =0; j < grid_ny; j++)
        {
            y = ypos[j];
            if(x*x + y*y < radius*radius)
            {
                epsilon[i][j] = epsR;
            }
            else
            {
                epsilon[i][j] = 1.0;
            }
            
            cout<<epsilon[i][j]<<" ";
            
        }
        cout<<endl;
    }*/
    
    
    return;
}


void InitializePML()
{
    
    //=========================================================
    // Auxiliary Variables
    //=========================================================
    //
    // ALOCATE MEMORY & INITIALIZE
    //---------------------------------------------------------

    q_e_xz = new double*[grid_nx];
    q_e_yz = new double*[grid_nx];
    q_h_xy = new double*[grid_nx];
    q_h_yx = new double*[grid_nx];

    
    
    
    kedx = new double [grid_nx];
    khdx = new double [grid_nx];
    kedy = new double [grid_ny];
    khdy = new double [grid_ny];

    
    be_x = new double [grid_nx];
    ce_x = new double [grid_nx];
    bh_x = new double [grid_nx];
    ch_x = new double [grid_nx];
    
    be_y = new double [grid_ny];
    ce_y = new double [grid_ny];
    bh_y = new double [grid_ny];
    ch_y = new double [grid_ny];

    
    
    
    for(int i = 0; i < grid_nx; ++i)
    {
        q_h_yx[i] = new double[grid_ny];
        q_h_xy[i] = new double[grid_ny];
        q_e_xz[i] = new double[grid_ny];
        q_e_yz[i] = new double[grid_ny];

        
        for(int j = 0; j < grid_ny; ++j)
        {

            q_h_yx[i][j] = 0.0;
            q_h_xy[i][j] = 0.0;
            q_e_xz[i][j] = 0.0;
            q_e_yz[i][j] = 0.0;

            
            be_y[j] = 0.0;
            ce_y[j] = 0.0;
            bh_y[j] = 0.0;
            ch_y[j] = 0.0;
            
            kedy[j] = 1.0;
            khdy[j] = 1.0;
            
        } // ++j
        
        be_x[i] = 0.0;
        ce_x[i] = 0.0;
        bh_x[i] = 0.0;
        ch_x[i] = 0.0;
        
        kedx[i] = 1.0;
        khdx[i] = 1.0;
        
    } // ++i
    
    CalculateCoeff(kedx, be_x, ce_x, khdx, bh_x,  ch_x, grid_nx, dl);
    CalculateCoeff(kedy, be_y, ce_y, khdy, bh_y,  ch_y, grid_ny, dl);
    
}
double J_inc(int i, int j, double n)
{
    double t;
    t = n*dt - (kx_hat*i*dl+ ky_hat*j*dl);
    //return gauss(t);
    return sin(2*M_PI*f*t);
}

double gauss(double t)
{
    double tw = 5*dt;
    double to = 4*tw;
    double x = (t-to)/tw;
    return -x*exp(-x*x);
}

void CalculateCoeff(double *kedl, double *be_l, double *ce_l, double *khdl, double *bh_l, double *ch_l, int grid_nl,double grid_dl)
{
    double lpos, temp;
    double depth;
    double cpml_m = 4.0;
    double cpml_ma = 4.0;
    double cpml_depth = (static_cast<double>(npml))*grid_dl;
    double cpml_sigmamax = (0.8*(static_cast<double>(cpml_m + 1)))/(grid_dl);
    double cpml_kappamax = 11.0;
    double cpml_alphamax = 0.1;
    double cpml_sigma;
    double cpml_kappa;
    double cpml_alpha;
    double grid_lsize = grid_dl*grid_nl;
    
    
    for(int k = 0; k < grid_nl; ++k)
    {
        
        //---------------------------------------------------------
        // E-FIELD
        //---------------------------------------------------------
        lpos = (static_cast<double>(k))*grid_dl;
        if(lpos < cpml_depth)
        {
            // get the depth into the pml layer
            depth = cpml_depth - lpos;
            
            cpml_sigma = cpml_sigmamax*pow((depth/cpml_depth), cpml_m);
            cpml_kappa = 1.0 + (cpml_kappamax - 1.0)*pow((depth/cpml_depth), cpml_m);
            cpml_alpha = cpml_alphamax*pow(((cpml_depth-depth)/cpml_depth), cpml_ma);
            
            temp = cpml_sigma + cpml_kappa*cpml_alpha;
            kedl[k] = cpml_kappa*grid_dl;
            
            be_l[k] = exp(0.0 - (cpml_sigma/cpml_kappa + cpml_alpha)*dt);
            ce_l[k] = cpml_sigma*(be_l[k] - 1.0)/(temp*kedl[k]);
            
        }
        else if(lpos >(grid_lsize - cpml_depth))
        {
            
            // get the depth into the pml layer
            depth = lpos - (grid_lsize - cpml_depth);
            
            cpml_sigma = cpml_sigmamax*pow((depth/cpml_depth), cpml_m);
            cpml_kappa = 1.0 + (cpml_kappamax - 1.0)*pow((depth/cpml_depth), cpml_m);
            cpml_alpha = cpml_alphamax*pow(((cpml_depth-depth)/cpml_depth), cpml_ma);
            
            kedl[k] = cpml_kappa*grid_dl;
            temp = cpml_sigma + cpml_kappa*cpml_alpha;
            
            be_l[k] = exp(0.0 - (cpml_sigma/cpml_kappa + cpml_alpha)*dt);
            ce_l[k] = cpml_sigma*(be_l[k] - 1.0)/(temp*kedl[k]);
        }
        else
        {
            cpml_kappa = 1.0;
            cpml_sigma = 0.0;
            cpml_alpha = 0.0;
            
            kedl[k] = cpml_kappa*grid_dl;
            be_l[k] = 0.0;
            ce_l[k] = 0.0;
        }
        
        //---------------------------------------------------------
        // H-FIELD
        //---------------------------------------------------------
        //lpos = grid_dl*(static_cast<double>(k)) + grid_dl/2.0;
        lpos = grid_dl*(static_cast<double>(k)) ;//PSTD
        
        if(lpos < cpml_depth)
        {
            // get the depth into the pml layer
            depth = cpml_depth - lpos;
            
            cpml_sigma = cpml_sigmamax*pow((depth/cpml_depth), cpml_m);
            cpml_kappa = 1.0 + (cpml_kappamax - 1.0)*pow((depth/cpml_depth), cpml_m);
            cpml_alpha = cpml_alphamax*pow(((cpml_depth-depth)/cpml_depth), cpml_ma);
            
            khdl[k] = cpml_kappa*grid_dl;
            temp = cpml_sigma + cpml_kappa*cpml_alpha;
            
            bh_l[k] = exp(0.0 - (cpml_sigma/cpml_kappa + cpml_alpha)*dt);
            ch_l[k] = cpml_sigma*(bh_l[k] - 1.0)/(temp*khdl[k]);
        }
        else if(lpos > (grid_lsize - cpml_depth))
        {
            // get the depth into the pml layer
            depth = lpos - (grid_lsize - cpml_depth);
            
            cpml_sigma = cpml_sigmamax*pow((depth/cpml_depth), cpml_m);
            cpml_kappa = 1.0 + (cpml_kappamax - 1.0)*pow((depth/cpml_depth), cpml_m);
            cpml_alpha = cpml_alphamax*pow(((cpml_depth-depth)/cpml_depth), cpml_ma);
            
            khdl[k] = cpml_kappa*grid_dl;
            temp = cpml_sigma + cpml_kappa*cpml_alpha;
            
            bh_l[k] = exp(0.0 - (cpml_sigma/cpml_kappa + cpml_alpha)*dt);
            ch_l[k] = cpml_sigma*(bh_l[k] - 1.0)/(temp*khdl[k]);
        }
        else
        {
            cpml_kappa = 1.0;
            cpml_sigma = 0.0;
            cpml_alpha = 0.0;
            
            khdl[k] = cpml_kappa*grid_dl;
            bh_l[k] = 0.0;
            ch_l[k] = 0.0;
        }
    }
    return;
}



void UpdateE()
{
    
    int i, j;
    double E_inc;
/*
    for (i = 1; i < grid_nx; i++ )
    {
        for (j = 1; j < grid_ny; j++)
        {

            ez[i][j] += dt*((hy[i][j] - hy[i-1][j])/khdx[i] -(hx[i][j] - hx[i][j-1])/khdy[j] + q_h_xy[i][j] - q_h_yx[i][j] + (1-epsilon[i][j])*J_inc(i,j,qtime))/epsilon[i][j];
            
            //cout<<(1-epsilon[i][j])*J_inc(i,j,qtime)/epsilon[i][j]<<endl;
        }
    }
    

    for (j = 0; j < grid_ny - 1; j++) {
        for (i = 0; i < grid_nz - 1; i++) {
            
            q_e_xz[i][j] = bh_x[i]*q_e_xz[i][j] + ch_x[i] * (ez[i + 1][j] - ez[i][j]);
            q_e_yz[i][j] = bh_y[j]*q_e_yz[i][j] + ch_y[j] * (ez[i][j + 1] - ez[i][j]);
            
            
        }
    }*/
    
    for (i = 0; i < grid_nx; i++ )
    {
        for (j = 0; j < grid_ny; j++)
        {
            //E_inc = J_inc(i, j, qtime);

            //ez[i][j] += dt*(hyx[i][j]/khdx[i] - hxy[i][j]/khdy[j] + q_h_xy[i][j] - q_h_yx[i][j] + (1-epsilon[i][j])*E_inc)/epsilon[i][j];
            ez[i][j] += dt*(hyx[i][j]/khdx[i] - hxy[i][j]/khdy[j] + q_h_xy[i][j] - q_h_yx[i][j]);
        }
    }
    
    i = 9;
    for (j = 0; j < grid_ny; j++)
    {
        E_inc = J_inc(i, j, qtime);
        ez[i][j] += E_inc;
    }
    CurlE();
    for (i = 0; i < grid_nx; i++) {
        for (j = 0; j < grid_ny; j++) {
            q_e_xz[i][j] = bh_x[i]*q_e_xz[i][j] + ch_x[i] * ezx[i][j];
            q_e_yz[i][j] = bh_y[j]*q_e_yz[i][j] + ch_y[j] * ezy[i][j];
        }
    }
    //abs(ez[i][j] - J_inc(i, j, qtime));
    printf("%11.10e\t%11.10e\t%11.10e\n",ez[grid_nx - 12][grid_ny/2],J_inc(grid_nx-12, grid_ny/2, qtime),abs(ez[grid_nx - 12][grid_ny/2] - J_inc(grid_nx-12, grid_ny/2, qtime)));
    


}



void UpdateH()
{
    
    
    int i, j;
    for (i = 0; i < grid_nx; i++ )
    {
        for (j = 0; j < grid_ny; j++)
        {
            
            hx[i][j] += dt*(- ezy[i][j]/kedy[j] - q_e_yz[i][j]);
            hy[i][j] += dt*(ezx[i][j]/kedx[i] + q_e_xz[i][j]);

        }
    }
    
    
    
    
    CurlH();
    for (i = 0; i < grid_nx; i++) {
        for (j = 0; j < grid_ny; j++) {
            q_h_yx[i][j] = be_y[j]*q_h_yx[i][j] + ce_y[j] * hxy[i][j];
            q_h_xy[i][j] = be_x[i]*q_h_xy[i][j] + ce_x[i] * hyx[i][j];
        }
    }
    
    
/*
    int i, j;
    for (i = 1; i < grid_nx - 1; i++ )
    {
        for (j = 1; j < grid_ny - 1; j++)
        {
            
            hx[i][j] += dt*(- (ez[i][j+1] - ez[i][j])/kedy[j] - q_e_yz[i][j]);
            hy[i][j] += dt*((ez[i+1][j] - ez[i][j])/kedx[i] + q_e_xz[i][j]);
        }
    }
    
    for (j = 1; j < grid_ny; j++) {
        for (i = 1; i < grid_nx; i++) {
            
            q_h_yx[i][j] = be_y[j]*q_h_yx[i][j] + ce_y[j] * (hx[i][j] - hx[i][j - 1]);
            q_h_xy[i][j] = be_x[i]*q_h_xy[i][j] + ce_x[i] * (hy[i][j] - hy[i - 1][j]);
            
            
        }
    }*/
}

void FastChebyshevTransform(double *f, int N, int s)
{
    
    for (int j = 0; j < N; j++) {
        f[j] = f[j]/(N-1);   //Fast Cosine Transform with fftw
    }
    if (s == -1) {
        f[0] = 2*f[0];
        f[N-1] = 2*f[N-1]/2;
    }
    fftw_plan p;
    p = fftw_plan_r2r_1d(N, f, f, FFTW_REDFT00, FFTW_ESTIMATE);
    fftw_execute(p);
    if (s == 1) {
        f[0] = f[0]/2.0;
        f[N-1] = f[N-1]/2.0;
    }
    fftw_destroy_plan(p);
}

void FastChebyshevDerivative(double* f, int N)
{
    
    FastChebyshevTransform(f, N,1);
    ChebyshevDerivativeCoefficients(f, N);
    FastChebyshevTransform(f, N,-1);
    
}

void ChebyshevDerivativeCoefficients(double* f, int N)
{
    double* g;
    g = new double[N];
    g[N-1] = 0.0;
    g[N-2] = (2*(N-1))*f[N-1];
    for (int k = N - 3; k >=1 ; k--) {
        g[k] = 2*(k+1)*f[k+1] + g[k+2];
    }
    g[0] = f[1] + g[2]/2;
    for (int i = 0; i < N; i++) {
        f[i] = -g[i];
    }
    delete [] g;
    
}

double* ChebyshevGaussLobattoNodes(int N)
{
    double* x;
    x = new double[N];
    for (int j = 0; j < N ; j++) {
        x[j] = -cos(j*M_PI/(N-1));
    }
    return x;
}


void PerformPSTD()
{
    FILE * outfile = fopen( "Edata.dat", "w" );
    fprintf(outfile,"#eFields\n");
    
    xdiff = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid_nx);
    ydiff = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid_ny);

    for(qtime = 0; qtime < maxTime; qtime++)
    {
        UpdateE();
        UpdateH();
        for (int i = 0; i < grid_nx; i++) {
            for (int j = 0; j < grid_ny; j++) {
                fprintf(outfile, "%11.10e\t", ez[i][j]);
            }
            fprintf(outfile, "\n");
        }
        
    }
    
    fclose(outfile);
    fftw_free(xdiff);
    fftw_free(ydiff);
    
}

void PseudoSpecDiff(fftw_complex *data, int N)
{
    double temp, imomentum;
    fftw_plan pb;
    fftw_plan pf;
    pb = fftw_plan_dft_1d(N, data, data, FFTW_BACKWARD, FFTW_ESTIMATE);
    pf = fftw_plan_dft_1d(N, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
    
    
    fftw_execute(pb);
    for (int  k = 0 ; k < N/2 ; k++) {
        imomentum = 2*M_PI* ((double)k)/((double)N);
        temp = data[k][0];
        data[k][0] = imomentum * data[k][1];
        data[k][1] = -imomentum * temp;
    }
    data[N/2][0] = 0.0;
    data[N/2][1] = 0.0;
    for (int  k =  N/2 + 1; k < N; k++) {
        imomentum = 2*M_PI* ((double)(k-N))/((double)N);
        temp = data[k][0];
        data[k][0] = imomentum*data[k][1];
        data[k][1] = -imomentum*temp;
    }
    fftw_execute(pf);
    
    for (int  j = 0;  j < N; j++) {
        data[j][0] = data[j][0]/((double)N);
        data[j][1] = data[j][1]/((double)N);
    }
    
    fftw_destroy_plan(pb);
    fftw_destroy_plan(pf);
    return;
}

void CurlH()
{

    for (int j = 0; j < grid_ny; j++) {
        for (int i = 0; i < grid_nx; i++) {
            xdiff[i][0] = hy[i][j];
            xdiff[i][1] = 0.0;
        }
        PseudoSpecDiff(xdiff, grid_nx);
        for (int i = 0; i < grid_nx; i++) {
            hyx[i][j] = xdiff[i][0];
        }
    }
    
    for (int i = 0; i < grid_nx; i++) {
        
        for (int j = 0; j < grid_ny; j++) {
            ydiff[j][0] = hx[i][j];
            ydiff[j][1] = 0.0;
            
        }
        PseudoSpecDiff(ydiff, grid_ny);
        for (int j = 0; j < grid_ny; j++) {
            hxy[i][j] = ydiff[j][0];
        }
        
    }

/*  double *hydx, *hxdy;
    hydx = new double [grid_nx];
    hxdy = new double [grid_ny];
    for (int j = 0; j < grid_ny; j++) {
        for (int i = 0; i < grid_nx; i++) {
            xdiff[i][0] = hy[i][j];
            xdiff[i][1] = 0.0;
            //hydx[i] = hy[i][j];
        }
        PseudoSpecDiff(xdiff, grid_nx);
        //FastChebyshevDerivative(hydx,grid_nx);
        for (int i = 0; i < grid_nx; i++) {
            hyx[i][j] = xdiff[i][0];
            //hyx[i][j] = hydx[i];
        }
    }
    
    for (int i = 0; i < grid_nx; i++) {

        for (int j = 0; j < grid_ny; j++) {
            //ydiff[j][0] = hx[i][j];
            //ydiff[j][1] = 0.0;
            hxdy[j] = hx[i][j];
            
        }
        //PseudoSpecDiff(ydiff, grid_ny);
        FastChebyshevDerivative(hxdy,grid_ny);
        for (int j = 0; j < grid_ny; j++) {
            //hxy[i][j] = ydiff[j][0];
            hxy[i][j] = hxdy[j];
        }
        
    }
    delete [] hxdy;
    delete [] hydx;*/
    return;
}


void CurlE()
{
    
    for (int j = 0; j < grid_ny; j++) {
        
        for (int i = 0; i < grid_nx; i++) {
            xdiff[i][0] = 0.0;
            xdiff[i][1] = ez[i][j];
        }
        PseudoSpecDiff(xdiff, grid_nx);
        for (int i = 0; i < grid_nx; i++) {
            ezx[i][j] = xdiff[i][1];
        }
    }
    for (int i = 0; i < grid_nx; i++) {

        for (int j = 0; j < grid_ny; j++) {
            ydiff[j][0] = 0.0;
            ydiff[j][1] = ez[i][j];
            
        }
        PseudoSpecDiff(ydiff, grid_ny);
        for (int j = 0; j < grid_ny; j++) {
            ezy[i][j] = ydiff[j][1];
        }
    }
    
/*    double *ezdx, *ezdy;
    
    ezdx = new double [grid_nx];
    ezdy = new double [grid_ny];
    
    for (int j = 0; j < grid_ny; j++) {
        for (int i = 0; i < grid_nx; i++) {
            //xdiff[i][0] = hy[i][j];
            //xdiff[i][1] = 0.0;
            ezdx[i] = ez[i][j];
        }
        PseudoSpecDiff(xdiff, grid_nx);
        //FastChebyshevDerivative(ezdx,grid_nx);
        for (int i = 0; i < grid_nx; i++) {
            hyx[i][j] = xdiff[i][0];
            //ezx[i][j] = ezdx[i];
        }
    }
    
    for (int i = 0; i < grid_nx; i++) {
        
        for (int j = 0; j < grid_ny; j++) {
            ydiff[j][0] = hx[i][j];
            ydiff[j][1] = 0.0;
            //ezdy[j] = ez[i][j];
            
        }
        PseudoSpecDiff(ydiff, grid_ny);
        //FastChebyshevDerivative(ezdy,grid_ny);
        for (int j = 0; j < grid_ny; j++) {
            hxy[i][j] = ydiff[j][0];
            //ezy[i][j] = ezdy[j];
        }
        
    }
    delete [] ezdy;
    delete [] ezdx;*/
    

    return;
}






void FreeMemory()
{
    
    for (int i = 0; i < grid_nx; i++)
    {

        delete [] ez[i];
        delete [] hx[i];
        delete [] hy[i];
        

        delete [] q_e_xz[i];
        delete [] q_e_yz[i];
        
        delete [] q_h_xy[i];
        delete [] q_h_yx[i];
        
        delete [] epsilon[i];
        
        
        
    }

    delete [] ez;
    
    delete [] hx;
    delete [] hy;
    

    delete [] q_e_xz;
    delete [] q_e_yz;

    
    delete [] q_h_xy;
    delete [] q_h_yx;

    
    delete [] epsilon;
    
    

    
    delete [] be_x;
    delete [] ce_x;
    delete [] bh_x;
    delete [] ch_x;
    
    delete [] kedx;
    delete [] khdx;
    delete [] be_y;
    delete [] ce_y;
    delete [] bh_y;
    delete [] ch_y;
    
    delete [] kedy;
    delete [] khdy;
    
    
    for (int i = 0; i < grid_nx; i++)
    {

        delete [] ezx[i];
        delete [] ezy[i];
        
        delete [] hxy[i];
        delete [] hyx[i];
        
        
    }
    

    delete [] ezx;
    delete [] ezy;
    
    delete [] hxy;
    delete [] hyx;
    
}



