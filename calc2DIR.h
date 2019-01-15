#include <string.h>
#include <complex.h>
#include <fstream>

typedef double vec3[3];
using namespace std;

#ifndef calc2DIR_H
#define calc2DIR_H
#define IR2DOK 0
#define PI M_PI
#define HBAR 5.308837367 // in cm-1*ps
class IR2D
{
    public:
        // variables
        string _efile_="e.dat";         // name for energy file
        string _dfile_="d.dat";         // name for dipole file
        string _ofile_="spectrum";      // name for output file
        double dt = 0.010 ;             // time step in ps
        double t1_max=10. ;             // t1 max in ps
        double t2    =0.  ;             // t2 in ps
        double t3_max=10. ;             // t3 max in ps
        double lifetime_T1=0.2;         // lifetime T2 in ps
        double lifetime_T2;             // lifetime T1 in ps
        double anharm=14.0;             // anharmonicity in cm-1
        int    nsamples=1 ;             // number of samples
        double sample_every=10.;        // how often to take a new sample in ps
        int    nchrom=1;                // number of uncoupled chromophores to consider
        int    t1_npoints;              // number of data points for t1
        int    t3_npoints;              // number of data points for t3
        double window0 = 1400;          // lower limit of spectral window in cm-1
        double window1 = 1700;          // upper limit if spectral window in cm-1
        double shift;                   // reference frequency in cm-1

        // arrays to hold energy and dipole
        double *energy_t1;              // energy at frame t1
        double *energy_t1_last;         // energy at frame t1 - 1
        double *energy_t3;              // energy at frame t3
        double *energy_t3_last;         // energy at frame t3 - 1
        vec3   *dipole_t0;              // dipole vector at frame t0
        vec3   *dipole_t1;              // dipole vector at frame t1
        vec3   *dipole_t2;              // dipole vector at frame t1 + t2
        vec3   *dipole_t3;              // dipole vector at frame t1 + t2 + t3

        // constants
        const complex<double> img          = {0.,1.};   
        const complex<double> complex_one  = {1.,0.};   
        const complex<double> complex_zero = {0.,0.};
      
        // response functions
        complex<double> *R1D;           // Linear response function
        complex<double> *R2D_R1;        // third order response function R1=R2
        complex<double> *R2D_R3;        // third order response function R3
        complex<double> *R2D_R4;        // third order response function R4=R5
        complex<double> *R2D_R6;        // third order response function R6
        complex<double> *eint_t1;       // integral from 0 to t1 of dw(tau)_01
        complex<double> *eint_t3;       // integral from t1+t2 to t1+t2+t3 of dw(tau)_01

        // file handles
        ifstream efile, dfile;

        // default constructor and destructor
        IR2D( string _inpf_ );
        ~IR2D();

        // functions
        void fileOpenErr( string _fn_ );
        int readParam( string _inpf_ );
        template<class T> void tellParam( string param, T value );
        int readEframe( int frame, string time );
        int readDframe( int frame, string time );
        int get_eint( int t1, string which );
        complex<double> getR1D( int t1 );
        complex<double> getR2D( int t1, int t2, string which );
        int writeR1D();
        int writeR2D();
        int write1Dfft();
        int write2DRabs();
        int write2Dout( complex<double> *data, string fn, string which, int length);
        double dot3( vec3 x, vec3 y );
};
#endif
