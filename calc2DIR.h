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
        string efn="e.dat";             // name for energy file
        string dfn="d.dat";             // name for dipole file
        string ofn="spectrum";          // name for output file
        double dt = 0.010 ;             // time step in ps
        double t1_max=10. ;             // t1 max in ps
        double t2    =0.  ;             // t2 in ps
        double t3_max=10. ;             // t3 max in ps
        double RT2=0.3;                 // T2 relaxation time
        int    nsamples=1 ;             // number of samples
        double sample_every=10.;        // how often to take a new sample in ps
        int    nchrom=1;                // number of uncoupled chromophores to consider
        int    R1D_npoints;             // number of data points in 1D response function
        double window0 = 1400;          // lower limit of spectral window in cm-1
        double window1 = 1700;          // upper limit if spectral window in cm-1
        double shift;                   // reference frequency in cm-1

        // arrays to hold energy and dipole
        double *energy;                  // energy
        double *energy_last;             // energy from previous frame
        vec3   *dipole;                  // dipole vector
        vec3   *dipole0;                 // dipole vector at t0

        // constants
        const complex<double> img          = {0.,1.};   
        const complex<double> complex_one  = {1.,0.};   
        const complex<double> complex_zero = {0.,0.};
      
        // response functions
        complex<double> *R1D;            // Linear response function
        complex<double> *egt;            // to keep track of the lineshape g(t)

        // file handles
        ifstream efile, dfile;

        // default constructor and destructor
        IR2D( string _inpf_ );
        ~IR2D();

        // functions
        void fileOpenErr( string _fn_ );
        int readParam( string _inpf_ );
        int readEframe( int frame );
        int readDframe( int frame );
        complex<double> getR1D( int t1 );
        int write1D();

        double dot3( vec3 x, vec3 y );
};
#endif
