#include <string.h>
#include <complex.h>
#include <fstream>

typedef float vec3[3];
using namespace std;

#ifndef calc2DIR_H
#define calc2DIR_H
#define IR2DOK 0
class IR2D
{
    public:
        // variables
        string efn="e.dat";             // name for energy file
        string dfn="d.dat";             // name for dipole file
        float  dt = 0.010 ;             // time step in ps
        float  t1_max=10. ;             // t1 max in ps
        float  t2    =0.  ;             // t2 in ps
        float  t3_max=10. ;             // t3 max in ps
        int    nsamples=1 ;             // number of samples
        float  sample_every=10.;        // how often to take a new sample in ps
        int    nframes=0;               // number of frames in energy and dipole files to consider
        int    nchrom=1;                // number of uncoupled chromophores to consider

        float *energy;                  // energy
        vec3  *dipole;                  // dipole vector

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

};
#endif
