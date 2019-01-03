#include <string.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <fftw3.h>
#include <math.h>
#include <algorithm>
#include "calc2DIR.h"

#define PI M_PI

using namespace std;

IR2D::IR2D( string _inpf_ )
// Default Constructor
{
    int err, lengthTmp;
    string line;

    // read input file
    if ( err = readParam( _inpf_ ) != IR2DOK ) { 
        cerr << "Warning:: readParam returned with error " << err << "." << endl; 
        exit( EXIT_FAILURE );
    }

    // allocate variables for energy and dipole
    energy = new float[nchrom]();
    dipole = new vec3[nchrom]();

    // open the energy and dipole files
    efile.open( efn, ios::binary );    
    if ( ! efile.is_open() ) { fileOpenErr( efn ); exit(EXIT_FAILURE);}
    dfile.open( dfn, ios::binary );    
    if ( ! dfile.is_open() ) { fileOpenErr( dfn ); exit(EXIT_FAILURE);}

}

IR2D::~IR2D()
// Default Destructor
{
    // destroy arrays
    delete [] energy;
    delete [] dipole;

    // close files
    efile.close();
    dfile.close();
}

int IR2D::readParam( string _inpf_ )
// read the input file to get parameters
{
    string line, para, value;

    ifstream inpf( _inpf_.c_str() );
    if ( ! inpf.is_open() ) { fileOpenErr( _inpf_ ); return 1; }
    cout << ">>> Reading simulation parameters from " << _inpf_ << endl;

    // parse input file
    while ( getline( inpf, line ) )
    {
        // create string stream and parse parameter name and value
        istringstream line2( line );
        line2 >> para >> value;

        // skip comments
        if ( para[0] == '#' ) continue;

        // transform parameter name to lower case
        transform( para.begin(), para.end(), para.begin(), ::tolower);
        //cout << para << endl;

        // save parameters as class variable
        if      ( para.compare("energy_file") == 0 )  efn          = value;
        else if ( para.compare("dipole_file") == 0 )  dfn          = value;
        else if ( para.compare("time_step") == 0 )    dt           = stof(value);
        else if ( para.compare("t1_max") == 0 )       t1_max       = stof(value);
        else if ( para.compare("t3_max") == 0 )       t3_max       = stof(value);
        else if ( para.compare("t2") == 0 )           t2           = stof(value);
        else if ( para.compare("nsamples") == 0 )     nsamples     = stoi(value);
        else if ( para.compare("sample_every") == 0 ) sample_every = stof(value);
        else if ( para.compare("nframes") == 0 )      nframes      = stoi(value);
        else if ( para.compare("nchrom") == 0 )       nchrom       = stoi(value);
        else cerr << "\tWARNING:: parameter " << para << " is not recognized." << endl;
    }
    inpf.close();

    // some output to confirm parameters
    cout << "\tSetting energy_file to: " << efn << endl;
    cout << "\tSetting dipole_file to: " << dfn << endl;
    cout << "\tSetting time_step to: " << dt << " ps" << endl;
    cout << "\tSetting t1_max to: " << t1_max << " ps" << endl;
    cout << "\tSetting t3_max to: " << t3_max << " ps" << endl;
    cout << "\tSetting t2 to: " << t2 << " ps" << endl;
    cout << "\tSetting nsamples to: " << nsamples << endl;
    cout << "\tSetting sample_every to: " << sample_every << " ps" << endl;
    cout << "\tSetting nframes to: " << nframes << endl;
    cout << "\tSetting nchrom to: " << nchrom << endl;
    cout << ">>> Done reading simulation parameters from " << _inpf_ << endl;

    return IR2DOK;
}

void IR2D::fileOpenErr( string _fn_ )
// give an error message when I cant open a file
{
    cerr << "ERROR:: Could not open " << _fn_ << "." << endl;
}

int IR2D::readEframe( int frame )
// Read the energy file
{
    int frameTmp; 
    float energyTmp[ nchrom*(nchrom+1)/2 ]; 
    int64_t file_offset;
                                            
    file_offset = frame*(sizeof(int)+sizeof(float)*nchrom*(nchrom+1)/2);
    efile.seekg( file_offset );

    efile.read( (char*)&frameTmp, sizeof(int) );
    efile.read( (char*)energyTmp, sizeof(float)*nchrom*(nchrom+1)/2 );
    return 1;

    // only keep the energies, ignore the couplings for now
    int col=0;
    for ( int i = 0; i < nchrom; i ++ ){
        energy[nchrom + i] = energyTmp[col];
        col += nchrom-i;
    }
}

int IR2D::readDframe( int frame )
// Read the dipole file
{
    int frameTmp;
    float dipoleTmp[ nchrom*3 ];
    int64_t file_offset;

    file_offset = frame*(sizeof(int)+sizeof(float)*nchrom*3);
    dfile.seekg( file_offset );

    dfile.read( (char*)&frameTmp, sizeof(int) );
    dfile.read( (char*)dipoleTmp, sizeof(float)*nchrom*3 );

    // put these into the dipole vector variable
    // note that in the bin file all x's come first, then y's, etc
    for ( int i = 0; i < 3; i ++ ){
        for ( int chrom = 0; chrom < nchrom; chrom ++ ){
            dipole[chrom][i] = dipoleTmp[i*nchrom + chrom];
        }
    }
}

int main( int argc, char* argv[] )
{

    if ( argc != 2 ){
        cout << "ERROR:: Program expects the name of the input" << 
                "file as the only argument. Aborting." << endl;
        exit( EXIT_FAILURE );
    }

    // get input file name
    string inpf( argv[1] );

    // Initialize the IR2D class
    IR2D spectrum( inpf );

    // Loop over the trajectory
    for ( int frame = 0; frame < spectrum.nframes; frame ++ )
    {
        // read in energy and dipole
        spectrum.readEframe(frame);
        spectrum.readDframe(frame);
    }

}
