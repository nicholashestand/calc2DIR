#include <string.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <algorithm>
#include "calc2DIR.h"

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
    energy          = new double[nchrom]();
    energy_last     = new double[nchrom]();
    dipole          = new vec3[nchrom]();
    dipole0         = new vec3[nchrom]();
    R1D             = new complex<double>[R1D_npoints]();
    egt             = new complex<double>[nchrom]();

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
    delete [] energy_last;
    delete [] dipole;
    delete [] dipole0;
    delete [] R1D;
    delete [] egt;

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

        // save parameters as class variable
        if      ( para.compare("energy_file") == 0 )  efn          = value;
        else if ( para.compare("dipole_file") == 0 )  dfn          = value;
        else if ( para.compare("output_file") == 0 )  ofn          = value;
        else if ( para.compare("time_step") == 0 )    dt           = stof(value);
        else if ( para.compare("t1_max") == 0 )       t1_max       = stof(value);
        else if ( para.compare("t3_max") == 0 )       t3_max       = stof(value);
        else if ( para.compare("t2") == 0 )           t2           = stof(value);
        else if ( para.compare("rt2") == 0 )          RT2          = stof(value);
        else if ( para.compare("nsamples") == 0 )     nsamples     = stoi(value);
        else if ( para.compare("sample_every") == 0 ) sample_every = stof(value);
        else if ( para.compare("nchrom") == 0 )       nchrom       = stoi(value);
        else if ( para.compare("window0") == 0 )      window0      = stof(value);
        else if ( para.compare("window1") == 0 )      window1      = stof(value);
        else cerr << "\tWARNING:: parameter " << para << " is not recognized." << endl;
    }
    inpf.close();

    // some output to confirm parameters
    cout << "\tSetting energy_file to: " << efn << endl;
    cout << "\tSetting dipole_file to: " << dfn << endl;
    cout << "\tSetting output_file to: " << ofn << endl;
    cout << "\tSetting time_step to: " << dt << " ps" << endl;
    cout << "\tSetting t1_max to: " << t1_max << " ps" << endl;
    cout << "\tSetting t3_max to: " << t3_max << " ps" << endl;
    cout << "\tSetting t2 to: " << t2 << " ps" << endl;
    cout << "\tSetting relaxation time T2 to: " << RT2 << " ps" << endl;
    cout << "\tSetting nsamples to: " << nsamples << endl;
    cout << "\tSetting sample_every to: " << sample_every << " ps" << endl;
    cout << "\tSetting nchrom to: " << nchrom << endl;
    cout << "\tSetting spectral limits to: " << window0 << " - " << window1 << " cm." << endl;
    cout << ">>> Done reading simulation parameters from " << _inpf_ << endl;

    // set number of points in response functions based on t1_max and dt
    R1D_npoints = static_cast<int>(t1_max/dt);

    // shift energies by mean of spectral window to avoid high frequency oscillations
    shift = (window1 + window0)/2.;
    cout << ">>> Shifting frequencies by: " << shift << " cm." << endl;


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

    // only keep the energies, ignore the couplings for now
    int col=0;
    for ( int i = 0; i < nchrom; i ++ ){
        energy[i] = energyTmp[col] - shift;
        col += nchrom-i;
    }

    return IR2DOK;
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

    return IR2DOK;
}

complex<double> IR2D::getR1D( int t1 )
// return the linear response function at a given t1 for a given chromophore
// See Eq 7.10 from Hamm and Zanni
{
    double mu[nchrom];
    complex<double> arg;
    complex<double> R1Dtmp;

    // reset all integrals to zero
    if ( t1 == 0 ){
        //   initialize egt to one
        for ( int chrom = 0; chrom < nchrom; chrom ++ ) egt[chrom] = complex_one;
    }
    else {
        for ( int chrom = 0; chrom < nchrom; chrom ++ ){
            // integrate the equation for each chromophore using the trapezoid rule
            arg = -img*dt*(energy_last[chrom]+energy[chrom])/(2.*HBAR);
            egt[chrom] *= exp(arg);
        }
    }
    
    // get dipole part
    for ( int chrom = 0; chrom < nchrom; chrom ++ ){
        mu[chrom] = dot3(dipole0[chrom],dipole[chrom]);
    }

    // return the response function
    // average the response function over all of the chromophores then return
    R1Dtmp = complex_zero;
    for ( int chrom = 0; chrom < nchrom; chrom ++ ){
        R1Dtmp += img*mu[chrom]*egt[chrom]/(1.*nchrom);
    }

    return R1Dtmp;
}

double IR2D::dot3( vec3 x, vec3 y )
// dot product of 3 vector
{
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2]; 
}

int IR2D::write1D()
// write response function and FTIR spectrum
{
    string fn;
    ofstream ofile;
    const int length=32768;
    complex<double> fftIn[length];
    complex<double> fftOut[length];
    fftw_plan plan;
    double convert, freq, intens, abs[length], dis[length];

    fn = ofn+"-R1D.dat";
    ofile.open(fn);
    if ( ! ofile.is_open() ) { fileOpenErr( fn ); exit(EXIT_FAILURE);}
    cout << ">>> Writing " << fn << "." << endl;

    ofile << "# time (ps) Real Imag" << endl;
    for ( int i = 0; i < R1D_npoints; i ++ ){
        ofile << i * dt << " " << R1D[i].real() << " " << R1D[i].imag() << endl;
    }
    ofile.close();

    // do fft
    plan = fftw_plan_dft_1d( length, reinterpret_cast<fftw_complex*>(fftIn), \
                                     reinterpret_cast<fftw_complex*>(fftOut),\
                                     FFTW_BACKWARD, FFTW_ESTIMATE );
    
    // See Eq 4.8 from Hamm and Zanni -- Absorptive part
    for ( int i = 0; i < length ; i ++ ) fftIn[i] = complex_zero;
    for ( int i = 0; i < R1D_npoints; i ++ ){
        fftIn[i]        = img*R1D[i]*exp(-i*dt/RT2);
        fftIn[length-i] = conj(img*R1D[i])*exp(-i*dt/RT2);
    }
    fftw_execute(plan);
    for ( int i = 0; i < length; i ++ ) abs[i] = fftOut[i].real();

    // See Eq 4.11 from Hamm and Zanni -- Dispersive part
    for ( int i = 0; i < length ; i ++ ) fftIn[i] = complex_zero;
    for ( int i = 0; i < R1D_npoints; i ++ ){
        fftIn[i]        = -R1D[i]*exp(-i*dt/RT2);
        fftIn[length-i] = conj(-R1D[i])*exp(-i*dt/RT2);
    }
    fftw_execute(plan);
    for ( int i = 0; i < length; i ++ ) dis[i] = fftOut[i].real();

    fn = ofn+"-FTIR.dat";
    ofile.open(fn);
    if ( ! ofile.is_open() ) { fileOpenErr( fn ); exit(EXIT_FAILURE);}
    cout << ">>> Writing " << fn << "." << endl;
    ofile << "# frequency (cm-1) Absorption Dispersion" << endl;


    intens = -dt*length/(2.*PI*HBAR*sqrt(length));
    // scale intensity to keep total area of the spectrum conserved and normalize
    // by root(length) since FFTW3 does not
    // multiply the signal by -1 so it looks like absorption and not signal intensity
        
    // negative frequencies are stored at the end of the fftOut array
    for ( int i = length/2; i < length; i ++ )
    {
        // get frequency in wavenumber
        freq   = 2.*PI*HBAR*(i-length)/(dt*length) + shift;
        if ( freq < window0 or freq > window1 ) continue;
        ofile << freq << " " << intens*abs[i] << " " << intens*dis[i] << endl;
    }
    for ( int i = 0; i < length/2; i ++ ){
        // get frequency in wavenumber
        freq   = 2.*PI*HBAR*i/(dt*length) + shift;
        if ( freq < window0 or freq > window1 ) continue;
        ofile << freq << " " << intens*abs[i] << " " << intens*dis[i] << endl;
    }
    ofile.close();

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
    for ( int sample = 0; sample < spectrum.nsamples; sample ++ ){
        int frame0 = static_cast<int>(spectrum.sample_every/spectrum.dt);
        // read and save dipole at t0
        spectrum.readDframe(frame0);
        memcpy( spectrum.dipole0, spectrum.dipole, spectrum.nchrom*3*sizeof(double) );
        
        // 1D spectrum -- t1 is in units of dt
        for ( int t1 = 0; t1 < spectrum.R1D_npoints; t1 ++ )
        {
            // get frame number for current t1
            int frame = frame0 + t1;

            // read in energy and dipole
            spectrum.readEframe(frame);
            spectrum.readDframe(frame);
            
            spectrum.R1D[t1] += spectrum.getR1D( t1 );

            // save energy for next
            memcpy( spectrum.energy_last, spectrum.energy, spectrum.nchrom*sizeof(double) );
        }
    }

    // normalize r1d by number of samples
    for ( int t1 = 0; t1 < spectrum.R1D_npoints; t1 ++ ){
        spectrum.R1D[t1]/=(1.*spectrum.nsamples);
    }

    // write the response functions and the spectrum
    spectrum.write1D();
}
