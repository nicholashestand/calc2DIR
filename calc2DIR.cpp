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

    // allocate variable arrays
    energy_t1       = new double[nchrom]();
    energy_t1_last  = new double[nchrom]();
    energy_t3       = new double[nchrom]();
    energy_t3_last  = new double[nchrom]();
    dipole_t0       = new vec3[nchrom]();
    dipole_t1       = new vec3[nchrom]();
    dipole_t2       = new vec3[nchrom]();
    dipole_t3       = new vec3[nchrom]();
    R1D             = new complex<double>[t1_npoints]();
    R2D_R1          = new complex<double>[t1_npoints*t3_npoints]();
    R2D_R3          = new complex<double>[t1_npoints*t3_npoints]();
    R2D_R4          = new complex<double>[t1_npoints*t3_npoints]();
    R2D_R6          = new complex<double>[t1_npoints*t3_npoints]();
    eint_t1         = new complex<double>[nchrom]();
    eint_t3         = new complex<double>[nchrom]();

    // open the energy and dipole files
    efile.open( _efile_, ios::binary );    
    if ( ! efile.is_open() ) { fileOpenErr( _efile_ ); exit(EXIT_FAILURE);}
    dfile.open( _dfile_, ios::binary );    
    if ( ! dfile.is_open() ) { fileOpenErr( _dfile_ ); exit(EXIT_FAILURE);}

}

IR2D::~IR2D()
// Default Destructor
{
    // destroy arrays
    delete [] energy_t1;
    delete [] energy_t1_last;
    delete [] energy_t3;
    delete [] energy_t3_last;
    delete [] dipole_t0;
    delete [] dipole_t1;
    delete [] dipole_t2;
    delete [] dipole_t3;
    delete [] R1D;
    delete [] R2D_R1;
    delete [] R2D_R3;
    delete [] R2D_R4;
    delete [] R2D_R6;
    delete [] eint_t1;
    delete [] eint_t3;

    // close files
    efile.close();
    dfile.close();
}

int IR2D::readParam( string _inpf_ )
// read the input file to get parameters
{
    string line, para, value;

    // open file for reading
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
        if      ( para.compare("energy_file") == 0 )  _efile_      = value;
        else if ( para.compare("dipole_file") == 0 )  _dfile_      = value;
        else if ( para.compare("output_file") == 0 )  _ofile_      = value;
        else if ( para.compare("time_step") == 0 )    dt           = stof(value);
        else if ( para.compare("t1_max") == 0 )       t1_max       = stof(value);
        else if ( para.compare("t3_max") == 0 )       t3_max       = stof(value);
        else if ( para.compare("t2") == 0 )           t2           = stof(value);
        else if ( para.compare("rt2") == 0 )          RT2          = stof(value);
        else if ( para.compare("rt1") == 0 )          RT1          = stof(value);
        else if ( para.compare("anharm") == 0 )       anharm       = stof(value);
        else if ( para.compare("nsamples") == 0 )     nsamples     = stoi(value);
        else if ( para.compare("sample_every") == 0 ) sample_every = stof(value);
        else if ( para.compare("nchrom") == 0 )       nchrom       = stoi(value);
        else if ( para.compare("window0") == 0 )      window0      = stof(value);
        else if ( para.compare("window1") == 0 )      window1      = stof(value);
        else cerr << "\tWARNING:: parameter " << para << " is not recognized." << endl;
    }
    inpf.close();

    // some output to confirm parameters
    tellParam<string>( "energy_file", _efile_ );
    tellParam<string>( "dipole_file", _dfile_ );
    tellParam<string>( "output_file", _ofile_ );
    tellParam<double>( "time_step", dt );
    tellParam<double>( "t1_max", t1_max );
    tellParam<double>( "t3_max", t3_max );
    tellParam<double>( "t2", t2 );
    tellParam<double>( "T2", RT2 );
    tellParam<double>( "T1", RT1 );
    tellParam<double>( "anharm", anharm );
    tellParam<int>( "nsamples", nsamples );
    tellParam<int>( "sample_every", sample_every );
    tellParam<int>( "nchrom", nchrom );
    tellParam<double>( "window0", window0 );
    tellParam<double>( "window1", window1 );
    cout << ">>> Done reading simulation parameters from " << _inpf_ << endl;


    // set number of points in response functions based on t1_max and dt
    t1_npoints = static_cast<int>(t1_max/dt);
    t3_npoints = static_cast<int>(t3_max/dt);

    // shift energies by mean of spectral window to avoid high frequency oscillations
    shift = (window1 + window0)/2.;
    cout << ">>> Shifting frequencies by: " << shift << " cm." << endl;


    return IR2DOK;
}

template<class T>
void IR2D::tellParam( string param, T value )
// write parameter and value to output
{
    cout << "\tSetting " << param << " to " << value << "." << endl; 
}

void IR2D::fileOpenErr( string _fn_ )
// give an error message when I cant open a file
{
    cerr << "ERROR:: Could not open " << _fn_ << "." << endl;
}

int IR2D::readEframe( int frame, string t )
// Read the energy file
{
    int frameTmp; 
    float energyTmp[ nchrom*(nchrom+1)/2 ]; 
    int64_t file_offset;

    file_offset = frame*(sizeof(int)+sizeof(float)*nchrom*(nchrom+1)/2);
    efile.seekg( file_offset );

    // for files to use with NISE, the hamiltonian is stored in upper tridiagonal format
    // on each line with an integer at the beginning of the line
    efile.read( (char*)&frameTmp, sizeof(int) );
    efile.read( (char*)energyTmp, sizeof(float)*nchrom*(nchrom+1)/2 );

    // only keep the energies, ignore the couplings for now
    int col=0;
    for ( int i = 0; i < nchrom; i ++ ){
        // set to variable depending on t value
        if ( t.compare("t1") == 0 )      energy_t1[i] = energyTmp[col] - shift;
        else if ( t.compare("t3") == 0 ) energy_t3[i] = energyTmp[col] - shift;
        else return 1;
        col += nchrom-i;
    }

    return IR2DOK;
}

int IR2D::readDframe( int frame, string t )
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
            if ( t.compare("t0") == 0 )      dipole_t0[chrom][i] = dipoleTmp[i*nchrom + chrom];
            else if ( t.compare("t1") == 0 ) dipole_t1[chrom][i] = dipoleTmp[i*nchrom + chrom];
            else if ( t.compare("t2") == 0 ) dipole_t2[chrom][i] = dipoleTmp[i*nchrom + chrom];
            else if ( t.compare("t3") == 0 ) dipole_t3[chrom][i] = dipoleTmp[i*nchrom + chrom];
            else return 1;
        }
    }

    return IR2DOK;
}

complex<double> IR2D::getR1D( int t1 )
// return the linear response function at a given t1 for a given chromophore
// See Eq 7.10 from Hamm and Zanni
{
    double mu;
    complex<double> R1Dtmp;

    R1Dtmp = complex_zero;
    for ( int chrom = 0; chrom < nchrom; chrom ++ ){
        // dipole part -- NOTE: I DO NOT MAKE THE CONDON APPROXIMATION
        // Also note that this is an isotropic averaged spectrum
        mu      = dot3(dipole_t0[chrom],dipole_t1[chrom]);
        // average the response function over all of the chromophores then return
        // note that we have included the T2 relaxation here
        R1Dtmp += img*mu*eint_t1[chrom]*exp(-t1*dt/RT2);
    }

    return R1Dtmp/(1.*nchrom);
}

complex<double> IR2D::getR2D( int t1, int t3, string which )
// return the third order response function at a given t3 for a given chromophore
// See Eq 7.35 from Hamm and Zanni
{
    double mu;
    complex<double> R2Dtmp;

    R2Dtmp = complex_zero;
    for ( int chrom = 0; chrom < nchrom; chrom ++ ){
        // dipole part -- NOTE: I DO NOT MAKE THE CONDON APPROXIMATION
        // Also note that this is an isotropic spectrum where all three
        // pulses have the same polarization 
        mu      = dot3(dipole_t0[chrom],dipole_t1[chrom])*
                  dot3(dipole_t2[chrom],dipole_t3[chrom]);

        // average the response function over all of the chromophores then return
        // note that we have included the T2 relaxation here
        if ( which.compare("R1") == 0 ){ // rephasing
            R2Dtmp += img*mu* 
                      conj(eint_t1[chrom])*exp(-t1*dt/RT2)  // first pulse oscillating in 01 coherence
                                          *exp(-t2*dt/RT1)  // second pulse population relaxation
                          *eint_t3[chrom] *exp(-t3*dt/RT2); // third pulse oscillating in 01 coherence
        }
        else if ( which.compare("R4") == 0 ){ // non-rephasing
            R2Dtmp += img*mu*  
                      eint_t1[chrom]*exp(-t1*dt/RT2)        // first pulse oscillating in 01 coherence
                                    *exp(-t2*dt/RT1)        // second pulse population relaxation
                    *eint_t3[chrom] *exp(-t3*dt/RT2);       // third pulse oscillating in 01 coherence
        }
        else if ( which.compare("R3") == 0 ){ // rephasing
            R2Dtmp += 2.*img*mu*                               // the factor of 2 assumes the transition dipoles scale like a harmonic oscillator (see p 68 of Hamm and Zanni)
                      conj(eint_t1[chrom])*exp(-t1*dt/RT2)  // first pulse oscillating in 01 coherence
                                          *exp(-t2*dt/RT1)  // second pulse population relaxation
                          *eint_t3[chrom] *exp(-t3*dt/RT2)  // third pulse oscillating in 01 coherence
                          *exp(img*(dt*t3)*anharm/HBAR);         // include anharmonicity term for the 12 transition
                            
        }
        else if ( which.compare("R6") == 0 ){ // non-rephasing
            R2Dtmp += 2.*img*mu*                               // the factor of 2 assumes the transition dipoles scale like a harmonic oscillator (see p 68 of Hamm and Zanni)
                      eint_t1[chrom]*exp(-t1*dt/RT2)        // first pulse oscillating in 01 coherence
                                    *exp(-t2*dt/RT1)        // second pulse population relaxation
                    *eint_t3[chrom] *exp(-t3*dt/RT2)        // third pulse oscillating in 01 coherence
                    *exp(img*(dt*t3)*anharm/HBAR);
        }
        else {
            cout << "ERROR:: IR2D::getR2D which= " << which << " is unknown. Aborting." << endl;
            exit(EXIT_FAILURE);
        }
    }

    return R2Dtmp/(1.*nchrom);
}



int IR2D::get_eint_t1( int t1 )
// integral from t0 to t1 for linear and third order response functions of 01 frequency
// See Eq 7.10 from Hamm and Zanni for linear response function
// See Eq 7.35 from Hamm and Zanni for third order response function
{
    complex<double> arg;

    for ( int chrom = 0; chrom < nchrom; chrom ++ ){
        // reset integral to zero at t1 = 0 (the exponential will be 1)
        if ( t1 == 0 ) eint_t1[chrom] = complex_one;
        else {
            // integrate the equation for each chromophore using the trapezoid rule
            arg = -img*dt*(energy_t1_last[chrom]+energy_t1[chrom])/(2.*HBAR);
            eint_t1[chrom] *= exp(arg);
        }
    }
    
    return IR2DOK;
}

int IR2D::get_eint_t3( int t3 )
// integral from t1+t2 to t1+t2+t3 for third order response functions of 01 frequency
// See Eq 7.35 from Hamm and Zanni for third order response function
{
    complex<double> arg;

    for ( int chrom = 0; chrom < nchrom; chrom ++ ){
        // reset integral to zero at t3 = 0 (the exponential will be 1)
        if ( t3 == 0 ) eint_t3[chrom] = complex_one;
        else {
            // integrate the equation for each chromophore using the trapezoid rule
            arg = -img*dt*(energy_t3_last[chrom]+energy_t3[chrom])/(2.*HBAR);
            eint_t3[chrom] *= exp(arg);
        }
    }
    
    return IR2DOK;
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
    double convert, freq, scale, abs[length], dis[length];

    fn = _ofile_+"-R1D.dat";
    ofile.open(fn);
    if ( ! ofile.is_open() ) { fileOpenErr( fn ); exit(EXIT_FAILURE);}
    cout << ">>> Writing " << fn << "." << endl;

    ofile << "# time (ps) Real Imag" << endl;
    for ( int i = 0; i < t1_npoints; i ++ ){
        ofile << i * dt << " " << R1D[i].real() << " " << R1D[i].imag() << endl;
    }
    ofile.close();

    // do fft
    plan = fftw_plan_dft_1d( length, reinterpret_cast<fftw_complex*>(fftIn), \
                                     reinterpret_cast<fftw_complex*>(fftOut),\
                                     FFTW_BACKWARD, FFTW_ESTIMATE );
    
    // See Eq 4.8 from Hamm and Zanni -- Absorptive part
    for ( int i = 0; i < length ; i ++ ) fftIn[i] = complex_zero;
    for ( int i = 0; i < t1_npoints; i ++ ){
        fftIn[i]        = img*R1D[i];
        fftIn[length-i] = conj(img*R1D[i]);
    }
    fftw_execute(plan); for ( int i = 0; i < length; i ++ ) abs[i] = fftOut[i].real();

    // See Eq 4.11 from Hamm and Zanni -- Dispersive part
    for ( int i = 0; i < length ; i ++ ) fftIn[i] = complex_zero;
    for ( int i = 0; i < t1_npoints; i ++ ){
        fftIn[i]        = -R1D[i];
        fftIn[length-i] = conj(-R1D[i]);
    }
    fftw_execute(plan); for ( int i = 0; i < length; i ++ ) dis[i] = fftOut[i].real();

    fn = _ofile_+"-FTIR.dat";
    ofile.open(fn);
    if ( ! ofile.is_open() ) { fileOpenErr( fn ); exit(EXIT_FAILURE);}
    cout << ">>> Writing " << fn << "." << endl;
    ofile << "# frequency (cm-1) Absorptive Dispersive" << endl;

    scale = -dt*length/(2.*PI*HBAR*sqrt(length));
    // scale intensity to keep total area of the spectrum conserved and normalize
    // by root(length) since FFTW3 does not
    // multiply the signal by -1 so it looks like absorption and not signal intensity
    scale = -dt*length/(2.*PI*HBAR*sqrt(length));
        
    // negative frequencies are stored at the end of the fftOut array
    for ( int i = length/2; i < length; i ++ )
    {
        // get frequency in wavenumber, add back the shift
        freq   = 2.*PI*HBAR*(i-length)/(dt*length) + shift;
        if ( freq < window0 or freq > window1 ) continue;
        ofile << freq << " " << scale*abs[i] << " " << scale*dis[i] << endl;
    }
    // positive frequencies are stored at the beginning
    for ( int i = 0; i < length/2; i ++ ){
        // get frequency in wavenumber, add back the shift
        freq   = 2.*PI*HBAR*i/(dt*length) + shift;
        if ( freq < window0 or freq > window1 ) continue;
        ofile << freq << " " << scale*abs[i] << " " << scale*dis[i] << endl;
    }
    ofile.close();

}

int main( int argc, char* argv[] )
{
    int sample;
    int frame0, frame_t1, frame_t2, frame_t3;
    int t1, t2, t3;
    int err;

    if ( argc != 2 ){
        cout << "ERROR:: Program expects the name of the input" << 
                "file as the only argument. Aborting." << endl;
        exit( EXIT_FAILURE );
    }

    // get input file name
    string inpf( argv[1] );

    // Initialize the IR2D class
    IR2D spectrum( inpf );

    // set t2 waiting time in units of frame
    t2 = static_cast<int>(spectrum.t2/spectrum.dt);

    // Loop over the trajectory
    for ( sample = 0; sample < spectrum.nsamples; sample ++ ){

        // get frame number, read and save dipole at t0
        frame0 = static_cast<int>(spectrum.sample_every/spectrum.dt);
        spectrum.readDframe(frame0, "t0" );
        cout << sample << "/" << spectrum.nsamples << endl;
        
        // loop over t1
        for ( int t1 = 0; t1 < spectrum.t1_npoints; t1 ++ ){
            // get frame number for current t1
            frame_t1 = frame0 + t1;

            // read in energy and dipole
            if ( err = spectrum.readEframe(frame_t1, "t1") != IR2DOK ){ 
                cout << "ERROR: IR2D::readEframe returned error " << err << endl;
            }
            if ( err = spectrum.readDframe(frame_t1, "t1") != IR2DOK ){ 
                 cout << "ERROR: IR2D::readDframe returned error " << err << endl;
            }
            
            // get exponential integral, calculate 1D response function, and save energy in t1_last
            spectrum.get_eint_t1( t1 );
            spectrum.R1D[t1] += spectrum.getR1D( t1 );
            memcpy( spectrum.energy_t1_last, spectrum.energy_t1, spectrum.nchrom*sizeof(double) );

            // get frame number and dipole for time t2
            frame_t2 = frame0 + t1 + t2;
            if ( err = spectrum.readDframe(frame_t2, "t2") != IR2DOK ){ 
                    cout << "ERROR: IR2D::readDframe returned error " << err << endl;
            }

            // loop over t3 at current t1
            // TODO:: Could parallelize with open mp if desired... may require restructuring
            // may be easier to parallelize over samples
            for ( t3 = 0; t3 < spectrum.t3_npoints; t3 ++ ){
                // get frame number for current t3
                frame_t3 = frame0 + t1 + t2 + t3;

                // read in energy and dipole
                if ( err = spectrum.readEframe(frame_t3, "t3") != IR2DOK ){ 
                    cout << "ERROR: IR2D::readEframe returned error " << err << endl;
                }
                if ( err = spectrum.readDframe(frame_t3, "t3") != IR2DOK ){ 
                    cout << "ERROR: IR2D::readDframe returned error " << err << endl;
                }

                // calculate 2D response function then save energy in t3_last
                spectrum.get_eint_t3(t3);
                spectrum.R2D_R1[ t1 * spectrum.t3_npoints + t3 ] += spectrum.getR2D(t1, t3, "R1" );
                spectrum.R2D_R3[ t1 * spectrum.t3_npoints + t3 ] += spectrum.getR2D(t1, t3, "R3" );
                spectrum.R2D_R4[ t1 * spectrum.t3_npoints + t3 ] += spectrum.getR2D(t1, t3, "R4" );
                spectrum.R2D_R6[ t1 * spectrum.t3_npoints + t3 ] += spectrum.getR2D(t1, t3, "R6" );
                memcpy( spectrum.energy_t3_last, spectrum.energy_t3, spectrum.nchrom*sizeof(double) );

            }
        }
    }

    // normalize response functions by number of samples
    for ( int t1 = 0; t1 < spectrum.t1_npoints; t1 ++ ){
        spectrum.R1D[t1]/=(1.*spectrum.nsamples);
        for ( t3 = 0; t3 < spectrum.t3_npoints; t3 ++ ){
            spectrum.R2D_R1[ t1 * spectrum.t3_npoints + t3 ] /= ( 1.*spectrum.nsamples );
            spectrum.R2D_R3[ t1 * spectrum.t3_npoints + t3 ] /= ( 1.*spectrum.nsamples );
            spectrum.R2D_R4[ t1 * spectrum.t3_npoints + t3 ] /= ( 1.*spectrum.nsamples );
            spectrum.R2D_R6[ t1 * spectrum.t3_npoints + t3 ] /= ( 1.*spectrum.nsamples );
        }
    }

    // write the response functions and the spectrum
    spectrum.write1D();
    // TODO: 2D response functions
}
