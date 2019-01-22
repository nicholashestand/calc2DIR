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
    string line;

    // read input file
    if ( readParam( _inpf_ ) != IR2DOK ) exit( EXIT_FAILURE );

    // allocate variable arrays
    R1D             = new complex<double>[t1t3_npoints]();
    R2D_R1          = new complex<double>[t1t3_npoints*t1t3_npoints]();
    R2D_R3          = new complex<double>[t1t3_npoints*t1t3_npoints]();
    R2D_R4          = new complex<double>[t1t3_npoints*t1t3_npoints]();
    R2D_R6          = new complex<double>[t1t3_npoints*t1t3_npoints]();

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
    delete [] R1D;
    delete [] R2D_R1;
    delete [] R2D_R3;
    delete [] R2D_R4;
    delete [] R2D_R6;

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
        else if ( para.compare("trjlen")      == 0 )  trjlen       = stoi(value);
        else if ( para.compare("output_file") == 0 )  _ofile_      = value;
        else if ( para.compare("time_step")   == 0 )  dt           = stof(value);
        else if ( para.compare("t1t3_max")    == 0 )  t1t3_max     = stof(value);
        else if ( para.compare("t2")          == 0 )  t2           = stof(value);
        else if ( para.compare("lifetimet1")  == 0 )  lifetime_T1  = stof(value);
        else if ( para.compare("anharm")      == 0 )  anharm       = stof(value);
        else if ( para.compare("nsamples")    == 0 )  nsamples     = stoi(value);
        else if ( para.compare("sample_every") == 0 ) sample_every = stof(value);
        else if ( para.compare("fftlen")      == 0 )  fftlen       = stoi(value);
        else if ( para.compare("window0")     == 0 )  window0      = stof(value);
        else if ( para.compare("window1")     == 0 )  window1      = stof(value);
        else cerr << "\tWARNING:: parameter " << para << " is not recognized." << endl;
    }
    inpf.close();

    // some output to confirm parameters
    tellParam<string>( "energy_file", _efile_ );
    tellParam<string>( "dipole_file", _dfile_ );
    tellParam<int>( "trjlen", trjlen );
    tellParam<string>( "output_file", _ofile_ );
    tellParam<double>( "time_step", dt );
    tellParam<double>( "t1t3_max", t1t3_max );
    tellParam<double>( "t2", t2 );
    tellParam<double>( "lifetimeT1", lifetime_T1 );
    tellParam<double>( "anharm", anharm );
    tellParam<int>( "nsamples", nsamples );
    tellParam<int>( "sample_every", sample_every );
    tellParam<int>( "fftlen", fftlen );
    tellParam<double>( "window0", window0 );
    tellParam<double>( "window1", window1 );
    cout << ">>> Done reading simulation parameters from " << _inpf_ << endl;

    if ( trjlen < static_cast<int>(((nsamples-1)*sample_every + (2*t1t3_max + t2))/dt) ){
        cout << "WARNING:: The given trajectory length is not long enough.\n" << 
                      "\t  Must be " <<  ((nsamples-1)*sample_every + (2*t1t3_max + t2))/dt << 
                      " frames long.\n\t  Check input file. Aborting." << endl;
        exit(EXIT_FAILURE);
    }

    // set number of points in response functions based on t1t3_max and dt
    t1t3_npoints = static_cast<int>(t1t3_max/dt);

    // shift energies by mean of spectral window to avoid high frequency oscillations
    shift = 0.5*(window1 + window0);
    cout << ">>> Shifting frequencies by: " << shift << " cm." << endl;

    // set T2 lifetime as 2T1 (see Hamm and Zanni p29 for discussion)
    // also see Liang and Jansen JCTC 2012 eq 14 and 16
    lifetime_T2 = 2*lifetime_T1;
    cout << ">>> Setting T2 to 2 T1: " << lifetime_T2 << " ps." << endl;

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

void IR2D::fileReadErr( string _fn_ )
// give an error message when I cant open a file
{
    cerr << "ERROR:: Reading " << _fn_ << " failed. Probably reached EOF...Aborting." << endl;
}

int IR2D::readEframe( int frame, string which )
// Read the energy file
{
    int     frameTmp;
    float   energyTmp;
    int64_t file_offset;

    file_offset = frame*(sizeof(int)+sizeof(float));
    efile.seekg( file_offset );

    // the energy for each frame is stored on a seperate line
    // with an integer frame number at the beginning of the line
    efile.read( (char*)&frameTmp , sizeof(int) );
    if ( not efile.good() ){ fileReadErr( _efile_ ); return 1;}
    efile.read( (char*)&energyTmp, sizeof(float) );
    if ( not efile.good() ){ fileReadErr( _efile_ ); return 1;}

    // only keep the energies, ignore the couplings for now
    if      ( which.compare("t1") == 0 ) energy_t1 = energyTmp - shift;
    else if ( which.compare("t3") == 0 ) energy_t3 = energyTmp - shift;
    else{
        cout << "ERROR:: IR2D::readEframe which= " << which << " unknown." << endl;
        return 1;
    }

    return IR2DOK;
}

int IR2D::readDframe( int frame, string which )
// Read the dipole file
{
    int     frameTmp, i;
    float   dipoleTmp[ 3 ];
    int64_t file_offset;

    file_offset = frame*(sizeof(int)+sizeof(float)*3);
    dfile.seekg( file_offset );

    dfile.read( (char*)&frameTmp, sizeof(int) );
    if ( not dfile.good() ){ fileReadErr( _dfile_ ); return 1;}
    dfile.read( (char*)dipoleTmp, sizeof(float)*3 );
    if ( not dfile.good() ){ fileReadErr( _dfile_ ); return 1;}

    // put these into the dipole vector variable
    // note that in the bin file all x's come first, then y's, etc
    for ( i = 0; i < 3; i ++ ){
        if      ( which.compare("t0") == 0 ) dipole_t0[i] = dipoleTmp[i];
        else if ( which.compare("t1") == 0 ) dipole_t1[i] = dipoleTmp[i];
        else if ( which.compare("t2") == 0 ) dipole_t2[i] = dipoleTmp[i];
        else if ( which.compare("t3") == 0 ) dipole_t3[i] = dipoleTmp[i];
        else{
            cout << "ERROR: IR2D::readDframe which= " << which << " unknown." << endl;
            return 1;
        }
    }

    return IR2DOK;
}

complex<double> IR2D::getR1D( int t1 )
// return the linear response function at a given t1
// See Eq 7.10 from Hamm and Zanni
{
    double mu;
    complex<double> R1D;
    vec3 polx={1.,0.,0.}, poly={0.,1.,0.}, polz={0.,0.,1.};

    // get dipole part for isotropically averaged spectrum
    // NOTE: The Condon approximation is NOT made here
    mu = 0.;
    mu+=dot3(dipole_t0,polx)*dot3(dipole_t1,polx);
    mu+=dot3(dipole_t0,poly)*dot3(dipole_t1,poly);
    mu+=dot3(dipole_t0,polz)*dot3(dipole_t1,polz);

    // get the response function
    R1D = img*mu*eint_t1*exp(-t1*dt/lifetime_T2);

    return R1D;
}

complex<double> IR2D::getR2D( int t1, int t3, string which )
// return the third order response function at a given t1, t3
// See Eq 7.35 from Hamm and Zanni
{
    double mu;
    complex<double> R2D;
    vec3   polx={1.,0.,0.}, poly={0.,1.,0.}, polz={0.,0.,1.};
    double lifetime_T2_12 = 2.*lifetime_T1/3.;// See Hamm and Zanni eq 4.21
    lifetime_T2_12 = lifetime_T2;             // set equal here, see Jansen 2012 

    R2D = complex_zero;
    // get dipole part for isotropically averaged ZZZZ spectrum
    // NOTE: The Condon approximation is NOT made here
    // all four pulses have the same polarization 
    // see Hamm and Zanni eq 5.35
    mu = 0;
    mu+=dot3(dipole_t0,polx)*dot3(dipole_t1,polx)*
        dot3(dipole_t2,polx)*dot3(dipole_t3,polx);
    mu+=dot3(dipole_t0,poly)*dot3(dipole_t1,poly)*
        dot3(dipole_t2,poly)*dot3(dipole_t3,poly);
    mu+=dot3(dipole_t0,polz)*dot3(dipole_t1,polz)*
        dot3(dipole_t2,polz)*dot3(dipole_t3,polz);

    // get the response function
    if ( which.compare("R1") == 0 ){ // rephasing
        // dipole
        // first pulse oscillating in 01 coherence
        // second pulse population relaxation (note t2 is in ps already)
        // third pulse oscillating in 01 coherence
        R2D += img*mu* 
               conj(eint_t1)*exp(-t1*dt/lifetime_T2)
                            *exp(-t2   /lifetime_T1)
                   *eint_t3 *exp(-t3*dt/lifetime_T2);
    }
    else if ( which.compare("R4") == 0 ){ // non-rephasing
        // dipole
        // first pulse oscillating in 01 coherence
        // second pulse population relaxation (note t2 is in ps already)
        // third pulse oscillating in 01 coherence
        R2D += img*mu* 
               eint_t1*exp(-t1*dt/lifetime_T2)
                      *exp(-t2   /lifetime_T1)
              *eint_t3*exp(-t3*dt/lifetime_T2);
    }
    else if ( which.compare("R3") == 0 ){ // rephasing
        // dipole, the factor of 2 assumes the transition dipoles scale 
        // like a harmonic oscillator (see p 68 of Hamm and Zanni)
        // first pulse oscillating in 01 coherence
        // second pulse population relaxation (note t2 is in ps already)
        // third pulse oscillating in 12 coherence -- see eq 4.21 for relaxation
        // include anharmonicity term for the 12 transition
        R2D -= 2.*img*mu*
               conj(eint_t1)*exp(-t1*dt/lifetime_T2)
                            *exp(-t2   /lifetime_T1)
                   *eint_t3 *exp(-t3*dt/lifetime_T2_12)
                   *exp(img*(dt*t3)*anharm/HBAR);              
    }
    else if ( which.compare("R6") == 0 ){ // non-rephasing
        // dipole, the factor of 2 assumes the transition dipoles scale 
        // like a harmonic oscillator (see p 68 of Hamm and Zanni)
        // first pulse oscillating in 01 coherence
        // second pulse population relaxation (note t2 is in ps already)
        // third pulse oscillating in 12 coherence -- see eq 4.21 for relaxation
        // include anharmonicity term for the 12 transition
        R2D -= 2.*img*mu*
               eint_t1*exp(-t1*dt/lifetime_T2)
                      *exp(-t2   /lifetime_T1)
              *eint_t3*exp(-t3*dt/lifetime_T2_12)
              *exp(img*(dt*t3)*anharm/HBAR);
    }
    else {
        cout << "ERROR:: IR2D::getR2D which= " << which << " is unknown. Aborting." << endl;
        exit(EXIT_FAILURE);
    }

    return R2D;
}

int IR2D::get_eint( int t, string which )
// integral from t0 to t1 for linear and third order response functions of 01 frequency
// See Eq 7.10 from Hamm and Zanni for linear response function
// See Eq 7.35 from Hamm and Zanni for third order response function
{
    complex<double> arg;

    if ( which.compare("t1") == 0 ){
        // reset integral to zero at t1 = 0 (the exponential will be 1)
        if ( t == 0 ) eint_t1 = complex_one;
        else {
            // integrate the equation using the trapezoid rule
            arg = -img*dt*(energy_t1_last+energy_t1)/(2.*HBAR);
            eint_t1 *= exp(arg);
        }
        // save current energies for next time step so can do integration
        energy_t1_last = energy_t1;
    }
    else if ( which.compare("t3") == 0 ){
        // reset integral to zero at t3 = 0 (the exponential will be 1)
        if ( t == 0 ) eint_t3 = complex_one;
        else {
            // integrate the equation using the trapezoid rule
            arg = -img*dt*(energy_t3_last+energy_t3)/(2.*HBAR);
            eint_t3 *= exp(arg);
        }
        // save current energies for next time step so can do integration
        energy_t3_last = energy_t3;
    }
    else{
        cout << "ERROR:: IR2D::get_eint which= " << which << " unknown." << endl;
        return 1;
    }
    
    return IR2DOK;
}

double IR2D::dot3( vec3 a, vec3 b )
// dot product of 3 vector
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; 
}

int IR2D::writeR1D()
// write R1D to file
{
    string   fn=_ofile_+"-R1D.dat";
    ofstream ofile;
    int t1;

    ofile.open( fn );
    if ( ! ofile.is_open() ) { fileOpenErr( fn ); return 1;}
    cout << ">>> Writing " << fn << "." << endl;

    ofile << "# time (ps) Real Imag" << endl;
    for ( t1 = 0; t1 < t1t3_npoints; t1 ++ ){
        ofile << t1 * dt << " " << R1D[t1].real() << " " << R1D[t1].imag() << endl;
    }
    ofile.close();

    return IR2DOK;
}

int IR2D::writeR2D()
// write R2D to file
{
    string fn;
    ofstream ofile;
    int t1, t3;
    complex<double> Rtmp;

    // write rephasing response function
    fn=_ofile_+"-RparI.dat";
    ofile.open( fn );
    if ( ! ofile.is_open() ) { fileOpenErr( fn ); return 1;}
    cout << ">>> Writing " << fn << "." << endl;

    ofile << "# Rephasing parallel ZZZZ polarized response function, t2 = " << t2 << endl;
    ofile << "# t1 (ps) t3 (ps) Real Imag" << endl;
    for ( t1 = 0; t1 < t1t3_npoints; t1 ++ ){
        for ( t3 = 0; t3 < t1t3_npoints; t3 ++ ){
            // see Hamm and Zanni eq 4.23 and note R1=R2, hence the factor of 2
            // The experiment can only see all rephasing or non-rephasing, not the individual
            // response functions, so add them here.
            Rtmp = 2.*R2D_R1[t1 * t1t3_npoints + t3] + R2D_R3[t1 * t1t3_npoints + t3];
            ofile << t1 * dt << " " << t3 * dt << " " << Rtmp.real() << " " << Rtmp.imag() << endl;
        }
    }
    ofile.close();

    // write non-rephasing response function
    fn=_ofile_+"-RparII.dat";
    ofile.open( fn );
    if ( ! ofile.is_open() ) { fileOpenErr( fn ); return 1;}
    cout << ">>> Writing " << fn << "." << endl;

    ofile << "# Non-rephasing parallel ZZZZ polarized response function, t2 = " << t2 << endl;
    ofile << "# t1 (ps) t3 (ps) Real Imag" << endl;
    for ( t1 = 0; t1 < t1t3_npoints; t1 ++ ){
        for ( t3 = 0; t3 < t1t3_npoints; t3 ++ ){
            // see Hamm and Zanni eq 4.23 and note R4=R5, hence the factor of 2
            // The experiment can only see all rephasing or non-rephasing, not the individual
            // response functions, so add them here.
            Rtmp = 2.*R2D_R4[t1 * t1t3_npoints + t3] + R2D_R6[t1 * t1t3_npoints + t3];
            ofile << t1 * dt << " " << t3 * dt << " " << Rtmp.real() << " " << Rtmp.imag() << endl;
        }
    }
    ofile.close();

    return IR2DOK;
}

int IR2D::write1Dfft()
// write fft of R1D
{
    string   fn;
    ofstream ofile;
    complex<double> *fftIn, *fftOut, *res;
    fftw_plan plan;
    double   freq, scale;
    int      t1, i;

    // allocate arrays
    fftIn  = new complex<double>[fftlen]();
    fftOut = new complex<double>[fftlen]();
    res    = new complex<double>[fftlen]();

    // scale fftOut to keep total area of the spectrum conserved and normalize
    // by root(fftlen) since FFTW3 does not
    // the -1 makes it look like absorption spectrum
    scale = -dt*fftlen/(2.*PI*HBAR*sqrt(fftlen));
     
    // do fft
    plan = fftw_plan_dft_1d( fftlen, reinterpret_cast<fftw_complex*>(fftIn), \
                                     reinterpret_cast<fftw_complex*>(fftOut),\
                                     FFTW_BACKWARD, FFTW_ESTIMATE );

    if ( 2*t1t3_npoints > fftlen ){
        cout << "ERROR:: fftlen = " << fftlen << " < " << "2*t1t3_max/dt = " 
             << 2*t1t3_npoints << endl;
        cout << "Specify longer fftlen in input file. Aborting." << endl;
        exit(EXIT_FAILURE);
    }
    
    // See Eq 4.8 from Hamm and Zanni -- Absorptive part
    for ( i = 0; i < fftlen ; i ++ ) fftIn[i] = complex_zero;
    for ( t1 = 0; t1 < t1t3_npoints; t1 ++ ){
        fftIn[t1]          = img*R1D[t1];
        if ( t1 == 0 ) continue; // dont take hermitian at t=0
        fftIn[fftlen-t1] = conj(img*R1D[t1]);
    }
    fftw_execute(plan); 
    for ( t1 = 0; t1 < fftlen; t1 ++ ) res[t1].real(fftOut[t1].real());

    // See Eq 4.11 from Hamm and Zanni -- Dispersive part
    for ( i = 0; i < fftlen ; i ++ ) fftIn[i] = complex_zero;
    for ( t1 = 0; t1 < t1t3_npoints; t1 ++ ){
        fftIn[t1]          = -R1D[t1];
        if ( t1 == 0 ) continue; // dont take hermitian at t=0
        fftIn[fftlen-t1] = conj(-R1D[t1]);
    }
    fftw_execute(plan); 
    fftw_destroy_plan(plan);
    for ( t1 = 0; t1 < fftlen; t1 ++ ) res[t1].imag(fftOut[t1].real());

    // write file
    fn = _ofile_+"-R1Dw.dat";
    ofile.open(fn);
    if ( ! ofile.is_open() ) { fileOpenErr( fn ); return 1;}
    cout << ">>> Writing " << fn << "." << endl;
    ofile << "# Frequency (cm-1) Absorptive Dispersive" << endl;

    // negative frequencies are stored at the end of the fftOut array
    for ( i = fftlen/2; i < fftlen; i ++ )
    {
        // get frequency in wavenumber, add back the shift
        freq   = 2.*PI*HBAR*(i-fftlen)/(dt*fftlen) + shift;
        if ( freq < window0 or freq > window1 ) continue;
        ofile << freq << " " << scale*res[i].real() << " " << scale*res[i].imag() << endl;
    }
    // positive frequencies are stored at the beginning of fftOut
    for ( i = 0; i < fftlen/2; i ++ ){
        // get frequency in wavenumber, add back the shift
        freq   = 2.*PI*HBAR*i/(dt*fftlen) + shift;
        if ( freq < window0 or freq > window1 ) continue;
        ofile << freq << " " << scale*res[i].real() << " " << scale*res[i].imag() << endl;
    }
    ofile.close();

    // delete arrays
    delete [] fftIn;
    delete [] fftOut;
    delete [] res;

    return IR2DOK;
}

int IR2D::write2DRabs()
// write the purely absorptive 2D IR spectrum
{
    string    fn;
    complex<double> *fftIn, *fftOut, *res;
    fftw_plan plan;
    int       t1, t3, i, j;

    // allocate arrays
    fftIn  = new complex<double>[fftlen*fftlen]();
    fftOut = new complex<double>[fftlen*fftlen]();
    res    = new complex<double>[fftlen*fftlen]();

    // do fft plan
    plan = fftw_plan_dft_2d( fftlen, fftlen, \
                             reinterpret_cast<fftw_complex*>(fftIn), \
                             reinterpret_cast<fftw_complex*>(fftOut),\
                             FFTW_BACKWARD, FFTW_ESTIMATE );
    
    // fourier transform rephasing response functions, see Hamm and Zanni eq 4.31
    // see Hamm and Zanni eq 4.23 and note R1=R2, hence the factor of 2
    for ( i = 0; i < fftlen*fftlen; i ++ ) fftIn[i] = complex_zero;
    for ( t1 = 0; t1 < t1t3_npoints; t1 ++ ){
        for ( t3 = 0; t3 < t1t3_npoints; t3 ++ ){
            fftIn[ t1*fftlen + t3 ] = 2.*img*R2D_R1[ t1*t1t3_npoints + t3 ]\
                                      +  img*R2D_R3[ t1*t1t3_npoints + t3 ];
            // divide t=0 point by 2 (see Hamm and Zanni sec 9.5.3)
            if ( t1 == 0 and t3 == 0 ) fftIn[ t1*fftlen + t3 ] /=2.;
        }
    }
    fftw_execute(plan);
    fn=_ofile_+"-RparIw.dat";
    if ( write2Dout( fftOut, fn, "rephasing", fftlen ) != IR2DOK ) return 1;

    // Save rephasing contribution to purly absorptive spectrum
    // here we map w1 to -w1
    // See Hamm and Zanni eq 4.36
    for ( i = 0; i < fftlen; i ++ ){
        for ( j = 0; j < fftlen; j ++ ){
            if ( i == 0 ) res[ i*fftlen + j ] = fftOut[ i*fftlen + j ]; // w1=0 goes to w1=0
            else          res[ i*fftlen + j ] = fftOut[ (fftlen - i)*fftlen + j ]; // w1 to -w1
        }
    }

    // fourier transform non-rephasing response functions, see Hamm and Zanni eq 4.31
    // see Hamm and Zanni eq 4.23 and note R4=R5, hence the factor of 2
    for ( i = 0; i < fftlen*fftlen; i ++ ) fftIn[i] = complex_zero;
    for ( t1 = 0; t1 < t1t3_npoints; t1 ++ ){
        for ( t3 = 0; t3 < t1t3_npoints; t3 ++ ){
            fftIn[ t1*fftlen + t3 ] = 2.*img*R2D_R4[ t1*t1t3_npoints + t3 ]\
                                      +  img*R2D_R6[ t1*t1t3_npoints + t3 ];
            // divide t=0 point by 2 (see Hamm and Zanni sec 9.5.3)
            if ( t1 == 0 and t3 == 0 ) fftIn[ t1*fftlen + t3 ] /=2.;
        }
    }
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fn=_ofile_+"-RparIIw.dat";
    if ( write2Dout( fftOut, fn, "non-rephasing", fftlen ) != IR2DOK ) return 1;
    
    // Save rephasing contribution to purly absorptive spectrum
    // here w1 goes to w1
    // See Hamm and Zanni eq 4.36
    for ( i = 0; i < fftlen; i ++ ){
        for ( j = 0; j < fftlen; j ++ ){
            res[ i*fftlen + j ] += fftOut[ i*fftlen + j ];
        }
    }

    // write purly absorptive spectrum 
    fn=_ofile_+"-RparAbs.dat";
    if ( write2Dout( res, fn, "rabs", fftlen ) != IR2DOK ) return 1;

    // delete arrays
    delete [] fftIn;
    delete [] fftOut;
    delete [] res;

    return IR2DOK;
}


int IR2D::write2Dout( complex<double> *data, string fn, string which, int n )
// write 2D fourier transformed output
{
    int i, j;
    double shift_w1, shift_w3, window0_w1, window0_w3, window1_w1, window1_w3;
    double scale, w1, w3;
    ofstream ofile;

    // scaling
    scale = dt*n/(2.*PI*HBAR*sqrt(n));
    scale*= scale;  // since 2D scaling

    ofile.open( fn );
    if ( ! ofile.is_open() ) { fileOpenErr( fn ); return 1;}
    cout << ">>> Writing " << fn << "." << endl;

    // assign shifts and spectral window limits
    shift_w1   = shift;
    shift_w3   = shift;
    window0_w1 = window0;
    window1_w1 = window1;
    window0_w3 = window0;
    window1_w3 = window1;

    if ( which.compare("rephasing") == 0 ){
        ofile << "# Rephasing parallel ZZZZ polarized response function, t2 = " << t2 << endl;
        shift_w1   = -shift;
        window0_w1 = -window1;
        window1_w1 = -window0;
    }
    else if ( which.compare("non-rephasing") == 0 ){
        ofile << "# Non-rephasing parallel ZZZZ polarized response function, t2 = " << t2 << endl;
    }
    else if ( which.compare("rabs") == 0 ){
        ofile << "# Pure parallel ZZZZ polarized absorption, t2 = " << t2 << endl;
    }
    else{
        cout << "ERROR:: write2Dout which= " << which << " unknown." << endl;
        return 1;
    }

    ofile << "# w1 (cm-1) w3 (cm-1) Real Imag" << endl;
    // negative frequencies are stored at the end of data array
    for ( i = n/2; i < n; i ++ ){
        w1 = 2.*PI*HBAR*(i-n)/(dt*n) + shift_w1;
        if ( w1 < window0_w1 or w1 > window1_w1 ) continue;
        for ( j = n/2; j < n; j ++ ){
            w3 = 2.*PI*HBAR*(j-n)/(dt*n) + shift_w3;
            if ( w3 < window0_w3 or w3 > window1_w3 ) continue;
            ofile << w1 << " " << w3 << " " << scale*data[i*n+j].real() 
                        << " " << scale*data[i*n+j].imag() << endl;
        }
        for ( j = 0; j < n/2; j ++ ){
            w3 = 2.*PI*HBAR*j/(dt*n) + shift_w3;
            if ( w3 < window0_w3 or w3 > window1_w3 ) continue;
            ofile << w1 << " " << w3 << " " << scale*data[i*n+j].real() 
                        << " " << scale*data[i*n+j].imag() << endl;
        }
    }
    // positive frequencies are stored at the beginning of data array
    for ( i = 0; i < n/2; i ++ ){
        w1 = 2.*PI*HBAR*i/(dt*n) + shift_w1;
        if ( w1 < window0_w1 or w1 > window1_w1 ) continue;
        for ( j = n/2; j < n; j ++ ){
            w3 = 2.*PI*HBAR*(j-n)/(dt*n) + shift_w3;
            if ( w3 < window0_w3 or w3 > window1_w3 ) continue;
            ofile << w1 << " " << w3 << " " << scale*data[i*n+j].real() 
                  << " " << scale*data[i*n+j].imag() << endl;
        }
        for ( j = 0; j < n/2; j ++ ){
            w3 = 2.*PI*HBAR*j/(dt*n) + shift_w3;
            if ( w3 < window0_w3 or w3 > window1_w3 ) continue;
            ofile << w1 << " " << w3 << " " << scale*data[i*n+j].real() 
                  << " " << scale*data[i*n+j].imag() << endl;
        }
    }
    ofile.close();

    return IR2DOK;
}

void printProgress( int currentStep, int totalSteps )
// print a progress bar to keep updated on usage
{
    float percentage = (float) currentStep / (float) totalSteps;
    int lpad = (int) (percentage*PWID);
    int rpad = PWID - lpad;
    fprintf(stderr, "\r [%.*s%*s]%3d%%", lpad, PSTR, rpad, "",(int) (percentage*100));
}

int main( int argc, char* argv[] )
{
    int sample;
    int frame_t0, frame_t1, frame_t2, frame_t3;
    int t1, t2, t3;

    if ( argc != 2 ){
        cout << "ERROR:: Program expects the name of the input" << 
                "file as the only argument. Aborting." << endl;
        exit( EXIT_FAILURE );
    }

    // get input file name and initialize IR2D class
    IR2D spectrum( argv[1] );

    // set t2 waiting time in units of frames
    t2 = static_cast<int>(spectrum.t2/spectrum.dt);

    // Loop over the trajectory
    for ( sample = 0; sample < spectrum.nsamples; sample ++ ){

        // get frame number, read and save dipole at t0
        frame_t0 = sample*static_cast<int>(spectrum.sample_every/spectrum.dt);
        fprintf(stderr, "    Now processing sample %d/%d starting at %.2f ps\n", \
                sample+1, spectrum.nsamples, frame_t0*spectrum.dt ); fflush(stderr);
        if ( spectrum.readDframe(frame_t0, "t0" ) != IR2DOK ) exit(EXIT_FAILURE);;

        // loop over t1
        for ( t1 = 0; t1 < spectrum.t1t3_npoints; t1 ++ ){
            // get frame number for current t1
            frame_t1 = frame_t0 + t1;

            // read in energy and dipole
            if ( spectrum.readEframe(frame_t1, "t1") != IR2DOK ) exit(EXIT_FAILURE);
            if ( spectrum.readDframe(frame_t1, "t1") != IR2DOK ) exit(EXIT_FAILURE);
            
            // get exponential integral and 1D response function
            if ( spectrum.get_eint( t1, "t1" ) != IR2DOK ) exit(EXIT_FAILURE);
            spectrum.R1D[t1] += spectrum.getR1D( t1 );

            // get frame number and dipole for time t2
            frame_t2 = frame_t1 + t2;
            if ( spectrum.readDframe(frame_t2, "t2") != IR2DOK ) exit(EXIT_FAILURE);

            // loop over t3 at current t1
            for ( t3 = 0; t3 < spectrum.t1t3_npoints; t3 ++ ){
                // get frame number for current t3
                frame_t3 = frame_t2 + t3;

                // read in energy and dipole
                if ( spectrum.readEframe(frame_t3, "t3") != IR2DOK ) exit(EXIT_FAILURE);
                if ( spectrum.readDframe(frame_t3, "t3") != IR2DOK ) exit(EXIT_FAILURE); 

                // get exponential integral and 2D response function 
                if ( spectrum.get_eint(t3, "t3") != IR2DOK ) exit(EXIT_FAILURE);
                spectrum.R2D_R1[ t1 * spectrum.t1t3_npoints + t3 ] += spectrum.getR2D(t1, t3, "R1" );
                spectrum.R2D_R3[ t1 * spectrum.t1t3_npoints + t3 ] += spectrum.getR2D(t1, t3, "R3" );
                spectrum.R2D_R4[ t1 * spectrum.t1t3_npoints + t3 ] += spectrum.getR2D(t1, t3, "R4" );
                spectrum.R2D_R6[ t1 * spectrum.t1t3_npoints + t3 ] += spectrum.getR2D(t1, t3, "R6" );
            }
            printProgress( t1+1, spectrum.t1t3_npoints );
        }
        cerr << endl;
    }

    // normalize response functions by number of samples
    for ( t1 = 0; t1 < spectrum.t1t3_npoints; t1 ++ ){
        spectrum.R1D[t1]/=(1.*spectrum.nsamples);
        for ( t3 = 0; t3 < spectrum.t1t3_npoints; t3 ++ ){
            spectrum.R2D_R1[ t1 * spectrum.t1t3_npoints + t3 ] /= ( 1.*spectrum.nsamples );
            spectrum.R2D_R3[ t1 * spectrum.t1t3_npoints + t3 ] /= ( 1.*spectrum.nsamples );
            spectrum.R2D_R4[ t1 * spectrum.t1t3_npoints + t3 ] /= ( 1.*spectrum.nsamples );
            spectrum.R2D_R6[ t1 * spectrum.t1t3_npoints + t3 ] /= ( 1.*spectrum.nsamples );
        }
    }

    // do fourier transforms and write them out
    if ( spectrum.writeR1D()    != IR2DOK ) exit(EXIT_FAILURE);
    if ( spectrum.writeR2D()    != IR2DOK ) exit(EXIT_FAILURE);
    if ( spectrum.write1Dfft()  != IR2DOK ) exit(EXIT_FAILURE);
    if ( spectrum.write2DRabs() != IR2DOK ) exit(EXIT_FAILURE);
    cout << ">>> Done!" << endl;
}
