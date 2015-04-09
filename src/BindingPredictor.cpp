#include <algorithm>
#include <ctime>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_siman.h>
#include <gsl/gsl_randist.h>

#include "BindingPredictor.h"

double FactorIntFuncBinary::compFactorInt( double normalInt, double dist, bool orientation ) const
{
	assert( dist >= 0 );

    if ( dist > distThr ) return 1.0;
    
	double spacingTerm = normalInt; 
	double orientationTerm = orientation ? 1.0 : orientationEffect;	
	return spacingTerm * orientationTerm;	
}

double FactorIntFuncLinear::compFactorInt( double normalInt, double dist, bool orientation ) const
{
	assert( dist >= 0 );

    if ( dist > distThr ) return 1.0;
    
	double spacingTerm = dist > range ? 1.0 + ( distThr - dist ) * ( normalInt - 1 ) / ( distThr - range ) : normalInt;
	double orientationTerm = orientation ? 1.0 : orientationEffect;	
	return spacingTerm * orientationTerm;	
}

double FactorIntFuncPeriodic::compFactorInt( double normalInt, double dist, bool orientation ) const
{
	assert( dist >= 0 );

    if ( dist > distThr ) return 1.0;
    
	double phase = 2 * M_PI * dist / period + ( orientation ? phi : phi0 );
    return normalInt * pow( spacingEffect, sin( phase ) ); 
}

double FactorIntFuncGeometric::compFactorInt( double normalInt, double dist, bool orientation ) const
{
	assert( dist >= 0 );
	
	double spacingTerm = max( 1.0, dist <= distThr ? normalInt : normalInt * pow( spacingEffect, dist - distThr ) ); 
	double orientationTerm = orientation ? 1.0 : orientationEffect;	
	return spacingTerm * orientationTerm;
}


BindingFunc::BindingFunc( const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< double >& _maxBindingWts, const Matrix& _factorIntMat ) : motifs( _motifs ), intFunc( _intFunc ), maxBindingWts( _maxBindingWts ), factorIntMat( _factorIntMat )
{
	int nFactors = motifs.size();
	assert( maxBindingWts.size() == nFactors );
	assert( factorIntMat.nRows() == nFactors );
	assert( factorIntMat.nCols() == nFactors );

//	cout << factorIntMat << endl;
	
	// the TF-TF matrix should be symmetric
	for ( int i = 0; i < nFactors; i++ ) {
		for ( int j = 0; j < nFactors; j++ ) assert( factorIntMat( i, j ) == factorIntMat( j, i ) );	
	}	
}

//qcheng75
void BindingFunc::releaseFactorIntMat(){//qcheng75
	if (this->factorIntMat.getData()!=NULL) this->factorIntMat.freeMatrix();
}

void BindingFunc::clear(){//qcheng75
	//if (! motifs.empty()) motifs.clear();
	if (! maxBindingWts.empty()) maxBindingWts.clear();
	releaseFactorIntMat();
	sites.clear();
	boundaries.clear();
	Z.clear();
	Zt.clear();

}

double BindingFunc::predictBinding( const SiteVec& targetSites, int factorIdx )
{
//     cout << "Sequence" << endl;
//     for ( int i = 0; i < targetSites.size(); i++ ) cout << targetSites[i] << endl;
	bindingWts.clear();
	boundaries.clear();
    int n = targetSites.size();
    
	// store the sequence and create the boundaries
	sites = targetSites;
    sites.insert( sites.begin(), Site( 0, true, 0 ) );  // start with a pseudo-site at position 0
    boundaries.push_back( 0 );
    for ( int i = 1; i <= n; i++ ) {
        int j; 
        for ( j = i - 1; j >= 1; j-- ) {
            if ( ( sites[i].start - sites[j].start ) > intFunc->getMaxDist() ) break; 
        }
        int boundary = j;
        boundaries.push_back( boundary );
    }
//     cout << "Boundaries";
//     for ( int i = 1; i <= n; i++ ) cout << "\t" << boundaries[i];
// 	cout << endl;
    
	// compute the Boltzman weights of binding for all sites
    bindingWts.push_back( 1.0 );
	for ( int i = 1; i <= n; i++ ) {
		bindingWts.push_back( compBindingWt( sites[ i ] ) );	
	}
		
	// initialization
	vector< double > Y( n + 1 );
    Y[0] = 0;
// 	for ( int i = 1; i <= n; i++ ) Y.push_back( 1 );
	vector< double > Yt( n + 1 );
    Yt[0] = 0;
    
	// compute the partition function of binding
	double Z_bind = compPartFunc();
// 	cout << "Z_bind = " << Z_bind << endl;
    
	// recurrence 
	for ( int i = 1; i <= n; i++ ) {
		int match = sites[ i ].factorIdx == factorIdx ? 1 : 0; 
		double sum = Yt[boundaries[i]] + match * Zt[boundaries[i]];
		for ( int j = boundaries[i] + 1; j < i; j++ ) {
			if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
			sum += compFactorInt( sites[ i ], sites[ j ] ) * ( Y[ j ] +  match * Z[ j ] );
		}
		Y[ i ] = bindingWts[ i ] * sum;
        Yt[i] = Y[i] + Yt[i - 1];
	}
	
	// the expected number of molecules bound
	double Y_total = Yt[n];
// 	for ( int i = 0; i < sites.size(); i++ ) {
// 		Y_total += Y[ i ];	
// 	}
	return Y_total / Z_bind;
}

double BindingFunc::compPartFunc()
{
    int n = sites.size() - 1;
    
	// initialization
    Z = vector< double >( n + 1 );
    Z[0] = 1.0;
    Zt = vector< double >( n + 1 );
    Zt[0] = 1.0;
	
	// recurrence 
	for ( int i = 1; i <= n; i++ ) {
		double sum = Zt[boundaries[i]];
		for ( int j = boundaries[i] + 1; j < i; j++ ) {
			if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
			sum += compFactorInt( sites[ i ], sites[ j ] ) * Z[ j ];	
		}
		Z[i] = bindingWts[ i ] * sum;
        Zt[i] = Z[i] + Zt[i - 1];
	}
	
	// the partition function 
// 	double Z_bind = 1;
// 	for ( int i = 0; i < sites.size(); i++ ) {
// 		Z_bind += Z[ i ];	
// 	}
    double Z_bind = Zt[n];
	return Z_bind;
}

double BindingFunc::compBindingWt( const Site& a ) const
{
	return ( maxBindingWts[ a.factorIdx ] * a.wtRatio );	
}

double BindingFunc::compFactorInt( const Site& a, const Site& b ) const
{
	assert( !siteOverlap( a, b, motifs ) );
	
	double normalInt = factorIntMat( a.factorIdx, b.factorIdx );
	double dist = abs( a.start - b.start );
	bool orientation = ( a.strand == b.strand ); 
	return intFunc->compFactorInt( normalInt, dist, orientation );	
}

string BindingPar::GLOBAL_INIT_PARAMETER_MODEL_FILE_NAME; // ## qcheng75
BindingPar::BindingPar( int _nInFactors, int _nOutFactors )
{	
	assert( _nInFactors > 0 && _nOutFactors >= 0 );
	
    //if file exists
    if (! GLOBAL_INIT_PARAMETER_MODEL_FILE_NAME.empty()){// ## qcheng75 -->
        readFreeParameterModelSource(GLOBAL_INIT_PARAMETER_MODEL_FILE_NAME, _nInFactors, _nOutFactors  );
        print();
        cout<<"\n"<<endl;
    }else{    // ## qcheng75 -->

		int nFactors = _nInFactors + _nOutFactors;
		for ( int i = 0; i < nFactors; i++ ) {
			maxBindingWts.push_back( BindingPar::initMaxBindingWt );
		}

		for ( int i = 0; i < _nInFactors; i++ ) {
			vector< double > ints( i + 1, BindingPar::initInt );
			inFactorIntMat.push_back( ints );
		}

		if ( _nOutFactors > 0 ) outFactorIntMat.setDimensions( _nOutFactors, _nInFactors );
		for ( int i = 0; i < _nOutFactors; i++ ) {
			for ( int j = 0; j < _nInFactors; j++ ) outFactorIntMat( i, j ) = BindingPar::initInt;
		}
		
		for ( int i = 0; i < _nInFactors; i++ ) expRatios.push_back( BindingPar::initExpRatio );

    }
}

/*
 * ## qcheng75
 */
BindingPar::BindingPar( int _nInFactors, int _nOutFactors, const vector< double >& _maxBindingWts )
{
	assert( _nInFactors > 0 && _nOutFactors >= 0 );

    //if file exists
    if (! GLOBAL_INIT_PARAMETER_MODEL_FILE_NAME.empty()){// ## qcheng75 -->
        readFreeParameterModelSource(GLOBAL_INIT_PARAMETER_MODEL_FILE_NAME, _nInFactors, _nOutFactors  );
        print();
        cout<<"\n"<<endl;
    }else{    // ## qcheng75 -->

		int nFactors = _nInFactors + _nOutFactors;
		for ( int i = 0; i < nFactors; i++ ) {
			maxBindingWts.push_back( _maxBindingWts[i] );
		}

		for ( int i = 0; i < _nInFactors; i++ ) {
			vector< double > ints( i + 1, BindingPar::initInt );
			inFactorIntMat.push_back( ints );
		}

		if ( _nOutFactors > 0 ) outFactorIntMat.setDimensions( _nOutFactors, _nInFactors );
		for ( int i = 0; i < _nOutFactors; i++ ) {
			for ( int j = 0; j < _nInFactors; j++ ) outFactorIntMat( i, j ) = BindingPar::initInt;
		}

		for ( int i = 0; i < _nInFactors; i++ ) expRatios.push_back( BindingPar::initExpRatio );

    }
}

BindingPar::BindingPar( const vector< double >& _maxBindingWts, const vector< vector< double > >& _inFactorIntMat, const Matrix& _outFactorIntMat, const vector< double >& _expRatios ) : maxBindingWts( _maxBindingWts ), inFactorIntMat( _inFactorIntMat ), outFactorIntMat( _outFactorIntMat ), expRatios( _expRatios )
{
	if ( !outFactorIntMat.isEmpty() ) assert( outFactorIntMat.nRows() == nOutFactors() ); 	
	if ( !outFactorIntMat.isEmpty() ) assert( outFactorIntMat.nCols() == nInFactors() ); 
	assert( expRatios.size() == nInFactors() );
}

BindingPar::BindingPar( const vector< double >& pars, const BindingPar& fixed, const vector< bool >& option )
{	
	assert( option.size() == BindingPar::nparts );
	int _nInFactors = fixed.nInFactors();
	int _nOutFactors = fixed.nOutFactors();
	int counter = 0;
	
	// set maxBindingWts of in-factors
	if ( option[ 0 ]  && (BindingPar::log10_min_wt<BindingPar::log10_max_wt)) { //qcheng75 add and condition
	   for ( int i = 0; i < _nInFactors; i++ ) maxBindingWts.push_back( pars[ counter++ ] );
	} else {
		for ( int i = 0; i < _nInFactors; i++ ) maxBindingWts.push_back( fixed.maxBindingWts[ i ] );
	}
	
	// set maxBindingWts of out-factors
	if ( option[ 1 ]  && (BindingPar::log10_min_wt<BindingPar::log10_max_wt)) { //qcheng75 add and condition
		for ( int i = 0; i < _nOutFactors; i++ ) maxBindingWts.push_back( pars[ counter++ ] );
	} else {
		for ( int i = 0; i < _nOutFactors; i++ ) maxBindingWts.push_back( fixed.maxBindingWts[ i + _nInFactors ] );
	}

    // set the dimensions of inFactorIntMat
    for ( int i = 0; i < _nInFactors; i++ ) {
        vector< double > ints( i + 1 );
        inFactorIntMat.push_back( ints );
    }

    // set the self-cooperativity terms
    if ( option[2] && (BindingPar::log10_min_int<BindingPar::log10_max_int)) { //qcheng75 add (BindingPar::log10_min_int<BindingPar::log10_max_int) for self cooperativity
        for ( int i = 0; i < _nInFactors; i++ ) {
            inFactorIntMat[i][i] = pars[counter++];
        }
    } else {
        for ( int i = 0; i < _nInFactors; i++ ) {
            inFactorIntMat[i][i] = fixed.inFactorIntMat[i][i];
        }       
    }
    
	// set the interaction terms between different factors
    if ( option[ 3 ]  && (BindingPar::log10_min_int<BindingPar::log10_max_int)) { //qcheng75 add (BindingPar::log10_min_int<BindingPar::log10_max_int) for self cooperativity
		for ( int i = 0; i < _nInFactors; i++ ) {
			for ( int j = 0; j < i; j++ ) { 
                inFactorIntMat[i][j] = pars[counter++];
            }
		}
		
		if ( _nOutFactors > 0 ) outFactorIntMat.setDimensions( _nOutFactors, _nInFactors );
		for ( int i = 0; i < _nOutFactors; i++ ) {
			for ( int j = 0; j < _nInFactors; j++ ) outFactorIntMat( i, j ) = pars[ counter++ ];	
		}		
	} else {
		for ( int i = 0; i < _nInFactors; i++ ) {
			for ( int j = 0; j < i; j++ ) { 
                inFactorIntMat[i][j] = fixed.inFactorIntMat[i][j];
            }
		}    
		
		if ( !fixed.outFactorIntMat.isEmpty() ) outFactorIntMat = fixed.outFactorIntMat;
	}
	
	// set the expRatios
    if ( option[ 4 ]  && (BindingPar::log10_min_ratio<BindingPar::log10_max_ratio)) { //qcheng75 add and condition
		for ( int i = 0; i < _nInFactors; i++ ) expRatios.push_back( pars[ counter++ ] );		
	} else {
		expRatios = fixed.expRatios; 
	}	
}

BindingPar::BindingPar( const vector< double >& pars, int nInFactors, int nOutFactors, const vector< bool >& option )
{
	assert( option.size() == BindingPar::nparts );
	int _nInFactors = nInFactors;
	int _nOutFactors = nOutFactors;
	int counter = 0;

	// set maxBindingWts of in-factors
	if ( option[ 0 ]  && (BindingPar::log10_min_wt<BindingPar::log10_max_wt)) { //qcheng75 add and condition
	   for ( int i = 0; i < _nInFactors; i++ ) maxBindingWts.push_back( pars[ counter++ ] );
	} else {
		for ( int i = 0; i < _nInFactors; i++ ) maxBindingWts.push_back( BindingPar::initMaxBindingWt );
	}

	// set maxBindingWts of out-factors
	if ( option[ 1 ]  && (BindingPar::log10_min_wt<BindingPar::log10_max_wt)) { //qcheng75 add and condition
		for ( int i = 0; i < _nOutFactors; i++ ) maxBindingWts.push_back( pars[ counter++ ] );
	} else {
		for ( int i = 0; i < _nOutFactors; i++ ) maxBindingWts.push_back( BindingPar::initMaxBindingWt );
	}

    // set the dimensions of inFactorIntMat
    for ( int i = 0; i < _nInFactors; i++ ) {
        vector< double > ints( i + 1 );
        inFactorIntMat.push_back( ints );
    }

    // set the self-cooperativity terms
    if ( option[2] && (BindingPar::log10_min_int<BindingPar::log10_max_int)) { //qcheng75 add (BindingPar::log10_min_int<BindingPar::log10_max_int) for self cooperativity
        for ( int i = 0; i < _nInFactors; i++ ) {
            inFactorIntMat[i][i] = pars[counter++];
        }
    } else {
        for ( int i = 0; i < _nInFactors; i++ ) {
            inFactorIntMat[i][i] = BindingPar::initInt;
        }
    }

	// set the interaction terms between different factors
    if ( option[ 3 ]  && (BindingPar::log10_min_int<BindingPar::log10_max_int)) { //qcheng75 add (BindingPar::log10_min_int<BindingPar::log10_max_int) for self cooperativity
		for ( int i = 0; i < _nInFactors; i++ ) {
			for ( int j = 0; j < i; j++ ) {
                inFactorIntMat[i][j] = pars[counter++];
            }
		}

		if ( _nOutFactors > 0 ) outFactorIntMat.setDimensions( _nOutFactors, _nInFactors );
		for ( int i = 0; i < _nOutFactors; i++ ) {
			for ( int j = 0; j < _nInFactors; j++ ) outFactorIntMat( i, j ) = pars[ counter++ ];
		}
	}

	// set the expRatios
    if ( option[ 4 ]  && (BindingPar::log10_min_ratio<BindingPar::log10_max_ratio)) { //qcheng75 add and condition
		for ( int i = 0; i < _nInFactors; i++ ) expRatios.push_back( pars[ counter++ ] );
	}
}
BindingPar::BindingPar( const vector< double >& pars, const BindingPar& fixed, const vector< bool >& option, bool isParBound  )
{
	assert( option.size() == BindingPar::nparts );
	int _nInFactors = fixed.nInFactors();
	int _nOutFactors = fixed.nOutFactors();
	int counter = 0;

	// set maxBindingWts of in-factors
	if ( option[ 0 ]  && (BindingPar::log10_min_wt<BindingPar::log10_max_wt)) { //qcheng75 add and condition
		double bw=0;
		for ( int i = 0; i < _nInFactors; i++ ) {
			if (pars[ counter ] < pow( 10.0,BindingPar::log10_min_wt) ){
				bw = pow( 10.0,BindingPar::log10_min_wt + BindingPar::min_float_value);
			}else if (pars[ counter ] > pow( 10.0,BindingPar::log10_max_wt) ){
				bw = pow( 10.0,BindingPar::log10_max_wt - BindingPar::min_float_value);
			}else{
				bw = pars[ counter ];
			}
			counter++;
			maxBindingWts.push_back( bw  );
		}
	} else {
		for ( int i = 0; i < _nInFactors; i++ ) maxBindingWts.push_back( fixed.maxBindingWts[ i ] );
	}

	// set maxBindingWts of out-factors
	if ( option[ 1 ]  && (BindingPar::log10_min_wt<BindingPar::log10_max_wt)) { //qcheng75 add and condition
		double bw=0;
		for ( int i = 0; i < _nOutFactors; i++ ) {
			if (pars[ counter ] < pow( 10.0,BindingPar::log10_min_wt) ){
				bw = pow( 10.0,BindingPar::log10_min_wt + BindingPar::min_float_value);
			}else if (pars[ counter ] > pow( 10.0,BindingPar::log10_max_wt ) ){
				bw = pow( 10.0,BindingPar::log10_max_wt - BindingPar::min_float_value);
			}else{
				bw = pars[ counter ];
			}
			counter++;
			maxBindingWts.push_back( bw  );
		}
	} else {
		for ( int i = 0; i < _nOutFactors; i++ ) maxBindingWts.push_back( fixed.maxBindingWts[ i + _nInFactors ] );
	}

    // set the dimensions of inFactorIntMat
    for ( int i = 0; i < _nInFactors; i++ ) {
        vector< double > ints( i + 1 );
        inFactorIntMat.push_back( ints );
    }

    // set the self-cooperativity terms
    if ( option[2] && (BindingPar::log10_min_int<BindingPar::log10_max_int)) { //qcheng75 add (BindingPar::log10_min_int<BindingPar::log10_max_int) for self cooperativity

		double cw=0;
		for ( int i = 0; i < _nInFactors; i++ ) {
			if (pars[ counter ] < pow( 10.0,BindingPar::log10_min_int) ){
				cw = pow( 10.0,BindingPar::log10_min_int + BindingPar::min_float_value);
			}else if (pars[ counter ] > pow( 10.0,BindingPar::log10_max_int) ){
				cw = pow( 10.0,BindingPar::log10_max_int - BindingPar::min_float_value);
			}else{
				cw = pars[ counter ];
			}
			counter++;
			inFactorIntMat[i][i] = cw;
		}
    } else {
        for ( int i = 0; i < _nInFactors; i++ ) {
            inFactorIntMat[i][i] = fixed.inFactorIntMat[i][i];
        }
    }

	// set the interaction terms between different factors
    if ( option[ 3 ]  && (BindingPar::log10_min_int<BindingPar::log10_max_int)) { //qcheng75 add (BindingPar::log10_min_int<BindingPar::log10_max_int) for self cooperativity
    	double cw=0;
    	for ( int i = 0; i < _nInFactors; i++ ) {
			for ( int j = 0; j < i; j++ ) {
				if (pars[ counter ] < pow( 10.0,BindingPar::log10_min_int) ){
					cw = pow( 10.0,BindingPar::log10_min_int + BindingPar::min_float_value);
				}else if (pars[ counter ] > pow( 10.0,BindingPar::log10_max_int) ){
					cw = pow( 10.0,BindingPar::log10_max_int - BindingPar::min_float_value);
				}else{
					cw = pars[ counter ];
				}
				counter++;
				inFactorIntMat[i][j] = cw;
			}
    	}

    	cout << "<_nOutFactors, _nInFactors>=<" << _nOutFactors << "," << _nInFactors << ">"<< endl;
		if ( _nOutFactors > 0 ) outFactorIntMat.setDimensions( _nOutFactors, _nInFactors );

		cw=0;
		for ( int i = 0; i < _nOutFactors; i++ ) {
			for ( int j = 0; j < _nInFactors; j++ ) {
				if (pars[ counter ] < pow( 10.0,BindingPar::log10_min_int) ){
					cw = pow( 10.0,BindingPar::log10_min_int + BindingPar::min_float_value);
				}else if (pars[ counter ] > pow( 10.0,BindingPar::log10_max_int) ){
					cw = pow( 10.0,BindingPar::log10_max_int - BindingPar::min_float_value);
				}else{
					cw = pars[ counter ];
				}
				counter++;
				outFactorIntMat( i, j ) = cw;
			}
		}
	} else {
		for ( int i = 0; i < _nInFactors; i++ ) {
			for ( int j = 0; j < i; j++ ) {
                inFactorIntMat[i][j] = fixed.inFactorIntMat[i][j];
            }
		}

		if ( !fixed.outFactorIntMat.isEmpty() ) outFactorIntMat = fixed.outFactorIntMat;
	}

	// set the expRatios
    //LEFT
    if ( option[ 4 ]  && (BindingPar::log10_min_ratio<BindingPar::log10_max_ratio)) { //qcheng75 add and condition
		for ( int i = 0; i < _nInFactors; i++ ) expRatios.push_back( pars[ counter++ ] );
	} else {
		expRatios = fixed.expRatios;
	}
}
// ## qcheng75
void BindingPar::readFreeParameterModelSource(string para_model_file, int _nInFactors, int _nOutFactors ) {

    ifstream fdata( para_model_file.c_str() );
    if (fdata){
    	cout<< "[BindingPar::readFreeParameterModelSource] _nInFactors=" << _nInFactors<< " _nInFactors=" << _nOutFactors<<endl;
        readFreeParameterModelSource(fdata, _nInFactors, _nOutFactors  );
        fdata.close();
    }
}

// ## qcheng75
void BindingPar::readFreeParameterModelSource(std::ifstream& _fdata, int _nInFactors, int _nOutFactors ){
	string lineS;
	int nFactors = _nInFactors + _nOutFactors;
	int oldInFactor=0,oldOutFactor=0;
	//cout<<"Read Interaction Factors: "<<_nInFactors << " " << _nOutFactors<<endl;
	// read parameters
	while ( getline( _fdata, lineS ) ) {
		//const char* line=lineS.data();
		size_t ptr = lineS.find('=');
		string tagName = lineS.substr(0, int(ptr));
		if (tagName.compare("maxBindingWts")==0){
			string factorBindingWtsList=lineS.substr(int(ptr)+1, lineS.length()-int(ptr)-1);
			char* factorWts = strtok ((char *)(factorBindingWtsList.data())," \t");
			int wtsCount=0;
			while (factorWts != NULL){
				maxBindingWts.push_back(atof(factorWts));
			    factorWts = strtok (NULL, " \t");
			    wtsCount++;
			}
			//cout<<"Read maxBindingWts: "<<wtsCount << " " << nFactors<<endl;
			while (wtsCount < nFactors) {
				maxBindingWts.push_back( BindingPar::initMaxBindingWt );
				wtsCount++;
			}
		}else if (tagName.compare("expRatios")==0){
			string expRatio=lineS.substr(int(ptr)+1, lineS.length()-int(ptr)-1);
			char* factorExpRatio = strtok ((char *)(expRatio.data())," \t");
			int expRCount=0;
			while (factorExpRatio != NULL){
				expRatios.push_back(atof(factorExpRatio));
				factorExpRatio = strtok (NULL, " \t");
			    expRCount++;
			}
			//cout<<"Read factorExpRatio: "<<expRCount << " " << _nInFactors<<endl;
			while (expRCount < _nInFactors) {
				expRatios.push_back( BindingPar::initExpRatio );
				expRCount++;
			}
			//for ( int i = 0; i < _nInFactors; i++ ) expRatios.push_back( BindingPar::initExpRatio );

		}else if (tagName.compare("numOfInFactors")==0){
			string numOfInFactors=lineS.substr(int(ptr)+1, lineS.length()-int(ptr)-1);
			oldInFactor=atoi(numOfInFactors.data());
		}else if (tagName.compare("inFactorIntMat")==0){
			int i = 0;
			while ((oldInFactor>0) && (i < oldInFactor )){
				 getline( _fdata, lineS );
				 vector< double > ints;
				 char* f_fIntWts = strtok ((char *)(lineS.data())," \t");
				 int j=0;
				 while (f_fIntWts != NULL){
					 ints.push_back(atof(f_fIntWts));
					 f_fIntWts = strtok (NULL, " \t");
					 j+=1;
				 }
				 while (j<i) {
					 ints.push_back(BindingPar::initInt);
					 j++;
				 }
				 inFactorIntMat.push_back(ints);
				 i+=1;
			}
			//cout<<"Read inFactorIntMat: " << i << " " << _nInFactors<<endl;
			while (i<_nInFactors){
				vector< double > ints( i + 1, BindingPar::initInt );
				inFactorIntMat.push_back( ints );
				i++;
			}

		}else if (tagName.compare("numOfOutFactors")==0){
			string numOfOutFactors=lineS.substr(int(ptr)+1, lineS.length()-int(ptr)-1);
			oldOutFactor=atoi(numOfOutFactors.data());
		}else if (tagName.compare("outFactorIntMat")==0){
			int i = 0;
			//cout<<"Read outFactorIntMat: " << i << " " << oldOutFactor<< " " << _nOutFactors<<endl;
			if ( _nOutFactors > 0 ){
				outFactorIntMat.setDimensions( _nOutFactors, _nInFactors );

				while ((oldOutFactor>0) && (i < oldOutFactor )){
					 getline( _fdata, lineS );
					 char* f_fIntWts = strtok ((char *)(lineS.data())," \t");
					 int j=0;
					 while (f_fIntWts != NULL){
						 outFactorIntMat(i,j)=atof(f_fIntWts);
						 f_fIntWts = strtok (NULL, " \t");
						 j+=1;
					 }
					 while (j<_nInFactors)  {
						 outFactorIntMat( i, j ) = BindingPar::initInt;
						 j++;
					 }
					 i+=1;
				}
				//cout<<"Read outFactorIntMat: " << i << " " << oldOutFactor<< " " << _nOutFactors<<endl;
				while (i<_nOutFactors){
					for ( int j = 0; j < _nInFactors; j++ ) outFactorIntMat( i, j ) = BindingPar::initInt;
					i++;
				}
			}
		}
	}
}

// BindingPar* BindingPar::copy_construct() const
// {
// 	return new BindingPar( maxBindingWts, inFactorIntMat, outFactorIntMat, expRatios );		
// }

// void BindingPar::copy( const BindingPar* src )
// {
// 	maxBindingWts = src->maxBindingWts;
// 	inFactorIntMat = src->inFactorIntMat;
// 	outFactorIntMat = src->outFactorIntMat;
// 	expRatios = src->expRatios;
// }

void BindingPar::clear(){//qcheng75
	if (! maxBindingWts.empty()) maxBindingWts.clear();
	if (! expRatios.empty()) expRatios.clear();
	outFactorIntMat.freeMatrix();
	if ( outFactorIntMat.getData()!= NULL) delete outFactorIntMat.getData();
	for (int i=0;i<inFactorIntMat.size();i++){
		inFactorIntMat[i].clear();
		//delete &(inFactorIntMat[i]);
	}
	inFactorIntMat.clear();
}

void BindingPar::createFactorIntMat( Matrix& factorIntMat ) const
{
	int nFactors = nInFactors() + nOutFactors();
	factorIntMat.setDimensions( nFactors, nFactors );
	
	for ( int i = 0; i < nFactors; i++ ) {
		for ( int j = 0; j < nFactors; j++ ) {
			if ( i < nInFactors() ) {
				// in-in
				if ( j < nInFactors() ) factorIntMat( i, j ) = ( i >= j ) ? inFactorIntMat[ i ][ j ] : inFactorIntMat[ j ][ i ];
				// in-out
				else factorIntMat( i, j ) = outFactorIntMat( j - nInFactors(), i );	
			} else {
				// out-in
				if ( j < nInFactors() ) factorIntMat( i, j ) = outFactorIntMat( i - nInFactors(), j );
				// out-out
				else factorIntMat( i, j ) = 1;
			}
		}	
	}	
}

void BindingPar::getFreePars( vector< double >& pars, const vector< bool >& option ) const
{
	assert( option.size() == BindingPar::nparts );
	pars.clear();
		
	// write maxBindingWts of in-factors
	if ( option[ 0 ] && (BindingPar::log10_min_wt<BindingPar::log10_max_wt)) { //qcheng75
		for ( int i = 0; i < nInFactors(); i++ ) pars.push_back( maxBindingWts[ i ] );
	} 
	
	// write maxBindingWts of out-factors
	if ( option[ 1 ] && (BindingPar::log10_min_wt<BindingPar::log10_max_wt)) {//qcheng75
		for ( int i = 0; i < nOutFactors(); i++ ) pars.push_back( maxBindingWts[ i + nInFactors() ] );
	} 

    // write self-cooperativity terms
	if ( option[2] && (BindingPar::log10_min_int<BindingPar::log10_max_int) ) {//qcheng75
        for ( int i = 0; i < nInFactors(); i++ ) pars.push_back( inFactorIntMat[i][i] );
    }
    
	// write the interaction terms between different factors
	if ( option[3] && (BindingPar::log10_min_int<BindingPar::log10_max_int) ) {//qcheng75
		for ( int i = 0; i < nInFactors(); i++ ) {
			for ( int j = 0; j < i; j++ ) pars.push_back( inFactorIntMat[ i ][ j ] );
		}
		
		for ( int i = 0; i < nOutFactors(); i++ ) {
			for ( int j = 0; j < nInFactors(); j++ ) pars.push_back( outFactorIntMat( i, j ) );
		}		
	} 
	
	// write the expRatios
	if ( option[ 4 ] && (BindingPar::log10_min_ratio<BindingPar::log10_max_ratio) ) {//qcheng75
		for ( int i = 0; i < nInFactors(); i++ ) pars.push_back( expRatios[ i ] );	 
	}	
}

// get the box-boundary parameter constraints qcheng75
void BindingPar::getParsLB( vector< double >& parsLB, const vector< bool >& option ) const
{
	assert( option.size() == BindingPar::nparts );
	parsLB.clear();

	// write maxBindingWts of in-factors
	if ( option[ 0 ] ) {
		for ( int i = 0; i < nInFactors(); i++ ) parsLB.push_back( BindingPar::log10_min_wt );
	}

	// write maxBindingWts of out-factors
	if ( option[ 1 ] ) {
		for ( int i = 0; i < nOutFactors(); i++ ) parsLB.push_back( BindingPar::log10_min_wt );
	}

    // write self-cooperativity terms
	if ( option[2]  ) {
        for ( int i = 0; i < nInFactors(); i++ ) parsLB.push_back( BindingPar::log10_min_int );
    }

	// write the interaction terms between different factors
	if ( option[3] ) {
		for ( int i = 0; i < nInFactors(); i++ ) {
			for ( int j = 0; j < i; j++ ) parsLB.push_back( BindingPar::log10_min_int  );
		}

		for ( int i = 0; i < nOutFactors(); i++ ) {
			for ( int j = 0; j < nInFactors(); j++ ) parsLB.push_back( BindingPar::log10_min_int  );
		}
	}

	// write the expRatios
	if ( option[ 4 ] ){
		for ( int i = 0; i < nInFactors(); i++ ) parsLB.push_back( BindingPar::log10_min_ratio );
	}
}

// get the box-boundary parameter constraints qcheng75
void BindingPar::getParsUB( vector< double >& parsUB, const vector< bool >& option ) const
{
	assert( option.size() == BindingPar::nparts );
	parsUB.clear();

	// write maxBindingWts of in-factors
	if ( option[ 0 ] ) {
		for ( int i = 0; i < nInFactors(); i++ ) parsUB.push_back( BindingPar::log10_max_wt );
	}

	// write maxBindingWts of out-factors
	if ( option[ 1 ] ) {
		for ( int i = 0; i < nOutFactors(); i++ ) parsUB.push_back( BindingPar::log10_max_wt );
	}

    // write self-cooperativity terms
	if ( option[2]  ) {
        for ( int i = 0; i < nInFactors(); i++ ) parsUB.push_back( BindingPar::log10_max_int );
    }

	// write the interaction terms between different factors
	if ( option[3] ) {
		for ( int i = 0; i < nInFactors(); i++ ) {
			for ( int j = 0; j < i; j++ ) parsUB.push_back( BindingPar::log10_max_int  );
		}

		for ( int i = 0; i < nOutFactors(); i++ ) {
			for ( int j = 0; j < nInFactors(); j++ ) parsUB.push_back( BindingPar::log10_max_int  );
		}
	}

	// write the expRatios
	if ( option[ 4 ] ){
		for ( int i = 0; i < nInFactors(); i++ ) parsUB.push_back( BindingPar::log10_max_ratio );
	}
}

void BindingPar::print() const
{
	cout << "maxBindingWts = " << maxBindingWts << endl;
// 	cout << "inFactorIntMat = " << endl;
// 	for ( int i = 0; i < inFactorIntMat.size(); i++ ) {
// 		for ( int j = 0; j <= i; j++ ) {
// 			cout << inFactorIntMat[ i ][ j ];
// 			if ( j != i ) cout << " ";
// 		}	
// 		cout << endl;
// 	}
// 	cout << "outFactorIntMat = " << endl << outFactorIntMat;	

	cout << "numOfInFactors=" <<nInFactors()<< endl;
	cout << "inFactorIntMat=" << endl;
	for ( int i = 0; i < inFactorIntMat.size(); i++ ) {
		for ( int j = 0; j <= i; j++ ) {
			cout << inFactorIntMat[ i ][ j ];
			if ( j != i ) cout << " ";
		}
		cout << endl;
	}
	cout << "numOfOutFactors=" <<nOutFactors()<< endl;
	cout << "outFactorIntMat=" <<endl << outFactorIntMat;
	cout << "expRatios = " << expRatios << endl; 
}

//qcheng75
void BindingPar::printParaModel2Stream(std::ofstream& fout) const
{
	fout << "maxBindingWts=" << maxBindingWts << endl;

	//createFactorIntMat( m );
	//fout << "factorIntMat=" << endl << m;
	fout << "numOfInFactors=" <<nInFactors()<< endl;
 	fout << "inFactorIntMat=" << endl;
 	for ( int i = 0; i < inFactorIntMat.size(); i++ ) {
 		//fout << "*";
 		for ( int j = 0; j <= i; j++ ) {
 			fout << inFactorIntMat[ i ][ j ];
 			if ( j != i ) fout << " ";
 		}
 		fout << endl;
 	}
	fout << "numOfOutFactors=" <<nOutFactors()<< endl;
 	fout << "outFactorIntMat=" <<endl << outFactorIntMat;
	fout << "expRatios=" << expRatios << endl;
}

int BindingPar::nparts = 5;
double BindingPar::initMaxBindingWt = 1.0;
double BindingPar::initInt = 1.0;
double BindingPar::initExpRatio = 1.0;
double BindingPar::log10_min_wt = -4;		
double BindingPar::log10_max_wt = 4;		
double BindingPar::log10_min_int = -3;
double BindingPar::log10_max_int = 3;
double BindingPar::log10_min_ratio = -3;
double BindingPar::log10_max_ratio = 3;
double BindingPar::log10_range_wt = 2;
double BindingPar::wt_step = 0.2;
double BindingPar::int_step = 0.2;
double BindingPar::ratio_step = 0.2;

double BindingPar::min_float_value = 1.0E-6;
double BindingPar::max_float_value = 1.0E+6;

int BindingPar::FIX_COOP_PROB = 0; //## qcheng75
int BindingPar::NELDER_MEAD_CONSTRAINED = 0;
int BindingPar::NEWTON_KTCONSTRAINED = 0;
int BindingPar::PENALTY_BASED_CONSTRAINED = 1;
double BindingPar::PENALTY_COEFFIENT = 1000.0;
int BindingPar::IS_LOG_SCALE_STEP = 0;

// void BindingPar::assignParams( const vector< double >& pars, int nInFactors, int nOutFactors )
// {
// 	int counter = 0;
// 	int nFactors = nInFactors + nOutFactors; 
// 	int nParams = nFactors + nInFactors * ( nInFactors + 1 ) / 2 + nOutFactors * nInFactors + nInFactors; 
// 	assert ( pars.size() == nParams ); 

// 	// clear the current parameters	
// 	maxBindingWts.clear();
// 	inFactorIntMat.clear();	
// 	expRatios.clear();
// 		
// 	// read maxBindingWts
// 	for ( int i = 0; i < nFactors; i++ ) maxBindingWts.push_back( pars[ counter++ ] );
// 	
// 	// read inFactorIntMat
// 	for ( int i = 0; i < nInFactors; i++ ) {
// 		vector< double > ints;
// 		for ( int j = 0; j <= i; j++ ) { ints.push_back( pars[ counter++ ] ); }
// 		inFactorIntMat.push_back( ints );	
// 	}
// 	
// 	// read outFactorIntMat
// 	if ( nOutFactors > 0 ) outFactorIntMat.setDimensions( nOutFactors, nInFactors );
// 	for ( int i = 0; i < nOutFactors; i++ ) {
// 		for ( int j = 0; j < nInFactors; j++ ) outFactorIntMat( i, j ) = pars[ counter++ ];	
// 	}
// 	
// 	// read expRatios
// 	for ( int i = 0; i < nInFactors; i++ ) expRatios.push_back( pars[ counter++ ] );	
// }

BasePredictor::BasePredictor( const vector< Motif >& _motifs, const vector< double >& _energyThrs, const vector< vector< Sequence > >& _seqs, const vector< vector< double > >& _bindingData ) : motifs( _motifs ), energyThrs( _energyThrs ), seqs( _seqs ), bindingData( _bindingData ), par_model( seqs.size(), motifs.size() - seqs.size() )
{
	assert( motifs.size() > 0 );
	assert( energyThrs.size() == motifs.size() );
	assert( bindingData.size() == seqs.size() );	
}

/*
 * ## qcheng75
 */
BasePredictor::BasePredictor( const vector< Motif >& _motifs, const vector< double >& _energyThrs, const vector< vector< Sequence > >& _seqs, const vector< vector< double > >& _bindingData, const vector< double >& _maxBindingWts  ) : motifs( _motifs ), energyThrs( _energyThrs ), seqs( _seqs ), bindingData( _bindingData ), par_model( seqs.size(), motifs.size() - seqs.size(), _maxBindingWts )
{
	assert( motifs.size() > 0 );
	assert( energyThrs.size() == motifs.size() );
	assert( bindingData.size() == seqs.size() );
}

/*
 * ## qcheng75
 */
BasePredictor::BasePredictor( const vector< Motif >& _motifs, const vector< double >& _energyThrs, const vector< vector< Sequence > >& _seqs, const vector< vector< double > >& _bindingData, const vector< double >& _allWts, int nInFactors, int nOutFactors, bool isAllParam  ) : motifs( _motifs ), energyThrs( _energyThrs ), seqs( _seqs ), bindingData( _bindingData ), par_model( _allWts,nInFactors, nOutFactors, BindingPredictor::parOption)
{
	assert( motifs.size() > 0 );
	assert( energyThrs.size() == motifs.size() );
	assert( bindingData.size() == seqs.size() );
}

BasePredictor::~BasePredictor(){//qcheng75

}
void BasePredictor::clear(){//qcheng75
	if (! motifs.empty()) {
		for (int i=0; i< motifs.size();i++){
			motifs[i].clear();
		}
	}
	motifs.clear();
	if (! energyThrs.empty()) energyThrs.clear();

	par_model.clear();

	/*
	if (! this->bAccData.empty()) {
		for (int i=0;i<bAccData.size();i++){
			bAccData[i].clear();
			//delete &(bAccData[i]);
		}
		bAccData.clear();
	}
	*/
}

double BasePredictor::test( const vector< vector< Sequence > >& testSeqs, const vector< vector< double > >& testBindingData, int perfOption ) const
{
	double ret=0;
	assert( perfOption == 0 );	
	int nExps = seqs.size();
	assert( testSeqs.size() == nExps && testBindingData.size() == nExps ); 
	
	// make predictions
	vector< vector< double > > predicted( nExps );
	for ( int i = 0; i < testSeqs.size(); i++ ) {
		for ( int j = 0; j < testSeqs[ i ].size(); j++ ) {
			assert( testBindingData[ i ].size() == testSeqs[ i ].size() );
			predicted[ i ].push_back( predict( testSeqs[ i ][ j ], i ) );
		} 	
	}
	
	// Pearson correlation between predictions and observations
	if ( perfOption == 0 ) {
		vector< double > corrs; 
		for ( int i = 0; i < nExps; i++ ) {
			corrs.push_back( correlation( predicted[ i ], testBindingData[ i ] ) );
		}
		
		ret=mean( corrs );
		corrs.clear();
	}

	for (int i=0; i<predicted.size(); i++){//qcheng75
		predicted[i].clear();
		//delete &(predicted[i]);
	}
	predicted.clear();//qcheng75
	return ret;
}

int BindingPredictor::GLOBAL_NRANDSTARTS; // ## SS
BindingPredictor::BindingPredictor( const vector< Motif >& _motifs, const vector< double >& _energyThrs, const vector< vector< Sequence > >& _seqs, const vector< vector< double > >& _bindingData, const FactorIntFunc* _intFunc ) : BasePredictor( _motifs, _energyThrs, _seqs, _bindingData ), intFunc( _intFunc ), seqSites( seqs.size() )
{
	// sequence annotation
	SeqAnnotator ann( motifs, energyThrs );
	for ( int i = 0; i < seqs.size(); i++ ) {
		for ( int j = 0; j < seqs[ i ].size(); j++ ) {
			SiteVec sites;
			ann.annot( seqs[ i ][ j ], sites );
			seqSites[ i ].push_back( sites );
			sites.clear();//qcheng75
		}
	}	
	ann.clear();//qcheng75
}

/*
 * ## qcheng75
 */
BindingPredictor::BindingPredictor( const vector< Motif >& _motifs, const vector< double >& _energyThrs, const vector< vector< Sequence > >& _seqs, const vector< vector< double > >& _bindingData, const FactorIntFunc* _intFunc, const vector< double >& _maxBindingWts  ) : BasePredictor( _motifs, _energyThrs, _seqs, _bindingData, _maxBindingWts  ), intFunc( _intFunc ), seqSites( seqs.size() )
{
	// sequence annotation
	SeqAnnotator ann( motifs, energyThrs );
	for ( int i = 0; i < seqs.size(); i++ ) {
		for ( int j = 0; j < seqs[ i ].size(); j++ ) {
			SiteVec sites;
			ann.annot( seqs[ i ][ j ], sites );

			seqSites[ i ].push_back( sites );
			sites.clear();//qcheng75
		}
	}
	ann.clear();//qcheng75
}

/*
 * ## qcheng75
 */
BindingPredictor::BindingPredictor( const vector< Motif >& _motifs, const vector< double >& _energyThrs, const vector< vector< Sequence > >& _seqs, const vector< vector< double > >& _bindingData, const FactorIntFunc* _intFunc, const vector< double >& _allWts, int nInFactors, int nOutFactors, bool isAllParam  ) : BasePredictor( _motifs, _energyThrs, _seqs, _bindingData, _allWts,nInFactors, nOutFactors, isAllParam ), intFunc( _intFunc ), seqSites( seqs.size() )
{
	// sequence annotation
	SeqAnnotator ann( motifs, energyThrs );
	for ( int i = 0; i < seqs.size(); i++ ) {
		for ( int j = 0; j < seqs[ i ].size(); j++ ) {
			SiteVec sites;
			ann.annot( seqs[ i ][ j ], sites );

			seqSites[ i ].push_back( sites );
			sites.clear();//qcheng75
		}
	}
	ann.clear();//qcheng75
}

void BindingPredictor::clear(){ //qcheng75
	 cout << "clean" <<endl;
	for (int i=0; i<seqSites.size(); i++){
		seqSites[i].clear();
	}
	seqSites.clear();		// the extracted sites for all sequences

	BasePredictor::clear();
}

BindingPredictor::~BindingPredictor(){
	//delete intFunc;
}

// get the box-boundary parameter constraints qcheng75
void BindingPredictor::getParBoundaryConstraints( vector< double >& parsLB,  vector< double >& parsUB ) const
{
	assert( parOption.size() == BindingPar::nparts );
	parsLB.clear();

	// write maxBindingWts of in-factors
	if ( parOption[ 0 ] ) {
		for ( int i = 0; i < nInFactors(); i++ ) {
			parsLB.push_back( BindingPar::IS_LOG_SCALE_STEP==1?BindingPar::log10_min_wt: pow(double(10.0),BindingPar::log10_min_wt));
			parsUB.push_back( BindingPar::IS_LOG_SCALE_STEP==1?BindingPar::log10_max_wt: pow(double(10.0),BindingPar::log10_max_wt));
		}
	}

	// write maxBindingWts of out-factors
	if ( parOption[ 1 ] ) {
		for ( int i = 0; i < nOutFactors(); i++ ) {
			parsLB.push_back( BindingPar::IS_LOG_SCALE_STEP==1?BindingPar::log10_min_wt: pow(double(10.0),BindingPar::log10_min_wt));
			parsUB.push_back( BindingPar::IS_LOG_SCALE_STEP==1?BindingPar::log10_max_wt: pow(double(10.0),BindingPar::log10_max_wt));
		}
	}

    // write self-cooperativity terms
	if ( parOption[2]  ) {
        for ( int i = 0; i < nInFactors(); i++ ) {
        	parsLB.push_back( BindingPar::IS_LOG_SCALE_STEP==1?BindingPar::log10_min_int: pow(double(10.0),BindingPar::log10_min_int));
        	parsUB.push_back( BindingPar::IS_LOG_SCALE_STEP==1?BindingPar::log10_max_int: pow(double(10.0),BindingPar::log10_max_int));
        }
    }

	// write the interaction terms between different factors
	if ( parOption[3] ) {
		for ( int i = 0; i < nInFactors(); i++ ) {
			for ( int j = 0; j < i; j++ ) {
				parsLB.push_back( BindingPar::IS_LOG_SCALE_STEP==1?BindingPar::log10_min_int: pow(double(10.0),BindingPar::log10_min_int));
				        	parsUB.push_back( BindingPar::IS_LOG_SCALE_STEP==1?BindingPar::log10_max_int: pow(double(10.0),BindingPar::log10_max_int));
			}
		}

		for ( int i = 0; i < nOutFactors(); i++ ) {
			for ( int j = 0; j < nInFactors(); j++ ) {
				parsLB.push_back( BindingPar::IS_LOG_SCALE_STEP==1?BindingPar::log10_min_int: pow(double(10.0),BindingPar::log10_min_int));
				        	parsUB.push_back( BindingPar::IS_LOG_SCALE_STEP==1?BindingPar::log10_max_int: pow(double(10.0),BindingPar::log10_max_int));
			}
		}
	}

	// write the expRatios
	if ( parOption[ 4 ] ){
		for ( int i = 0; i < nInFactors(); i++ ) {
			parsLB.push_back( BindingPar::IS_LOG_SCALE_STEP==1?BindingPar::log10_min_ratio: pow(double(10.0),BindingPar::log10_min_ratio));
			parsUB.push_back( BindingPar::IS_LOG_SCALE_STEP==1?BindingPar::log10_max_ratio: pow(double(10.0),BindingPar::log10_max_ratio));
		}
	}
}

double BindingPredictor::objFunc( const BindingPar& par ) const
{
	return compCorr( par );
}

int BindingPredictor::train( const BindingPar& par_init, int nAlternations )
{
	if ( nAlternations == 0 ) return 0;
	
	// parameter fitting
	par_model = par_init;
	double obj_result;
	BindingPar par_result;
// 	simplex_minimize( par_result, obj_result );
// 	gradient_minimize( par_result, obj_result );
// 	SA_minimize( par_result, obj_result );
	
	// alternate between two different methods
	for ( int i = 0; i < nAlternations; i++ ) {
		if ((BindingPar::FIX_COOP_PROB>0) && (BindingPar::NEWTON_KTCONSTRAINED)){
			gradient_minimize_KTconstrained( par_result, obj_result );
			par_model = par_result;
		}else{
			gradient_minimize( par_result, obj_result );
			par_model = par_result;

			simplex_minimize( par_result, obj_result );
			par_model = par_result;
		}
	}
	
	// commit the parameters and the value of the objective function
	par_model = par_result; 
	obj_model = obj_result;
		
	par_result.clear(); //qcheng75
	return 0;	
}

int BindingPredictor::train()
{	
	// training using the default initial values
	int nAlternationsShallow = 1;	
	BindingPar par_default( nInFactors(), nOutFactors() );
	train( par_default, nAlternationsShallow );		
	
	// training with random starts
	gsl_rng* rng;
	gsl_rng_env_setup();
	const gsl_rng_type * T = gsl_rng_default;	// create rng type
	rng = gsl_rng_alloc( T );
	gsl_rng_set( rng, time( 0 ) );		// set the seed equal to simulTime(0)
		
	int nRandStarts = GLOBAL_NRANDSTARTS; // 0; ## SS
	BindingPar par_init = par_default, par_best = par_model;
	double obj_best = obj_model;
	for ( int i = 0; i < nRandStarts; i++ ) {
		cout << endl << "Random start " << i + 1 << endl;
		randSamplePar( rng, par_init ); 
		train( par_init, nAlternationsShallow );
		if ( obj_model > obj_best ) {
			par_best = par_model;
			obj_best = obj_model;	
		}
	}	
	
// 	par_model = par_best;
// 	obj_model = obj_best;
				
	// deep search
	int nAlternationsDeep = 0;
// 	cout << endl << "Deep search " << endl;
	train( par_best, nAlternationsDeep );
	if ((nAlternationsDeep>0) && (obj_model > obj_best )) {//qcheng75
		par_best = par_model;
		obj_best = obj_model;
		cout<<"--get better obj_model="<< obj_model <<endl;
	}
	gsl_rng_free(rng);//qcheng75
	cout<<"finish training "<<endl;

	par_init.clear(); //qcheng75
	par_best.clear();
	par_default.clear();
	return 0;	
}

int BindingPredictor::analyzeFeatures( vector< vector< double > >& featureEffects ) const
{
    featureEffects = vector< vector< double > >( nOutFactors(), nInFactors() );
    
// 	if ( nInFactors() != 1 ) { cerr << "analyzeFeatures() is only applicable when nInFactors = 1" << endl; exit( 1 ); }

    // parameters of in-factors
    BindingPar par_in( nInFactors(), nOutFactors() );
    for ( int k = 0; k < nInFactors(); k++ ) par_in.maxBindingWts[ k ] = par_model.maxBindingWts[ k ];
    par_in.inFactorIntMat = par_model.inFactorIntMat;
    par_in.expRatios = par_model.expRatios;
    BindingFunc* func_in = createBindingFunc( par_in );

    // parameters of out-factors
    vector< BindingPar > par_outs( nOutFactors(), par_in );
    for ( int k = 0; k < nOutFactors(); k++ ) {
        par_outs[ k ].maxBindingWts[ k + nInFactors() ] = par_model.maxBindingWts[ k + nInFactors() ];
        par_outs[ k ].outFactorIntMat.setRow( k, par_model.outFactorIntMat.getRow( k ) );
    }
    
    for ( int i = 0; i < nInFactors(); i++ ) {
    // 	cout << "Seq\tReal\tPredicted\tPredicted_in\tPredicted_out\tEffect" << endl;
        vector< vector< double > > effects( nOutFactors() );
        for ( int j = 0; j < seqs[ i ].size(); j++ ) {
    // 		cout << j << "\t" << bindingData[ 0 ][ j ] << "\t" << predict( seqs[ 0 ][ j ], 0 ) << "\t";
            
            // binding if only the in-factor is involved
            double predicted_in = func_in->predictBinding( seqSites[ i ][ j ], i );
    // 		cout << predicted_in; 
                
            // binding when each of the out-factor is involved
            for ( int k = 0; k < nOutFactors(); k++ ) {
                BindingFunc* func_out = createBindingFunc( par_outs[ k ] );	
                double predicted_out = func_out->predictBinding( seqSites[ i ][ j ], i );
                double effect = predicted_in == 0 ? 0 : ( predicted_out - predicted_in ) / predicted_in;
    // 			cout << "\t" << predicted_out << "\t" << effect;
                effects[ k ].push_back( effect );
            }
        }
// 		cout << endl;	

        // mean effect of each out-factor in this experiment
        for ( int k = 0; k < nOutFactors(); k++ ) {
            featureEffects[ k ][ i ] = mean( effects[ k ] );	
        }
    }
    
	return 0;
}

double BindingPredictor::predict( const Sequence& targetSeq, int factorIdx ) const
{	
	assert( factorIdx >= 0 && factorIdx < nInFactors() );
	
	// create site representation of the target sequence
	SiteVec targetSites;
	SeqAnnotator ann( motifs, energyThrs );	
	ann.annot( targetSeq, targetSites );
		
	// predict the binding
	BindingFunc* func = createBindingFunc( par_model );	
	double predicted = func->predictBinding( targetSites, factorIdx );

	func->clear();//qcheng75
	func->releaseFactorIntMat();//qcheng75
	delete func;//qcheng75
	ann.clear();//qcheng75

	return predicted;
}

const bool parOptionArray[] = { 1, 1, 1, 1, 0 }; 
vector< bool > BindingPredictor::parOption( parOptionArray, parOptionArray + sizeof( parOptionArray ) / sizeof( *parOptionArray ) );	

BindingFunc* BindingPredictor::createBindingFunc( const BindingPar& par ) const
{
	int nFactors = nInFactors() + nOutFactors();
	Matrix factorIntMat( nFactors, nFactors );
	par.createFactorIntMat( factorIntMat );
	
	//qcheng75 changed
	//return new BindingFunc( motifs, intFunc, par.maxBindingWts, factorIntMat );
	BindingFunc *p= new BindingFunc( motifs, intFunc, par.maxBindingWts, factorIntMat );
	factorIntMat.freeMatrix(); //qcheng75 added
	return p;
}

// double BindingPredictor::compRMSE( const BindingPar* par ) const
// {
// 	// create the binding function
// 	BindingFunc* func = createBindingFunc( par );
// 		
// 	// error of each sequence and RMSE
// 	int N = 0;
// 	double err = 0;	
// 	for ( int i = 0; i < seqs.size(); i++ ) {
// 		for ( int j = 0; j < seqs[ i ].size(); j++ ) {			
// 			// predicted binding of the i-th factor to the j-th sequence in the i-th experiment
// 			double predicted = func->predictBinding( seqSites[ i ][ j ], i );
// 			
// 			// observed binding of the i-th factor to the j-th sequence in the i-th experiment
// 			double observed = bindingData[ i ][ j ];
// 			
// 			err += ( predicted - observed ) * ( predicted - observed );
// 		}
// 		N += seqs[ i ].size();
// 	}	
// 	
// 	return sqrt( err / N );
// }	

double BindingPredictor::compCorr( const BindingPar& par ) const
{	
	// create the binding function
	BindingFunc* func = createBindingFunc( par );
		
	// correlation of each TF experiment
	int N = 0;
	vector< double > corrs; 	
	for ( int i = 0; i < seqs.size(); i++ ) {
		vector< double > predicted; 
		for ( int j = 0; j < seqs[ i ].size(); j++ ) {			
			// predicted binding of the i-th factor to the j-th sequence in the i-th experiment
			predicted.push_back( func->predictBinding( seqSites[ i ][ j ], i ) );			
		}
		corrs.push_back( correlation( predicted, bindingData[ i ] ) );
		predicted.clear();//qcheng75
	}	
	
	func->clear();//qcheng75
	func->releaseFactorIntMat(); //qcheng75
	delete func;//qcheng75

	double CC = mean( corrs );
	corrs.clear();//qcheng75
	return CC;
	//return mean( corrs );
}	

bool BindingPredictor::testPar( const BindingPar& par ) const
{
    vector< bool > option = BindingPredictor::parOption;

    // test the factor interaction matrix
    if (BindingPar::FIX_COOP_PROB==0){ //qcheng75 added the conditional statement
        // test binding weights
        for ( int i = 0; i < nInFactors() + nOutFactors(); i++ ) {
            if ( par.maxBindingWts[i] < BindingPar::min_float_value || par.maxBindingWts[i] > BindingPar::max_float_value ) return false;
        }
    	// test the factor interaction matrix
		if ( option[2] ) {
			for ( int i = 0; i < nInFactors(); i++ ) {
				for ( int j = 0; j <= i; j++ ) {
					if ( par.inFactorIntMat[i][j] < BindingPar::min_float_value || par.inFactorIntMat[i][j] > BindingPar::max_float_value ) return false;
				}
			}

			for ( int i = 0; i < nOutFactors(); i++ ) {
				for ( int j = 0; j < nInFactors(); j++ ) {
					if ( par.outFactorIntMat( i, j ) < BindingPar::min_float_value || par.outFactorIntMat( i, j ) > BindingPar::max_float_value )
						return false;
				}
			}
		}


		// test the exp. ratios
		if ( option[ 3 ] ) {
			for ( int i = 0; i < nInFactors(); i++ ) {
				if ( par.expRatios[ i ] < BindingPar::min_float_value || par.expRatios[i] > BindingPar::max_float_value ) return false;
			}
		}
    }else {
    	for ( int i = 0; i < nInFactors() + nOutFactors(); i++ ) {
    		//cout << "--testPar: " << i << " "<< par.maxBindingWts[i] << " " << pow( 10.0,BindingPar::log10_min_wt) << " " << pow( 10.0,BindingPar::log10_max_wt) << endl;
			if ( par.maxBindingWts[i] < pow( 10.0,BindingPar::log10_min_wt) || par.maxBindingWts[i] > pow( 10.0,BindingPar::log10_max_wt)) {
				return false;
			}
		}

		if ( option[2] ) {  //qcheng75 added
			for ( int i = 0; i < nInFactors(); i++ ) {
				//cout << "--testPar sint: " << i << par.inFactorIntMat[i][i] <<endl;
					if ( par.inFactorIntMat[i][i] < pow( 10.0,BindingPar::log10_min_int) || par.inFactorIntMat[i][i] > pow( 10.0,BindingPar::log10_max_int ) ) return false;
			}
		}

		if ( option[3] ) {  //qcheng75 changed  [2]
			for ( int i = 0; i < nInFactors(); i++ ) {
				for ( int j = 0; j < i; j++ ) {
					//cout << "--testPar mint-: " << i << "," <<j << par.inFactorIntMat[i][j] <<endl;
					if ( par.inFactorIntMat[i][j] < pow( 10.0,BindingPar::log10_min_int) || par.inFactorIntMat[i][j] > pow( 10.0,BindingPar::log10_max_int ) ) return false;
				}
			}
			//cout << "~~~1" << endl;
			for ( int i = 0; i < nOutFactors(); i++ ) {
				for ( int j = 0; j < nInFactors(); j++ ) {
					//cout << "--testPar mint: " << i << "," <<j << " "<< par.outFactorIntMat( i, j ) <<endl;
					if ( par.outFactorIntMat( i, j ) < pow( 10.0,BindingPar::log10_min_int) || par.outFactorIntMat( i, j ) > pow( 10.0,BindingPar::log10_max_int ) )
						return false;
				}
			}
			//cout << "~~~2" << endl;
		}


		// test the exp. ratios
		/*if ( option[ 4 ] ) {  //qcheng75 changed  [3]
			for ( int i = 0; i < nInFactors(); i++ ) {
				if ( par.expRatios[ i ] < BindingPar::min_float_value || par.expRatios[i] > BindingPar::max_float_value ) return false;
			}	
		}
		*/
    }
    return true;
}

void BindingPredictor::randSamplePar( const gsl_rng* rng, BindingPar& par ) const 
{
	vector< bool > option = BindingPredictor::parOption;
	
    if (BindingPar::FIX_COOP_PROB==0){ //qcheng75 added the conditional statement
    	// sample maxBindingWts
    	for ( int i = 0; i < nInFactors() + nOutFactors(); i++ ) {
    		double min_wt = BindingPar::log10_min_wt;
    		double max_wt = BindingPar::log10_max_wt;
    		double rand_wt = gsl_ran_flat( rng, min_wt, max_wt );
    		if ( ( i < nInFactors() && option[ 0 ] ) || ( i >= nInFactors() && option[ 1 ] ) )
    			par.maxBindingWts[ i ] = pow( 10.0, rand_wt );
    	}

    	// sample the interaction matrices
    	if ( option[ 2 ] ) {
    		for ( int i = 0; i < nInFactors(); i++ ) {
    			for ( int j = 0; j <= i; j++ ) {
    				double min_int = BindingPar::log10_min_int;
    				double max_int = BindingPar::log10_max_int;
    				double rand_int = gsl_ran_flat( rng, min_int, max_int );
    				par.inFactorIntMat[ i ][ j ] = pow( 10.0, rand_int );
    			}
    		}

    		for ( int i = 0; i < nOutFactors(); i++ ) {
    			for ( int j = 0; j < nInFactors(); j++ ) {
    				double min_int = BindingPar::log10_min_int;
    				double max_int = BindingPar::log10_max_int;
    				double rand_int = gsl_ran_flat( rng, min_int, max_int );
    				par.outFactorIntMat( i, j ) = pow( 10.0, rand_int );
    			}
    		}
    	}

    	// sample expRatios
    	if ( option[ 3 ] ) {
    		for ( int i = 0; i < nInFactors(); i++ ) {
    			double min_ratio = BindingPar::log10_min_ratio;
    			double max_ratio = BindingPar::log10_max_ratio;
    			double rand_ratio = gsl_ran_flat( rng, min_ratio, max_ratio );
    			par.expRatios[ i ] = pow( 10.0, rand_ratio );
    		}
    	}
    }else{
		// sample maxBindingWts
		if (BindingPar::log10_min_wt< BindingPar::log10_max_wt){ //qcheng75 add the condition control
			for ( int i = 0; i < nInFactors() + nOutFactors(); i++ ) {
				double min_wt = BindingPar::log10_min_wt;
				double max_wt = BindingPar::log10_max_wt;
				double rand_wt = gsl_ran_flat( rng, min_wt, max_wt );
				if ( ( i < nInFactors() && option[ 0 ] ) || ( i >= nInFactors() && option[ 1 ] ) )
					par.maxBindingWts[ i ] = pow( 10.0, rand_wt );
			}

			// sample the interaction matrices
			// sample the interaction matrices
			if (BindingPar::log10_min_int < BindingPar::log10_max_int){
				if ( option[ 2 ] ) {   //qcheng75
					for ( int i = 0; i < nInFactors(); i++ ) {
						double min_int = BindingPar::log10_min_int;
						double max_int = BindingPar::log10_max_int;
						double rand_int = gsl_ran_flat( rng, min_int, max_int );
						par.inFactorIntMat[ i ][ i ] = pow( 10.0, rand_int );
					}
				}

				if ( option[ 3 ] ) { //qcheng75 changed [2]
					for ( int i = 0; i < nInFactors(); i++ ) {
						for ( int j = 0; j < i; j++ ) {//qcheng75<=
							double min_int = BindingPar::log10_min_int;
							double max_int = BindingPar::log10_max_int;
							double rand_int = gsl_ran_flat( rng, min_int, max_int );
							par.inFactorIntMat[ i ][ j ] = pow( 10.0, rand_int );
						}
					}

					for ( int i = 0; i < nOutFactors(); i++ ) {
						for ( int j = 0; j < nInFactors(); j++ ) {
							double min_int = BindingPar::log10_min_int;
							double max_int = BindingPar::log10_max_int;
							double rand_int = gsl_ran_flat( rng, min_int, max_int );
							par.outFactorIntMat( i, j ) = pow( 10.0, rand_int );
						}
					}
				}
			}

			// sample expRatios
			if (BindingPar::log10_min_ratio < BindingPar::log10_max_ratio){
				if ( option[ 4 ] ) {  //qcheng75 changed [3]
					for ( int i = 0; i < nInFactors(); i++ ) {
						double min_ratio = BindingPar::log10_min_ratio;
						double max_ratio = BindingPar::log10_max_ratio;
						double rand_ratio = gsl_ran_flat( rng, min_ratio, max_ratio );
						par.expRatios[ i ] = pow( 10.0, rand_ratio );
					}
				}
			}
		}
    }
	option.clear();//qcheng75
}
		
int BindingPredictor::simplex_minimize( BindingPar& par_result, double& obj_result ) const
{
// 	cout << "Start minimization" << endl;
	int nFactors = nInFactors() + nOutFactors();

	// extract initial parameters
	vector< double > v;
	BindingPar par_init = par_model;
	par_init.getFreePars( v, parOption ); 
		
	// set the objective function
	gsl_multimin_function my_func;
	my_func.f = &gsl_obj_f;
	my_func.n = v.size();
 	my_func.params = (void*)this;
	
	// set the initial values to be searched
 	gsl_vector* x;
	vector< double > lbv,ubv;
 	if ((BindingPar::FIX_COOP_PROB>0) && (BindingPar::NELDER_MEAD_CONSTRAINED)){ //qcheng75
 		par_init.getParsLB( lbv, parOption );
 		par_init.getParsUB( ubv, parOption );
 		x = bvector2gslInf(v, lbv, ubv,BindingPar::IS_LOG_SCALE_STEP);
 	}else{
 		x = vector2gsl(v,BindingPar::IS_LOG_SCALE_STEP);
 	}

	// CHECK POINT: evaluate gsl_obj_f() function
 	// cout << "binding at the initial value of parameters = " << gsl_obj_f( x, (void*)this ) << endl;
		
	// choose the method of optimization and set its parameters
	const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex;
	gsl_vector* ss = gsl_vector_alloc( my_func.n );
	gsl_vector_set_all( ss, 1.0 );
		
	// create the minimizer
	gsl_multimin_fminimizer* s = gsl_multimin_fminimizer_alloc( T, my_func.n );
	gsl_multimin_fminimizer_set( s, &my_func, x, ss );
	
	// iteration
	size_t iter = 0;
	int status;
	double size;	
	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate( s );
	     
		// check for error
		if ( status ) break;

        // check if the current values of parameters are valid
		vector< double > expv;
		if ((BindingPar::FIX_COOP_PROB>0) && (BindingPar::NELDER_MEAD_CONSTRAINED)){
			expv=infgsl2bvector( s->x, lbv, ubv,BindingPar::IS_LOG_SCALE_STEP );
		}else{
			expv=gsl2vector( s->x,BindingPar::IS_LOG_SCALE_STEP );
		}
        BindingPar par_curr( expv, par_init, parOption );
        expv.clear(); //qcheng75
        if ( !testPar( par_curr ) ) {
        	par_curr.clear();//qcheng75
        	break;
        }

		// check for stopping condition
		size = gsl_multimin_fminimizer_size( s );
		status = gsl_multimin_test_size( size, 1e-1 ); 
		par_curr.clear();//qcheng75

		//if ( status == GSL_SUCCESS ) { cout << "converged to minimum at " << iter << endl; }
	     
		//open the comment--<
		//cout << iter << "\t";
		//for ( int i = 0; i < my_func.n; i++ ) {
		//	printf( "%7.3f ", (BindingPar::IS_LOG_SCALE_STEP)?pow((double)10, gsl_vector_get( s->x, i ) ):gsl_vector_get( s->x, i ) );
		//}
		//printf( "f() = %7.3f size = %.3f\n", s->fval, size );
		//open the comment-->
	} while ( status == GSL_CONTINUE && iter < 200 );

	// get the results
	vector< double > expv;
	if ((BindingPar::FIX_COOP_PROB>0) && (BindingPar::NELDER_MEAD_CONSTRAINED)){
		expv=infgsl2bvector( s->x, lbv, ubv,BindingPar::IS_LOG_SCALE_STEP );
	}else{
		expv=gsl2vector( s->x,BindingPar::IS_LOG_SCALE_STEP );
	}
	par_result = BindingPar( expv, par_init, parOption );
	obj_result = -s->fval;
	
	if ((BindingPar::FIX_COOP_PROB>0) && (BindingPar::NELDER_MEAD_CONSTRAINED)){
		lbv.clear();
		ubv.clear();
	}
	// free the minimizer
	gsl_vector_free( x );
	gsl_vector_free( ss );
	gsl_multimin_fminimizer_free( s );	
	
	expv.clear();//qcheng75
	v.clear();//qcheng75
	par_init.clear();//qcheng75

	return 0;
}

int BindingPredictor::gradient_minimize( BindingPar& par_result, double& obj_result ) const
{
 	//cout << "Start gradient_minimize" << endl;
	int nFactors = nInFactors() + nOutFactors();

	// extract initial parameters
	vector< double > v;
	BindingPar par_init = par_model;
	par_init.getFreePars( v, parOption ); 
		
	// set the objective function and its gradient
	gsl_multimin_function_fdf my_func;
	my_func.f = &gsl_obj_f;
	my_func.df = &gsl_obj_df;
	my_func.fdf = &gsl_obj_fdf;
	my_func.n = v.size();
 	my_func.params = (void*)this;
	
	// set the initial values to be searched
	gsl_vector* x = vector2gsl(v,BindingPar::IS_LOG_SCALE_STEP);

	// CHECK POINT: evaluate gsl_obj_f() function
// 	cout << "binding at the initial value of parameters = " << gsl_obj_f( x, (void*)this ) << endl; 
		
	// choose the method of optimization and set its parameters
// 	const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_conjugate_pr;	
	const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_vector_bfgs;
		
	// create the minimizer
	gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc( T, my_func.n );
	double init_step = 0.02, tol = 0.1;
	gsl_multimin_fdfminimizer_set( s, &my_func, x, init_step, tol );
		
	// iteration
	size_t iter = 0;
	int status;
	do {
		iter++;
		double prev_f = s->f;			
		status = gsl_multimin_fdfminimizer_iterate( s );
	    double curr_f = s->f;
// 	    if ( prev_f - curr_f < 0.001 ) break;
		// check for error
		if ( status ) break;

        // check if the current values of parameters are valid
        vector< double > expv = gsl2vector(s->x,BindingPar::IS_LOG_SCALE_STEP);
        BindingPar par_curr( expv, par_init, parOption );
        if ( !testPar( par_curr ) ) {
        	expv.clear();//qcheng75
        	par_curr.clear();//qcheng75
        	break;
        }
        
		// check for stopping condition
		status = gsl_multimin_test_gradient( s->gradient, 1e-3 );
// 		if ( status == GSL_SUCCESS ) { cout << "converged to minimum at " << iter << endl; }

		expv.clear();//qcheng75
		par_curr.clear();//qcheng75

		//open the comment--<
		//cout << iter << "\t";
		//for ( int i = 0; i < my_func.n; i++ ) {
		//	printf( "%7.3f ", (BindingPar::IS_LOG_SCALE_STEP)?pow((double)10, gsl_vector_get( s->x, i ) ):gsl_vector_get( s->x, i ) );
		//}
		//printf( "f() = %7.3f\n", s->f );
		//open the comment-->
	} while ( status == GSL_CONTINUE && iter < 50 );

	// get the results
	vector< double > expv=gsl2vector( s->x,BindingPar::IS_LOG_SCALE_STEP );
	par_result = BindingPar( expv, par_init, parOption );

	obj_result = -s->f;
	
	// free the minimizer
	//qcheng75 changed
	//gsl_multimin_fdfminimizer_free( s );
	//gsl_vector_free( x );
	expv.clear();//qcheng75
	// free the minimizer
	gsl_multimin_fdfminimizer_free( s );
	gsl_vector_free( x ); //qcheng75
	v.clear();//qcheng75
	par_init.clear();//qcheng75
	
	return 0;
}

// box-boundary constrained newton method by adding the KUHN-TUCKER conditions
double gsl_obj_f_KTconstained( const gsl_vector* v, void* params )
{
	// the BindingPredictor object
	BindingPredictor* predictor = (BindingPredictor*)params;

	// parse the variables (parameters to be optimized)
	vector< double > expv;
	for ( int i = 0; i < v->size/3; i++ ) {
		//cout<<gsl_vector_get( v, i )<<endl;
		expv.push_back( pow((double)10, gsl_vector_get( v, i ) ) );
	}
	BindingPar par( expv, predictor->getPar(), BindingPredictor::parOption );

	// call the BindingPredictor object to evaluate the objective function
	double obj = predictor->objFunc( par );
	cout << obj << " ";

	// handling with lower boundary conditions ## qcheng75
	vector<double > parLB,parUB;
	predictor->getParBoundaryConstraints(parLB,parUB);

	double extraobj=0.0, diff=0.0;
	int j=0;
	for ( int i = v->size/3; i < 2* v->size/3; i++ ) {
		j = i- v->size/3;
		diff=1/(pow((double)10, parUB[ j ])-pow((double)10, parLB[ j ]));
		cout << diff << " " << pow((double)10, parLB[ j ]) << " " << gsl_vector_get( v, j ) << " " << pow((double)10, gsl_vector_get( v, j ) ) << " - ";
		extraobj+= diff * ( gsl_vector_get( v, i )) * (gsl_vector_get( v, i )) * (pow((double)10, parLB[ j ]) - pow((double)10, gsl_vector_get( v, j ) ));
	}

	// handling with upper boundary conditions ## qcheng75
	for ( int i = 2*v->size/3; i < v->size; i++ ) {
		j = i- 2* v->size/3;
		diff=1/(pow((double)10, parUB[ j ])-pow((double)10, parLB[ j ]));
		cout << diff << " " << pow((double)10, parUB[j]) << " " << gsl_vector_get( v, j ) << " " << pow((double)10, gsl_vector_get( v, j ) ) << endl;
		extraobj+= diff * ( gsl_vector_get( v, i )) * (gsl_vector_get( v, i )) * (-pow((double)10, parUB[j]) + pow((double)10, gsl_vector_get( v, j ) ));
	}

	//## qcheng75
	expv.clear();
	par.clear();
	parLB.clear();
	parUB.clear();

	return -obj+extraobj;
}

void gsl_obj_df_KTconstained( const gsl_vector* v, void* params, gsl_vector* grad )
{
	double step = 1.0E-3;
	numeric_deriv( grad, gsl_obj_f_KTconstained, v, params, step );
}

void gsl_obj_fdf_KTconstained( const gsl_vector* v, void* params, double* result, gsl_vector* grad )
{
	*result = gsl_obj_f_KTconstained( v, params );
	gsl_obj_df_KTconstained( v, params, grad );
}

int BindingPredictor::gradient_minimize_KTconstrained( BindingPar& par_result, double& obj_result ) const
{
// 	cout << "Start minimization" << endl;
	int nFactors = nInFactors() + nOutFactors();

	// extract initial parameters
	vector< double > v;
	BindingPar par_init = par_model;
	par_init.getFreePars( v, parOption );

	// set the objective function and its gradient

	gsl_multimin_function_fdf my_func;
	my_func.f = &gsl_obj_f_KTconstained;
	my_func.df = &gsl_obj_df_KTconstained;
	my_func.fdf = &gsl_obj_fdf_KTconstained;
	my_func.n = v.size() * 3;
 	my_func.params = (void*)this;

	// set the initial values to be searched
	gsl_vector* x = gsl_vector_alloc( v.size()*3 );
	for ( int i = 0; i < v.size(); i++ ) gsl_vector_set( x, i, log10( v[ i ] ) );
	if ((BindingPar::FIX_COOP_PROB>0) && (BindingPar::NEWTON_KTCONSTRAINED)){ //set up initial values for other parameters ##qcheng75
		for ( int i = v.size(); i < 3 * v.size(); i++ ) {
			gsl_vector_set( x, i, 1 );
		}
	}
	// CHECK POINT: evaluate gsl_obj_f() function
	// cout << "binding at the initial value of parameters = " << gsl_obj_f( x, (void*)this ) << endl;

	// choose the method of optimization and set its parameters
	// const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_conjugate_pr;
	const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_vector_bfgs;

	// create the minimizer
	gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc( T, my_func.n );
	double init_step = 0.02, tol = 0.1;
	gsl_multimin_fdfminimizer_set( s, &my_func, x, init_step, tol );

	// iteration
	size_t iter = 0;
	int status;
	do {
		iter++;
		double prev_f = s->f;
		status = gsl_multimin_fdfminimizer_iterate( s );
	    double curr_f = s->f;
// 	    if ( prev_f - curr_f < 0.001 ) break;

		// check for error
		if ( status ) break;

        // check if the current values of parameters are valid
        vector< double > expv = gsl2vector(s->x,BindingPar::IS_LOG_SCALE_STEP);
        BindingPar par_curr( expv, par_init, parOption );
        if ( !testPar( par_curr ) ) {
        	expv.clear();//qcheng75
        	par_curr.clear();//qcheng75
        	break;
        }

		// check for stopping condition
		status = gsl_multimin_test_gradient( s->gradient, 1e-3 );
// 		if ( status == GSL_SUCCESS ) { cout << "converged to minimum at " << iter << endl; }

		expv.clear();//qcheng75
		par_curr.clear();//qcheng75

		//open the comment--<
		cout << iter << "\t";
		for ( int i = 0; i < my_func.n; i++ ) {
			printf( "%7.3f ", (BindingPar::IS_LOG_SCALE_STEP)?pow((double)10, gsl_vector_get( s->x, i ) ):gsl_vector_get( s->x, i ) );
		}
		printf( "f() = %7.3f\n", s->f );
		//open the comment-->
	} while ( status == GSL_CONTINUE && iter < 50 );

	// get the results
	vector< double > expv = gsl2vector(s->x,BindingPar::IS_LOG_SCALE_STEP);
	par_result = BindingPar( expv, par_init, parOption );

	obj_result = -s->f;

	// free the minimizer
	//qcheng75 changed
	//gsl_multimin_fdfminimizer_free( s );
	//gsl_vector_free( x );
	expv.clear();//qcheng75
	// free the minimizer
	gsl_multimin_fdfminimizer_free( s );
	gsl_vector_free( x ); //qcheng75
	v.clear();//qcheng75
	par_init.clear();//qcheng75

	return 0;
}

int BindingPredictor::SA_n_tries = 0;
int BindingPredictor::SA_iters_fixed_T = 40;
double BindingPredictor::SA_step_size = 1.0;
double BindingPredictor::SA_K = 1.0;
double BindingPredictor::SA_T_init = 2.0;
double BindingPredictor::SA_mu_T = 1.01;
double BindingPredictor::SA_T_min = 0.1;		

int BindingPredictor::SA_minimize( BindingPar& par_result, double& obj_result ) const
{
	// parameters for SA
	gsl_siman_params_t params = { SA_n_tries, SA_iters_fixed_T, SA_step_size, SA_K, SA_T_init, SA_mu_T, SA_T_min };

	// random number generator
	gsl_rng* rng;
	gsl_rng_env_setup();
	const gsl_rng_type * T = gsl_rng_default;	// create rng type
	rng = gsl_rng_alloc( T );
	gsl_rng_set( rng, time( 0 ) );		// set the seed equal to simulTime(0)
	
	// initialization
	SAState* state = new SAState( BindingPar( nInFactors(), nOutFactors() ), this );
	
	// SA
	gsl_siman_solve( rng, state, SA_Efunc, SA_step, 0, SA_print, SA_copy, SA_copy_construct, SA_destroy, 0, params );
	
	// results	
	par_result = state->par;
	obj_result = objFunc( par_result );
		
	// free the memory
	gsl_rng_free( rng );
	
	return 0;
}

double gsl_obj_f( const gsl_vector* v, void* params )
{
	//cout<< "enter gsl_obj_f" << endl;
	// the BindingPredictor object
	BindingPredictor* predictor = (BindingPredictor*)params;
		
	// parse the variables (parameters to be optimized)	
	vector< double > expv = gsl2vector(v,BindingPar::IS_LOG_SCALE_STEP);
	BindingPar par( expv, predictor->getPar(), BindingPredictor::parOption );
		
	// call the BindingPredictor object to evaluate the objective function 
	double obj = predictor->objFunc( par );	

	double extraobj=0;
	if ((BindingPar::FIX_COOP_PROB>0) && (BindingPar::PENALTY_BASED_CONSTRAINED)){
		vector<double > parLB,parUB;
		predictor->getParBoundaryConstraints(parLB,parUB);

		double maxJ=0.0,valJ=0.0;
		for ( int i = 0; i < v->size; i++ ) {
			valJ = (gsl_vector_get( v, i )-parUB[ i ])*(gsl_vector_get( v, i )-parLB[ i ]);
			maxJ = valJ <=0? 0 : valJ;
			if ((maxJ<1) && (maxJ>0)) BindingPar::PENALTY_COEFFIENT=10000000000.0;
			extraobj+=BindingPar::PENALTY_COEFFIENT * maxJ;
			//cout <<"extraobj=" << extraobj << " maxJ=" << maxJ << endl;
		}
	}

	expv.clear();//qcheng75
	par.clear();//qcheng75

	return -obj + extraobj;
}

void gsl_obj_df( const gsl_vector* v, void* params, gsl_vector* grad )
{
	//cout<< "enter gsl_obj_df" << endl;
	double step = 1.0E-3;
	numeric_deriv( grad, gsl_obj_f, v, params, step );	
}

void gsl_obj_fdf( const gsl_vector* v, void* params, double* result, gsl_vector* grad )
{
	//cout<< "enter gsl_obj_fdf" << endl;
	*result = gsl_obj_f( v, params ); 
	gsl_obj_df( v, params, grad );		
}

double SA_Efunc( void* xp )
{
	SAState* state = (SAState*)xp;
	double r = state->predictor->objFunc( state->par );	
	//double e = -log( ( 1 + r ) / ( 1 - r ) );
	double e = -r;
	
	return e;
}

double SA_metric( void* xp, void* yp )
{
// 	SAState* state_x = (SAState*)xp;
// 	SAState* state_y = (SAState*)yp;
// 	int option = BindingPredictor::parOption;
// 	
// 	vector< double > pars_x;
// 	vector< double > pars_y;
// 	state_x->getFreePars( pars_x, option );
// 	state_y->getFreePars( pars_y, option );
// 	return Eucledian_dist( pars_x, pars_y );
}

void SA_step( const gsl_rng* r, void* xp, double step_size )
{
	SAState* state = (SAState*)xp;
	vector< bool > option = BindingPredictor::parOption;
	BindingPar newPar = state->par;
	
	// sample maxBindingWts of the in-factors
	if ( option[ 0 ] ) {
		for ( int i = 0; i < state->predictor->nInFactors(); i++ ) {
			double init_wt = state->predictor->getPar().maxBindingWts[ i ];
			double min_wt = init_wt == BindingPar::initMaxBindingWt ? BindingPar::log10_min_wt : log10( state->predictor->getPar().maxBindingWts[ i ] ) - BindingPar::log10_range_wt;  
			double max_wt = init_wt == BindingPar::initMaxBindingWt ? BindingPar::log10_max_wt : log10( state->predictor->getPar().maxBindingWts[ i ] ) + BindingPar::log10_range_wt;  
			double u = gsl_rng_uniform( r );
			double step = step_size * BindingPar::wt_step; 
			double newWeight = log10( state->par.maxBindingWts[ i ] ) + 2 * u * step - step;
			if ( newWeight > max_wt ) newWeight = max_wt - ( newWeight - max_wt );
			if ( newWeight < min_wt ) newWeight = min_wt + ( min_wt - newWeight );
			newPar.maxBindingWts[ i ] = pow( 10.0, newWeight );
		}	
	}
	
	// sample maxBindingWts of the out-factors
	if ( option[ 1 ] ) {
		for ( int i = 0; i < state->predictor->nOutFactors(); i++ ) {
			double min_wt = BindingPar::log10_min_wt;  
			double max_wt = BindingPar::log10_max_wt;  
			double u = gsl_rng_uniform( r );
			double step = step_size * BindingPar::wt_step; 
			double newWeight = log10( state->par.maxBindingWts[ i + state->predictor->nInFactors() ] ) + 2 * u * step - step;
			if ( newWeight > max_wt ) newWeight = max_wt - ( newWeight - max_wt );
			if ( newWeight < min_wt ) newWeight = min_wt + ( min_wt - newWeight );
			newPar.maxBindingWts[ i + state->predictor->nInFactors() ] = pow( 10.0, newWeight );				
		}	
	}
	
	// sample the interaction matrices
	if ( option[ 2 ] ) {
		for ( int i = 0; i < state->predictor->nInFactors(); i++ ) {
			for ( int j = 0; j <= i; j++ ) {
				double min_int = BindingPar::log10_min_int;
				double max_int = BindingPar::log10_max_int;				
				double u = gsl_rng_uniform( r );
				double step = step_size * BindingPar::int_step;
				double newInt = log10( state->par.inFactorIntMat[ i ][ j ] ) + 2 * u * step - step;
				if ( newInt > max_int ) newInt = max_int - ( newInt - max_int );
				if ( newInt < min_int ) newInt = min_int + ( min_int - newInt );
				newPar.inFactorIntMat[ i ][ j ] = pow( 10.0, newInt );
			}	
		}	
	
		for ( int i = 0; i < state->predictor->nOutFactors(); i++ ) {
			for ( int j = 0; j < state->predictor->nInFactors(); j++ ) {
				double min_int = BindingPar::log10_min_int;
				double max_int = BindingPar::log10_max_int;				
				double u = gsl_rng_uniform( r );
				double step = step_size * BindingPar::int_step;
				double newInt = log10( state->par.outFactorIntMat( i, j ) ) + 2 * u * step - step;
				if ( newInt > max_int ) newInt = max_int - ( newInt - max_int );
				if ( newInt < min_int ) newInt = min_int + ( min_int - newInt );
				newPar.outFactorIntMat( i, j ) = pow( 10.0, newInt );
			}	
		}				
		
	}
	
	// sample expRatios
	if ( option[ 3 ] ) {
		for ( int i = 0; i < state->predictor->nInFactors(); i++ ) {
			double u = gsl_rng_uniform( r );
			double step = step_size * BindingPar::ratio_step;
			double newRatio	= log10( state->par.expRatios[ i ] ) + 2 * u * step - step;
			newPar.expRatios[ i ] = pow( 10.0, newRatio );
		}	
	}

	state->par = newPar;
}

void SA_print( void *xp )
{
	SAState* state = (SAState*)xp;
	vector< double > pars;
	state->par.getFreePars( pars, BindingPredictor::parOption );	
	cout << "\t" << pars << "\t";
}

void SA_copy( void* src, void* dest )
{
	SAState* state_src = (SAState*)src;
	SAState* state_dest = (SAState*)dest;
	state_dest->par = state_src->par;
	state_dest->predictor = state_src->predictor;
}

void* SA_copy_construct( void* xp )
{
	SAState* state = (SAState*)xp;
	SAState* new_state = new SAState( state->par, state->predictor );	
	return new_state;
}

void SA_destroy( void* xp )
{
	if ( !xp ) return;
	SAState* state = (SAState*)xp;
	delete state;
}

StepwisePredictor::StepwisePredictor( const vector< Motif >& _motifs, const vector< double >& _energyThrs, const vector< vector< Sequence > >& _seqs, const vector< vector< double > >& _bindingData, const FactorIntFunc* _intFunc ) : BasePredictor( _motifs, _energyThrs, _seqs, _bindingData ), intFunc( _intFunc ), model( 0 )
{
}

int StepwisePredictor::train()
{
	// initialize the factor set
	vector< bool > currSet( motifs.size(), 0 );
	for ( int i = 0; i < nInFactors(); i++ ) currSet[ i ] = true;
	BindingPredictor* currPredictor = evalSet( currSet );
	double currPerf = currPredictor->getObj(); 

	// stepwise procedure
	while ( true ) {		
		// forward test
		int bestFactor, bestPredictorIdx; 
		double bestPerf = currPerf; 		
		vector< BindingPredictor* > bps;
		for ( int i = nInFactors(); i < motifs.size(); i++ ) {
			if ( currSet[ i ] ) continue;
			bps.push_back( forwardTest( currSet, i ) );  
			if ( bps.back()->getObj() > bestPerf ) {
				bestFactor = i;
				bestPredictorIdx = bps.size() - 1;
				bestPerf = bps.back()->getObj();
			}
		}
		
		for ( int i = 0; i < bps.size(); i++ ) {
			if ( i == bestPredictorIdx ) continue;
			bps[ i ]->clear();//qcheng75
			delete bps[ i ];	
		}
		
		if ( ( bestPerf - currPerf ) / currPerf > testThr ) {
			currSet[ bestFactor ] = 1;
			if ( currPredictor ) {
				currPredictor->clear();//qcheng75
				delete currPredictor;
			}
			currPredictor = bps[ bestPredictorIdx ];
			currPerf = bestPerf;	
		} else break; 
		
		// backward test
		for ( int i = nInFactors(); i < motifs.size(); i++ ) {
			if ( !currSet[ i ] || i == bestFactor ) continue;
			BindingPredictor* bp = backwardTest( currSet, i );
			if ( ( currPerf - bp->getObj() ) / currPerf < testThr ) {
				currSet[ i ] = 0;
                if ( currPredictor ) {
                	currPredictor->clear();//qcheng75
                	delete currPredictor;
                }
                currPredictor = bp;
                currPerf = currPredictor->getObj();
			} else {
				bp->clear();//qcheng75
                delete bp;
            }
		}
	}

	// the result
	indicators = currSet;	
	model = currPredictor; 
	mapPar( model, currSet, par_model );
	obj_model = model->getObj();
	
	return 0;		
}

int StepwisePredictor::analyzeFeatures( vector< vector< double > >& featureEffects ) const
{
	featureEffects.clear();

//     cout << "Model learned:" << endl;
//     model->getPar().print();
//     cout << "nInFactors = " << model->nInFactors() << endl << "nOutFactors = " << model->nOutFactors() << endl;
	vector< vector< double > > subEffects;
	model->analyzeFeatures( subEffects );
	int counter = 0;
	for ( int i = 0; i < nOutFactors(); i++ ) {
        vector< double > effects = indicators[ nInFactors() + i ] ? subEffects[ counter++ ] : vector< double >( nInFactors(), 0.0 );
		featureEffects.push_back( effects );
	}
	
	return 0;
}

double StepwisePredictor::predict( const Sequence& targetSeq, int factorIdx ) const
{
	return model->predict( targetSeq, factorIdx );
} 

double StepwisePredictor::testThr = 0.02;

BindingPredictor* StepwisePredictor::evalSet( const vector< bool >& indSet )
{
	assert( indSet.size() == motifs.size() );
	cout << "Evaluating motif set: " << indSet << endl;
	vector< Motif > subMotifs;
	vector< double > subEnergyThrs;
	
	for ( int i = 0; i < indSet.size(); i++ ) {
		if ( indSet[ i ] ) { 
			subMotifs.push_back( motifs[ i ] );
			subEnergyThrs.push_back( energyThrs[ i ] );
		}		
	}
	
	BindingPredictor* bp = new BindingPredictor( subMotifs, subEnergyThrs, seqs, bindingData, intFunc );
	bp->train();
	return bp; 
}
	
BindingPredictor* StepwisePredictor::forwardTest( const vector< bool >& currSet, int testFactor )
{
	assert( currSet[ testFactor ] == 0 );
	vector< bool > set = currSet;
	set[ testFactor ] = 1;
	return evalSet( set );
}

BindingPredictor* StepwisePredictor::backwardTest( const vector< bool >& currSet, int testFactor )
{
	assert( currSet[ testFactor ] == 1 );
	vector< bool > set = currSet;
	set[ testFactor ] = 0;
	return evalSet( set );	
}

void StepwisePredictor::mapPar( const BindingPredictor* bp, const vector< bool >& indSet, BindingPar& par )
{
	BindingPar par_sub = bp->getPar();
	par = BindingPar( nInFactors(), nOutFactors() );
	int counter = 0;
	
	// set the maxBindingWts
	counter = 0;
	for ( int i = 0; i < nInFactors() + nOutFactors(); i++ ) {
		if ( indSet[ i ] ) par.maxBindingWts[ i ] = par_sub.maxBindingWts[ counter++ ];
	}
	
	// set the interaction matrix
	par.inFactorIntMat = par_sub.inFactorIntMat;
	counter = 0;
	for ( int i = 0; i < nOutFactors(); i++ ) {
		if ( indSet[ nInFactors() + i ] ) par.outFactorIntMat.setRow( i, par_sub.outFactorIntMat.getRow( counter++ ) );
	}
	
	// set the expRatios
	par.expRatios = par_sub.expRatios;	
	par_sub.clear();//qcheng75
}

FeatureSelector::FeatureSelector( const vector< Motif >& _motifs, const vector< double >& _energyThrs, const vector< vector< Sequence > >& _seqs, const vector< vector< double > >& _bindingData, const FactorIntFunc* _intFunc ) : motifs( _motifs ), energyThrs( _energyThrs ), seqs( _seqs ), bindingData( _bindingData ), intFunc( _intFunc )
{
	assert( motifs.size() > 0 );
	assert( energyThrs.size() == motifs.size() );
	assert( bindingData.size() == seqs.size() );	
}

int FeatureSelector::select()
{
    // initialization
    indicators = vector< bool >( nInFactors() + nOutFactors(), false );
    for ( int i = 0; i < nInFactors(); i++ ) indicators[i] = true;
    vector< Motif > inFactors;
    vector< double > inEnergyThrs;
    for ( int i = 0; i < nInFactors(); i++ ) {
        inFactors.push_back( motifs[ i ] );
        inEnergyThrs.push_back( energyThrs[ i ] );
    }
    
    // testing each out-factor
    for ( int i = 0; i < nOutFactors(); i++ ) {
        // run the analysis with in-factors and this out-factor
        vector< Motif > testFactors = inFactors;
        vector< double > testEnergyThrs = inEnergyThrs;
        testFactors.push_back( motifs[ nInFactors() + i ] );
        testEnergyThrs.push_back( energyThrs[ nInFactors() + i ] );
        BindingPredictor* bp = new BindingPredictor( testFactors, testEnergyThrs, seqs, bindingData, intFunc );
        bp->train();

        // the effects of the out-factor in all experiments
        vector< vector< double > > effectsAll;
        bp->analyzeFeatures( effectsAll );
        vector< double > effects;   // effects of this out-factor in each experiment
        for ( int j = 0; j < nInFactors(); j++ ) effects.push_back( effectsAll[ 0 ][ j ] );

        // choose an out-factor if its effect in any experiment is greater than the threshold
        for ( int j = 0; j < nInFactors(); j++ ) {
            if ( abs( effects[ j ] ) >= effectThr ) {
                indicators[ nInFactors() + i ] = true; break;
            }
        }
    }

    return 0;
}

double FeatureSelector::effectThr = 0.5;

BindingCrossValidator::BindingCrossValidator( const vector< vector< Sequence > >& _seqs, const vector< vector< string > >& _names, const vector< vector< double > >& _bindingData ) : seqs( _seqs ), names( _names ), bindingData( _bindingData )
{
	int nInFactors = seqs.size();
	assert( names.size() == nInFactors );
	assert( bindingData.size() == nInFactors );
}

BindingCrossValidator::BindingCrossValidator( const vector< vector< Sequence > >& _seqs, const vector< vector< double > >& _bindingData ) : seqs( _seqs ), bindingData( _bindingData )
{
	int nInFactors = seqs.size();
	assert( bindingData.size() == nInFactors );
	
	names.resize( nInFactors );
	for ( int i = 0; i < nInFactors; i++ ) {
		for ( int j = 0; j < seqs[ i ].size(); j++ ) {
			char buffer[ 20 ];
			sprintf( buffer, "Exp%iSeq%i", i, j );
			names[ i ].push_back( buffer );	
		}	
	}
}

void BindingCrossValidator::partition( int nFolds )
{
	assert( nFolds > 0 );
	int nExps = seqs.size();
// 	srand( static_cast<unsigned>( time( 0 ) ) ); 
	
	// partition the indices of data points
	vector< vector< vector< int > > > indices( nFolds );
	for ( int k = 0; k < nFolds; k++ ) {
		indices[ k ] = vector< vector< int > >( nExps ); 	
	}
	for ( int i = 0; i < nExps; i++ ) {
		int nDataPoints = seqs[ i ].size();
		
		// reshuffle the indices for this experiment	
		vector< int > shuffledIndices;
		for ( int j = 0; j < nDataPoints; j++ ) shuffledIndices.push_back( j );
		random_shuffle( shuffledIndices.begin(), shuffledIndices.end() );		
		
		// add the index sets
		int sampleSize = nDataPoints / nFolds;
		for ( int k = 0; k < nFolds; k++ ) {
			vector< int > indexSet;
			for ( int j = k * sampleSize; j < min( nDataPoints, ( k + 1 ) * sampleSize ); j++ ) {
				indexSet.push_back( shuffledIndices[ j ] );	
			}	
			indices[ k ][ i ] = indexSet;
		}
	}

	// partition the data
	trainSeqs.resize( nFolds );
	trainNames.resize( nFolds );
	trainBindingData.resize( nFolds );
	testSeqs.resize( nFolds );
	testNames.resize( nFolds );
	testBindingData.resize( nFolds );
	for ( int k = 0; k < nFolds; k++ ) {
		trainSeqs[ k ] = vector< vector< Sequence > >( nExps );
		trainNames[ k ] = vector< vector< string > >( nExps );		
		trainBindingData[ k ] = vector< vector< double > >( nExps );
		testSeqs[ k ] = vector< vector< Sequence > >( nExps );
		testNames[ k ] = vector< vector< string > >( nExps );		
		testBindingData[ k ] = vector< vector< double > >( nExps );			
	}
	for ( int i = 0; i < nExps; i++ ) {
		// create training and testing data for each experiment
		for ( int k = 0; k < nFolds; k++ ) {
			int nDataPoints = seqs[ i ].size();
			
			// testing data
			vectSubset( seqs[ i ], indices[ k ][ i ], testSeqs[ k ][ i ] );
			vectSubset( names[ i ], indices[ k ][ i ], testNames[ k ][ i ] );			
			vectSubset( bindingData[ i ], indices[ k ][ i ], testBindingData[ k ][ i ] );
			
			// training data
			vector< int > complIndexSet;
			indexCompl( indices[ k ][ i ], nDataPoints, complIndexSet );
			vectSubset( seqs[ i ], complIndexSet, trainSeqs[ k ][ i ] );
			vectSubset( names[ i ], complIndexSet, trainNames[ k ][ i ] );			
			vectSubset( bindingData[ i ], complIndexSet, trainBindingData[ k ][ i ] );			
		}
	}					
}

void BindingCrossValidator::runThermoModel( vector< double >& results, const vector< Motif >& motifs, const vector< double >& energyThrs, const FactorIntFunc* intFunc, bool featureSelection, int modelOption )
{
	int nFactors = motifs.size();
	int nInFactors = seqs.size();
	assert( energyThrs.size() == nFactors );
	
	results.clear();
	int nFolds = trainSeqs.size();
	
	for ( int k = 0; k < nFolds; k++ ) {
		// training
		cout << "Fold " << k + 1 << endl;
        vector< Motif > activeMotifs;
        vector< double > activeEnergyThrs;        
        if ( featureSelection ) {
            FeatureSelector selector( motifs, energyThrs, trainSeqs[ k ], trainBindingData[k], intFunc );
            selector.select();
            vector< bool > indicators = selector.getIndicators(); 
            for ( int i = 0; i < motifs.size(); i++ ) {
                if ( indicators[ i ] ) {
                    activeMotifs.push_back( motifs[ i ] );
                    activeEnergyThrs.push_back( energyThrs[ i ] );
                }
            }
        } else {
            activeMotifs = motifs;
            activeEnergyThrs = energyThrs;
        }
        
		BasePredictor* predictor;
		if ( modelOption == 0 ) predictor = new BindingPredictor( activeMotifs, activeEnergyThrs, trainSeqs[ k ], trainBindingData[ k ], intFunc );
		else predictor = new StepwisePredictor( activeMotifs, activeEnergyThrs, trainSeqs[ k ], trainBindingData[ k ], intFunc );
		predictor->train();
		predictor->getPar().print();
			
		// testing	
		double result = predictor->test( testSeqs[ k ], testBindingData[ k ], perfOption );
		cout << "Performance = " << result << endl << endl;

		delete predictor; //qcheng75
		results.push_back( result );				 	
	}	
}

void BindingCrossValidator::print() const
{
	int nFolds = trainSeqs.size();
	int nExps = trainSeqs[ 0 ].size();
	for ( int k = 0; k < nFolds; k++ ) {
		for ( int i = 0; i < nExps; i++ ) {
			cout << "Fold " << k << "\tExp " << i << endl;
			// training data
			cout << "Training data" << endl;
			for ( int j = 0; j < trainSeqs[ k ][ i ].size(); j++ ) {
				cout << trainNames[ k ][ i ][ j ] << "\t" << trainBindingData[ k ][ i ][ j ] << endl;
			}
			
			// testing data	
			cout << "Testing data" << endl;
			for ( int j = 0; j < testSeqs[ k ][ i ].size(); j++ ) {
				cout << testNames[ k ][ i ][ j ] << "\t" << testBindingData[ k ][ i ][ j ] << endl;
			}			
		}
	}	
}

int BindingCrossValidator::perfOption = 0;

