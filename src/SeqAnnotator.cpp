#include <sstream>
#include "SeqAnnotator.h"

bool isNt( int a )
{
	if ( a < 0 || a > 3 ) return false;
	else return true;	
}

int complement( int a )
{
	assert( a >= 0 && a < ALPHABET_SIZE );
		
	if ( a == 0 ) return 3;
	if ( a == 1 ) return 2;
	if ( a == 2 ) return 1;
	if ( a == 3 ) return 0;	
	if ( a == MISSING ) return MISSING;
	if ( a == GAP ) return GAP;	
}

int symbolToInt( char c )
{
	char upper = toupper( c );
	for ( int i = 0; i < ALPHABET_SIZE; i++ ) {
		if ( ALPHABET[ i ] == upper ) return i;	
	}
	
	return -1;
}

char strand2char( bool strand )
{
	if ( strand ) return '+';
	else return '-';	
}

Sequence::Sequence( const Sequence& other, int start, int length, bool strand )
{
	assert( start >= 0 && length >= 0 && ( start + length ) <= other.size() );	

	for ( int i = 0; i < length; i++ ) {
		if ( strand ) {	nts.push_back( other[ start + i ] ); }
		else { nts.push_back( complement( other[ start + length - 1 - i ] ) ); }
	}	
}

void Sequence::clear() {
	if (! nts.empty()) nts.clear();
}

Sequence::~Sequence(){
	clear();
}

int Sequence::push_back( int nt )
{	
	assert( nt >= 0 && nt < ALPHABET_SIZE );
	nts.push_back( nt );
	
	return 0;
}

int Sequence::push_back( const Sequence& elem )
{
	for ( int i = 0; i < elem.size(); i++ ) push_back( elem[ i ] );	
	return 0;
}

Sequence Sequence::compRevCompl() const
{
	return Sequence( *this, 0, size(), false );	
}

void Sequence::getNtCounts( vector< int >& counts ) const
{
	counts.clear();
	for ( int i = 0; i < NBASES; i++ ) {
		counts.push_back( 0 );
	}
	
	for ( int i = 0; i < nts.size(); i++ ) {
		if ( nts[ i ] != GAP ) counts[ nts[ i ] ]++;	
	}
}

bool Sequence::containsMissing() const
{
	for ( int i = 0; i < nts.size(); i++ ) {
		if ( nts[ i ] == MISSING ) return true;	
	}	
	
	return false;
}

int Sequence::load( const string& file, string& name, int format )
{
	vector< Sequence > seqs;
	vector< string > names;
	int rval = readSequences( file, seqs, names, format );
	if ( rval == RET_ERROR ) return RET_ERROR;
	
	copy( seqs[ 0 ] );
	name = names[ 0 ];
	return rval;
}

int Sequence::load( const string& file, int format )
{
	string name;
	int rval = load( file, name, format );
	
	return rval;	
}

ostream& operator<<( ostream& os, const Sequence& seq )
{
	// output the nts
	for ( int i = 0; i < seq.size(); i++ ) {
		os << ALPHABET[ seq[ i ] ];	
	}	
			
	return os;
}

int readSequences( const string& file, vector< Sequence >& seqs, vector< string >& names, int format )
{
	// check if the format character is legal
	if ( format != FASTA ) { return RET_ERROR; }
	 
	// 	open the file
	ifstream fin( file.c_str() );
	if ( !fin ) { cerr << "Cannot open" << file << endl; exit( 1 ); }

	string line;
	Sequence seq;
	
	// read sequences: FASTA format
	if ( format == FASTA ) {
		while ( getline( fin, line ) ) {
			// add the sequence and start a new sequence if the line starts with >
			//cout << line << endl;
			if ( line[ 0 ] == '>' ) { 	
				if ( seq.size() ) {
					seqs.push_back( seq );
					seq.clear();	
				}
				
				stringstream ss( line.substr( 1 ) );
				string name; 
				ss >> name;
				names.push_back( name );
			} else { 
				// check if the line contains content
				int start = line.find_first_not_of( " \t\r" );
				int last = line.find_last_not_of( " \t\r" );
				if ( start == string::npos || last == string::npos ) continue;
				
				// append the sequence	
				for ( int i = start; i <= last; i++ ) {
					int nt = symbolToInt( line[ i ] );	// could be a NNN or gap
					if ( nt >= 0 && nt < ALPHABET_SIZE ) {
						seq.push_back( nt );
					} else {
						//cerr << "Illegal symbol: " << nt << " in " << file << endl;
						return RET_ERROR;	
					} 
				}
			}			
		}
		
		// add the last sequence
		if( seq.size() ) seqs.push_back( seq );
		seq.clear(); //qcheng75

		return 0;
	}	
}

int readSequences( const string& file, vector< Sequence >& seqs, int format )
{
	vector< string > names;
	int rval = readSequences( file, seqs, names, format );	
	return rval;
}

int writeSequences( const string& file, const vector< Sequence >& seqs, const vector< string >& names, int format )
{
	assert( seqs.size() == names.size() );
	
	// check if the format character is legal
	if ( format != FASTA ) { return RET_ERROR; }
		
	ofstream fout( file.c_str() );
	
	if ( format == FASTA ) {
		for ( int i = 0; i < seqs.size(); i++ ) {
			fout << ">" << names[ i ] << endl;
			fout << seqs[ i ] << endl;
		}
	}
	
	return 0;
}

int writeSequences( const string& file, const vector< Sequence >& seqs, int format )
{
	// default name: integer starting from 1
	vector< string > names;
	for ( int i = 0; i < seqs.size(); i++ ) {
		char buffer[ 10 ];
		sprintf( buffer, "%i", i );
		names.push_back( string( buffer ) );	
	}	
	
	// print
	return writeSequences( file, seqs, names, format );
}

void compWtmx( const Matrix& countMatrix, double pseudoCount, Matrix& pwm )//qcheng75 added pwm
{
    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
    
    int l = countMatrix.nRows();		// l: the length of motif		
    pwm.setDimensions(l, 4); //qcheng75
    //Matrix pwm( l, 4 );

//     // the sum of each position/column should be a const. (number of sequences)
//     double n = 0;		// number of sites used in the count matrix
//     for ( int j = 0; j < 4; j++ ) {
//         n += countMatrix( 0, j );
//     }
//     for ( int i = 1; i < l; i++ ) {
//         double count = 0;
//         for ( int j = 0; j < 4; j++ ) {
//             count += countMatrix( i, j );
//         }
//         if ( count != n ) { cout << "count matrix incorrect" << endl; exit( 1 ); }
//     }
    
    // the multinomial distribution at each column
    for ( int i = 0; i < l; i++ ) {
        double n = 0;       // total counts at this position 
        for ( int j = 0; j < 4; j++ ) {
            n += countMatrix( i, j );
        }
        for ( int j = 0; j < 4; j++ ) {
            pwm( i, j ) = ( countMatrix( i, j ) + pseudoCount ) / ( n + 4.0 * pseudoCount );
        }	
    }

    //return pwm;		//qcheng75
}

Motif::Motif( const Matrix& _pwm, const vector< double >& _background ) : pwm( _pwm ), background( _background ), LLRMat( pwm.nRows(), 4 )
{
	assert( background.size() == 4 );	
	
	init();
}

Motif::Motif( const Matrix& countMatrix, double pseudoCount, const vector< double >& _background ) : background( _background ), LLRMat( countMatrix.nRows(), 4 )
{
	assert( background.size() == 4 );
	
	compWtmx( countMatrix, pseudoCount, pwm);//qcheng75
	//pwm=compWtmx( countMatrix, pseudoCount);

	init();
}

double Motif::LLR( const Sequence& elem ) const
{
	int l = pwm.nRows();
	if ( elem.size() != l ) return GSL_NEGINF;
	if ( elem.containsMissing() ) return GSL_NEGINF;
	
	double result = 0;
	for ( int i = 0; i < l; i++ ) {
		result += LLRMat( i, elem[ i ] ); 	
	}
	
	return result;
}

double Motif::energy( const Sequence& elem ) const
{
	return ( -LLR( elem ) + maxLLR );	
}

void Motif::sample( const gsl_rng* rng, Sequence& elem, bool strand ) const
{
	assert( rng != NULL );
	
	int l = pwm.nRows();
	Sequence sampleElem;
	for ( int i = 0; i < l; i++ ) {
		// nt. distribution at position i
		vector< double > distr = pwm.getRow( i );
		
		// sample nt. from this distribution	
		int nt = sampleMul( rng, distr );
		sampleElem.push_back( nt );

		distr.clear();//qcheng75
	}		
	
	if ( strand == 0 ) elem = sampleElem.compRevCompl();
	else elem = sampleElem;
}

int Motif::load( const string& file, const vector< double >& background, string& name )
{
	vector< Motif > motifs;
	vector< string > names;
	int rval = readMotifs( file, background, motifs, names );
	if ( rval == RET_ERROR ) return RET_ERROR;
	
	copy( motifs[ 0 ] );
	name = names[ 0 ];
	return rval;				
}

int Motif::load( const string& file, const vector< double >& background )
{
	string name;
	int rval = load( file, background, name );
	
	return rval;	
}

ostream& operator<<( ostream& os, const Motif& motif )
{
	os << motif.pwm;
	
	return os;
}

int Motif::GLOBAL_isFeatureStat=0; // ## qcheng75
int Motif::SITE_THRESHOLD_BY_LLRthres=0; // ## qcheng75
void Motif::init()
{
	int l = pwm.nRows();
	
	// compute the LLR matrix
	for ( int i = 0; i < l; i++ ) {
		for ( int j = 0; j < 4; j++ ) {			
			LLRMat( i, j ) = log( pwm( i, j ) / background[ j ] );
		}
	}
	
	// the strongest site
	for ( int i = 0; i < l; i++ ) {
		int b_max;
		max( pwm.getRow( i ), b_max );
		maxSite.push_back( b_max );	
	}
	
	// compute the LLR of the strongest site
	maxLLR = 0;
	for ( int i = 0; i < l; i++ ) {
		maxLLR += LLRMat( i, maxSite[ i ] );	
	}

    // ## qcheng75 motif specificity measured by information content sum_{k} sum_{b \in \{A,C,G,T\}} F(b,k)log (F(b,k)/background(b))
	if (GLOBAL_isFeatureStat==1){
		infoContent = 0;
		//cout << ">> checking motif IC: " << endl;
		for ( int i = 0; i < l; i++ ) {
			for ( int j = 0; j < 4; j++ ) {
				infoContent += pwm( i, j )*log( pwm( i, j ) / background[ j ] );
			}
		}
		//cout << infoContent << "\n" << endl;
	}
}

//qcheng75
void Motif::clear(){

	if ( ! background.empty()) background.clear();
	if (! maxSite.getNTS().empty()) maxSite.clear();
	if (pwm.getData()) pwm.freeMatrix();
	if (LLRMat.getData()) LLRMat.freeMatrix();

}

int readMotifs( const string& file, const vector< double >& background, vector< Motif >& motifs, vector< string >& names )
{
	// 	open the file
	ifstream fin( file.c_str() );
	if ( !fin ) { cerr << "Cannot open" << file << endl; exit( 1 ); }

	string line;
	
	// read the motifs
	do {
		getline( fin, line );
		
		if ( line[ 0 ] != '>' ) continue;
		
		// read the names, length and pseudocount
		int MAX_SIZE = 100;
		char lineStr[ MAX_SIZE ];
		strcpy( lineStr, ( line.substr( 1 ) ).c_str() );
		char *name, *lengthStr, *pseudoCountStr;
		name = strtok( lineStr, " \t" );
		lengthStr = strtok( NULL, " \t" );
		pseudoCountStr = strtok( NULL, " \t" );
		int length;
		double pseudoCount;
		if ( lengthStr ) length = atoi( lengthStr );
		else { return RET_ERROR; }
		if ( pseudoCountStr ) pseudoCount = atof( pseudoCountStr );
		else pseudoCount = PSEUDO_COUNT;
		
		// read the count matrix
		Matrix countMat( length, NBASES );
		for ( int i = 0; i < length; ++i ) {
			for ( int j = 0; j < NBASES; ++j ) {
				fin >> countMat( i, j );
			}	
		}
		
		// create the motif
		names.push_back( string( name ) );
//         cout << "Creating motif " << name << endl;
		motifs.push_back( Motif( countMat, pseudoCount, background ) );	

		countMat.freeMatrix(); //qcheng75
	} while ( !fin.eof() );
					
	return 0;
}

int readMotifs( const string& file, const vector< double >& background, vector< Motif >& motifs )
{
	vector< string > names;
	return readMotifs( file, background, motifs, names );	
}

ostream& operator<<( ostream& os, const Site& site )
{
	char strandChar = site.strand ? '+' : '-';
	os << site.start << "\t" << strandChar << "\t" << site.factorIdx << "\t" << site.energy << "\t" << site.wtRatio;
	
	return os;
}

bool siteOverlap( const Site& a, const Site& b, const vector< Motif >& motifs )
{
	if ( a.start + motifs[ a.factorIdx ].length() <= b.start ) return false;
	if ( b.start + motifs[ b.factorIdx ].length() <= a.start ) return false;
	
	return true;	
}

//qcheng75
void SeqAnnotator::clear(){
	motifs.clear();	// all the TF binding motifs
	energyThrs.clear();	// energy thresholds for all the motifs
}

// ## qcheng75
bool SeqAnnotator::isSite(const Sequence& elem, const Motif& motif, double threshold)  const{
	if (Motif::SITE_THRESHOLD_BY_LLRthres==1){
		if ( motif.LLR(elem) >= threshold) {
			return true;
		}else{
			return false;
		}
	}else{
		if ( motif.energy( elem ) <= threshold ) {
			return true;
		}else{
			return false;
		}
	}
}

int SeqAnnotator::annot( const Sequence& seq, SiteVec& sites ) const
{
	// cout << "start annotation:" << endl;
	//qcheng75 changed
	//sites.clear();
	if (! sites.empty()){
		for (int i=0; i< sites.size(); i++){
			delete &(sites[i]);
		}
	}
	
	// scan the sequence for the sites of all motifs
	for ( int i = 0; i < seq.size(); i++ ) {
		// test for each motif
		for ( int k = 0; k < motifs.size(); k++ ) {
			int l = motifs[ k ].length();
			if ( i + l > seq.size() ) continue;
			double energy;
			
			// positive strand
			Sequence elem( seq, i, l, 1 );
			energy = motifs[ k ].energy( elem );
			// cout << elem << Site( i, 1, k, exp( energy ) ) << endl;
			if (isSite(elem, motifs[ k ], energyThrs[ k ])){
			//if ( energy <= energyThrs[ k ] ) {
				sites.push_back( Site( i, 1, k, energy ) );
			}	
			
			// negative strand
			Sequence rcElem( seq, i, l, 0 );
			energy = motifs[ k ].energy( rcElem );
			//	cout << rcElem << Site( i, 0, k, exp( energy ) ) << endl;
			if (isSite(rcElem, motifs[ k ], energyThrs[ k ])){
			//if ( energy <= energyThrs[ k ] ) {
				sites.push_back( Site( i, 0, k, energy ) );
			}				
		}	
	}
	
	// cout << "end annotation" << endl;
	return sites.size();
}

void SeqAnnotator::doFeatureStatistics(std::ofstream& fsout, string seqname, const Sequence& seq, int motif, double chipscore) { // ## qcheng75 added it
	int numP=0,numN=0, numSites=0;
	double LLRSum_P=0; // absolute LLR sum at positive strand
	double mismatchEnergy_P=0; //relative LLR Sum at positive strand
	double thermoK_P=1;
	double LLRSum_N=0; // absolute LLR sum  at negative strand
	double mismatchEnergy_N=0; //relative LLR Sum at negative strand
	double thermoK_N=1;
	double LLRSum=0; // absolute LLR sum  in total
	double mismatchEnergy=0; //relative LLR Sum in total
	double thermoK=1;
	int k=motif;
	int l = motifs[ k ].length();
	double TRAP_lamda = 0.7, TRAP_R_0 = exp(0.585*l - 5.66);  // TRAP affinity = sum_{l in both strands} TRAP_R_0 * e^{-energy(site) * 1/TRAP_lamda }/(1+ TRAP_R_0 * e^{-energy(site) * 1/TRAP_lamda })
	double TRAP_Affinity_P =0, TRAP_Affinity_N =0, TRAP_Affinity =0;
	// scan the sequence for the sites of all motifs
	for ( int i = 0; i < seq.size(); i++ ) {
			if ( i + l > seq.size() ) continue;
			double energy;

			// positive strand
			Sequence elem( seq, i, l, 1 );
			energy = motifs[ k ].energy( elem );
			if ( energy <= energyThrs[ k ] ) {
				//fsout << "P1 " << i << " " << k << " "<< motifs[ k ].getMaxLLR() << " " << motifs[ k ].LLR(elem) << " " << energy << " " << exp( -energy ) << endl;
				numP++;
				LLRSum_P += motifs[ k ].LLR(elem);
				mismatchEnergy_P += energy;
				thermoK_P = thermoK_P * exp(-energy);
				TRAP_Affinity_P += (TRAP_R_0 * exp( -energy * (1/TRAP_lamda)))/(1+ (TRAP_R_0 * exp( -energy * (1/TRAP_lamda))));
			}

			// negative strand
			Sequence rcElem( seq, i, l, 0 );
			energy = motifs[ k ].energy( rcElem );
			if ( energy <= energyThrs[ k ] ) {
				//fsout << "P-1 " << i << " " << k << " "<< motifs[ k ].getMaxLLR() << " " << motifs[ k ].LLR(rcElem) << " " << energy << " " << exp( -energy )<< endl;
				numN++;
				LLRSum_N += motifs[ k ].LLR(rcElem);
				mismatchEnergy_N += energy;
				thermoK_N = thermoK_N * exp(-energy);
				TRAP_Affinity_N += (TRAP_R_0 * exp( -energy * (1/TRAP_lamda)))/(1+ (TRAP_R_0 * exp( -energy * (1/TRAP_lamda))));

			}
			numSites = numP + numN;
			LLRSum = LLRSum_P + LLRSum_N;
			mismatchEnergy = mismatchEnergy_P + mismatchEnergy_N;
			thermoK = thermoK_P * thermoK_N;
			TRAP_Affinity = TRAP_Affinity_P + TRAP_Affinity_N;
	}
	fsout<< k << "\t" << seqname << "\t" << chipscore << "\t" << numSites << "\t" << LLRSum << "\t" << mismatchEnergy << "\t" << thermoK << "\t" << TRAP_Affinity << endl;
	fsout.flush();
	return;
}
