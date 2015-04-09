/*
* the library of mathematics
*/
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

#include "Tools.h"

const log_add_table table( -10.0, 0, 500 );		// global variable

// create a vector from gsl_vector
vector< double > gsl2vector( const gsl_vector* v, int IS_GSL_MIN_OPT_STEP_LOG_SCALE )
{
	//cout << "enter gsl2vector....." << v->size << endl;
	vector< double > result;
	double val=0.0;
	for ( int i = 0; i < v->size; i++ ) {
		val = IS_GSL_MIN_OPT_STEP_LOG_SCALE==1? pow((double)10,gsl_vector_get( v, i )): gsl_vector_get( v, i );
		result.push_back( val );
		//cout << "para" << i << " " << val << endl;
	}
	return result;
}

gsl_vector* vector2gsl( const vector< double >& v, int IS_GSL_MIN_OPT_STEP_LOG_SCALE )
{
	gsl_vector* result = gsl_vector_alloc( v.size() );
	for ( int i = 0; i < v.size(); i++ ) gsl_vector_set( result, i, IS_GSL_MIN_OPT_STEP_LOG_SCALE==1? log10(v[ i ]): v[ i ] );
	return result;
}


// create a vector from gsl_vector with box-boundary constraints  //qcheng75 add
vector< double > infgsl2bvector( const gsl_vector* v, const vector<double >& lbV, const vector<double >& ubV, int IS_GSL_MIN_OPT_STEP_LOG_SCALE )
{
	vector< double > result;
	for ( int i = 0; i < v->size; i++ ) {
		result.push_back(  IS_GSL_MIN_OPT_STEP_LOG_SCALE==1? pow((double)10,inverse_infty_transform(gsl_vector_get( v, i ), lbV[i], ubV[i]) ): inverse_infty_transform(gsl_vector_get( v, i ), lbV[i], ubV[i]));
	}
	return result;
}

//qcheng75 add
gsl_vector* bvector2gslInf( const vector< double >& v, const vector<double >& lbV, const vector<double >& ubV, int IS_GSL_MIN_OPT_STEP_LOG_SCALE )
{
	gsl_vector* result = gsl_vector_alloc( v.size() );
	for ( int i = 0; i < v.size(); i++ ) {
		gsl_vector_set( result, i, infty_transform( IS_GSL_MIN_OPT_STEP_LOG_SCALE==1? log10(v[ i ]):v[i], lbV[i], ubV[i] ));
	}
	return result;
}

double infty_transform( double x, double a, double b )
{
    // assert x \in [a,b]
	cout << x << " " << a << " " << b ;
    assert( x > a && x < b );

    // transformation
    double y = ( x - a ) / ( b - a );
    double z = log( y / ( 1.0 - y ) );

    cout << "~" << z << endl;
    return z;
}

double inverse_infty_transform( double z, double a, double b )
{
	cout << z << " " << a << " " << b;
    assert( a <= b );

    double y = exp( z ) / ( 1.0 + exp( z ) );
    double x = a + y * ( b - a );

    cout << "-" << x << endl;
    return x;
}

// Matrix constructor
Matrix::Matrix( int nRows, int nCols ) 
{
	assert( nRows > 0 && nCols > 0 );
	data = gsl_matrix_alloc( nRows, nCols );
}

// Matrix constructor
Matrix::Matrix( const gsl_matrix* _data )
{
	assert( _data );
	setDimensions( _data->size1, _data->size2 ); 
	gsl_matrix_memcpy( data, _data );
}

// Matrix: constructor
Matrix::Matrix( const double **_data, int nRows, int nCols )
{
	assert( _data && nRows > 0 && nCols > 0 );
	setDimensions( nRows, nCols );
	for ( int i = 0; i < nRows; ++i ) {
		for ( int j = 0; j < nCols; ++j ) {
			gsl_matrix_set( data, i, j, _data[ i ][ j ] ); 	
		}	
	}
}

// destructor
Matrix::~Matrix()
{
	//if ( data ) gsl_matrix_free( data );	
}

//qcheng75
void Matrix::freeMatrix()
{
	//cout << "Release Matrix ...." << endl;
	if ( data != NULL ) gsl_matrix_free( data );
	data = NULL;
}

// get the row vector
vector< double > Matrix::getRow( int row ) const
{
	gsl_vector* v = gsl_vector_alloc( nCols() );
	gsl_matrix_get_row( v, data, row );

	//qcheng75
	vector<double> ret=gsl2vector( v, 0);
	gsl_vector_free(v);
	return ret;
	//return gsl2vector( v );
}

// get the column vector
vector< double > Matrix::getCol( int col ) const
{
	gsl_vector* v = gsl_vector_alloc( nRows() );
	gsl_matrix_get_col( v, data, col );
	return gsl2vector( v,0 );
}

// set the matrix dimensions
void Matrix::setDimensions( int nRows, int nCols)
{
	assert( nRows > 0 && nCols > 0 );

	/*// qcheng75 deleted
	if ( data && ( data->size1 != nRows || data->size2 != nCols ) ) gsl_matrix_free( data );
	//if ( data ) gsl_matrix_free( data );
	data = gsl_matrix_alloc( nRows, nCols );
	*/

	//qcheng75 added
	if ( data && ( data->size1 != nRows || data->size2 != nCols ) ) {
		gsl_matrix_free( data );
		data = gsl_matrix_alloc( nRows, nCols );
	}else if (! data){
		data = gsl_matrix_alloc( nRows, nCols );
	}
}

// equality test
bool Matrix::operator==( const Matrix& other ) const
{
	// compare dimensions
	if ( nRows() != other.nRows() || nCols() != other.nCols() )
		return false;
	
	// compare elements
	for ( int i = 0; i < nRows(); i++ ) {
		for ( int j = 0; j < nCols(); j++ ) {
			if ( getElement( i, j ) != other.getElement( i, j ) ) return false;	
		}	
	}
	
	return true;
}

// check if square
bool Matrix::isSquare() const
{
	if ( nRows() == nCols() ) 
		return true;
	else 
		return false;	
}

// Matrix: set row
void Matrix::setRow( int row, const vector< double >& v )
{
	assert( row >= 0 && row < data->size1 );
	
	for ( int j = 0; j < data->size2; j++ ) 
		setElement( row, j, v[ j ] );	
}

// Matrix: set column
void Matrix::setCol( int col, const vector< double >& v )
{
	assert( col >= 0 && col < data->size2 );
	
	for ( int i = 0; i < data->size1; i++ ) 
		setElement( i, col, v[ i ] );	
		
}

// Matrix: set all rows
void Matrix::setRows( const vector< double >& v )
{
	for ( int i = 0; i < data->size1; i++ )
		setRow( i, v );	
}

// Matrix: set all columns
void Matrix::setCols( const vector< double >& v )
{
	for ( int j = 0; j < data->size2; j++ ) 
		setCol( j, v );	
}

// set the identity matrix of dimension n
void Matrix::setIdentityMatrix( int n )		
{
	assert( n > 0 );
	setDimensions( n, n );	
	setZero();
	for ( int i = 0; i < n; i++ ) setElement( i, i, 1.0 );
}

// set the diagonal matrix
void Matrix::setDiagonalMatrix( const vector< double >& diag )
{
	int n = diag.size();
	
	setDimensions( n, n );
	setZero();
	for ( int i = 0; i < n; i++ ) setElement( i, i, diag[ i ] );
}

// Matrix transpose
Matrix Matrix::transpose() const
{
	Matrix result( nCols(), nRows() );
	
	for ( int i = 0; i < nRows(); i++ ) {
		for ( int j = 0; j < nCols(); j++ ) {
			result( j, i ) = getElement( i, j );
		}
	}
	return( result );
}

// matrix addition
Matrix Matrix::operator+( const Matrix& other ) const
{
	// check dimensions
	assert( nRows() == other.nRows() && nCols() == other.nCols() );
	
	// addition
	Matrix result( *this );
	gsl_matrix_add( result.data, other.data );	

	return result;
}

// matrix addition
Matrix& Matrix::operator+=( const Matrix& other )
{
	// check dimensions
	assert( nRows() == other.nRows() && nCols() == other.nCols() );

	// addition
	gsl_matrix_add( data, other.data );
	
	return ( *this );
}

// matrix multiplication by a constant
Matrix Matrix::operator*( const double c ) const
{
	Matrix result( *this );
	gsl_matrix_scale( result.data, c );	

	return result;
}	

// matrix multiplication by a constant
Matrix& Matrix::operator*=( const double c )
{
	gsl_matrix_scale( data, c );	
	
	return ( *this );	
}

// matrix multiplication
Matrix Matrix::operator*( const Matrix& other ) const
{
	// check dimensions
	assert( nCols() == other.nRows() );
	
	// matrix multiplication
	Matrix result( nRows(), other.nCols() );
	gsl_linalg_matmult( data, other.data, result.data );
	
	return result;
}

// output the Matrix
ostream& operator<<( ostream& os, const Matrix& m )
{
	os.setf( ios::fixed );
	
	for ( int i = 0; i < m.nRows(); i++ ) {
		for ( int j = 0; j < m.nCols(); j++ ) {
			os << m( i, j );
			if ( j != m.nCols() - 1 ) os << " ";
		}
		os << endl;
	}

	return os;
}

// load Matrix from a file
// If readDims = true, read dimensions (first row of the file); o/w, only read data, assuming the correct dimensions are known
int Matrix::load( const string& fileName, bool readDims )
{
	// open the file
	ifstream fin( fileName.c_str() );
	if ( !fin ){ cerr << "Cannot open " << fileName << endl;	exit( 1 ); } 
	
	// read the matrix dimensions
	int nRows, nColumns;
	if ( readDims ) {
		fin >> nRows >> nColumns;
		assert( nRows > 0 && nColumns > 0 );
		setDimensions( nRows, nColumns );		
	} else {
		nRows = data->size1;
		nColumns = data->size2;
	}

	// read the matrix 
	for ( int i = 0; i < nRows; i++ ) {
		for ( int j = 0; j < nColumns; j++ ) {
			double elem;
			fin >> elem;
			gsl_matrix_set( data, i, j, elem );
		}
	}
	
	return 0;
}

// save Matrix to a file
void Matrix::save( const string& fileName ) const
{
	// open the file
	ofstream fout( fileName.c_str(), ofstream::out | ofstream::app );
	
	// write the matrix
	fout << *this;	
}

// log_add_table constructor: compute the table used for log_add algorithm
log_add_table::log_add_table( double a, double b, int N ) : x_array( N + 1 ), T( N + 1 )
{
	delta = ( b - a ) / N;
	x_array[ 0 ] = a;
	T[ 0 ] = log( 1 + exp( a ) );
	for ( int i = 1; i < N; i++ ) {
		double x = x_array[ i - 1 ] + delta;
		x_array[ i ] = x;
		T[ i ] = log( 1 + exp( x ) );	
	}
	x_array[ N ] = b; 
	T[ N ] = log( 1 + exp( b ) );
}

// log_add algorithm
double log_add( double p, double q )
{
	// check if p or q is -inf
	if ( gsl_isinf( p ) == -1 ) return q;
	if ( gsl_isinf( q ) == -1 ) return p;
	
	// log_add
	int N = table.x_array.size();
	double x, val;
	if ( p > q ) {
		x = q - p; val = p;
	}
	else {
		x = p - q; val = q;
	}
	if ( x < table.x_array[ 0 ] ) {
		val += 0;
	} else {
		int k = (int)floor( ( x - table.x_array[ 0 ] ) / table.delta );
		val += table.T[ k ];
		//if ( k == N ) 
		//	val += table.T[ N ];
		//else 
		//	val += table.T[ k ] + ( x - table.x_array[ k ] ) * ( table.T[ k + 1 ] - table.T[ k ] ) / table.delta;	
	}	
	
	return val;
}

// log_add algorithm
double log_add( const vector< double >& p )
{
	double result = p[ 0 ];
	
	for ( int i = 1; i < p.size(); i++ ) {
		if ( gsl_isinf( p[ i ] ) == -1 ) continue;
		result = log_add( result, p[ i ] );	
	}	
	
	return result;
}

// log_add algorithm
double log_add( double p, double q, double r )
{
	double result;
	result = log_add( p, q );
	result = log_add( result, r );
	
	return result;	
}

vector< double > log( const vector< double >& v )
{
	vector< double > result;
	for ( int i = 0; i < v.size(); i++ ) result.push_back( log( v[ i ] ) );
	return result;	
}

vector< double > exp( const vector< double >& v )
{
	vector< double > result;
	for ( int i = 0; i < v.size(); i++ ) result.push_back( exp( v[ i ] ) );
	return result;			
}

// log. of a matrix
Matrix log( const Matrix& M )
{
	int m = M.nRows(), n = M.nCols();
	
	Matrix result( m, n );	
	for ( int i = 0; i < m; i++ ) {
		for ( int j = 0; j < n; j++ ) {
				result( i, j ) = log( M( i, j ) );
		}
	}
		
	return result;
}	

// maximum of a set of numbers
double max( const vector< double >& v, int &arg )
{
	double v_max = GSL_NEGINF;
	for ( int i = 0; i < v.size(); ++i ) {
		if ( v[ i ] > v_max ) {
			arg = i;
			v_max = v[ i ];	
		}	
	}
	
	return v_max;
}

double mean( const vector< double >& x )
{
	if ( x.size() == 0 ) return 0;
	
	double sum = 0; 
	for ( int i = 0; i < x.size(); i++ ) sum += x[ i ];
	return ( sum / x.size() );
}

double median( const vector< double > &x )
{
	if ( x.size() == 0 ) { return 0; }
	
	vector< double > xCopy( x );
	sort( xCopy.begin(), xCopy.end() );
	int middle = (int)( ( xCopy.size() - 1 ) / 2.0 );
	return xCopy[ middle ];	
}

double std_dev( const vector< double >& x )
{
	double* data = new double[ x.size() ];
	for ( int i = 0; i < x.size(); i++ ) data[ i ] = x[ i ];
	return gsl_stats_sd( data, 1, x.size() );
}

double correlation( const vector< double >& v1, const vector< double >& v2 )
{
	assert( v1.size() == v2.size() );
	int n = v1.size();
	
	double* data1 = new double[ n ]; 
	double* data2 = new double[ n ];
	for ( int i = 0; i < n; i++ ) {
		data1[ i ] = v1[ i ];
		data2[ i ] = v2[ i ];	
	}
	
	//qcheng75
	//return gsl_stats_covariance( data1, 1, data2, 1, n ) / ( gsl_stats_sd( data1, 1, n ) * gsl_stats_sd( data2, 1, n) );
	double sd2=gsl_stats_sd( data1, 1, n ) * gsl_stats_sd( data2, 1, n) ;
	if (sd2<=1e-5) sd2=1e-5;
	double ret = gsl_stats_covariance( data1, 1, data2, 1, n ) / (sd2 );
	delete data1;
	delete data2;
	return ret;
}

// check if a vector of real numbers is probablity mass function: return false if not
bool isPmf( const vector< double > &p )
{
	double sum = 0; 
	for ( int i = 0; i < p.size(); ++i ) {
		sum += p[ i ];	
	}
	if ( abs( sum - 1.0 ) > numeric_limits< double >::epsilon() ) return false;
	else return true;
}

// sample one outcome from a multinomial distribution
int sampleMul( const gsl_rng * r, const vector< double > &p )
{
	int k = p.size();	// # of different outcomes
	
	// sample
	double u = gsl_rng_uniform( r );
	double c = 0;
	for ( int i = 0; i < k; ++i ) {
		c += p[ i ];
		if ( u < c ) { return i; }	
	}			
	
	return k - 1;
}

// sample a positive integer from a truncated geometric series: 1, r, r^2, ... r^(n - 1)
int sampleTruncatedGeometric( const gsl_rng* rng, double r, int n )
{
	double factor = ( 1.0 - r ) / ( 1.0 - pow( r, n ) );
	vector< double > probs;
	double curr = 1.0;
	for ( int i = 0; i < n; i++ ) {
		probs.push_back( curr * factor );
		curr *= r;
	}
	
	return 1 + sampleMul( rng, probs );
}

// mixture of two multinomial distributions, where the weight of the first one is "weight"
vector< double > multinomialMixture( const vector< double >& distr1, const vector< double >& distr2, double weight )
{
	int size = distr1.size();
	assert( distr2.size() == size && isPmf( distr1 ) && isPmf( distr2 ) );
	
	vector< double > mix;
	for ( int i = 0; i < size; i++ ) {
		mix.push_back( weight * distr1[ i ] + ( 1 - weight ) * distr2[ i ] );	
	}	
	
	return mix;
}

// the threshold value at a given percentage in a vector of numbers
// if highest is true, then the value at the highest p; otherwise, the value at the lowest p
double elementAt( const vector< double >& v, double p, bool highest )
{
	assert( p >= 0 && p <= 1 );
	
	// create a copy of data because we don't want to modify the original data
	vector< double > vCopy( v );

	// get the threshold
	sort( vCopy.begin(), vCopy.end() );
	int index;
	if ( highest == true ) {
		index = (int)ceil( vCopy.size() * ( 1 - p ) ) - 1;
	} else {
		index = (int)floor( vCopy.size() * p ) - 1;
	}
	if ( index < 0 ) index = 0;
	
	return vCopy[ index ];		
}

double trunc( double x )
{
	double lower = floor( x );
	double upper = ceil( x );
	if ( lower == upper ) return lower;
	
	if ( x - lower < upper - x ) return lower;
	else if ( x - lower > upper - x ) return upper;
	else return x >= 0 ? lower : upper;
}

double logit( double x )
{
	assert( x >= 0 && x <= 1 );
	
	return log( x / ( 1.0 - x ) );	
}

double inv_logit( double x )
{
	return ( exp( x ) / ( 1.0 + exp( x ) ) );
}

double Eucledian_dist( const vector< double >& x, const vector< double >& y )
{
	assert( x.size() == y.size() );
	
	double sum = 0;
	for ( int i = 0; i < x.size(); i++ ) {
		sum += ( x[ i ] - y[ i ] ) * ( x[ i ] - y[ i ] );	
	}	
	
	return sqrt( sum );
}

vector< double > transform( const vector< double > w )
{
	double w0 = 1.0 - sum( w );
	vector< double > results;
	for ( int k = 0; k < w.size(); k++ ) results.push_back( log( w0 ) - log( w[ k ] ) );	
	
	return results;
}

vector< double > inverse_transform( const vector< double > u )
{
	vector< double > neg_expo;
	double neg_expo_tot = 0;
	for ( int k = 0; k < u.size(); k++ ) {
		double temp = exp( -u[ k ] );
		neg_expo.push_back( temp );	
		neg_expo_tot += temp;
	}
	
	vector< double > results;
	for ( int k = 0; k < u.size(); k++ ) 
		results.push_back( neg_expo[ k ] / ( 1.0 + neg_expo_tot ) );
		
	return results;
}

void numeric_deriv( gsl_vector* grad, double (*f)( const gsl_vector*, void* ), const gsl_vector* v, void* params, double step )
{
	int n = v->size;
	double f_val = (*f)( v, params );

	gsl_vector* dv = gsl_vector_alloc( n );
	for ( int i = 0; i < n; i++ ) {
		gsl_vector_memcpy( dv, v );
		gsl_vector_set( dv, i, gsl_vector_get( v, i ) + step );
		double partial_deriv = ( (*f)( dv, params ) - f_val ) / step;
		gsl_vector_set( grad, i , partial_deriv );
	}
}

/*
void numeric_deriv( gsl_vector* grad, double (*f)( const gsl_vector*, void* ), const gsl_vector* v, void* params, double step )
{
	int n = v->size;
	double f_val = (*f)( v, params );
	double val=0.00,nexval=0.00,deltx=0.001;
	
	gsl_vector* dv = gsl_vector_alloc( n );	
	for ( int i = 0; i < n; i++ ) {
		gsl_vector_memcpy( dv, v );
		val = gsl_vector_get( v, i );
		gsl_vector_set( dv, i, val + step );
		nexval =  (*f)( dv, params ) - f_val ;
		deltx=step;//(exp(val+step)-exp(val));
		double partial_deriv = nexval / deltx;
		gsl_vector_set( grad, i , partial_deriv );
		cout << val << " " << step << " " << deltx << " " << nexval << " "<< partial_deriv << endl;
	}		
	gsl_vector_free(dv); //qcheng75
}
*/
// read parameter file into a <field, value> table. Return 0 if successful, -1 otherwise
// Format: field = value [# ...] where # indicates the start of comments (ignored)
// Ex. kappa = 2.0 # kappa: the transition/transversion bias
int readParams( const string& file, map< string, string >& params, const char* comment )
{
	params.clear();
	
	// open the parameter file
	ifstream fin( file.c_str() );
	if ( !fin ){ cerr << "Cannot open " << file << endl; exit( 1 );}
	
	// read parameters
	string line;
	while ( getline( fin, line ) ) {
		// skip a line if it starts with a comment or non-alphnumerical charcter 
		int first = line.find_first_not_of( " \t\r" );
		if ( !isalnum( line[ first ] ) || strchr( comment, line[ first ] ) ) continue;
				
		// check if the line contains '=', if not error
		int posEq = line.find( '=' );
		if ( posEq == string::npos ) return RET_ERROR;
		
		// read the parameter name
		string name;
		int pos1 = line.find_first_not_of( " \t\r" );
		int pos2 = line.find_last_not_of( " \t\r", posEq - 1);
		if ( pos1 <= pos2 ) name = line.substr( pos1, pos2 - pos1 + 1 );
		
		// read the value, skipping the part of sentence starting with a comment symbol if any
		string value;
		pos1 = line.find_first_not_of( " \t\r", posEq + 1 );
		int posComment = line.find_first_of( comment );
		if ( posComment != string::npos ) {
			pos2 = line.find_last_not_of( " \t\r", posComment - 1 );
		} else {
			pos2 = line.find_last_not_of( " \t\r" );
		}
		if ( pos1 <= pos2 ) value = line.substr( pos1, pos2 - pos1 + 1 );
		
		// add to the table of parameters
		if ( name.size() && value.size() ) 
			params[name] = value;
		else 
			return RET_ERROR;	
	}
	
	return 0;
}

void split( vector< string >& result, const string& str, const string& sep )
{
	result.clear();
	
	int MAX_STRING = 1000;
	char buffer[ MAX_STRING ];
	strcpy( buffer, str.c_str() );
	char* token;
	token = strtok( buffer, sep.c_str() );
	while ( token ) {
		result.push_back( token );
		token = strtok( NULL, sep.c_str() );
	}
}

void stringToVector( vector< string >& result, const string& str, const string& leftBoundary, const string& rightBoundary, const string& sep )
{
	// find the substring containing the values
	int start = str.find_first_of( leftBoundary );
	int end = str.find_last_of( rightBoundary );
	string sub = str.substr( start + 1, end - start - 1 );

	// split the substring
	split( result, sub, sep );		
}

void stringToVector( vector< double >& result, const string& str, const string& leftBoundary, const string& rightBoundary, const string& sep )
{
	result.clear();
	vector< string > resultStr;
	stringToVector( resultStr, str, leftBoundary, rightBoundary, sep );
	
	for ( int i = 0; i < resultStr.size(); i++ ) {
		result.push_back( atof( resultStr[ i ].c_str() ) );		
	}
}

string IO_Ctrl::vector_delimiter = "\t";
string IO_Ctrl::map_delimiter = "\n";
string IO_Ctrl::field_separator = "\t";

int indexCompl( const vector< int >& indices, int n, vector< int >& complIndices )
{
	complIndices.clear();
	
	// bit vector of 1 to n
	vector< int > bitVect( n );
	for ( int i = 0; i < n; i++ ) bitVect[ i ] = 0;
	
	// mark all indices that are in the input set
	for ( int i = 0; i < indices.size(); i++ ) {
		if ( indices[ i ] >= n ) return RET_ERROR;
		bitVect[ indices[ i ] ] = 1;
	}
	
	// results
	for ( int i = 0; i < n; i++ ) {
		if ( bitVect[ i ] == 0 ) complIndices.push_back( i ); 		
	}
	
	return 0;
}
