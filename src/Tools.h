/*
* the library of general tools
*/
#ifndef TOOLS_H
#define TOOLS_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string.h>
#include <string>
#include <vector>
#include <limits>
#include <numeric>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>

using namespace std;

/* constants */
const int RET_ERROR = -1;		// return error of a function
const int INT_INF = numeric_limits< int >::max() / 2;

/*****************************************************
* Vectors and Matrices (GSL)
******************************************************/

/* functions for gsl_vector */
vector< double > gsl2vector( const gsl_vector* v, int IS_GSL_MIN_OPT_STEP_LOG_SCALE );
gsl_vector* vector2gsl( const vector< double >& v, int IS_GSL_MIN_OPT_STEP_LOG_SCALE );
vector< double > infgsl2bvector( const gsl_vector* v, const vector<double >& lbV, const vector<double >& ubV, int IS_GSL_MIN_OPT_STEP_LOG_SCALE );
gsl_vector* bvector2gslInf( const vector< double >& v, const vector<double >& lbV, const vector<double >& ubV, int IS_GSL_MIN_OPT_STEP_LOG_SCALE );

/* Matrix class */
class Matrix {
	gsl_matrix* data;	// data
public:
	// constructors
	Matrix() : data( NULL ) {}
	Matrix( int nRows, int nCols );
	Matrix( const gsl_matrix* _data );
	Matrix( const double **_data, int nRows, int nCols );
	void copy( const Matrix& other ) 
	{ 
		if ( other.isEmpty() ) { data = NULL; return; }
		setDimensions( other.nRows(), other.nCols() ); 
		gsl_matrix_memcpy( data, other.data );
	}
	Matrix( const Matrix& other ) : data( NULL ) { copy( other ); }

	// destructor
	~Matrix();
	void freeMatrix(); //qcheng75

	// assignment
	Matrix& operator=( const Matrix& other ) { copy(other); return *this; }
			
	// access methods
	int nRows() const { if ( !isEmpty() ) return data->size1; else return 0; }
	int nCols() const { if ( !isEmpty() ) return data->size2; else return 0; }
	double getElement( int row, int col ) const { return gsl_matrix_get( data, row, col ); } 
	void setElement( int row, int col, double v ) { gsl_matrix_set( data, row, col, v ); }
	const double& operator()( int row, int col) const {
		assert( row >= 0 && row < data->size1  && col >= 0 && col < data->size2 );
		return *gsl_matrix_ptr( data, row, col );		
	}
	double& operator()( int row, int col ) {
		assert( row >= 0 && row < data->size1  && col >= 0 && col < data->size2 );
		return *gsl_matrix_ptr( data, row, col );
	}
	gsl_matrix* getData() const { return data; }
	void setDimensions( int nRows, int nCols);
	vector< double > getRow( int row ) const;
	vector< double > getCol( int col ) const;
	
	// information
	bool operator==( const Matrix& other ) const;
	bool isSquare() const;
	bool isEmpty() const { if ( data ) return false; else return true; }
	
	// set matrix elements
	void setZero() { gsl_matrix_set_zero( data ); }
	void setElements( double v ) { gsl_matrix_set_all( data, v ); }
	void setRow( int row, const vector< double >& v );
	void setCol( int col, const vector< double >& v );
	void setRows( const vector< double >& v );
	void setCols( const vector< double >& v );
	
	// set special matrices 
	void setIdentityMatrix( int n );		// identity matrix of dimension n 
	void setDiagonalMatrix( const vector< double >& diag );	// diagonal matrix

	// transpose
	Matrix transpose() const; 
	
	// addition
	Matrix operator+( const Matrix& other ) const;
	Matrix& operator+=( const Matrix& other );
	
	// multiplication of a constant
	Matrix operator*( const double c ) const;
	Matrix& operator*=( const double c ); 
	
	// matrix multiplication
	Matrix operator*( const Matrix& other ) const;	
 	
	// output operator
	friend ostream& operator<<( ostream& os, const Matrix& m );
	
	// load Matrix from a file
	// If readDims = true, read dimensions (first row of the file); o/w, only read data, assuming the correct dimensions are known
	int load( const string& fileName, bool readDims = false );
	
	// save Matrix to a file
	void save( const string& fileName ) const;
};

/*****************************************************
* Mathematical Functions of Vectors and Matrices
******************************************************/

/* table used for log_add */
class log_add_table {
public:
	double delta;				// interval size
	vector< double > x_array;	// x[0], ..., x[N]
	vector< double > T;			// T[0], ..., T[N] where T[k] = log( 1 + exp( x[k] ) )
	
	// constructor
	log_add_table( double a, double b, int N );
};

/* log-add algorithm 
 * log_add( x, y ) = log( exp(x) + exp(y) ), etc. */
double log_add( double p, double q );
double log_add( const vector< double >& p );
double log_add( double p, double q, double r );
 
// log. of a vector 
vector< double > log( const vector< double >& v );

// exp. of a vector
vector< double > exp( const vector< double >& v );

// log. of a matrix
Matrix log( const Matrix& M );

// Eucledian distance bewteen two vectors
double Eucledian_dist( const vector< double >& x, const vector< double >& y );

double infty_transform( double x, double a, double b );
double inverse_infty_transform( double z, double a, double b );
/*****************************************************
* Probablity and Statistics 
******************************************************/

// sum
template< class T >
inline T sum( const vector< T >& v )
{
	T result = 0;
	for ( int i = 0; i < v.size(); i++ ) 
		result += v[ i ];
		
	return result;	
}

// max and min
double max( const vector< double >& v, int &arg );

// mean of a vector
double mean( const vector< double >& x );

// median of a vector
double median( const vector< double >& x );

// standar deviation 
double std_dev( const vector< double >& x );

// ## qcheng75 scaling  = normalization +
//void norm_rescaling( vector< double >& x );

// Pearson correlation of v1 and v2, they must have equal sizes
double correlation( const vector< double >& v1, const vector< double >& v2 );

// check if a vector of real numbers is probablity mass function
bool isPmf( const vector< double > &p );

// sample one outcome from a multinomial distribution
int sampleMul( const gsl_rng * r, const vector< double > &p );

// sample a positive integer from a truncated geometric series: 1, r, r^2, ... r^(n - 1)
int sampleTruncatedGeometric( const gsl_rng* rng, double r, int n );

// mixture of two multinomial distributions, where the weight of the first one is "weight"
vector< double > multinomialMixture( const vector< double >& distr1, const vector< double >& distr2, double weight );

// the threshold value at a given percentage in a vector of numbers
// if highest is true, then the value at the highest p; otherwise, the value at the lowest p
double elementAt( const vector< double >& v, double p, bool highest = true );

// rounding function: round x towards zero to the nearest integer (e.g. trunc (1.5) is 1.0 and trunc (-1.5) is -1.0 )
double trunc( double x );

// logit function & inverse logit function
double logit( double x );
double inv_logit( double x );

/*****************************************************
* Mathematics Miscellaneous 
******************************************************/

// transformation involving weight parameters (w: 0 <= w_i <= 1) s.t. the new parameters in (-inf,+inf)
vector< double > transform( const vector< double > w );		// w -> u
vector< double > inverse_transform( const vector< double > u );		// u -> w

// finite difference method of gradient of a function
void numeric_deriv( gsl_vector* grad, double (*f)( const gsl_vector*, void* ), const gsl_vector* v, void* params, double step );

/*****************************************************
* String & I/O Functions 
******************************************************/

// read parameter file into a <field, value> table. Return 0 if successful, -1 otherwise
// Format: field = value [# ...] where # indicates the start of comments (ignored)
// Ex. kappa = 2.0 # kappa: the transition/transversion bias
int readParams( const string& file, map< string, string >& params, const char* comment = "#" );

// split a string
void split( vector< string >& result, const string& str, const string& sep );

// convert a string in the format of "{bcd,Kr,hb}" (separated by , or space) into a vector of strings
void stringToVector( vector< string >& result, const string& str, const string& leftBoundary = "{", const string& rightBoundary = "}", const string& sep = " ," );

// convert a string in the format of "{0.3,0.2,0.2,0.3}" (separated by , or space) into a vector of numbers
void stringToVector( vector< double >& result, const string& str, const string& leftBoundary = "{", const string& rightBoundary = "}", const string& sep = " ," );

// control the output format
class IO_Ctrl {
public:
	static string vector_delimiter;
	static string map_delimiter;	
	static string field_separator;
};

// output map< T1, T2 >
template< class T1, class T2 >
ostream& operator<<( ostream& os, const map< T1, T2 >& table )
{
	class map< T1, T2 >::const_iterator it;
	for ( it = table.begin(); it != table.end(); it++ ) { 
		os << it->first << IO_Ctrl::field_separator << it->second << IO_Ctrl::map_delimiter;    
	}

	return os;
}

// output vector< T >
template< class T >
ostream& operator<<( ostream& os, const vector< T >& v )
{
	if ( !v.empty() ) os << v[ 0 ];
	for ( int i = 1; i < v.size(); i++ ) {
		os << IO_Ctrl::vector_delimiter << v[ i ];
	}
	
	return os;	
}

// print matrix
template< class T >
void printMatrix( ostream& os, T** &x, int m, int n )
{	
	for ( int i = 0; i < m; i++ ) {
		for ( int j = 0; j < n; j++ ) {
			os << x[ i ][ j ] << " ";	
		}	
		os << endl;
	}	
}

// read vector
template< class T >
int readVector( const string& file, vector< T >& v )
{
	v.clear();
	
	ifstream fin( file.c_str() );
	if ( !fin ) {
		cerr << "Cannot open " << file << endl;
		exit( 1 );
	}	
	T elem;
	while ( fin >> elem ) v.push_back( elem );
	
	return 0;
} 	


/*****************************************************
* Data Structure
******************************************************/

// transformation of an array to a vector
template< class T >
void array2vector( const T* arr, int n, vector< T >& v )
{
	v.clear();
	
	for ( int i = 0; i < n; i++ ) v.push_back( arr[ i ] );
}

// subset of a vector
template< class T >
int vectSubset( const vector< T >& v, const vector< int >& indices, vector< T >& results )
{
	results.clear();
	
	for ( int i = 0; i < indices.size(); i++ ) {
		if ( indices[ i ] >= v.size() ) return RET_ERROR;
		results.push_back( v[ indices[ i ] ] );	
	}	
	
	return 0;
}

// complement of an index set
int indexCompl( const vector< int >& indices, int n, vector< int >& complIndices );

/* TreeNode class: node of a tree */
template< class T >
class TreeNode {
public:
	typedef class map< int, TreeNode< T >* >::iterator iterator;
	typedef class map< int, TreeNode< T >* >::const_iterator const_iterator;
private:
	T data;
	map< int, TreeNode< T >* > branches;
public:
	// constructor
	TreeNode() {}
	TreeNode( const T& x ) { data = x; }
	void copy( const TreeNode< T >& other ) { data = other.data; branches = other.branches; }
	TreeNode( const TreeNode< T >& other ) { copy( other ); }
	
	// destructor
	~TreeNode() {}
	
	// assignment
	TreeNode< T >& operator=( const TreeNode< T >& other ) { copy( other ); return *this; }
	
	// information of the node
	int nChildren() const { return branches.size(); }
	bool isLeaf() const { return branches.empty(); }
	
	// access the data
	const T& get() const { return data; }
	T& get() { return data; }
	void set( const T& x ) { data = x; }
	
	// access a child node
	TreeNode< T >* &operator[]( int label ) { return branches[ label ]; }
	TreeNode< T >* find( int label ) const { return branches.find( label )->second; }
	
	// begin and end of the node
	iterator begin() { return branches.begin(); }
	const_iterator begin() const { return branches.begin(); } 
	iterator end() { return branches.end(); }
	const_iterator end() const { return branches.end(); }	
	
	// clear the node
	void clear() { branches.clear(); }
	
	// Block
	/*struct Block {
		T field1;
		map< int, TreeNode< T >* > field2;
		Block* next;
	};
	
	// free list
	static Block* freeList;
	
	// new 
	static void* operator new( size_t bytes );
	
	// delete
	static void operator delete( void* ptr );*/
};

// TreeNode: initialization of freeList
/*template< class T >
class TreeNode< T >::Block* TreeNode< T >::freeList = 0;

// TreeNode: new operator
template< class T >
void* TreeNode< T >::operator new( size_t bytes )
{
	//cout << "Overloaded new operator" << endl;
	if ( !freeList ) return malloc( sizeof( Block ) );
	
	Block* res = freeList;
	freeList = freeList->next;
	return res;
}

// TreeNode: delete operator
template< class T >
void TreeNode< T >::operator delete( void* ptr )
{
	//cout << "Overloaded delete" << endl;
	( (Block*)ptr )->next = freeList;
	freeList = (Block*)ptr;	
}*/

/* Tree class: a tree */
template< class T >
class Tree {
protected:
	typedef class TreeNode< T > Node_type;
	
	Node_type* root;
public:
	// constructors
	Tree() { root = NULL; }
	Tree( Node_type* _root ) : root( _root ) { assert( root ); }
	Tree( const T& x ) { root = new Node_type( x ); }
	
	// destructor
	//~Tree() { clear(); }
			
	// number of nodes
	int size() const { 
		if ( root ) return sizeAux( root ); 
		else return 0;
	}
		
	//int depth() const;	// max depth of the tree
	//int nLeaves() const;	// number of leaves of the tree
	
	// access methods
	Node_type* getRoot() const { return root; }
	void setRoot( Node_type* _root ) { root = _root; }
	
	// access the node, given the path to the node, return NULL if not found
	Node_type* getNode( const vector< int >& path );
	
	// access the data at the node, given the path to the node	
	const T& operator[]( const vector< int >& path ) const;
	T& operator[]( const vector< int >& path );
		
	// clear the tree
	//void clear() { if ( root ) clearAux( root ); }
	
	// search an element in the tree by DFS, return NULL if not found
	Node_type* DFS( const T& x ) const { 
		if ( root ) return DFSAux( x, root ); 
		else return NULL;
	}
	
	// print the tree
	void print( ostream& os ) const { if ( !root ) return; os << "(-,"; printAux( os, root, 0 ); }
private:
	// recursive routine of size()
	int sizeAux( Node_type* node ) const;
	
	// recursive routine of clear()
	//void clearAux( Node_type* node );
	
	// recursive routine of DFS()
	Node_type* DFSAux( const T& x, Node_type* node ) const;
	
	// recursive routine of print()
	void printAux( ostream& os, Node_type* node, int level ) const;
};

// Tree: recursive routine of size()
template< class T >
int Tree< T >::sizeAux( Node_type* node ) const
{
	int result = 1;	
	class TreeNode< T >::const_iterator it;
	for ( it = node->begin(); it != node->end(); it++ ) {
		result += sizeAux( it->second );	
	}
	
	return result;
}

// Tree class: access the node, given the path to the node
template< class T >
class Tree< T  >::Node_type* Tree< T >::getNode( const vector< int >& path )
{
	Node_type* curr = root;
	for ( int i = 0; i < path.size(); i++ ) {
		Node_type* child = (*curr)[ path[ i ] ];
		if ( child ) 
			curr = child;
		else 
			return NULL;
	}	
	
	return curr;
}

// Tree class: access the data at the node, given the path to the node	
template< class T >
const T& Tree< T >::operator[]( const vector< int >& path ) const
{
	Node_type* node = getNode( path );
	assert( node );
	return node->get();
}

// Tree class: access the data at the node, given the path to the node	
template< class T >
T& Tree< T >::operator[]( const vector< int >& path ) 
{
	Node_type* node = getNode( path );
	assert( node );
	return node->get();
}

// Tree class: recursive version of clear()
/*template< class T >
void Tree< T >::clearAux( Node_type* node )
{
	if ( node == NULL ) return;
	class TreeNode< T >::iterator it;
	for ( it = node->begin(); it != node->end(); it++ ) {
		clearAux( it->second );	
	}	
	delete node;
}*/

// Tree class: resursive routine of DFS()
template< class T >
class Tree< T >::Node_type* Tree< T >::DFSAux( const T& x, Node_type* node ) const
{
	if ( node->get() == x ) 
		return node;
	else {
		class TreeNode< T >::iterator it;
		for ( it = node->begin(); it != node->end(); it++ ) {
			DFSAux( x, it->second );	
		}		
	}	
}

// Tree: recursive routine of print()
template< class T >
void Tree< T >::printAux( ostream& os, Node_type* node, int level ) const
{
	os << node->get() << ")" << endl;
	class TreeNode< T >::iterator it;		
	for ( it = node->begin(); it != node->end(); it++ ) {
		for ( int i = 0; i < level + 1; i++ ) { os << "\t"; }
		os << "(" << it->first << ",";	
		printAux( os, it->second, level + 1 );
	}
}

#endif
