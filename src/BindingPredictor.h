#ifndef BINDING_PREDICTOR_H
#define BINDING_PREDICTOR_H

#include "SeqAnnotator.h"

/*****************************************************
* Predicting Binding of a TF to DNA Sequences
******************************************************/

/* FactorIntFunc class: distance-dependent function of TF-TF interaction  */
class FactorIntFunc {
public:
	// compute the factor interaction, given the normal interaction (when they are close enough)
	virtual double compFactorInt( double normalInt, double dist, bool orientation ) const = 0;

    // the maximum distance beyond which there is no interaction
    virtual double getMaxDist() const = 0;
};

/* FactorIntFuncBinary class: binary distance function */
class FactorIntFuncBinary : public FactorIntFunc {
public:
	// constructors
	FactorIntFuncBinary( double _distThr, double _orientationEffect = 1.0 ) : distThr( _distThr ), orientationEffect( _orientationEffect ) { assert( distThr > 0 ); }

	// compute the factor interaction
	double compFactorInt( double normalInt, double dist, bool orientation ) const;

    // the maximum distance beyond which there is no interaction
    double getMaxDist() const {
    	return distThr;
    }    
private:
	double distThr;		// if distance < thr, the "normal" value; otherwise 1 (no interaction)
	double orientationEffect;	// the effect of orientation: if at different strands, the effect should be multiplied this value	
};

/* FactorIntFuncLinear class: linear distance function */
class FactorIntFuncLinear : public FactorIntFunc {
public:
	// constructors
	FactorIntFuncLinear( double _distThr, double _range, double _orientationEffect = 1.0 ) : distThr( _distThr ), range( _range ), orientationEffect( _orientationEffect ) { assert( distThr > 0 && range > 0 ); }

	// compute the factor interaction
	double compFactorInt( double normalInt, double dist, bool orientation ) const;

    // the maximum distance beyond which there is no interaction
    double getMaxDist() const {
    	return distThr;
    }   
private:
    double distThr;    // distance threshold
    double range;   // range of maximum interaction, outside which the interaction starts to decline
	double orientationEffect;	// the effect of orientation: if at different strands, the effect should be multiplied this value	
};

/* FactorIntFuncPeriodic class: periodic distance function */
class FactorIntFuncPeriodic : public FactorIntFunc {
public: 
    // constructors
    FactorIntFuncPeriodic( double _distThr, double _period, double _spacingEffect, double _phi, double _phi0 ) : distThr( _distThr ), period( _period ), spacingEffect( _spacingEffect ), phi( _phi ), phi0( _phi0 ) { assert( distThr > 0 && period > 0 && spacingEffect > 0 );  } 
    
    // compute the factor interaction 
    double compFactorInt( double normalInt, double dist, bool orientation ) const;
    
    // the maximum distance beyond which there is no interaction
    double getMaxDist() const {
        return distThr; 
    }
private: 
    double distThr;    // distance threshold
    double period;  // period in bps
    double spacingEffect;    // the constant term due to DNA looping (exponential of the average free energy of DNA looping)
    double phi;     // phase of DNA looping when both sites are in the same strand
    double phi0;    // phase of DNA looping when sites are in the different strand
}; 

/* FactorIntFuncGeometric class: distance function decays geometrically (but never less than 1) */
class FactorIntFuncGeometric : public FactorIntFunc {
public:
	// constructors
	FactorIntFuncGeometric( double _distThr, double _spacingEffect, double _orientationEffect ) : distThr( _distThr ), spacingEffect( _spacingEffect ), orientationEffect( _orientationEffect ) { assert( distThr > 0 ); }

	// compute the factor interaction
	double compFactorInt( double normalInt, double dist, bool orientation ) const;

    // the maximum distance beyond which there is no interaction
    double getMaxDist() const {
        return distThr;
    } 
private:
	double distThr;		// if distance < thr, the "normal" value; otherwise decay with distance (by parameter spacingEffect)
	double spacingEffect;		// the effect of spacing
	double orientationEffect;	// the effect of orientation: if at different strands, the effect should be multiplied this value
};

/* BindingFunc class: predict the binding of a TF to a DNA sequence */
class BindingFunc {
public:
	// constructors
	BindingFunc( const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< double >& _maxBindingWts, const Matrix& _factorIntMat );
	
	// predict the binding of a given TF to a given sequence (its site representation, sorted by the start positions and the association constants computed) 
	double predictBinding( const SiteVec& targetSites, int factorIdx );

	void releaseFactorIntMat(); //qcheng75
	void clear();//qcheng75
private:
	// TF binding motifs
	const vector< Motif >& motifs; 	

	// function to compute distance-dependent TF-TF interactions 
	const FactorIntFunc* intFunc;		

	// TFs and their interaction
	vector< double > maxBindingWts;		// binding weight of the strongest site for each TF, q(S_max) = K(S_max) [TF]
	Matrix factorIntMat;	// TF-TF interaction matrix
		
	// the sequence whose binding is to be predicted
	SiteVec sites;
	vector< int > boundaries;   // left boundary of each site beyond which there is no interaction
	
	// intermediate computational results
	vector< double > bindingWts; 
	vector< double > Z;
    vector< double > Zt;
					
	// compute the partition function 
	double compPartFunc();
	
	// compute the Boltzmann weight of binding of a site
	double compBindingWt( const Site& a ) const;
	
	// compute the TF-TF interaction between two occupied sites
	double compFactorInt( const Site& a, const Site& b ) const;		
};

/* BindingPar class: the binding parameters */
// Whether the parameters are free or fixed is specified by a 5-bit vector: maxBindingWts of the in-factor(s); maxBindingWts of the out-factors; self-cooperativity (of in-factors); interactions between different factors; exp. ratios
class BindingPar {
public:
	// constructors 
	BindingPar() : outFactorIntMat() {}
	BindingPar( int _nInFactors, int _nOutFactors );		
	BindingPar( const vector< double >& _maxBindingWts, const vector< vector< double > >& _inFactorIntMat, const Matrix& _outFactorIntMat, const vector< double >& _expRatios );
	BindingPar( const vector< double >& pars, const BindingPar& fixed, const vector< bool >& option );	// construct from a "flat" vector of free parameters and a BindingPar object holding the fixed values. The meaning of option is the same as the one used in getFreePars()
	BindingPar( const vector< double >& pars, int nInFactors, int nOutFactors, const vector< bool >& option );	// ## qcheng75
	BindingPar( const vector< double >& pars, const BindingPar& fixed, const vector< bool >& option, bool isParBound );
	void copy( const BindingPar& other ) { maxBindingWts = other.maxBindingWts; inFactorIntMat = other.inFactorIntMat; outFactorIntMat = other.outFactorIntMat; expRatios = other.expRatios; }
	BindingPar( const BindingPar& other ) { copy( other ); }
	BindingPar( int _nInFactors, int _nOutFactors, const vector< double >& _maxBindingWts ); // ## qcheng75
	
	// assignment
	BindingPar& operator=( const BindingPar& other ) { copy( other ); return *this; }	
	
	// access methods
	int nInFactors() const { return inFactorIntMat.size(); }
	int nOutFactors() const { return maxBindingWts.size() - nInFactors(); }
	
	// create a copy of itself and return the pointer to the memory
// 	BindingPar* copy_construct() const;

	// copy from a source object, i.e. update the content
// 	void copy( const BindingPar* src ); 
	
	// construct the factor interaction matrix 
	void createFactorIntMat( Matrix& factorIntMat ) const;
	
	// get the free parameters. The argument option specifies which parameters are free.  
	void getFreePars( vector< double >& pars, const vector< bool >& option ) const; 

	// get the box-boundary parameter constraints ##qcheng75
	void getParsLB( vector< double >& parsLB, const vector< bool >& option ) const;
	void getParsUB( vector< double >& parsUB, const vector< bool >& option ) const;
	
	// print the parameters
	void print() const;

	void printParaModel2Stream(std::ofstream& fout) const;// ## qcheng75
	static string GLOBAL_INIT_PARAMETER_MODEL_FILE_NAME; // ## qcheng75
	void readFreeParameterModelSource(string para_model_file, int _nInFactors, int _nOutFactors ); // ## qcheng75
	void readFreeParameterModelSource(std::ifstream& _fdata, int _nInFactors, int _nOutFactors ); // ## qcheng75
	void clear(); //qcheng75

	// parameters
	vector< double > maxBindingWts;			// binding weight of the strongest site for each TF, q(S_max) = K(S_max) [TF]
	vector< vector< double > > inFactorIntMat;	// interactions among all in-factors (the factors in experiments): {(0,0); (1,0), (1,1); (2,0), (2,1), (2,2); ...}
	Matrix outFactorIntMat; 		// interactions between one out-factor with each in-factor, one row corresponds to one out-factor
	vector< double > expRatios; 		// constant factor of measurement to prediction for each experiment 
		
	static int nparts;		// number of different parts (types of parameters, binding term, interaction term, etc.)
	static double initMaxBindingWt;	// initial value of maximum binding weight (default)
	static double initInt;		// initial interaction (default)
	static double initExpRatio;		// initial ratio of an experiment (default)
	static double log10_min_wt;		// log. of the min. maxBindingWt
	static double log10_max_wt;		// log. of the max. maxBindingWt
	static double log10_min_int;	// log. of the min. interaction
	static double log10_max_int;	// log. of the max. interaction
	static double log10_min_ratio;		// log. of the min expRatio
	static double log10_max_ratio;		// log. of the max expRatio
	static double log10_range_wt;		// log10 of the range of maxBindingWt 
	static double wt_step;		// step of maxBindingWt (log10)
	static double int_step;		// step of interaction (log10)
	static double ratio_step;	// step of expRatio
    
	static int FIX_COOP_PROB; //cooperativity range change + option[2] setting qcheng75
	static int NELDER_MEAD_CONSTRAINED;
	static int NEWTON_KTCONSTRAINED;
	static int PENALTY_BASED_CONSTRAINED;
	static double PENALTY_COEFFIENT;
	static int IS_LOG_SCALE_STEP;

    static double min_float_value;      // minimum floating point value
    static double max_float_value;      // maximum floating point value
private: 
	// assign the parameters whose values are stored in a single vector (all parameters)
// 	void assignParams( const vector< double >& pars, int nInFactors, int nOutFactors ); 
};

/* BasePredictor class: the public interface of sequence-to-binding predictor */
class BasePredictor {
public: 
	// constructors
	BasePredictor( const vector< Motif >& _motifs, const vector< double >& _energyThrs, const vector< vector< Sequence > >& _seqs, const vector< vector< double > >& _bindingData );
	BasePredictor( const vector< Motif >& _motifs, const vector< double >& _energyThrs, const vector< vector< Sequence > >& _seqs, const vector< vector< double > >& _bindingData, const vector< double >& _maxBindingWts  ); // ## qcheng75
	BasePredictor( const vector< Motif >& _motifs, const vector< double >& _energyThrs, const vector< vector< Sequence > >& _seqs, const vector< vector< double > >& _bindingData,  const vector< double >& _allWts, int nInFactors, int nOutFactors, bool isAllParam  ); // ## qcheng75

	// access methods
	const vector< Motif >& getMotifs() const { return motifs; }
	const vector< vector< double > >& getBndingData() const { return bindingData; }	
	int nInFactors() const { return seqs.size(); }
	int nOutFactors() const { return motifs.size() - seqs.size(); }
	const BindingPar& getPar() const { return par_model; }
	double getObj() const { return obj_model; }
	
	// train the model
	virtual int train() = 0;	

	// analyze all features (out-factors) by their effects on in-factor binding. Only relevant if nInFactors = 1. Effect is measured by the average percent increase of binding over all sequences relative to the in-factor alone. Return the total effect of all out-factors. 
	virtual int analyzeFeatures( vector< vector< double > >& featureEffects ) const = 0;
		
	// predict binding of a given factor to a given sequence 
	virtual double predict( const Sequence& targetSeq, int factorIdx ) const = 0; 
				
	// test the model, perfOption = 0: Pearson correlation
	double test( const vector< vector< Sequence > >& testSeqs, const vector< vector< double > >& testBindingData, int perfOption = 0 ) const; 

	void clear();//qcheng75
	~BasePredictor();//qcheng75
protected:
	// fixed parameters
	vector< Motif > motifs;		// TF binding motifs
	vector< double > energyThrs;	// energy thresholds

	// training data
	const vector< vector< Sequence > >& seqs;	// sequences: one row per experiment (in-factor)
	const vector< vector< double > >& bindingData;		// binding intensities of sequences: one row per experiment (in-factor)
	
	// parameters of the model and the value of the objective function
	BindingPar par_model;
	double obj_model;				
};

/* BindingPredictor class: the thermodynamic model of sequence-to-binding */
class BindingPredictor : public BasePredictor {
public:
	// constructors
	BindingPredictor( const vector< Motif >& _motifs, const vector< double >& _energyThrs, const vector< vector< Sequence > >& _seqs, const vector< vector< double > >& _bindingData, const FactorIntFunc* _intFunc );
	BindingPredictor( const vector< Motif >& _motifs, const vector< double >& _energyThrs, const vector< vector< Sequence > >& _seqs, const vector< vector< double > >& _bindingData, const FactorIntFunc* _intFunc, const vector< double >& _maxBindingWts  ); // ## qcheng75
	BindingPredictor( const vector< Motif >& _motifs, const vector< double >& _energyThrs, const vector< vector< Sequence > >& _seqs, const vector< vector< double > >& _bindingData, const FactorIntFunc* _intFunc, const vector< double >& _allWts, int nInFactors, int nOutFactors, bool isAllParam ); // ## qcheng75

	// the objective function to be minimized
	double objFunc( const BindingPar& par ) const;
	
	// training the model
	int train( const BindingPar& par_init, int nAlternations = 1 ); 	// training with the initial values given
	int train();	// automatic training: first estimate the initial values, then train
	
	// analyze all features (out-factors) by their effects in each experiment. 
    // Effect is measured by the average percent increase of binding over all sequences relative to the in-factor alone.  
	int analyzeFeatures( vector< vector< double > >& featureEffects ) const;
		
	// predict binding of a given factor to a given sequence 
	double predict( const Sequence& targetSeq, int factorIdx ) const; 
		
	void clear();//qcheng75
	~BindingPredictor();//qcheng75

	// obtain the box-boundary constraints ##qcheng75
	void getParBoundaryConstraints( vector< double >& parsLB,  vector< double >& parsUB ) const;
	
	// option for model training: which parameters are free (to be estimated), the meaning is the same as the option of BindingPar.
	static vector< bool > parOption;

	// SA parameters
	static int SA_n_tries;
	static int SA_iters_fixed_T;
	static double SA_step_size;
	static double SA_K;
	static double SA_T_init;
	static double SA_mu_T;
	static double SA_T_min;			

	// optimization parameter ## SS
	static int GLOBAL_NRANDSTARTS; // ## SS
private:	
	// function to compute TF-TF interactions 
	const FactorIntFunc* intFunc;					
				
	// sequence annotation
	vector< vector< SiteVec > > seqSites;		// the extracted sites for all sequences
		
	// create the binding function
	BindingFunc* createBindingFunc( const BindingPar& par ) const;
		
	// different objective functions
// 	double compRMSE( const BindingPar& par ) const;		// root mean square error between predicted and observed expressions
 	double compCorr( const BindingPar& par ) const;		// Pearson correlation between predicted and observed expressions

    // check if parameter values are valid
    bool testPar( const BindingPar& par ) const;
    
 	// sample random parameter values (only those parameters that are not free)
 	void randSamplePar( const gsl_rng* rng, BindingPar& par ) const; 
 	
	// minimize the objective function, using the current model parameters as initial values
	int simplex_minimize( BindingPar& par_result, double& obj_result ) const;	// simplex	
	int gradient_minimize( BindingPar& par_result, double& obj_result ) const;	// methods using gradient (BFGS, conjugate gradient)
 	int SA_minimize( BindingPar& par_result, double& obj_result ) const;	// simulated annealing 	

 	// ## qcheng75
 	int gradient_minimize_KTconstrained( BindingPar& par_result, double& obj_result ) const;
};

// the objective function and its gradient of BindingPredictor::simplex_minimize or gradient_minimize
double gsl_obj_f( const gsl_vector* v, void* params );
void gsl_obj_df( const gsl_vector* v, void* params, gsl_vector* grad ); 
void gsl_obj_fdf( const gsl_vector* v, void* params, double* result, gsl_vector* grad ); 

/* SAState class: the state of simulated annealing search space (including both the actual variable search and the additional parameter) */
class SAState {
public:
	// constructor
	SAState( const BindingPar& _par, const BindingPredictor* _predictor ) : par( _par ), predictor( _predictor ) {}
	
	BindingPar par;	// binding parameters
	const BindingPredictor* predictor;	// hold the data and options
};
 
// functions used by gsl_siman_solve()
double SA_Efunc( void* xp );
double SA_metric( void* xp, void* yp );
void SA_step( const gsl_rng* r, void* xp, double step_size );
void SA_print( void *xp );
void SA_copy( void* src, void* dest );
void* SA_copy_construct( void* xp );
void SA_destroy( void* xp );

/* StepwisePredictor: stepwise method (model selection) for predicting binding */
class StepwisePredictor : public BasePredictor {
public:
	// constructor
	StepwisePredictor( const vector< Motif >& _motifs, const vector< double >& _energyThrs, const vector< vector< Sequence > >& _seqs, const vector< vector< double > >& _bindingData, const FactorIntFunc* _intFunc );
	
	// access methods
	const vector< bool >& getIndicators() const { return indicators; }
	
	// stepwise training
	int train();	

	// analyze all features (out-factors) by their effects on in-factor binding. 
	int analyzeFeatures( vector< vector< double > >& featureEffects ) const;
		
	// predict binding of a given factor to a given sequence 
	double predict( const Sequence& targetSeq, int factorIdx ) const; 
	
	// threshold for entry and exit test
	static double testThr;	
private:		
	// function to compute TF-TF interactions 
	const FactorIntFunc* intFunc;			
		
	// the best model and the set of motifs chosen (indicators)
	vector< bool > indicators;	
	BindingPredictor* model;
	
	// evaluate a motif set: return the model of this set
	BindingPredictor* evalSet( const vector< bool >& indSet ); 
	
	// forward test: adding a new factor to the current factor set
	BindingPredictor* forwardTest( const vector< bool >& currSet, int testFactor );
	
	// backward test: removing an existing factor from the current factor set
	BindingPredictor* backwardTest( const vector< bool >& currSet, int testFactor );	
	
	// map parameters: from the model and its associate indicator set to the parameters in the complete set of motifs
	void mapPar( const BindingPredictor* bp, const vector< bool >& indSet, BindingPar& par );
};

/*****************************************************
* Selecting Features (Motifs of Out-factors)
******************************************************/

class FeatureSelector {
public:
    // constructor
    FeatureSelector( const vector< Motif >& _motifs, const vector< double >& _energyThrs, const vector< vector< Sequence > >& _seqs, const vector< vector< double > >& _bindingData, const FactorIntFunc* _intFunc );

    // access methods
 	int nInFactors() const { return seqs.size(); }
	int nOutFactors() const { return motifs.size() - seqs.size(); }
    const vector< bool >& getIndicators() const {
         return indicators;
    }
    
    // select features
    int select();

    static double effectThr;    // the effect threshold for determining if a feature should be chosen
private:
	// the set of motifs to be chosen from
	vector< Motif > motifs;		// TF binding motifs
	vector< double > energyThrs;	// energy thresholds

	// data
	const vector< vector< Sequence > >& seqs;	// sequences: one row per experiment (in-factor)
	const vector< vector< double > >& bindingData;	// binding intensities of sequences: one row per experiment (in-factor)

	// function to compute TF-TF interactions 
	const FactorIntFunc* intFunc;	
    
    // result: whether a motif/feature is chosen
    vector< bool > indicators;
};

/*****************************************************
* Cross Validation of the Binding Prediction Methods
******************************************************/

/* BindingCrossValidator class: cross validation of the binding prediction methods */
class BindingCrossValidator {
public: 
	// constructor
	BindingCrossValidator( const vector< vector< Sequence > >& _seqs, const vector< vector< string > >& _names, const vector< vector< double > >& _bindingData ); 	
	BindingCrossValidator( const vector< vector< Sequence > >& _seqs, const vector< vector< double > >& _bindingData ); 
		
	// access methods
	const vector< vector< vector< Sequence > > >& getTrainSeqs() const { return trainSeqs; }
	const vector< vector< vector< double > > >& getTrainBindingData() const { return trainBindingData; }
	const vector< vector< vector< Sequence > > >& getTestSeqs() const { return testSeqs; }
	const vector< vector< vector< double > > >& getTestBindingData() const { return testBindingData; }
	
	// prepare for the training and testing data
	void partition( int nFolds ); 
		
	// cross validation of theromodyanmic model, modelOption = 0: basic model; = 1: stepwise model
	void runThermoModel( vector< double >& results, const vector< Motif >& motifs, const vector< double >& energyThrs, const FactorIntFunc* intFunc, bool featureSelection = false, int modelOption = 0 );
	
	void print() const;
	
	// performance option: 0 - Pearson correlation
	static int perfOption;	
private: 	
	// data
	vector< vector< Sequence > > seqs;	// sequences: one row per experiment (in-factor)
	vector< vector< string > > names;	// names: one row per experiment	
	vector< vector< double > > bindingData;		// binding intensities of sequences: one row per experiment (in-factor)
			
	
	// partitioned training and testing data: [ k ][ i ][ j ] - k-th fold, i-th experiment, j-th data point
	vector< vector< vector< Sequence > > > trainSeqs;	
	vector< vector< vector< string > > > trainNames;		
	vector< vector< vector< double > > > trainBindingData;	
	vector< vector< vector< Sequence > > > testSeqs;	
	vector< vector< vector< string > > > testNames;					
	vector< vector< vector< double > > > testBindingData;	
};

#endif
