#include "BindingPredictor.h"

int main( int argc, char* argv[] ) 
{
	// command line processing
	string motifFile, seqFile, dataFile;
	string testSeqFile, testDataFile;	// test data [optional]
	int nExps = 1;	// number of experiments
    bool featureSelection = false;  // whether do the feature selection
    double effectThr = 0.5;        // feature effect threshold
	int K = 0;		// K-fold cross-validation
	int modelOption = 0; 	// model option: 0 - basic model; 1 - stepwise model
    int coopOption = 1;     // 0 - no cooperativity at all; 1 - no self-cooperativity, but hetero-cooperativity; 2 - allow all cooperativities
    int perfOption = 0;		// option for measuring performance: 0 - Pearson correlation
	int interactionOption = 0;	// option for modeling factor-factor interaction: 0 - binary; 1 - linear; 2 - periodic; 3 - geometric
	double distThr = 150.0;		// distance threshold for interaction	
	double spacingEffect = M_E;		// effect of space/distance on interaction
	double orientationEffect = 1.0;		// orientation bias: multiply this constant if the factors bound in different orientations
    double range = 50.0;    // Linear model: range of maximum interaction, outside which the interaction starts to decline
    double period = 10.0;       // Periodic model: period 
    double phi = 0;     // Periodic model: phase angle if two sites are in the same strand
    double phi0 = 0;    // Periodic model: phase angle if two sites are in different strands
	double energyThr = 7.0;		// energy of a site should not be larger than this value
	string trainPredictionFile, testPredictionFile; // the predictions for training and testing sequences respectively	
	string paraModelFile, featureStatiscsTrFile, featureStatiscsTeFile, pminputfile; // ## qcheng75 parameter model will be saved to the file
	int nrandstarts = 1; //## SS
	int fixcooperror = 0, chiprescaling=0, isFeatureStat=0, isFeatureCollectionOnly=0, islogscalestep=0; //## qcheng75
	double fFreePMaxBinding=0, penaltyVal4BoundedOpt=100000.0;
	vector<double> allParaValues;				// Input weightedBinding Parameter and test  qcheng75
	bool isTestOnly=false, isFreePGiven=false, areAllPGiven=false;			// for decision of counting binding sites,...  qcheng75
	vector<double> motifSiteThr; // ##qcheng75
	for ( int i = 1; i < argc; ++i ) {
		if ( !strcmp( "-m", argv[ i ] ) )
			motifFile = argv[ ++i ];
		else if ( !strcmp( "-s", argv[ i ] ) )
			seqFile = argv[ ++i ];
		else if ( !strcmp( "-d", argv[ i ] ) )
			dataFile = argv[ ++i ];				
		else if ( !strcmp( "-n", argv[ i ] ) )
			nExps = atoi( argv[ ++i ] );			
		else if ( !strcmp( "-ts", argv[ i ] ) )
			testSeqFile = argv[ ++i ];	
		else if ( !strcmp( "-td", argv[ i ] ) )
			testDataFile = argv[ ++i ];	
		else if ( !strcmp( "-fs", argv[ i ] ) )
			featureSelection = true;		
        else if ( !strcmp( "-fet", argv[ i ] ) )
			effectThr = atof( argv[++i] );	            
		else if ( !strcmp( "-cv", argv[ i ] ) )
			K = atoi( argv[ ++i ] );
		else if ( !strcmp( "-mo", argv[ i ] ) )
			modelOption = atoi( argv[ ++i ] );
		else if ( !strcmp( "-co", argv[ i ] ) )
			coopOption = atoi( argv[ ++i ] );			
		else if ( !strcmp( "-po", argv[ i ] ) )
			perfOption = atoi( argv[ ++i ] );
		else if ( !strcmp( "-io", argv[ i ] ) )
			interactionOption = atoi( argv[ ++i ] );			
		else if ( !strcmp( "-dt", argv[ i ] ) )
			distThr = atof( argv[ ++i ] );
		else if ( !strcmp( "-se", argv[ i ] ) )
			spacingEffect = atof( argv[ ++i ] );
		else if ( !strcmp( "-oe", argv[ i ] ) )
			orientationEffect = atof( argv[ ++i ] );
        else if ( !strcmp( "-r", argv[ i ] ) )
			range = atof( argv[ ++i ] );
		else if ( !strcmp( "-T", argv[ i ] ) )
			period = atof( argv[ ++i ] );	
        else if ( !strcmp( "-phi", argv[ i ] ) )
			phi = atof( argv[ ++i ] );  
        else if ( !strcmp( "-phi0", argv[ i ] ) )
			phi0 = atof( argv[ ++i ] );    
		else if ( !strcmp( "-et", argv[ i ] ) )
			energyThr = atof( argv[ ++i ] );
		else if ( !strcmp( "-llrt", argv[ i ] ) ){ // ## qcheng75
			char* strMotifSiteThr = argv[ ++i ];
			Motif::SITE_THRESHOLD_BY_LLRthres=1;
			//split it by ';' => to obtain the motifs' llr score at at least pval=0.05
			char* pch=strtok (strMotifSiteThr,";");
			while (pch != NULL){
				motifSiteThr.push_back(atof(pch));
				pch = strtok (NULL, ";");
			}
		}else if ( !strcmp( "-p", argv[ i ] ) )
			trainPredictionFile = argv[++i];
		else if ( !strcmp( "-tp", argv[ i ] ) )
			testPredictionFile = argv[++i];            
		else if ( !strcmp( "-rr", argv[ i ] ) ) // ## SS
		  nrandstarts = atoi(argv[++i]);        // ## SS
		else if ( !strcmp( "-cc", argv[ i ] ) ) // ## qcheng75  Correct the Cooperativity range and parameter setting
			fixcooperror = atoi(argv[++i]);
		else if ( !strcmp( "-pm", argv[ i ] ) ) // ## qcheng75 save Parameter Model
			paraModelFile = argv[++i];
		else if ( !strcmp( "-cs", argv[ i ] ) ) // ## qcheng75 Chip data reScaling
			chiprescaling = atoi(argv[++i]);
		else if ( !strcmp( "-fstr", argv[ i ] ) ){ // ## qcheng75 Feature Statistics on TRaining data f(seq,motif)=<information content, num of bound sites, LLRenergy, TRAP score>
			featureStatiscsTrFile = argv[++i];
			isFeatureStat = 1;
		}
		else if ( !strcmp( "-fste", argv[ i ] ) ){ // ## qcheng75 Feature Statistics on TEsting data  f(seq,motif)=<information content, num of bound sites, LLRenergy, TRAP score>
			featureStatiscsTeFile = argv[++i];
			isFeatureStat = 1;
		}
		else if ( !strcmp( "-pin", argv[ i ] ) ){ // ## qcheng75 taking Parameter model file as the parameter INputs
			pminputfile = argv[++i];
		}else if ( !strcmp( "-bw", argv[ i ] ) ) {
			// Parameter: max binding weight ##qcheng75
			isFreePGiven=true;
			fFreePMaxBinding = atof( argv[ ++i ] );
		}else if ( !strcmp( "-ot", argv[ i ] ) || !strcmp( "-onlyTest", argv[ i ] )) {//qcheng75
			isTestOnly = true;
		}else if (!strcmp( "-pw", argv[ i ] ) ){
			penaltyVal4BoundedOpt = atof( argv[ ++i ] );
		}else if (!strcmp( "-ls", argv[ i ] ) ){
			islogscalestep = atoi( argv[ ++i ] );
		}else if (!strcmp( "-aw", argv[ i ] ) ){
			char* pch=strtok (argv[ i ],"_");
			while (pch != NULL){
			    allParaValues.push_back(atof(pch));
			    pch = strtok (NULL, "_");
			}
			areAllPGiven=true;
		}
	}
	if ( motifFile.empty() || seqFile.empty() || dataFile.empty() ) {
		cerr << "Usage: " << argv[ 0 ] << " -m motifFile -s seqFile -d dataFile [-n nExps -ts testSeqFile -td testDataFile -fs -fet effectThr -cv K -mo modelOption -co coopOption -po perfOption -io interactionOption -dt distThr -se spacingEffect -oe orientationEffect -r range -T period -phi phi -phi0 phi0 -et energyThr -p trainPredictionFile -tp testPredictionFile] -rr 20 -cc [1/0] -pm file2SaveParameterModel -cs isChipRescaleInNeed -fstr featureCollectionOnTraining -fste featureCollectionOnTesting -pin parameterModelInputFile" << endl;
		exit( 1 );
	}
	
    // set the number of random restarts
    BindingPredictor::GLOBAL_NRANDSTARTS = nrandstarts; // ## SS
    Motif::GLOBAL_isFeatureStat = isFeatureStat; // ## qcheng75

    BindingPar::GLOBAL_INIT_PARAMETER_MODEL_FILE_NAME =  pminputfile; // ## qcheng75

    // fix the cooperativity problem
    if (fixcooperror!=0){
		BindingPar::log10_min_wt=0;
		BindingPar::log10_max_wt=4;
		BindingPar::initMaxBindingWt=2;
		BindingPar::FIX_COOP_PROB=1;
		BindingPar::PENALTY_COEFFIENT = penaltyVal4BoundedOpt;
		BindingPar::IS_LOG_SCALE_STEP = islogscalestep;
		BindingPar::FIX_COOP_PROB=1;
		if (fixcooperror==11){
			BindingPar::log10_min_int=-6;
			BindingPar::log10_max_int=6;
		}else
		if (fixcooperror==1){
			BindingPar::log10_min_int=0;
			BindingPar::log10_max_int=6;
			BindingPar::FIX_COOP_PROB=1;
		}else if (fixcooperror==-1){
			BindingPar::log10_min_int=-6;
			BindingPar::log10_max_int=0;
			BindingPar::FIX_COOP_PROB=1;
		}if (fixcooperror==2){
			BindingPar::log10_min_int=0;
			BindingPar::log10_max_int=32000;
			BindingPar::FIX_COOP_PROB=2;
		}else if (fixcooperror==-2){
			BindingPar::log10_min_int=-32000;
			BindingPar::log10_max_int=0;
			BindingPar::FIX_COOP_PROB=2;
		}
    }
	//cout << fixcooperror << " " << BindingPar::log10_min_int << " " << BindingPar::log10_max_int  << " " << BindingPar::initInt << " " << BindingPar::FIX_COOP_PROB << endl;

	// set the options
	assert( perfOption >= 0 && perfOption <= 1 );	
	
	// read motifs
	vector< Motif > motifs;
	vector< string > motifNames;
	vector< double > background( 4, 0.25 );
 	background[ 0 ] = 0.3; background[ 1 ] = 0.2; background[ 2 ] = 0.2; background[ 3 ] = 0.3; // ## SS uncommented this
	int rval = readMotifs( motifFile, background, motifs, motifNames ); 
	assert( rval != RET_ERROR );
	int nFactors = motifs.size();
	
	// read the sequences
	vector< Sequence > allSeqs;
	vector< string > allNames;
	rval = readSequences( seqFile, allSeqs, allNames );
//assert( rval != RET_ERROR );
	int nSeqs = allSeqs.size() / nExps;
	vector< vector< Sequence > > seqs( nExps );
	vector< vector< string > > names( nExps );	
	int counter = 0;
	for ( int i = 0; i < nExps; i++ ) {
		for ( int j = 0; j < nSeqs; j++ ) {
			seqs[ i ].push_back( allSeqs[ counter ] );	
			names[ i ].push_back( allNames[ counter ] );
			counter++;
		}
	}

	// read the binding data: one row per experiment (binding of all sequences in that experiment)
	ifstream fdata( dataFile.c_str() );
    if ( !fdata ) {
        cerr << "Cannot find the binding data file " << dataFile << endl;
        exit( 1 );
    }
	vector< vector< double > > bindingData( nExps, vector< double>(nSeqs) );
	for ( int i = 0; i < nExps; i++ ) {
		for ( int j = 0; j < nSeqs; j++ ) {
			string name;
			fdata >> name;			
			if( name != names[ i ][ j ] ) { 
				cerr << "Error: " << names[ i ][ j ] << "\t" << name << endl;
				exit( 1 );
			}
			fdata >> bindingData[ i ][ j ];	
			if (chiprescaling!=0){ // ## qcheng75 rescale chip data, especially chip score is not corresponding to expected binding number
				bindingData[ i ][ j ] = (1/(1+exp(-bindingData[ i ][ j ]))-0.5)*200; //linear model
			}
		}
	}	
	
	// energy threshold
	vector< double > energyThrs( nFactors, energyThr );
	if (Motif::SITE_THRESHOLD_BY_LLRthres==1){
		for ( int i=0; i< motifs.size(); i++){
			energyThrs[i] = motifSiteThr[i];
		}
	}

    // create the interaction function
	FactorIntFunc* intFunc;
	if ( interactionOption == 0 ) {
		intFunc = new FactorIntFuncBinary( distThr, orientationEffect );
	} else if ( interactionOption == 1 ) {
		intFunc = new FactorIntFuncLinear( distThr, range, orientationEffect );		
	} else if ( interactionOption == 2 ) {
        intFunc = new FactorIntFuncPeriodic( distThr, period, spacingEffect, phi, phi0 );
    } else if ( interactionOption == 3 ) {
        intFunc = new FactorIntFuncGeometric( distThr, spacingEffect, orientationEffect );
    } else {
        cerr << "Invalid option of the interaction function" << endl; 
        exit( 1 );
    }
    
	// set the cooperativity option for training
    bool array[] = { 1, 1, 1, 1, 0 }; 
    switch ( coopOption ) {
        case 0: 
            array[2] = 0; array[3] = 0;
            break;
        case 1: 
            array[2] = 0; array[3] = 1;
            break;
        case 2: 
            array[2] = 1; array[3] = 1;
            break;
        default: 
            cerr << "Cooperativity option is invalid" << endl; exit( 1 );
    }
    BindingPredictor::parOption = vector< bool >( array, array + sizeof(array) / sizeof(*array) );
		
	//## qcheng75
	if ( !featureStatiscsTrFile.empty() ) { //## qcheng75
		ofstream featureStat( featureStatiscsTrFile.c_str() );
		if ( !featureStat ) {
			cerr << "Cannot open file " << featureStatiscsTrFile << endl;
			exit( 1 );
		}
		SeqAnnotator ann( motifs, energyThrs);
		for ( int i = 0; i < nExps; i++ ) {
			cout << i << endl;
			for ( int k = 0; k < motifs.size(); k++ ) {
				featureStat << ">>\t" << k << "\t" << motifNames[k] << "\t" << motifs[k].infoContent << "\t" << motifs[ k ].getMaxLLR() << endl;
				featureStat.flush();
				for ( int j = 0; j < nSeqs; j++ ) {
					ann.doFeatureStatistics(featureStat, names[ i ][ j ], seqs[ i ][ j ], k, bindingData[ i ][ j ]);
				}
			}
		}
		ann.clear();
		featureStat.close();
	}

	// cross validation
	if ( K > 0 ) {
		BindingCrossValidator::perfOption = perfOption;
		BindingCrossValidator cv( seqs, bindingData );		
		cv.partition( K );
// 		cv.print();
		vector< double > results; 
		cv.runThermoModel( results, motifs, energyThrs, intFunc, featureSelection, modelOption );
		for ( int k = 0; k < K; k++ ) cout << k << "\t" << results[ k ] << endl;
		cout << "Average performance = " << mean( results ) << endl;
		cout << "Standard deviation of performance = " << std_dev( results ) << endl;
		
		return 0;
	} 

    // feature selection
    if ( featureSelection ) {
        FeatureSelector::effectThr = effectThr;
        FeatureSelector selector( motifs, energyThrs, seqs, bindingData, intFunc );
        selector.select();
        vector< bool > indicators = selector.getIndicators();        
        vector< Motif > subMotifs;
        vector< double > subEnergyThrs;
        vector< string > subMotifNames;
        cout << "Feature selection:"; 
        for ( int i = 0; i < motifs.size(); i++ ) {
            if ( indicators[ i ] ) {
                subMotifs.push_back( motifs[ i ] );
                subEnergyThrs.push_back( energyThrs[ i ] );
                subMotifNames.push_back( motifNames[ i ] );
                cout << " " << subMotifNames.back();
            }
        }
        cout << endl;
        motifs = subMotifs;
        energyThrs = subEnergyThrs;
        motifNames = subMotifNames;
    }
	
	// training the model
    cout << "1) training the model " <<endl;
	BasePredictor* predictor;
	if (areAllPGiven){
		predictor = new BindingPredictor( motifs, energyThrs, seqs, bindingData, intFunc, allParaValues, seqs.size(),nFactors-seqs.size(), true );

	}else{
		if (isFreePGiven){
			vector< double > pBindingWts(seqs.size());
			for ( int i = 0; i < seqs.size(); i++ ) {
				pBindingWts[i]=fFreePMaxBinding;
				cout << fFreePMaxBinding << endl;
			}

			predictor = new BindingPredictor( motifs, energyThrs, seqs, bindingData, intFunc,pBindingWts );

		}else{
			if ( modelOption == 0 ) predictor = new BindingPredictor( motifs, energyThrs, seqs, bindingData, intFunc );
			else predictor = new StepwisePredictor( motifs, energyThrs, seqs, bindingData, intFunc );
		}
	}
	if (! isTestOnly ){
		predictor->train();
		if ( modelOption == 1 ) {
			cout << "Features chosen:";
			vector< bool > indicators = ( (StepwisePredictor*)predictor )->getIndicators();
			for ( int i = 0; i < indicators.size(); i++ ) {
				if ( indicators[ i ] ) cout << "\t" << motifNames[ i ];
			}
			cout << endl;
		}
		cout << "Estimated values of parameters:" << endl;
		predictor->getPar().print();
		if ( !paraModelFile.empty() ) {
			ofstream pmodel( paraModelFile.c_str() );
			predictor->getPar().printParaModel2Stream(pmodel);
			pmodel.close();
		}
		cout << "Training performance = " << predictor->getObj() << endl;
		if ( predictor->nOutFactors() > 0 ) {
			cout << "Analyze features:" << endl;
			vector< vector< double > > featureEffects;
			predictor->analyzeFeatures( featureEffects );
			for ( int k = 0; k < predictor->nOutFactors(); k++ )
				cout << motifNames[ k + predictor->nInFactors() ] << "\t" << featureEffects[ k ] << endl;
		}

		// print the training results: predictions for training data
		if ( !trainPredictionFile.empty() ) {
			ofstream ftrain( trainPredictionFile.c_str() );
			if ( !ftrain ) {
				cerr << "Cannot open file " << trainPredictionFile << endl;
				exit( 1 );
			}
			for ( int i = 0; i < nExps; i++ ) {
				for ( int j = 0; j < seqs[ i ].size(); j++ ) {
					double chip = bindingData[ i ][ j ];
					double predicted=predictor->predict( seqs[ i ][ j ], i );
					if (chiprescaling!=0){ // ## qcheng75 rescale chip data, especially chip score is not corresponding to expected binding number
						chip =log( (chip/200 + 0.5) / ((1 - (chip/200 + 0.5) )==0?0.00001: (1 - (chip/200 + 0.5) ))); //(1/(1+exp(-bindingData[ i ][ j ]))-0.5)*200; //linear model
						predicted =log( (predicted/200 + 0.5) / ((1 - (predicted/200 + 0.5) )==0?0.00001: (1 - (predicted/200 + 0.5) ))); //(1/(1+exp(-bindingData[ i ][ j ]))-0.5)*200; //linear model
						ftrain << names[ i ][ j ] << "\t" << chip << "\t" << predicted << "\t" << bindingData[ i ][ j ] << "\t" << predictor->predict( seqs[ i ][ j ], i ) << endl;
					}else{
						ftrain << names[ i ][ j ] << "\t" << bindingData[ i ][ j ] << "\t" << predictor->predict( seqs[ i ][ j ], i ) << endl;
					}
				}
			}
		}
	}
	// testing the model
	if ( !testSeqFile.empty() ) {
		cout << "2) testing the model " <<endl;
		 // reading test sequences
		vector< Sequence > allTestSeqs;
		vector< string > allTestNames;
		rval = readSequences( testSeqFile, allTestSeqs, allTestNames );
		int nTestSeqs = allTestSeqs.size() / nExps;
		vector< vector< Sequence > > testSeqs( nExps );
		vector< vector< string > > testNames( nExps );	
		int counter = 0;
		for ( int i = 0; i < nExps; i++ ) {
			for ( int j = 0; j < nTestSeqs; j++ ) {
				testSeqs[ i ].push_back( allTestSeqs[ counter ] );	
				testNames[ i ].push_back( allTestNames[ counter ] );
				counter++;
			}
		}		 
			 
		 // reading test data
		ifstream fTestData( testDataFile.c_str() );
		vector< vector< double > > testBindingData( nExps, vector< double >(nTestSeqs) );
		for ( int i = 0; i < nExps; i++ ) {
			for ( int j = 0; j < nTestSeqs; j++ ) {
				string name;
				fTestData >> name;			
				//cout << "assert(" << name << " == " << testNames[ i ][ j ]<<" );" << endl;
				assert( name == testNames[ i ][ j ] );
				fTestData >> testBindingData[ i ][ j ];	
				if (chiprescaling!=0){ // ## qcheng75 rescale chip data, especially chip score is not corresponding to expected binding number
					testBindingData[ i ][ j ] = (1/(1+exp(-testBindingData[ i ][ j ]))-0.5)*200; //linear model
				}
			}
		}	
			 
		//## qcheng75
		if ( !featureStatiscsTeFile.empty() ) { //## qcheng75
			ofstream featureStat( featureStatiscsTeFile.c_str() );
			if ( !featureStat ) {
				cerr << "Cannot open file " << featureStatiscsTeFile << endl;
				exit( 1 );
			}
			SeqAnnotator ann( motifs, energyThrs);
			for ( int i = 0; i < nExps; i++ ) {
				for ( int k = 0; k < motifs.size(); k++ ) {
					featureStat << ">>\t" << k << "\t" << motifNames[k] << "\t" << motifs[k].infoContent << "\t" << motifs[ k ].getMaxLLR() << endl;
					featureStat.flush();

					for ( int j = 0; j < nTestSeqs; j++ ) {
						ann.doFeatureStatistics(featureStat, testNames[ i ][ j ], testSeqs[ i ][ j ], k, testBindingData[ i ][ j ]);
					}
				}
			}
			ann.clear();
			featureStat.close();
		}

		 // testing
		double testResult;
		testResult = predictor->test( testSeqs, testBindingData, perfOption );	
		cout << "Testing performance = " << testResult << endl; 

        // print the predictions on the testing data
        if ( !testPredictionFile.empty() ) {
            ofstream ftest( testPredictionFile.c_str() );
            if ( !ftest ) {
                cerr << "Cannot open file " << testPredictionFile << endl;
                exit( 1 );
            }
            for ( int i = 0; i < nExps; i++ ) {
                for ( int j = 0; j < testSeqs[ i ].size(); j++ ) {
                	double testchip = testBindingData[ i ][ j ];
					double testpredicted=predictor->predict( testSeqs[ i ][ j ], i );
					if (chiprescaling!=0){ // ## qcheng75 rescale chip data, especially chip score is not corresponding to expected binding number
						testchip =log( (testchip/200 + 0.5) / ((1 - (testchip/200 + 0.5) )==0?0.00001: (1 - (testchip/200 + 0.5) ))); //(1/(1+exp(-bindingData[ i ][ j ]))-0.5)*200; //linear model
						testpredicted =log( (testpredicted/200 + 0.5) / ((1 - (testpredicted/200 + 0.5) )==0?0.00001: (1 - (testpredicted/200 + 0.5) ))); //(1/(1+exp(-bindingData[ i ][ j ]))-0.5)*200; //linear model
						ftest << testNames[ i ][ j ] << "\t" << testchip << "\t" << testpredicted << "\t" << testBindingData[ i ][ j ] << "\t" << predictor->predict( testSeqs[ i ][ j ], i ) << endl;
					}else{
						ftest << testNames[ i ][ j ] << "\t" << testBindingData[ i ][ j ] << "\t" << predictor->predict( testSeqs[ i ][ j ], i ) << endl;
					}
                }	
            }
        }

        //qcheng75
        for (int i=0; i<testBindingData.size(); i++){
			testBindingData[i].clear();
			//delete &(testBindingData[i]);
		}
		testBindingData.clear();
		allTestNames.clear();
		//allTestNames.clear();
		for (int i=0; i<allTestSeqs.size();i++){
			allTestSeqs[i].clear();
		}
		allTestSeqs.clear();
		for (int i=0; i<testSeqs.size();i++){
			testSeqs[i].clear();
		}
		testSeqs.clear();
		for (int i=0; i<testNames.size();i++){
			testNames[i].clear();
			//delete &(testNames[i]);
		}
		testNames.clear();
	}
	
	//qcheng75
	if ( modelOption == 0 ) {
		((BindingPredictor *)predictor)->clear();
		delete (BindingPredictor *)predictor;
	}else{
		predictor->clear();
		delete predictor;
	}

	//qcheng75
	if ( featureSelection ) {
		//to clean up memories for FeatureSelector
	}
	if ( K > 0 ) {
		//to clean up memories for BindingCrossValidator
	}
	BindingPredictor::parOption.clear();
	delete intFunc;
	for (int i=0; i<bindingData.size(); i++){
		bindingData[i].clear();
		//delete &(bindingData[i]);
	}
	bindingData.clear();
	allNames.clear();
	for (int i=0; i<allSeqs.size();i++){
		allSeqs[i].clear();
	}
	allSeqs.clear();
	for (int i=0; i<seqs.size();i++){
		seqs[i].clear();
	}
	seqs.clear();
	for (int i=0; i<names.size();i++){
		names[i].clear();
		//delete &(names[i]);
	}
	names.clear();
	energyThrs.clear();
	background.clear();
	motifNames.clear();
	for (int i=0; i<motifs.size();i++){
		motifs[i].clear();
	}
	motifs.clear();
	return 0;
}

