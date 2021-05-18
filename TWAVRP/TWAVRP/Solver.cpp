#include "Solver.h"

Solver::Solver(PbData *pbdata)
{
	this->pbdata = pbdata;
}
Solver::~Solver()
{
}
void Solver::generateEvaluateAllClusters() {
	generateAllCluster();
	evaluateAllCluster();
	fillClusterStructures();
}

double Solver::permute(vector<int> a, vector<int> &aBest, int l, int r, double cost, double &bestCost)
{
	// Base case 
	double minCost = cost, delta = 0;
	if (l == r) {
		//for (int i = 0; i < a.size(); ++i)
		//	cout << a.at(i) << "\t";
		//cout << "\t\t" << cost << "\t" << verifyCostSequence(a);
		//cout << endl;
	}

	else
	{
		// Permutations made 
		for (int i = l; i < r; i++)
		{

			// Swapping done 
			delta = permuteCost(a, l, i);
			cost += delta;
			swap(a[l], a[i]);
			// Recursion called 
			minCost = min(minCost, permute(a, aBest, l + 1, r, cost, bestCost));
			if (minCost < bestCost - 0.00001) {
				bestCost = minCost;
				//cout << minCost << "\t";
				for (int h = 0; h < aBest.size(); ++h) {
					//cout << a.at(h) << "\t";
					aBest.at(h) = a.at(h);
				}
				//cout << endl;
			}

			//backtrack 
			delta = permuteCost(a, l, i);
			cost += delta; 
			swap(a[l], a[i]);			
		}
	}
	return minCost;
}
double Solver::permuteCost(vector<int> a, int l, int r) {
	if (l == r) return 0;
	int prec, succ;
	if (l == 0) prec = 0;
	else prec = a.at(l - 1);
	if (r == a.size() - 1)
		succ = 0;
	else succ = a.at(r + 1);

	double delta = 0;

	if (l + 1 == r) {
		delta += this->pbdata->travelCosts.at(prec).at(a.at(r));
		delta += this->pbdata->travelCosts.at(a.at(r)).at(a.at(l));
		delta += this->pbdata->travelCosts.at(a.at(l)).at(succ);

		delta -= this->pbdata->travelCosts.at(prec).at(a.at(l));
		delta -= this->pbdata->travelCosts.at(a.at(l)).at(a.at(r));
		delta -= this->pbdata->travelCosts.at(a.at(r)).at(succ);
	}
	else {
		delta += this->pbdata->travelCosts.at(prec).at(a.at(r));
		delta += this->pbdata->travelCosts.at(a.at(r)).at(a.at(l + 1));

		delta -= this->pbdata->travelCosts.at(prec).at(a.at(l));
		delta -= this->pbdata->travelCosts.at(a.at(l)).at(a.at(l + 1));

		delta += this->pbdata->travelCosts.at(a.at(r - 1)).at(a.at(l));
		delta += this->pbdata->travelCosts.at(a.at(l)).at(succ);

		delta -= this->pbdata->travelCosts.at(a.at(r - 1)).at(a.at(r));
		delta -= this->pbdata->travelCosts.at(a.at(r)).at(succ);

	}

	return delta;
}
double Solver::verifyCostSequence(vector<int> sequence) {
	double cost = this->pbdata->travelCosts.at(0).at(sequence.at(0));
	for (int i = 1; i < sequence.size(); ++i) {
		cost += this->pbdata->travelCosts.at(sequence.at(i-1)).at(sequence.at(i));
	}
	cost += this->pbdata->travelCosts.at(sequence.at(sequence.size() - 1)).at(0);
	return cost;
}
void Solver::evaluateAllCluster() {
	/* calcualte the cost of the TSP in each cluster */
	vector<Cluster>::iterator cl = clusters.begin();
	int pos = 0;
	while (cl != clusters.end()) {
		double cost = this->pbdata->travelCosts.at(0).at(cl->cluster.at(0));
		for (int i = 1; i < cl->cluster.size(); ++i) {
			cost += this->pbdata->travelCosts.at(cl->cluster.at(i - 1)).at(cl->cluster.at(i));
		}
		cost += this->pbdata->travelCosts.at(cl->cluster.at(cl->cluster.size() - 1)).at(0);
		if (cl->cluster.size() > 1) {
			double bestCost = cost;
			cl->cost = permute(cl->cluster, cl->tsp, 0, cl->cluster.size(), cost, bestCost);
		}
		else cl->cost = cost;
		if (cl->cost > this->pbdata->ubExogenousTW.at(0) - this->pbdata->lbExogenousTW.at(0))
			cl = clusters.erase(cl);
		else {	
			bool twfeas = true;
			float t = this->pbdata->lbExogenousTW.at(0);
			t += this->pbdata->travelTimes.at(0).at(cl->tsp.at(0));
			if (t > this->pbdata->ubExogenousTW.at(cl->tsp.at(0))) {
				twfeas = false;
			}

			for (int i = 1; i < cl->cluster.size(); ++i) {
				if (t > this->pbdata->ubExogenousTW.at(cl->tsp.at(i))) {
					twfeas = false;
					break;
				}
				t = max(t, this->pbdata->lbExogenousTW.at(cl->tsp.at(i)));
				t += this->pbdata->travelTimes.at(cl->tsp.at(i - 1)).at(cl->tsp.at(i));
			}

			if (!twfeas) {
				cout << "tsp inf: " << cl->cost << endl;
				unordered_map<int, vector<int>> tempSol;
				tempSol.insert(make_pair(0, vector<int>()));
				tempSol.at(0).push_back(pos);
				double tempCost = solveSeparationUCVRP(tempSol);
				if (tempCost > this->pbdata->bigM - 1) {
					cl = clusters.erase(cl);
				}
				else {
					tempCost *= 1 / this->pbdata->scenarioProbability.at(0);
					cl->cost = tempCost;
					cout << "tsp inf: " << cl->cost << endl;
					cout << cl->print() << endl;
					cout << "*****\n";
					++cl; ++pos;
				}
			}
			else {
				++cl;
				++pos;
			}
		}
	}
}
void Solver::generateAllCluster() {
	vector<int> customers(this->pbdata->id.size() - 1);
	for (int i = 1; i < this->pbdata->id.size(); ++i)
		customers[i - 1] = this->pbdata->id.at(i);
	vector<vector<int>> baseClusters = subsets(customers);

	/* generate clusters that are feasible wrt capacity for each scenario */
	for (auto cl = baseClusters.begin(); cl != baseClusters.end(); ++cl) {
		if (cl->size() != 0) {
			Cluster cluster(this->pbdata);
			bool isfeas = true;
			for (int i = 0; i < cl->size(); ++i) {
				if (!cluster.addCustomer(cl->at(i))) {
					isfeas = false;
					break;
				}
			}
			if (isfeas)
				this->clusters.push_back(cluster);
		}
	}
}
void Solver::fillClusterStructures() {
	for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
		this->clusterFeasThatContain.push_back(unordered_map<int, vector<int>>());
	}
	for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
		for (auto it = this->pbdata->id.begin(); it != this->pbdata->id.end(); ++it) {
			this->clusterFeasThatContain.at(s).insert(make_pair(*it, vector<int>()));
		}
	}

	for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
		for (int cl = 0; cl < clusters.size(); ++cl) {
			if (clusters.at(cl).isFeasibleForScenario(s)) {
				for (int i = 0; i < clusters.at(cl).cluster.size(); ++i) {
					this->clusterFeasThatContain.at(s).at(clusters.at(cl).cluster.at(i)).push_back(cl);
				}
			}
		}
	}

	for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
		clusterFeasForScenario.push_back(vector<int>());
	}
	for (int cl = 0; cl < clusters.size(); ++cl) {
		for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
			if (clusters.at(cl).isFeasibleForScenario(s)) {
				this->clusterFeasForScenario.at(s).push_back(cl);
			}
		}
	}
}
void Solver::enumeration() {
	generateEvaluateAllClusters();

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	vector< unordered_map<int, IloIntVar > > u; //at position i contains the variables associated to the clusters
												//stored in position i in this->cluster
												// in position i there is a map which stores a 
												//pair key val where key is the scenario and val the variables
	vector< unordered_map<int, IloRange > > coveringCnst; //at position s contains the constraints associated to the scenario s
												// in position s there is a map which stores a 
												//pair key val where key is the customer and val the associated covering constraint

	vector<vector<int>> sol;					//in position s contains the vector of clusters associated to scenario s
												//the vector contains at position i contains the index of the related cluster

	vector<vector<int>> solBest;
	int nbCuts = 0;
	double lb, ub = DBL_MAX;
	double lb0, ub0;

	try {
		
		initializeVariablesUCVRP(env, u);
		initializeObjFunctionUCVRP(env, model, u);
		initializeCoveringCnstUCVRP(env, model, u, coveringCnst);
		
		cplex.solve();
		lb = cplex.getObjValue();
		lb0 = lb;
		determineClustersInSolUCVRP(cplex, u, sol);
		printSolUCVRP(sol);
		double cost = solveSeparationUCVRP(sol);
		if (cost < ub) {
			ub = cost;
			solBest = sol;
		}
		ub0 = ub;

		int iter = 0;
		while (lb < ub) {
			IloExpr expr(env);
			int rhs = 0;
			for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
				rhs += sol.at(s).size();
				for (int c = 0; c < sol.at(s).size(); ++c) {
					expr += u.at(sol.at(s).at(c)).at(s);
				}
			}
			model.add(expr <= rhs - 1);
			nbCuts++;

			cplex.solve();
			lb = cplex.getObjValue();

			determineClustersInSolUCVRP(cplex, u, sol);
			
			++iter;
			double cost = solveSeparationUCVRP(sol);
			if (cost < ub) {
				ub = cost;
				solBest = sol;
			}

			ofstream stream; stream.open(this->pbdata->instanceName + "_infos.dat", ios::app);
			stream << "******* " << iter << " ********\n";
			stream << "lb: " << lb << "\tub: " << ub << endl;
			stream << printSolUCVRP(sol) << endl;
			stream.close();


			cout << "lb: " << lb << endl;
			cout << "ub: " << ub << endl;
		}
	}
	catch (IloException &ex) {
		cerr << "enumeration error: " << ex << endl;
		cplex.end(); model.end(); env.end();
		return;
	}
	cplex.end(); model.end(); env.end();
	ofstream stream; stream.open("results.dat", ios::app);
	stream << this->pbdata->instanceName << "\t" << this->pbdata->nbCustomers << "\t" << this->pbdata->nbScenarios
		<< "\t" << calculateNbClusters()
		<< "\tlb0:" << lb0 << "\tub0:" << ub0 << "\tcuts: " << nbCuts
		<< "\tlb:" << lb << "\tub:" << ub << "\tcuts: " << nbCuts << endl;
	stream.close();
}
void Solver::enumerationCutSymmetries() {
	int tt = clock();
	cout << "gen started" << endl;
	generateEvaluateAllClusters();
	int tGen = clock() - tt;
	cout << "generation done" << endl;

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	vector< unordered_map<int, IloIntVar > > u; //at position i contains the variables associated to the clusters
												//stored in position i in this->cluster
												// in position i there is a map which stores a 
												//pair key val where key is the scenario and val the variables
												// u.at(c).at(s) equal 1 if cluster c is selected for scenario s

	vector<IloIntVar> y;						//at position i contains the variables associated to the clusters
												//stored in position i in this->cluster

	vector< unordered_map<int, IloIntVar > > z; //at position i contains the variables associated to the clusters
												//stored in position i in this->cluster
												// in position i there is a map which stores a 
												//pair key val where key is the scenario and val the variables


	vector< unordered_map<int, IloRange > > coveringCnst; //at position s contains the constraints associated to the scenario s
												// in position s there is a map which stores a 
												//pair key val where key is the customer and val the associated covering constraint

	vector< unordered_map<int, IloRange >> linkUandYCnst;

	unordered_map<int, vector<int>> sol;		//in position s contains the vector of clusters associated to scenario s
												//the vector contains at position i contains the index of the related cluster

	IloNumVar delta(env, 0.0, IloInfinity, "delta");
	unordered_map<int, unordered_map<int, IloNumVar>> deltaSS;

	unordered_map<int, unordered_map<int, IloNumVar>> deltaCC;

	unordered_map<int, unordered_map<int, IloNumVar>> gammaCC;

	unordered_map<int, vector<int>> solBest;
	int nbCutsF = 0, nbCutsI = 0, nbCutsP = 0;
	double lb, ub = this->pbdata->bigM;
	double lb0, ub0;


	UnorderedMapKeyPairOrderedInt memoryCuts(-2);
	UnorderedMapKeyPairOrderedInt memoryCutsY(-2);

	try {

		initializeVariablesUCVRP(env, u);
		initializeVariablesUCVRPCutSymmetries(env, y, z);
		//initializeVariablesUCVRPValidCuts(env, deltaSS);
		initializeObjFunctionUCVRP(env, model, u, delta);
		initializeCoveringCnstUCVRP(env, model, u, coveringCnst);
		initializeLinkUandYCnstUCVRPCutSymmetries(env, model, u, y, linkUandYCnst);

		cplex.setParam(IloCplex::Param::Threads, 1);
		cplex.solve();
		lb = cplex.getObjValue();
		lb0 = lb;
		

		int iter = 0;
		while (lb < ub - myEpsilon) {

			determineClustersInSolUCVRP(cplex, u, sol);
			double cost = solveSeparationUCVRP(sol);
		
			if (cost < ub) {
				ub = cost;
				solBest = sol;
			}
			ofstream stream; stream.open(this->pbdata->instanceName + "_" + toString(this->pbdata->nbScenarios) + "_infos.dat", ios::app);
			stream << "******* " << iter << " ********\n";
			stream << "lb: " << lb << "\tcost: " << cost << "\tub: " << ub << endl;
			stream << printSolUCVRP(sol) << endl;

			printRedMessage("lb: " + toString(lb) + "\nub: " + toString(ub) + "\n");
			
			ub0 = ub;
		
			if (cost < this->pbdata->bigM) {
				unordered_map<int, vector<int>> invSol;
				invertSolution(sol, invSol);
				nbCutsF += addFeasibilityCuts(env, model, u, z, invSol);				
			}
			else {
				nbCutsI += addInfeasibilityCuts(env, model, y, sol);
			}

			nbCutsP += addValidCutsOnY(env, model, delta, deltaCC, gammaCC, u, y, sol, iter, memoryCutsY);
			
			nbCutsP += addValidCuts(env, model, delta, u, y, sol, iter, memoryCuts);

			cplex.exportModel("model.lp");
			cplex.solve();
			lb = cplex.getObjValue();			
		}
	}
	catch (IloException &ex) {
		cerr << "enumeration error: " << ex << endl;
		cplex.end(); model.end(); env.end();
		return;
	}
	cplex.end(); model.end(); env.end();

	string message = this->pbdata->instanceName + "\t" + toString(this->pbdata->nbCustomers) +  "\t" 
		+ toString(this->pbdata->nbScenarios) + "\t" + toString(calculateNbClusters()) + "\tlb0: " + toString(lb0)
		+ "\tub0: " + toString(ub0) +  "\tlb: " + toString(lb) + "\tub: " + toString(ub) 
		+ "\tcutsI: " + toString(nbCutsI) + "\tcutsF: " + toString(nbCutsF) + "\tcustP: " + toString(nbCutsP)
		+ "\ttGen: " + toString(tGen) + "\ttTot: " + toString(clock() - tt);

	writeInStream("result.dat", message);

}
template<typename T>
void Solver::initializeVariablesUCVRP(IloEnv &env, vector< unordered_map<int, T > > &u) {
	for (int c = 0; c < this->clusters.size(); ++c) {
		u.push_back(unordered_map<int, T>());
		for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
			if (this->clusters.at(c).isFeasibleForScenario(s)) {
				string str = "u_s" + toString(s) + "_c" + toString(c);
				char* variableName = const_cast<char*>(str.c_str());
				T var;
				if(typeid(T) == typeid(IloIntVar))
					var = T(env, 0.0, 1.0, variableName);
				if (typeid(T) == typeid(IloNumVar))
					var = T(env, 0.0, IloInfinity, variableName);
				u.at(c).insert(make_pair(s, var));
			}
		}
	}
}
template<typename T>
void Solver::initializeVariablesUCVRPCutSymmetries(IloEnv &env, vector<T> &y, vector< unordered_map<int, T > > &z) {
	for (int c = 0; c < this->clusters.size(); ++c) {
		string str = "y_" + toString(c);
		char* variableName = const_cast<char*>(str.c_str());
		T var;
		if (typeid(T) == typeid(IloIntVar))
			var = T(env, 0.0, 1.0, variableName);
		if (typeid(T) == typeid(IloNumVar))
			var = T(env, 0.0, IloInfinity, variableName);
		y.push_back(var);

		z.push_back(unordered_map<int, T>());
		int k = 1;
		for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
			if (this->clusters.at(c).isFeasibleForScenario(s)) {
				string str = "z_" + toString(c) + "_" + toString(k);
				char* variableName = const_cast<char*>(str.c_str());
				T var;
				if (typeid(T) == typeid(IloIntVar))
					var = T(env, 0.0, 1.0, variableName);
				if (typeid(T) == typeid(IloNumVar))
					var = T(env, 0.0, IloInfinity, variableName);
				z.at(c).insert(make_pair(k, var));
				++k;
			}
		}
	}
}
void Solver::initializeVariablesUCVRPValidCuts(IloEnv &env, unordered_map<int, unordered_map<int, IloNumVar > > &lambdaSS) {
	for (int i = 0; i < this->pbdata->nbScenarios; ++i) {
		lambdaSS.insert(make_pair(i, unordered_map<int, IloNumVar >()));
		for (int j = i+1; j < this->pbdata->nbScenarios; ++j) {
			string str = "lambdaSS_" + toString(i) + "_" + toString(j);
			char* variableName = const_cast<char*>(str.c_str());
			IloNumVar var(env, 0.0, IloInfinity, variableName);
			lambdaSS.at(i).insert(make_pair(j, var));
		}
	}
}
template<typename T>
void Solver::initializeLinkUandYCnstUCVRPCutSymmetries(IloEnv &env, IloModel &model, vector< unordered_map<int, T > > &u, 
	vector<IloIntVar> &y, vector< unordered_map<int, IloRange > > &linkUandYCnst) {
	for (int c = 0; c < this->clusters.size(); ++c) {
		linkUandYCnst.push_back(unordered_map<int, IloRange >());
		for (auto s = u.at(c).begin(); s != u.at(c).end(); ++s) {
			IloExpr expr(env);
			string str = "link_" + toString(c) + "_" + toString(s->first);
			char* constraintName = const_cast<char*>(str.c_str());
			expr += s->second - y.at(c);
			IloRange constraint(env, -IloInfinity, expr, 0, constraintName);
			linkUandYCnst.at(c).insert(make_pair(s->first, constraint));
			model.add(constraint);
		}
	}
}
template<typename T>
void Solver::initializeObjFunctionUCVRP(IloEnv &env, IloModel &model, vector< unordered_map<int, T > > &u) {
	IloExpr expr(env);
	for (int c = 0; c < this->clusters.size(); ++c) {
		for (auto var = u.at(c).begin(); var != u.at(c).end(); ++var) {
			expr += this->pbdata->scenarioProbability.at(var->first) * this->clusters.at(c).cost * var->second;
		}
	}
	model.add(IloMinimize(env, expr));
}
template<typename T>
void Solver::initializeObjFunctionUCVRP(IloEnv &env, IloModel &model, vector< unordered_map<int, T > > &u, IloNumVar &delta) {
	IloExpr expr(env);
	for (int c = 0; c < this->clusters.size(); ++c) {
		for (auto var = u.at(c).begin(); var != u.at(c).end(); ++var) {
			expr += this->pbdata->scenarioProbability.at(var->first) * this->clusters.at(c).cost * var->second;
		}
	}
	expr += delta;
	model.add(IloMinimize(env, expr));
}
template<typename T>
void Solver::initializeCoveringCnstUCVRP(IloEnv &env, IloModel &model, vector< unordered_map<int, T > > &u, vector< unordered_map<int, IloRange > > &coveringCnst) {
	for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
		coveringCnst.push_back(unordered_map<int, IloRange >());
		for (int i = 0; i < this->pbdata->idCustomers.size(); ++i) {
			IloExpr expr(env);
			for (auto it = this->clusterFeasThatContain.at(s).at(this->pbdata->idCustomers.at(i)).begin(); it != this->clusterFeasThatContain.at(s).at(this->pbdata->idCustomers.at(i)).end(); ++it) {
				expr += u.at(*it).at(s);
			}
			string str = "covering_s" + toString(s) + "_c" + toString(this->pbdata->idCustomers.at(i));
			char* constraintName = const_cast<char*>(str.c_str());
			IloRange constraint(env, 1, expr, IloInfinity, constraintName);
			coveringCnst.at(s).insert(make_pair(this->pbdata->idCustomers.at(i), constraint));
			model.add(constraint);
		}
	}
}
template<typename T>
void Solver::determineClustersInSolUCVRP(IloCplex &cplex, vector< unordered_map<int, T > > &u, vector<vector<int>> &sol) {
	sol.clear();
	sol = vector<vector<int>>(this->pbdata->nbScenarios, vector<int>());
	/*for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
		sol.push_back(vector<int>());
	}*/
	for (int c = 0; c < u.size(); ++c) {
		for (auto var = u.at(c).begin(); var != u.at(c).end(); ++var) {
			if (cplex.getValue(var->second) > 0.99) {
				sol.at(var->first).push_back(c);
			}
		}
	}
}
template<typename T>
void Solver::determineClustersInSolUCVRP(IloCplex &cplex, vector< unordered_map<int, T > > &u, unordered_map<int, vector<int>> &sol) {
	sol.clear();
	for (int c = 0; c < u.size(); ++c) {
		for (auto var = u.at(c).begin(); var != u.at(c).end(); ++var) {
			if (cplex.getValue(var->second) > 0.99) {
				try {
					sol.at(var->first).push_back(c);
				}
				catch (out_of_range ex) {
					sol.insert(make_pair(var->first, vector<int>(1, c)));
				}
			}
		}
	}
}
template<typename T>
void Solver::printSolUCVRP(IloCplex &cplex, vector< unordered_map<int, T > > &u) {
	for (int c = 0; c < u.size(); ++c) {
		for (auto var = u.at(c).begin(); var != u.at(c).end(); ++var) {
			if (cplex.getValue(var->second) > 0.99) {
				cout << var->second.getName() << endl;
				cout << this->clusters.at(c).print(var->first) << endl;
			}
		}
	}
}
string Solver::printSolUCVRP(vector<vector<int>> &sol) {
	string str = "";
	for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
		for (int c = 0; c < sol.at(s).size(); ++c) {
			str += "s" + toString(s) + "_c" + toString(sol.at(s).at(c)) + "\n";
			str += this->clusters.at(sol.at(s).at(c)).print(s) + "\n";
		}
	}
	return str;
}
string Solver::printSolUCVRP(unordered_map<int, vector<int>> &sol) {

	map<int, vector<int>> ordered(sol.begin(), sol.end());
	
	string str = "";
	for (auto s = ordered.begin(); s != ordered.end(); ++s) {
		for (int c = 0; c < s->second.size(); ++c) {
			str += "s" + toString(s->first) + "_c" + toString(s->second.at(c)) + "\n";
			str += this->clusters.at(s->second.at(c)).print(s->first) + "\n";
		}
	}
	return str;
}
double Solver::solveSeparationUCVRP(vector<vector<int>> &sol) {

	//the first index of sol is the scenario, the second the position of the cluster in this->clusters

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	try {

		vector < unordered_map<int, unordered_map<int, IloIntVar>>> x;	// the first index is the scenario, then the two point indexes
		vector < unordered_map<int, IloNumVar>> t;						// the first index is the scenario, then the point index
		unordered_map<int, IloNumVar> y;								// the first index is the point index

		for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
			x.push_back(unordered_map<int, unordered_map<int, IloIntVar>>());
			t.push_back(unordered_map<int, IloNumVar>());
			for (int i = 0; i < this->pbdata->id.size(); ++i)
				x.at(s).insert(make_pair(this->pbdata->id.at(i), unordered_map<int, IloIntVar>()));
		}

		for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
			for (int i = 0; i < this->pbdata->idCustomers.size(); ++i) {
				string str = "x_" + toString(s) + "_" + toString(0) + "_" + toString(this->pbdata->idCustomers.at(i));
				char* variableName = const_cast<char*>(str.c_str());
				IloIntVar var(env, 0.0, 1.0, variableName);
				x.at(s).at(0).insert(make_pair(this->pbdata->idCustomers.at(i), var));

				str = "x_" + toString(s) + "_" + toString(this->pbdata->idCustomers.at(i)) + "_" + toString(this->pbdata->nbPoints);
				variableName = const_cast<char*>(str.c_str());
				var = IloIntVar(env, 0.0, 1.0, variableName);
				x.at(s).at(this->pbdata->idCustomers.at(i)).insert(make_pair(this->pbdata->nbPoints, var));
			}
			for (int i = 0; i < this->pbdata->id.size(); ++i) {
				string str = "t_" + toString(s) + "_" + toString(this->pbdata->id.at(i));
				char* variableName = const_cast<char*>(str.c_str());
				IloNumVar var(env, this->pbdata->lbExogenousTW.at(this->pbdata->id.at(i)), this->pbdata->ubExogenousTW.at(this->pbdata->id.at(i)), variableName);
				t.at(s).insert(make_pair(this->pbdata->id.at(i), var));
			}
			string str = "t_" + toString(s) + "_" + toString(this->pbdata->nbPoints);
			char* variableName = const_cast<char*>(str.c_str());
			IloNumVar var(env, 0.0, IloInfinity, variableName);
			t.at(s).insert(make_pair(this->pbdata->nbPoints, var));
		}

		for (int i = 0; i < this->pbdata->id.size(); ++i) {
			string str = "y_" + toString(this->pbdata->id.at(i));
			char* variableName = const_cast<char*>(str.c_str());
			IloNumVar var(env, this->pbdata->lbExogenousTW.at(this->pbdata->id.at(i)), this->pbdata->ubExogenousTW.at(this->pbdata->id.at(i)), variableName);
			y.insert(make_pair(this->pbdata->id.at(i), var));
		}
		string str = "y_" + toString(this->pbdata->nbPoints);
		char* variableName = const_cast<char*>(str.c_str());
		IloNumVar var(env, this->pbdata->lbExogenousTW.at(0), this->pbdata->ubExogenousTW.at(0), variableName);
		y.insert(make_pair(this->pbdata->nbPoints, var));

		for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
			for (int c = 0; c < sol.at(s).size(); ++c) {
				for (int i = 0; i < this->clusters.at(sol.at(s).at(c)).cluster.size(); ++i) {
					for (int j = i + 1; j < this->clusters.at(sol.at(s).at(c)).cluster.size(); ++j) {
						string str = "x_" + toString(s) + "_" + toString(this->clusters.at(sol.at(s).at(c)).cluster.at(i));
						str += "_" + toString(this->clusters.at(sol.at(s).at(c)).cluster.at(j));
						char* variableName = const_cast<char*>(str.c_str());
						IloIntVar var(env, 0.0, 1.0, variableName);
						x.at(s).at(this->clusters.at(sol.at(s).at(c)).cluster.at(i))
							.insert(make_pair(this->clusters.at(sol.at(s).at(c)).cluster.at(j), var));

						str = "x_" + toString(s) + "_" + toString(this->clusters.at(sol.at(s).at(c)).cluster.at(j));
						str += "_" + toString(this->clusters.at(sol.at(s).at(c)).cluster.at(i));
						variableName = const_cast<char*>(str.c_str());
						var = IloIntVar(env, 0.0, 1.0, variableName);
						x.at(s).at(this->clusters.at(sol.at(s).at(c)).cluster.at(j))
							.insert(make_pair(this->clusters.at(sol.at(s).at(c)).cluster.at(i), var));
					}
				}
			}
		}

		IloExpr expr(env);
		for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
			for (auto it = x.at(s).begin(); it != x.at(s).end(); ++it) {
				for (auto itj = it->second.begin(); itj != it->second.end(); ++itj) {
					int arrival = itj->first == this->pbdata->nbPoints ? 0 : itj->first;
					expr += this->pbdata->scenarioProbability.at(s) * this->pbdata->travelCosts.at(it->first).at(arrival) * itj->second;
				}
			}
		}
		model.add(IloMinimize(env, expr));

		/* flow out of the depot */
		for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
			IloExpr expr(env);
			for (auto it = x.at(s).at(0).begin(); it != x.at(s).at(0).end(); ++it) {
				expr += it->second;
			}
			model.add(expr == sol.at(s).size());
		}

		/* flow back into the depot */
		for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
			IloExpr expr(env);
			for (int i = 0; i < this->pbdata->idCustomers.size(); ++i) {
				expr += x.at(s).at(this->pbdata->idCustomers.at(i)).at(this->pbdata->nbPoints);
			}
			model.add(expr == sol.at(s).size());
		}


		/* flow conservation at nodes */
		for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
			for (int i = 0; i < this->pbdata->idCustomers.size(); ++i) {
				IloExpr expr(env);
				for (auto it = x.at(s).at(this->pbdata->idCustomers.at(i)).begin(); it != x.at(s).at(this->pbdata->idCustomers.at(i)).end(); ++it) {
					expr += it->second;
				}
				IloExpr out(env);
				for (int j = 0; j < this->pbdata->id.size(); ++j) {
					try {
						out += x.at(s).at(this->pbdata->id.at(j)).at(this->pbdata->idCustomers.at(i));
					}
					catch (out_of_range ex) {

					}
				}
				model.add(expr == 1);
				model.add(out == 1);
			}
		}

		/* scheduling of visits */
		int bigM = this->pbdata->ubExogenousTW.at(0);
		for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
			for (int i = 0; i < this->pbdata->id.size(); ++i) {
				for (auto it = x.at(s).at(this->pbdata->id.at(i)).begin(); it != x.at(s).at(this->pbdata->id.at(i)).end(); ++it) {
					int arrival = it->first == this->pbdata->nbPoints ? 0 : it->first;
					model.add(t.at(s).at(this->pbdata->id.at(i)) + it->second * this->pbdata->travelTimes.at(this->pbdata->id.at(i)).at(arrival)
						- t.at(s).at(it->first) - bigM * (1 - it->second) <= 0);
				}
			}
		}

		/* endogenous tw determination */
		for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
			for (int i = 0; i < this->pbdata->idCustomers.size(); ++i) {
				model.add(y.at(this->pbdata->idCustomers.at(i)) <= t.at(s).at(this->pbdata->idCustomers.at(i)));
				model.add(t.at(s).at(this->pbdata->idCustomers.at(i)) <= y.at(this->pbdata->idCustomers.at(i)) + this->pbdata->widthEndogenousTW.at(this->pbdata->idCustomers.at(i)));
			}
		}

		IloConstraintArray lazyCnst(env);
		int nbCnst = 0;
		for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
			for (int c = 0; c < sol.at(s).size(); ++c) {
				vector<vector<int>> subSets = subsets(this->clusters.at(sol.at(s).at(c)).cluster);
				for (int ss = 0; ss < subSets.size(); ++ss) {
					if (subSets.at(ss).size() > 1) {
						IloExpr expr(env);
						for (int i = 0; i < subSets.at(ss).size(); ++i) {
							for (int j = i + 1; j < subSets.at(ss).size(); ++j) {
								expr += x.at(s).at(subSets.at(ss).at(i)).at(subSets.at(ss).at(j));
								expr += x.at(s).at(subSets.at(ss).at(j)).at(subSets.at(ss).at(i));
							}
						}
						string str = "subtour_elim_" + toString(nbCnst);
						char* constraintName = const_cast<char*>(str.c_str());
						IloRange constraint(env, 0, expr, subSets.at(ss).size() - 1, constraintName);
						//lazyCnst.add(constraint);
						model.add(constraint);
						++nbCnst;
					}
				}
			}
		}

		//cplex.addLazyConstraints(lazyCnst);
		//cplex.addUserCuts(lazyCnst);

		//cplex.exportModel("sep.lp");
		cplex.solve();
		return cplex.getObjValue();
		cout << cplex.getObjValue() << endl;

		for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
			for (auto it = x.at(s).begin(); it != x.at(s).end(); ++it) {
				for (auto it1 = it->second.begin(); it1 != it->second.end(); ++it1) {
					if (cplex.getValue(it1->second) > 0.99) {
						cout << it1->second.getName() << endl;
						cout << cplex.getValue(t.at(s).at(it->first)) << endl;
						cout << cplex.getValue(t.at(s).at(it1->first)) << endl;
						cout << "----" << endl;
					}
				}
			}
		}
	}
	catch (IloException &ex) {
		cerr << "separation error: " << ex << endl;
		cplex.end(); model.end(); env.end();
		return this->pbdata->bigM;
	}
}
double Solver::solveSeparationUCVRP(unordered_map<int, vector<int>> &sol){

	//the first index of sol is the scenario, the second the position of the cluster in this->clusters

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	try {

		unordered_map <int, unordered_map<int, unordered_map<int, IloIntVar>>> x;
		unordered_map <int, unordered_map<int, IloNumVar>> t;						// the first index is the scenario, then the point index
		unordered_map<int, IloNumVar> y;											// the first index is the point index

		unordered_map<int, vector<IloNumVar>> c;									// the first index is the scenario

		for (auto s = sol.begin(); s != sol.end(); ++s) {
			x.insert(make_pair(s->first, unordered_map<int, unordered_map<int, IloIntVar>>()));
			t.insert(make_pair(s->first, unordered_map<int, IloNumVar>()));
			c.insert(make_pair(s->first, vector<IloNumVar>()));
		}


		for (auto s = sol.begin(); s != sol.end(); ++s) {
			x.at(s->first).insert(make_pair(0, unordered_map<int, IloIntVar>())); //vars x^s_0i
			for (int cl = 0; cl < s->second.size(); ++cl) {
				string str = "c_" + toString(s->first) + "_" + toString(s->second.at(cl));
				char* variableName = const_cast<char*>(str.c_str());
				IloNumVar varc(env, 0.0, IloInfinity, variableName);
				c.at(s->first).push_back(varc);

				IloExpr expr(env);
				for (int i = 0; i < this->clusters.at(s->second.at(cl)).cluster.size(); ++i) {
					x.at(s->first).insert(make_pair(this->clusters.at(s->second.at(cl)).cluster.at(i), unordered_map<int, IloIntVar>())); //vars x^s_ij
					
					string str = "x_" + toString(s->first) + "_" + toString(0) + "_" + toString(this->clusters.at(s->second.at(cl)).cluster.at(i));
					char* variableName = const_cast<char*>(str.c_str());
					IloIntVar var(env, 0.0, 1.0, variableName);
					x.at(s->first).at(0).insert(make_pair(this->clusters.at(s->second.at(cl)).cluster.at(i), var));
					expr += this->pbdata->travelCosts.at(0).at(this->clusters.at(s->second.at(cl)).cluster.at(i)) * var;

					for (int j = 0; j < this->clusters.at(s->second.at(cl)).cluster.size(); ++j) {
						if (i != j) {
							str = "x_" + toString(s->first) + "_" + toString(this->clusters.at(s->second.at(cl)).cluster.at(i)) + "_" + toString(this->clusters.at(s->second.at(cl)).cluster.at(j));
							variableName = const_cast<char*>(str.c_str());
							var = IloIntVar(env, 0.0, 1.0, variableName);
							x.at(s->first).at(this->clusters.at(s->second.at(cl)).cluster.at(i)).insert(make_pair(this->clusters.at(s->second.at(cl)).cluster.at(j), var));
							expr += this->pbdata->travelCosts.at(this->clusters.at(s->second.at(cl)).cluster.at(i)).at(this->clusters.at(s->second.at(cl)).cluster.at(j)) * var;
						}
					}

					str = "x_" + toString(s->first) + "_" + toString(this->clusters.at(s->second.at(cl)).cluster.at(i)) + "_" + toString(this->pbdata->nbPoints);
					variableName = const_cast<char*>(str.c_str());
					var = IloIntVar(env, 0.0, 1.0, variableName);
					x.at(s->first).at(this->clusters.at(s->second.at(cl)).cluster.at(i)).insert(make_pair(this->pbdata->nbPoints, var));
					expr += this->pbdata->travelCosts.at(this->clusters.at(s->second.at(cl)).cluster.at(i)).at(0) * var;
				}

				model.add(expr <= varc);
			}
		}

		for (auto s = sol.begin(); s != sol.end(); ++s) {			
			t.insert(make_pair(s->first, unordered_map<int, IloNumVar>()));
			string str = "t_" + toString(s->first) + "_" + toString(0);
			char* variableName = const_cast<char*>(str.c_str());
			IloNumVar var(env, this->pbdata->lbExogenousTW.at(this->pbdata->id.at(0)), this->pbdata->ubExogenousTW.at(this->pbdata->id.at(0)), variableName);
			t.at(s->first).insert(make_pair(0, var));
			
			for (int cl = 0; cl < s->second.size(); ++cl) {
				for (int i = 0; i < this->clusters.at(s->second.at(cl)).cluster.size(); ++i) {
					str = "t_" + toString(s->first) + "_" + toString(this->clusters.at(s->second.at(cl)).cluster.at(i));
					variableName = const_cast<char*>(str.c_str());
					var = IloNumVar(env, this->pbdata->lbExogenousTW.at(this->pbdata->id.at(this->clusters.at(s->second.at(cl)).cluster.at(i))), this->pbdata->ubExogenousTW.at(this->pbdata->id.at(this->clusters.at(s->second.at(cl)).cluster.at(i))), variableName);
					t.at(s->first).insert(make_pair(this->clusters.at(s->second.at(cl)).cluster.at(i), var));
				}
			}

			str = "t_" + toString(s->first) + "_" + toString(this->pbdata->nbPoints);
			variableName = const_cast<char*>(str.c_str());
			var = IloNumVar(env, this->pbdata->lbExogenousTW.at(this->pbdata->id.at(0)), this->pbdata->ubExogenousTW.at(this->pbdata->id.at(0)), variableName);
			t.at(s->first).insert(make_pair(this->pbdata->nbPoints, var));
		}

		for (auto s = sol.begin(); s != sol.end(); ++s) {
			for (int cl = 0; cl < s->second.size(); ++cl) {
				for (int i = 0; i < this->clusters.at(s->second.at(cl)).cluster.size(); ++i) {
					string str = "y_" + toString(this->clusters.at(s->second.at(cl)).cluster.at(i));
					char* variableName = const_cast<char*>(str.c_str());
					IloNumVar var(env, this->pbdata->lbExogenousTW.at(this->clusters.at(s->second.at(cl)).cluster.at(i)), this->pbdata->ubExogenousTW.at(this->clusters.at(s->second.at(cl)).cluster.at(i)), variableName);
					y.insert(make_pair(this->clusters.at(s->second.at(cl)).cluster.at(i), var));
				}
			}
		}

		IloExpr expr(env);
		/*for (auto s = x.begin(); s != x.end(); ++s) {
			for (auto it = s->second.begin(); it != s->second.end(); ++it) {
				for (auto it1 = it->second.begin(); it1 != it->second.end(); ++it1) {
					expr += this->pbdata->scenarioProbability.at(s->first) 
						* this->pbdata->travelCosts.at(it->first != this->pbdata->nbPoints ? it->first : 0)
						.at(it1->first != this->pbdata->nbPoints ? it1->first : 0) * it1->second;
				}
			}
		}*/
		for (auto s = c.begin(); s != c.end(); ++s) {
			for (int i = 0; i < s->second.size(); ++i) {				
				expr += this->pbdata->scenarioProbability.at(s->first) * this->clusters.at(sol.at(s->first).at(i)).coeffObjFunction * s->second.at(i);
			}
		}
		model.add(IloMinimize(env, expr));

		/* flow out of the depot */
		for (auto s = x.begin(); s != x.end(); ++s) {
			IloExpr expr(env);
			for (auto it = s->second.at(0).begin(); it != s->second.at(0).end(); ++it) {
				expr += it->second;
			}
			model.add(expr == sol.at(s->first).size());
		}

		/* flow back into the depot */
		for (auto s = x.begin(); s != x.end(); ++s) {
			IloExpr expr(env);
			for (auto it = s->second.begin(); it != s->second.end(); ++it) {
				try {
					expr += it->second.at(this->pbdata->nbPoints);
				} catch(out_of_range ex) {

				}
			}
			model.add(expr == sol.at(s->first).size());
		}


		/* flow conservation at nodes */
		for (auto s = x.begin(); s != x.end(); ++s) {
			for(auto it = s->second.begin(); it != s->second.end(); ++it) {
				if (it->first != 0) {
					IloExpr out(env);
					for (auto it1 = it->second.begin(); it1 != it->second.end(); ++it1) {
						out += it1->second;
					}
					IloExpr in(env);
					for (auto it2 = s->second.begin(); it2 != s->second.end(); ++it2) {
						try {
							in += it2->second.at(it->first);
						}
						catch (out_of_range ex) {

						}
					}
					model.add(out == 1);
					model.add(in == 1);
				}
				
			}
		}

		/* scheduling of visits */
		int bigM = this->pbdata->ubExogenousTW.at(0);
		for (auto s = x.begin(); s != x.end(); ++s) {
			for (auto it = s->second.begin(); it != s->second.end(); ++it) {
				for (auto it1 = it->second.begin(); it1 != it->second.end(); ++it1) {
					int arrival = it1->first == this->pbdata->nbPoints ? 0 : it1->first;
					model.add(t.at(s->first).at(it->first) + it1->second * this->pbdata->travelTimes.at(it->first).at(arrival)
						- t.at(s->first).at(it1->first) - bigM * (1 - it1->second) <= 0);
				}
			}
		}

		/* endogenous tw determination */
		for (auto s = sol.begin(); s != sol.end(); ++s) {
			for (int cl = 0; cl < s->second.size(); ++cl) {
				for (int i = 0; i < this->clusters.at(s->second.at(cl)).cluster.size(); ++i) {
					model.add(y.at(this->clusters.at(s->second.at(cl)).cluster.at(i)) <= t.at(s->first).at(this->clusters.at(s->second.at(cl)).cluster.at(i)));
					model.add(t.at(s->first).at(this->clusters.at(s->second.at(cl)).cluster.at(i)) <= 
						y.at(this->clusters.at(s->second.at(cl)).cluster.at(i)) + this->pbdata->widthEndogenousTW.at(this->clusters.at(s->second.at(cl)).cluster.at(i)));
				}
			}
		}

		IloConstraintArray lazyCnst(env);
		int nbCnst = 0;
		for (auto s = sol.begin(); s != sol.end(); ++s) {
			for (int c = 0; c < sol.at(s->first).size(); ++c) {
				vector<vector<int>> subSets = subsets(this->clusters.at(sol.at(s->first).at(c)).cluster);
				for (int ss = 0; ss < subSets.size(); ++ss) {
					if (subSets.at(ss).size() > 1) {
						IloExpr expr(env);
						for (int i = 0; i < subSets.at(ss).size(); ++i) {
							for (int j = i + 1; j < subSets.at(ss).size(); ++j) {
								expr += x.at(s->first).at(subSets.at(ss).at(i)).at(subSets.at(ss).at(j));
								expr += x.at(s->first).at(subSets.at(ss).at(j)).at(subSets.at(ss).at(i));
							}
						}
						string str = "subtour_elim_" + toString(nbCnst);
						char* constraintName = const_cast<char*>(str.c_str());
						IloRange constraint(env, 0, expr, subSets.at(ss).size() - 1, constraintName);
						//lazyCnst.add(constraint);
						model.add(constraint);
						++nbCnst;
					}
				}
			}
		}

		cplex.setParam(IloCplex::Param::Threads, 1);
		cplex.solve();
		double obj = cplex.getObjValue();
		cout << "obj: " <<  obj << endl;

		for (auto s = c.begin(); s != c.end(); ++s) {
			for (auto i = 0; i < s->second.size(); ++i) {
				this->clusters.at(sol.at(s->first).at(i)).costInSol = cplex.getValue(c.at(s->first).at(i));
			}
		}

		for (auto s = x.begin(); s != x.end(); ++s) {
			for (auto it = s->second.at(0).begin(); it != s->second.at(0).end(); ++it) {
				if (cplex.getValue(it->second) > 1 - myEpsilon) {
					
					int clusterIndex;
					for (int c = 0; c < sol.at(s->first).size(); ++c) {
						if (find(this->clusters.at(sol.at(s->first).at(c)).cluster.begin(), this->clusters.at(sol.at(s->first).at(c)).cluster.end(), it->first) != this->clusters.at(sol.at(s->first).at(c)).cluster.end()) {
							clusterIndex = sol.at(s->first).at(c);
							break;
						}
					}

					this->clusters.at(clusterIndex).tspInSol.clear();
					this->clusters.at(clusterIndex).tspInSol.push_back(it->first);


					int curr = it->first;
					do {

						auto it1 = s->second.at(curr).begin();
						while (cplex.getValue(it1->second) < myEpsilon) {
							it1++;
						}

						curr = it1->first;
						if (curr != this->pbdata->nbPoints) this->clusters.at(clusterIndex).tspInSol.push_back(curr);

					} while (curr != this->pbdata->nbPoints);
				}
			}
		}


		cplex.end(); model.end(); env.end();
		return obj;

	}
	catch (IloException &ex) {
		cerr << "separation error: " << ex << endl;
		cplex.end(); model.end(); env.end();
		return this->pbdata->bigM;
	}
}
int Solver::calculateNbClusters() {
	int nbClusters = 0;
	for (int c = 0; c < this->clusters.size(); ++c) {
		for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
			if (this->clusters.at(c).isFeasibleForScenario(s))
				nbClusters++;
		}
	}
	return nbClusters;
}
template<typename T>
int Solver::addFeasibilityCuts(IloEnv &env, IloModel &model, vector< unordered_map<int, T > > &u, vector< unordered_map<int, T > > &z, unordered_map<int, vector<int>> &invSol){
	for (auto itc = invSol.begin(); itc != invSol.end(); ++itc) {
		IloExpr expr(env);
		for (auto its = u.at(itc->first).begin(); its != u.at(itc->first).end(); ++its) {
			expr += its->second;
		}

		IloExpr expr1(env);
		for (auto its = z.at(itc->first).begin(); its != z.at(itc->first).end(); ++its) {
			expr -= its->first * z.at(itc->first).at(its->first);
			expr1 += z.at(itc->first).at(its->first);
		}
		model.add(expr == 0);
		model.add(expr1 <= 1);
	}

	IloExpr expr(env);
	for (auto itc = invSol.begin(); itc != invSol.end(); ++itc) {
		expr += z.at(itc->first).at(itc->second.size());
	}
	model.add(expr <= (int)invSol.size() - 1);
	return 1;
}
int Solver::addInfeasibilityCuts(IloEnv &env, IloModel &model, vector<IloIntVar> &y, unordered_map<int, vector<int>> &sol) {
	unordered_map<int, IloIntVar> tempy;
	for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
		for (int c = 0; c < sol.at(s).size(); ++c) {
			tempy.insert(make_pair(sol.at(s).at(c), y.at(sol.at(s).at(c))));
		}
	}
	IloExpr expr(env);
	for (auto it = tempy.begin(); it != tempy.end(); ++it) {
		expr += it->second;
	}
	int rhs = tempy.size();
	model.add(expr <= (int)tempy.size() - 1);
	return 1;
}
int Solver::addValidCutsOnY(IloEnv &env, IloModel &model, IloNumVar &delta, vector<IloIntVar> &y, unordered_map<int, vector<int>> &sol, int iter) {
	int nbCutsP = 0;
	ofstream cuts(this->pbdata->instanceName + "_" + toString(this->pbdata->nbScenarios) + "_cuts.dat", ios::app);
	cuts << "*****(y) " << iter << " (y)****\n";
	for (auto s = sol.begin(); s != sol.end(); ++s) {
		for (int c = 0; c < s->second.size(); ++c) {
			auto s1 = s;
			for (s1++; s1 != sol.end(); ++s1) {
				for (int c1 = 0; c1 < s1->second.size(); ++c1) {
					if (s->second.at(c) != s1->second.at(c1)) {
						double tspCost = (this->clusters.at(s->second.at(c)).cost + this->clusters.at(s1->second.at(c1)).cost) / (this->pbdata->nbScenarios);
						unordered_map<int, vector<int>> tempSol;
						tempSol.insert(make_pair(s->first, vector<int>(1, s->second.at(c))));
						if(s->first != s1->first)
							tempSol.insert(make_pair(s1->first, vector<int>(1, s1->second.at(c1))));
						else tempSol.at(s->first).push_back(s1->second.at(c1));
						double twaCost = solveSeparationUCVRP(tempSol);
						if (twaCost - tspCost > myEpsilon) {
							IloExpr expr(env);
							expr += (twaCost - tspCost) * y.at(s->second.at(c)) + (twaCost - tspCost) * y.at(s1->second.at(c1));
							expr -= delta + (twaCost - tspCost);
							model.add(expr <= 0);
							++nbCutsP;

							cuts << "----\n";
							cuts << this->clusters.at(s->second.at(c)).print() << endl;
							cuts << this->clusters.at(s1->second.at(c1)).print() << endl;
							cuts << "twacost: " << twaCost << "\ttspcost: " << tspCost << endl;
							
						}
					}
				}
			}
		}
	}

	cuts.close();
	return nbCutsP;
}
int Solver::addValidCutsOnY(IloEnv &env, IloModel &model, IloNumVar &delta, unordered_map<int, unordered_map<int, IloNumVar>> &deltaCC, unordered_map<int, unordered_map<int, IloNumVar>> &gammaCC, vector< unordered_map<int, IloIntVar > > &u, vector<IloIntVar> &y, unordered_map<int, vector<int>> &sol, int iter, UnorderedMapKeyPairOrderedInt& memoryCuts) {
	int nbCutsP = 0;
	
	vector<vector<ClustersPair>> setpairs;
	
	for (auto s = sol.begin(); s != sol.end(); ++s) {
		for (int c = 0; c < s->second.size(); ++c) {
			auto s1 = s;
			for (s1++; s1 != sol.end(); ++s1) {
				for (int c1 = 0; c1 < s1->second.size(); ++c1) {
					if (s->second.at(c) != s1->second.at(c1)) {
						if (memoryCuts.get(s->second.at(c), s1->second.at(c1)) == -2) { // cut not found
							double tspCost = (this->clusters.at(s->second.at(c)).cost + this->clusters.at(s1->second.at(c1)).cost) / (this->pbdata->nbScenarios);
							unordered_map<int, vector<int>> tempSol;
							tempSol.insert(make_pair(s->first, vector<int>(1, s->second.at(c))));
							if (s->first != s1->first)
								tempSol.insert(make_pair(s1->first, vector<int>(1, s1->second.at(c1))));
							else tempSol.at(s->first).push_back(s1->second.at(c1));
							double twaCost = solveSeparationUCVRP(tempSol);
							memoryCuts.insert(s->second.at(c), s1->second.at(c1), twaCost - tspCost);
							if (twaCost - tspCost > myEpsilon) {
								IloExpr expr(env);
								expr += (twaCost - tspCost) * y.at(s->second.at(c)) + (twaCost - tspCost) * y.at(s1->second.at(c1));
								expr -= delta + (twaCost - tspCost);
								model.add(expr <= 0);
								++nbCutsP;

								ofstream cuts(this->pbdata->instanceName + "_" + toString(this->pbdata->nbScenarios) + "_cuts.dat", ios::app);
								cuts << "----\n";
								cuts << this->clusters.at(s->second.at(c)).print() << endl;
								cuts << this->clusters.at(s1->second.at(c1)).print() << endl;
								cuts << "twacost: " << twaCost << "\ttspcost: " << tspCost << endl;
								cuts.close();

								determineDisjointPairsOfClusters(env, model, tspCost, twaCost, deltaCC, y, setpairs, s->second, s1->second, c, c1);

								determineMultipleOccurencesOfPairs(env, model, delta, tspCost, twaCost, gammaCC, u, s->second, s1->second, c, c1);

							}

							auto s2 = s1;
							for (s2++; s2 != sol.end(); ++s2) {
								for (int c2 = 0; c2 < s2->second.size(); ++c2) {
									if (s->second.at(c) != s2->second.at(c2) && s1->second.at(c1) != s2->second.at(c2)) {
										tspCost += (this->clusters.at(s2->second.at(c2)).cost / this->pbdata->nbScenarios);
										if (s->first != s2->first && s1->first != s2->first)
											tempSol.insert(make_pair(s2->first, vector<int>(1, s2->second.at(c2))));
										else tempSol.at(s2->first).push_back(s2->second.at(c2));
										double twaCost = solveSeparationUCVRP(tempSol);

										if (twaCost - tspCost > myEpsilon) {

											IloExpr expr(env);
											expr += (twaCost - tspCost) * (y.at(s->second.at(c)) + y.at(s1->second.at(c1)) + y.at(s2->second.at(c2)));
											expr -= delta + 2 * (twaCost - tspCost);
											model.add(expr <= 0);
											++nbCutsP;

											ofstream cuts(this->pbdata->instanceName + "_" + toString(this->pbdata->nbScenarios) + "_cuts.dat", ios::app);
											cuts << "---- triplets ---- \n";
											cuts << this->clusters.at(s->second.at(c)).print() << endl;
											cuts << this->clusters.at(s1->second.at(c1)).print() << endl;
											cuts << this->clusters.at(s2->second.at(c2)).print() << endl;
											cuts << "twacost: " << twaCost << "\ttspcost: " << tspCost << endl;
											cuts.close();
										}
									}
								}
							}
						}
						else if (memoryCuts.get(s->second.at(c), s1->second.at(c1)) > myEpsilon) {
							ClustersPair cp(&this->clusters.at(s->second.at(c)), &this->clusters.at(s1->second.at(c1)), s->second.at(c), s1->second.at(c1));
							insertClusterPairInSets(cp, setpairs);
						}
					}
				}
			}
		}
	}

	ofstream cuts(this->pbdata->instanceName + "_" + toString(this->pbdata->nbScenarios) + "_cuts.dat", ios::app);
	cuts << "*****(y) " << iter << " (y)****\n";
	cuts << "pairs cut\n";

	for (auto it = setpairs.begin(); it != setpairs.end(); ++it) {
		IloExpr expr(env);
		cuts << "new set\n";
		for (auto itpair = it->begin(); itpair != it->end(); ++itpair) {
			cuts << "\t" << this->clusters.at(itpair->pos1).print() << endl;
			cuts << "\t" << this->clusters.at(itpair->pos2).print() << endl;
			cuts << "---\n";
			expr += deltaCC.at(itpair->pos1).at(itpair->pos2);
		}
		model.add(expr <= delta);
		++nbCutsP;
	}


	cuts.close();
	return nbCutsP;
}
void Solver::determineMultipleOccurencesOfPairs(IloEnv &env, IloModel &model, IloNumVar &delta, double tspCost, double twaCost, unordered_map<int, unordered_map<int, IloNumVar>> &gammaCC, vector< unordered_map<int, IloIntVar > > &u, vector<int> &s, vector<int> s1, int c, int c1) {
	int id = s.at(c);
	int id1 = s1.at(c1);

	if (id > id1) {
		id = s1.at(c1);
		id1 = s.at(c);
	}


	if (gammaCC.find(id) == gammaCC.end()) {
		gammaCC.insert(make_pair(id, unordered_map<int, IloNumVar>()));
	}

	if (gammaCC.find(id) != gammaCC.end() && gammaCC.at(id).find(id1) == gammaCC.at(id).end()) {
		string str = "gamma_" + toString(id) + "_" + toString(id1);
		char* variableName = const_cast<char*>(str.c_str());
		IloIntVar vargammaCC(env, 0, 1, variableName);
		gammaCC.at(id).insert(make_pair(id1, vargammaCC));
	}

	IloExpr v(env);
	for (auto it = u.at(id).begin(); it != u.at(id).end(); ++it) {
		v += it->second;
	}
	IloExpr v1(env);
	for (auto it = u.at(id1).begin(); it != u.at(id1).end(); ++it) {
		v1 += it->second;
	}

	model.add(v - v1 <= this->pbdata->nbScenarios * (1 - gammaCC.at(id).at(id1)));
	model.add(v1 - v <= this->pbdata->nbScenarios * gammaCC.at(id).at(id1));
	model.add((twaCost - tspCost) * v1 <= delta + this->pbdata->nbScenarios * gammaCC.at(id).at(id1));
	model.add((twaCost - tspCost) * v <= delta + this->pbdata->nbScenarios * (1 - gammaCC.at(id).at(id1)));

	ofstream cuts(this->pbdata->instanceName + "_" + toString(this->pbdata->nbScenarios) + "_cuts.dat", ios::app);
	cuts << "multiple counts\n";
	cuts << this->clusters.at(s.at(c)).print() << endl;
	cuts << this->clusters.at(s1.at(c1)).print() << endl;
	cuts << "twacost: " << twaCost << "\ttspcost: " << tspCost << endl;
	cout << "----\n";
	cuts.close();

}
void Solver::determineDisjointPairsOfClusters(IloEnv &env, IloModel &model, double tspCost, double twaCost, unordered_map<int, unordered_map<int, IloNumVar>> &deltaCC, vector<IloIntVar> &y, vector<vector<ClustersPair>> &setpairs, vector<int> &s, vector<int> s1, int c, int c1){
	ClustersPair cp(&this->clusters.at(s.at(c)), &this->clusters.at(s1.at(c1)), s.at(c), s1.at(c1));

	int id = s.at(c);
	int id1 = s1.at(c1);

	if (id > id1) {
		id = s1.at(c1);
		id1 = s.at(c);
	}

	if (deltaCC.find(id) == deltaCC.end()) {
		deltaCC.insert(make_pair(id, unordered_map<int, IloNumVar>()));
	}

	if (deltaCC.find(id) != deltaCC.end() && deltaCC.at(id).find(id1) == deltaCC.at(id).end()) {
		string str = "delta_" + toString(id) + "_" + toString(id1);
		char* variableName = const_cast<char*>(str.c_str());
		IloNumVar vardeltaCC(env, 0.0, IloInfinity, variableName);
		deltaCC.at(id).insert(make_pair(id1, vardeltaCC));
	}

	IloExpr expr(env);
	expr += (twaCost - tspCost) * y.at(s.at(c)) + (twaCost - tspCost) * y.at(s1.at(c1));
	expr -= deltaCC.at(id).at(id1) + (twaCost - tspCost);
	model.add(expr <= 0);
	
	ofstream cuts(this->pbdata->instanceName + "_" + toString(this->pbdata->nbScenarios) + "_cuts.dat", ios::app);
	cuts << "----\n";
	cuts << this->clusters.at(s.at(c)).print() << endl;
	cuts << this->clusters.at(s1.at(c1)).print() << endl;
	cuts << "twacost: " << twaCost << "\ttspcost: " << tspCost << endl;
	cuts.close();

	insertClusterPairInSets(cp, setpairs);
}
int Solver::addValidCutsOnU(IloEnv &env, IloModel &model, IloNumVar &delta, unordered_map<int, unordered_map<int, IloNumVar > > &deltaSS, vector< unordered_map<int, IloIntVar > > &u, unordered_map<int, vector<int>> &sol) {
	int nbCutsP = 0;
	for (auto s = sol.begin(); s != sol.end(); ++s) {
		for (int c = 0; c < s->second.size(); ++c) {
			for (auto s1 = sol.begin(); s1 != sol.end(); ++s1) {
				for (int c1 = 0; c1 < s1->second.size(); ++c1) {					
					double tspCost = (this->clusters.at(s->second.at(c)).cost + this->clusters.at(s1->second.at(c1)).cost) / (this->pbdata->nbScenarios);
					unordered_map<int, vector<int>> tempSol;
					tempSol.insert(make_pair(s->first, vector<int>(1, s->second.at(c))));
					tempSol.insert(make_pair(s1->first, vector<int>(1, s1->second.at(c1))));
					double twaCost = solveSeparationUCVRP(tempSol);
					if (twaCost - tspCost > myEpsilon) {
						IloExpr expr(env);
						expr += (twaCost - tspCost) * u.at(s->second.at(c)).at(s->first) 
							+ (twaCost - tspCost) * u.at(s1->second.at(c1)).at(s1->first);
						(s->first < s1->first) ? expr -= deltaSS.at(s->first).at(s1->first) + (twaCost - tspCost) 
							: expr -= deltaSS.at(s1->first).at(s->first) + (twaCost - tspCost);
						model.add(expr <= 0);
						++nbCutsP;
					}					
				}
			}
		}
	}
	return nbCutsP;
}
int Solver::addValidCuts(IloEnv &env, IloModel &model, IloNumVar &delta, vector< unordered_map<int, IloIntVar > > &u, vector<IloIntVar> &y, unordered_map<int, vector<int>> &sol, int iter, UnorderedMapKeyPairOrderedInt& memoryCuts) {
	HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(hConsole, 10);
	cout << "********** " << "cuts in " << " **********\n";
	SetConsoleTextAttribute(hConsole, 15);
	int nbCutsP = 0;

	ofstream cuts(this->pbdata->instanceName + "_" + toString(this->pbdata->nbScenarios) + "_cuts.dat", ios::app);
	cuts << "*******(g) " << iter << " (g)********\n";

	for (auto s = sol.begin(); s != sol.end(); ++s) {
		for (int c = 0; c < s->second.size(); ++c) {
			auto s1 = s;
			for (++s1; s1 != sol.end(); ++s1) {
				for (int c1 = 0; c1 < s1->second.size(); ++c1) {

					HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
					SetConsoleTextAttribute(hConsole, 10);
					cout << "********** " << nbCutsP << " **********\n";
					cout << s->first << "\t" << s1->first << "\t" << c << "\t" << c1 << endl;
					cout << "********************\n";
					SetConsoleTextAttribute(hConsole, 15);						

					if (s->second.at(c) != s1->second.at(c1)) {
						if (memoryCuts.get(s->second.at(c), s1->second.at(c1)) == -2) {
							unordered_map<int, vector<int>> tempSol;
							tempSol.insert(make_pair(s->first, vector<int>(1, s->second.at(c))));
							if (s->first != s1->first)
								tempSol.insert(make_pair(s1->first, vector<int>(1, s1->second.at(c1))));
							else tempSol.at(s->first).push_back(s1->second.at(c1));

							this->clusters.at(s->second.at(c)).coeffObjFunction = 1;
							this->clusters.at(s1->second.at(c1)).coeffObjFunction = 0;

							double tspCost = this->pbdata->scenarioProbability.at(s->first) * this->clusters.at(s->second.at(c)).cost;
							double twaCost = solveSeparationUCVRP(tempSol);

							double R0 = twaCost - tspCost;
							this->clusters.at(s->second.at(c)).coeffObjFunction = 0;
							this->clusters.at(s1->second.at(c1)).coeffObjFunction = 1;

							tspCost = this->pbdata->scenarioProbability.at(s1->first) * this->clusters.at(s1->second.at(c1)).cost;
							twaCost = solveSeparationUCVRP(tempSol);
							double R1 = twaCost - tspCost;
							memoryCuts.insert(s->second.at(c), s1->second.at(c1), R0 + R1);
							if (R0 + R1 > myEpsilon) {
								IloExpr expr(env);
								for (auto scen = u.at(s->second.at(c)).begin(); scen != u.at(s->second.at(c)).end(); ++scen) {
									expr += R0 * scen->second;
								}
								for (auto scen = u.at(s1->second.at(c1)).begin(); scen != u.at(s1->second.at(c1)).end(); ++scen) {
									expr += R1 * scen->second;
								}
								expr -= delta;
								double M = max(R0, R1);
								model.add(expr <= M * (2 - y.at(s->second.at(c)) - y.at(s1->second.at(c1))));
								++nbCutsP;
								cuts << "---\n";
								cuts << this->clusters.at(s->second.at(c)).print() << "\n" << this->clusters.at(s1->second.at(c1)).print() << endl;
								cuts << "M: " << M << "\t(" << R0 << ", " << R1 << ")" << endl;
							}
							this->clusters.at(s->second.at(c)).coeffObjFunction = 1;
						}
					}
				}
			}
		}
	}
	
	cuts.close();


	hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(hConsole, 10);
	cout << "********** " << "cuts out" << " **********\n";
	SetConsoleTextAttribute(hConsole, 15);

	return nbCutsP;
}
void Solver::invertSolution(unordered_map<int, vector<int>> &sol, unordered_map<int, vector<int>> &invSol) {
	for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
		for (int c = 0; c < sol.at(s).size(); ++c) {
			try {
				invSol.at(sol.at(s).at(c)).push_back(s);
			}
			catch (out_of_range ex) {
				vector<int> temp; temp.push_back(s);
				invSol.insert(make_pair(sol.at(s).at(c), temp));
			}
		}
	}
}
void Solver::test() {
	
	//c42
	Cluster c(this->pbdata);
	//c.cluster = { 1,	3,	5,	6,	7 };
	//c.tsp = { 6,	3,	1,	5,	7 };
	//c.cost = 9.68;
	
	c.addCustomer(1);
	c.addCustomer(3);
	c.addCustomer(5);
	c.addCustomer(6);
	c.addCustomer(7);

	//c9
	Cluster c1(this->pbdata);
	//c1.cluster = { 1,	2,	3,	4,	5,	6,	7,	8,	10 };
	//c1.tsp = { 7,	5,	1,	3,	6,	2,	4,	10,	8 };
	//c1.cost = 15.83;

	c1.addCustomer(1);
	c1.addCustomer(2);
	c1.addCustomer(3);
	c1.addCustomer(4);
	c1.addCustomer(5);
	c1.addCustomer(6);
	c1.addCustomer(7);
	c1.addCustomer(10);
	c1.addCustomer(8);
	//c965
	Cluster c2(this->pbdata);
	//c2.cluster = { 1,	2,	3,	5,	6,	7,	9 };
	//c2.tsp = { 7,	5,	1,	3,	6,	2,	9 };
	//c2.cost = 10.73;

	c2.addCustomer(2);
	c2.addCustomer(4);
	c2.addCustomer(8);
	c2.addCustomer(9);

	int r = c.intersection(c1);
	int r1 = c.intersection(c2);
	int r2 = c1.intersection(c2);


	this->clusters.push_back(c);
	this->clusters.push_back(c1);
	this->clusters.push_back(c2);

	unordered_map<int, vector<int>> sol;
	sol.insert(make_pair(0, vector<int>()));
	sol.insert(make_pair(1, vector<int>()));
	sol.insert(make_pair(2, vector<int>()));
	sol.at(0).push_back(0);
	sol.at(1).push_back(1);
	sol.at(2).push_back(2);


	cout << "twa cost: " << solveSeparationUCVRP(sol) << endl;
	//cout << "tsp cost: " << (c1.cost + c.cost) / this->pbdata->nbScenarios << endl;

	for (auto s = sol.begin(); s != sol.end(); ++s) {
		for (int c = 0; c < s->second.size(); ++c) {
			cout << "twa cost: " << this->clusters.at(s->second.at(c)).costInSol 
				<< "\ttsp cost: "<< this->clusters.at(s->second.at(c)).cost << endl;
		}
	}
}
void Solver::insertClusterPairInSets(ClustersPair &cp, vector<vector<ClustersPair>> &setpairs) {
	bool pairInserted = false;
	for (auto itset = setpairs.begin(); itset != setpairs.end(); ++itset) {
		bool intersection = false;
		for (auto itpair = itset->begin(); itpair != itset->end(); ++itpair) {
			ClustersPair ccp = *itpair;
			if (ccp.intersection(cp)) {
				intersection = true;
				break;
			}
		}
		if (!intersection) {
			itset->push_back(cp);
			pairInserted = true;
		}
	}
	if (!pairInserted) {
		vector<ClustersPair> v; v.push_back(cp);
		setpairs.push_back(v);
	}
}
