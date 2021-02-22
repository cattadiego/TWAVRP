#include "Solver.h"

Solver::Solver(PbData *pbdata)
{
	this->pbdata = pbdata;
}
Solver::~Solver()
{
}
void Solver::generateAllClusters() {
	vector<int> customers;
	for (int i = 1; i < this->pbdata->id.size(); ++i)
		customers.push_back(this->pbdata->id.at(i));

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
				clusters.push_back(cluster);
		}
	}

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

double Solver::permute(vector<int> a, int l, int r, double cost)
{
	// Base case 
	double minCost = cost, delta = 0;
	if (l == r) {
		//for (int i = 0; i < a.size(); ++i)
		//	cout << a.at(i) << "\t";
		//cout << "\t\t" << cost; // << "\t" << verifyCostSequence(a);
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
			minCost = min(minCost, permute(a, l + 1, r, cost));

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
	for (auto cl = clusters.begin(); cl != clusters.end(); ++cl) {
		double cost = this->pbdata->travelCosts.at(0).at(cl->cluster.at(0));
		for (int i = 1; i < cl->cluster.size(); ++i) {
			cost += this->pbdata->travelCosts.at(cl->cluster.at(i - 1)).at(cl->cluster.at(i));
		}
		cost += this->pbdata->travelCosts.at(cl->cluster.at(cl->cluster.size() - 1)).at(0);
		if (cl->cluster.size() > 1) {
			cl->cost = permute(cl->cluster, 0, cl->cluster.size(), cost);
		}
		else cl->cost = cost;
	}
}
void Solver::enumeration() {
	generateAllClusters();
	evaluateAllCluster();

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

	try {
		
		initializeVariablesUCVRP(env, u);
		initializeObjFunctionUCVRP(env, model, u);
		initializeCoveringCnstUCVRP(env, model, u, coveringCnst);
		
		cplex.solve();
		lb = cplex.getObjValue();
		determineClustersInSolUCVRP(cplex, u, sol);
		printSolUCVRP(sol);
		double cost = solveSeparationUCVRP(sol);
		if (cost < ub) {
			ub = cost;
			solBest = sol;
		}

		while (lb < ub) {
			IloExpr expr(env);
			int rhs = 0;
			for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
				rhs += sol.at(s).size();
				for (int c = 0; c < sol.at(s).size(); ++c) {
					expr += u.at(sol.at(s).at(c)).at(s);
				}
			}
			model.add(expr <= rhs - 2);
			nbCuts++;

			cplex.solve();
			lb = cplex.getObjValue();

			determineClustersInSolUCVRP(cplex, u, sol);
			printSolUCVRP(sol);
			double cost = solveSeparationUCVRP(sol);
			if (cost < ub) {
				ub = cost;
				solBest = sol;
			}
		}			
	}
	catch (IloException &ex) {
		cerr << "enumeration error: " << ex << endl;
		cplex.end(); model.end(); env.end();
		return;
	}
	ofstream stream; stream.open("results.dat", ios::app);
	stream << this->pbdata->instanceName << "\t" << this->pbdata->nbCustomers << "\t" << this->pbdata->nbScenarios 
		<<  "lb:" << lb << "\tub:" << ub << "\tcuts: " << nbCuts << endl;
	stream.close();
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
	for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
		sol.push_back(vector<int>());
	}
	for (int c = 0; c < u.size(); ++c) {
		for (auto var = u.at(c).begin(); var != u.at(c).end(); ++var) {
			if (cplex.getValue(var->second) > 0.99) {
				sol.at(var->first).push_back(c);
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
void Solver::printSolUCVRP(vector<vector<int>> &sol) {
	for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
		for (int c = 0; c < sol.at(s).size(); ++c) {
			cout << "s" << s << "_c" << sol.at(s).at(c) << endl;
			cout << this->clusters.at(sol.at(s).at(c)).print(s) << endl;
		}
	}
}
double Solver::solveSeparationUCVRP(vector<vector<int>> &sol) {
	
	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	try {

		vector < unordered_map<int, unordered_map<int, IloIntVar>>> x;
		vector < unordered_map<int, IloNumVar>> t;
		unordered_map<int, IloNumVar> y;

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

		for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
			for (int i = 0; i < this->pbdata->idCustomers.size(); ++i) {
				model.add(y.at(this->pbdata->idCustomers.at(i)) <= t.at(s).at(this->pbdata->idCustomers.at(i)));
				model.add(t.at(s).at(this->pbdata->idCustomers.at(i)) <= y.at(this->pbdata->idCustomers.at(i)) + this->pbdata->widthEndogenousTW.at(this->pbdata->idCustomers.at(i)));
			}
		}

		cplex.exportModel("sep.lp");
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
		return DBL_MAX;
	}
}