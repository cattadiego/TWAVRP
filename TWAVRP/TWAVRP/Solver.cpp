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


	try {
		for (int c = 0; c < this->clusters.size(); ++c) {
			u.push_back(unordered_map<int, IloIntVar>());
			for (int s = 0; s < this->pbdata->nbScenarios; ++s) {
				if (this->clusters.at(c).isFeasibleForScenario(s)) {
					string str = "u_s" + toString(s) + "_c" + toString(c);
					char* variableName = const_cast<char*>(str.c_str());
					IloIntVar var(env, 0.0, 1.0, variableName);
					u.at(c).insert(make_pair(s, var));
				}
			}
		}

		IloExpr expr(env);
		for (int c = 0; c < this->clusters.size(); ++c) {
			for (auto var = u.at(c).begin(); var != u.at(c).end(); ++var) {
				expr += this->pbdata->scenarioProbability.at(var->first) * this->clusters.at(c).cost * var->second;
			}
		}
		model.add(IloMinimize(env, expr));

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

		cplex.exportModel("model.lp");
		cplex.solve();

		cout << cplex.getObjValue() << endl;


		for (int c = 0; c < u.size(); ++c) {
			for (auto var = u.at(c).begin(); var != u.at(c).end(); ++var) {
				if (cplex.getValue(var->second) > 0.99) {
					cout << var->second.getName() << endl;
					cout << this->clusters.at(c).print(var->first) << endl;
				}
			}
		}

			
	}
	catch (IloException &ex) {
		cerr << "enumeration error: " << ex << endl;
		cplex.end(); model.end(); env.end();
		return;
	}
}