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
	vector<Cluster> clusters;

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
	/* calcualte the cost of the TSP in each cluster */
	for (auto cl = clusters.begin(); cl != clusters.end(); ++cl) {
		double cost = this->pbdata->travelCosts.at(0).at(cl->cluster.at(0));
		for (int i = 1; i < cl->cluster.size(); ++i) {
			cost += this->pbdata->travelCosts.at(cl->cluster.at(i-1)).at(cl->cluster.at(i));
		}
		cost += this->pbdata->travelCosts.at(cl->cluster.at(cl->cluster.size() - 1)).at(0);
		if (cl->cluster.size() > 1) {
			cl->cost = permute(cl->cluster, 0, cl->cluster.size(), cost);
		}
	}
}
double Solver::permute(vector<int> a, int l, int r, double cost)
{
	// Base case 
	double minCost = cost, delta = 0;
	if (l == r) {
		for (int i = 0; i < a.size(); ++i)
			cout << a.at(i) << "\t";
		cout << "\t\t" << cost; // << "\t" << verifyCostSequence(a);
		cout << endl;
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
