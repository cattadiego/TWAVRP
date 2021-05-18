#include "ClustersPair.h"

ClustersPair::ClustersPair(Cluster *c1, Cluster *c2, int pos1, int pos2) {
	if (pos1 < pos2) {
		this->c1 = c1;
		this->c2 = c2;
		this->pos1 = pos1;
		this->pos2 = pos2;
	}
	else {
		this->c1 = c2;
		this->c2 = c1;
		this->pos1 = pos2;
		this->pos2 = pos1;
	}
}
bool ClustersPair::intersection(ClustersPair &cp) {
	// if returns false the intersection is empty

	if (this->pos1 == cp.pos1) return 1;
	if (this->pos1 == cp.pos2) return 1;
	if (this->pos2 == cp.pos1) return 1;
	if (this->pos2 == cp.pos2) return 1;

	return 0;
}	