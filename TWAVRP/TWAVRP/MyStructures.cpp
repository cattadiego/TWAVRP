#include "MyStructures.h"

SerialMatrix::SerialMatrix(int rows, int columns) : rows(rows), columns(columns){

}




SerialIntMatrix::SerialIntMatrix(int rows, int columns) : SerialMatrix(rows, columns){
	this->matrix = vector<int>(rows * columns);
}

int SerialIntMatrix::getValue(int row, int col) {
	return this->matrix[col * SerialMatrix::columns + row];
}


UnorderedMapKeyPairOrderedInt::UnorderedMapKeyPairOrderedInt(float not_found_value) {
	this->not_found_value = not_found_value;
}

bool UnorderedMapKeyPairOrderedInt::insert(int k1, int k2, float val) {
	if (k1 <= k2) {
		if (this->my_map.find(k1) == my_map.end()) {
			my_map.insert(make_pair(k1, unordered_map<int, float>()));
			my_map.at(k1).insert(make_pair(k2, val));
		}
		else {
			if (this->my_map.at(k1).find(k2) == my_map.at(k1).end()) {
				my_map.at(k1).insert(make_pair(k2, val));
				return true;
			}
			return false;
		}
	}
	else insert(k2, k1, val);
	return false;
}
float UnorderedMapKeyPairOrderedInt::get(int k1, int k2) {
	if (k1 <= k2) {
		if (this->my_map.find(k1) == my_map.end())
			return this->not_found_value;
		if(this->my_map.at(k1).find(k2) == my_map.at(k1).end())
			return this->not_found_value;
		return this->my_map.at(k1).at(k2);
	}
	else {
		return get(k2, k1);
	}
}
float UnorderedMapKeyPairOrderedInt::getNotFoundValue() {
	return this->getNotFoundValue;
}