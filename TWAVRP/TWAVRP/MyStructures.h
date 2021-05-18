#pragma once

#include <vector>
#include <unordered_map>

using namespace std;

class SerialMatrix {
protected:
	int columns;
	int rows;
public:
	
	SerialMatrix(int rows, int columns);
};

class SerialIntMatrix : public SerialMatrix {
private:
	vector<int> matrix;
public:
	SerialIntMatrix(int rows, int columns);
	int getValue(int row, int col);
};

class UnorderedMapKeyPairOrderedInt {
private:
	unordered_map<int, unordered_map<int, float>> my_map;
	float not_found_value;
public:

	UnorderedMapKeyPairOrderedInt(float not_found_value);

	bool insert(int k1, int k2, float val);
	float get(int k1, int k2);
	float getNotFoundValue();
};